/* file gbmaker.c: main program for the grain boundary maker program, 
* for creating models of grain boundaries
********************************************************************/

#include <stdio.h>	/*  ANSI-C libraries */
#include <time.h>

#include <vector>
#include <algorithm>

#include "floatdef.h"
#include "stemtypes_fftw3.h"
#include "data_containers.h"
//#include "memory_fftw3.h"	/* memory allocation routines */
#include "readparams.h"
#include "matrixlib.h"
#include "fileio_fftw3.h"

#define NAME_BUF_LEN 64
#define CRYSTALLINE 0
#define AMORPHOUS 1
#define SPECIAL_GRAIN 2

MULS *muls = new MULS();
grainBox *grains = NULL;
int nGrains = 0;
superCellBox superCell;

/* table of atomic radii:                                                      
* Al,P,S?
* Taken from http://www.physlink.com/Reference/PeriodicTable.cfm
* Needs to be completed one day ...,
*/
float_tt atRadf[]  = {0.79f,0.49f,
2.05f,1.40f,1.17f,0.91f,0.75f,0.65f,0.57f,0.51f,
2.23f,1.72f,0.00f,1.46f,0.00f,0.00f,0.97f,0.88f,
2.77f,2.23f,2.09f,2.00f,1.92f,1.85f,1.79f,1.72f,1.67f,1.62f,1.57f,1.53f,1.81f,1.52f,1.33f,1.22f,1.12f,1.03f};
float_tt covRadf[] = {0.32f,0.93f,
1.23f,0.90f,0.82f,0.77f,0.75f,0.73f,0.72f,0.71f,
1.54f,1.36f,0.00f,1.90f,0.00f,0.00f,0.99f,0.98f,
2.03f,1.74f,1.44f,1.32f,1.63f,1.18f,1.17f,1.17f,1.16f,1.15f,1.17f,1.25f,1.26f,1.22f,1.20f,1.16f,1.14f,1.12f};
float_tt *covRad,*atRad;

int readParams(char *datFileName);
void showData();
void makeSuperCell();
void makeAmorphous();
void makeSpecial(int distPlotFlag);
int removeVacancies(std::vector<atom> atoms,int natoms);
void computeCenterofMass();
void makeDistrPlot(std::vector<atom> atoms,float_tt ax);
float_tt xDistrFun1(float_tt xcenter,float_tt width);
float_tt xDistrFun2(float_tt xcenter,float_tt width1,float_tt width2);


int main(int argc, char *argv[]) {
	char outFileName[64],datFileName[64],moldyName[128];
	char *str;
	// atom *atomPtr;
	int g,nstart,j,ak,count;
	FILE *fp;
	float_tt charge;
	int moldyFlag = 0;     // suppress creation of ..._moldy.in file
	int distPlotFlag = 0;  // suppress creation of disList.dat

	/* Let's set the radii of certain elements by hand:
	*/
	//	atRad  = (float_tt *)realloc(atRadf, 109*sizeof(float_tt));
	// covRad = (float_tt *)realloc(covRadf,109*sizeof(float_tt));
	atRad  = (float_tt *)malloc(109*sizeof(float_tt));
	covRad = (float_tt *)malloc(109*sizeof(float_tt));
	memcpy(atRad,  atRadf,36*sizeof(float_tt));
	memcpy(covRad,covRadf,36*sizeof(float_tt));
	atRad[39-1]=2.27f; covRad[39-1]=2.62f;  // Y
	atRad[57-1]=2.74f; covRad[57-1]=1.69f;  // La
	atRad[71-1]=2.25f; covRad[71-1]=1.56f;  // Lu

	if (argc < 2)
		sprintf(datFileName,"gb.gbm");
	else
		strcpy(datFileName,argv[1]);
	// read a flag:
	if (argc > 2) {
		if (strncmp(argv[2],"-m",2) == 0) {
			moldyFlag = 1;	// also write a moldy file!
			printf("Saving moldy input file!\n");
		}
	}

	muls->nCellX = 1;
	muls->nCellY = 1;
	muls->nCellZ = 1;
	muls->ctiltx = 0;
	muls->ctilty = 0;
	muls->ctiltz = 0;
	superCell.natoms = 0;

	if (!readParams(datFileName))
		exit(0);
	if (nGrains > 0) {
		// loop, until we find crystalline grains:
		for (g=0; g<nGrains; g++) if (grains[g].amorphFlag == CRYSTALLINE) break;
		/* make the crystalline part of the super cell */
		if (g<nGrains) makeSuperCell();


		/* if there is also an amorphous part, then add it now, and also write a
		* MD input file, so that the amorphous phase atoms can be relaxed
		*/
		for (g=0; g<nGrains; g++) if (grains[g].amorphFlag != CRYSTALLINE) break;
		if (g<nGrains) {
			if (moldyFlag) {
				sprintf(moldyName,"%s",datFileName);
				moldyName[strlen(datFileName)-4] = '\0';
				strcat(moldyName,"_moldy.in");
				if ((fp=fopen(moldyName,"w")) == NULL) {
					printf("Could not open moldy input file %s!\n",moldyName);
					exit(0);
				}
			}
			else {
				fp = NULL;
			}
			if (nGrains > 1) {
				writeFrameWork(fp,superCell);
				computeCenterofMass();
				nstart = superCell.natoms;
				switch (grains[g].amorphFlag) {
				case 1: makeAmorphous();
					break;
				case 2: makeSpecial(distPlotFlag);
					break;
				}	
				writeAmorphous(fp,superCell,nstart,superCell.natoms);
			}
			else {
				switch (grains[g].amorphFlag) {
				case 1: makeAmorphous();
					break;
				case 2: makeSpecial(distPlotFlag);
					break;
				}	
				writeAmorphous(fp,superCell,0,superCell.natoms);
			}
			if (moldyFlag)	fclose(fp);


		}	 
	}
	if (0) { // (moldyFlag) {
		///////////////////////////////////////////////////////////////
		// write Moldy input file, without presence of amorphous phase:
		sprintf(moldyName,"%s",datFileName);
		moldyName[strlen(datFileName)-4] = '\0';
		strcat(moldyName,"_moldy.in");
		if ((fp=fopen(moldyName,"w")) == NULL) {
			printf("Could not open moldy input file %s!\n",moldyName);
			exit(0);
		}
		// writeFrameWork(fp,superCell);
		// computeCenterofMass();
		// superCell2Moldy(fp,superCell);
		fclose(fp);
	} // end of: if moldyFlag ...
	strcpy(outFileName,datFileName);

	// atomPtr = readUnitCell(&natoms,fileName,&muls);
	// writePDB(atomPtr,nat  /* reset the input file and advance to the next crystal row */

	str = strchr(outFileName,'.');
	if (str == NULL) str=outFileName+strlen(outFileName);
	sprintf(str,".cfg");
	muls->ax = (float)superCell.ax;
	muls->by = (float)superCell.by;
	muls->c	= (float)superCell.cz;

	superCell.natoms = removeVacancies(superCell.atoms,superCell.natoms);

	printf("will write cfg file to %s\n",outFileName);
	writeCFG(superCell.atoms, superCell.natoms, outFileName, *muls);
	printf("wrote cfg file to %s\n",outFileName);

	/**************************************************************
	* find the charge for the Y-atoms, in order to remain neutral:
	*/
	charge = 0.0;
	if (0) {
		for (ak=0;ak<muls->atomKinds;ak++) {
			count =0;
			for (j=0;j<superCell.natoms;j++) {
				if (muls->Znums[ak] == superCell.atoms[j].Znum) count++;
			}
			printf("Z=%3d: %d\n",muls->Znums[ak],count);
			switch (muls->Znums[ak]) {
			case  7: charge += static_cast<float_tt>(count*(-3.0)); break;
			case  8: charge += static_cast<float_tt>(count*(-2.0));  break;
			case  38: charge += static_cast<float_tt>(count*(2.0));  break;
			case  22: charge += static_cast<float_tt>(count*(4.0));  break;
			case 14: charge += static_cast<float_tt>(count*  4.0);  break;
			}	 
		}
	}
	else {
		for (j=0;j<superCell.natoms;j++) {
			charge += superCell.atoms[j].q*superCell.atoms[j].occ;
		}
	}
	// printf("Total charge: %g, i.e. %g %s\n",charge,charge,(charge > 0) ? "holes" : "electrons");
	printf("Total charge: %g",charge);
	if (charge > 0) printf(", i.e. %g holes\n",charge);
	if (charge < 0) printf(", i.e. %g electrons\n",-charge);

	delete(muls);
	return 0;
}


int removeVacancies(std::vector<atom> atoms,int natoms) {
	int natomsFinal;
	int i,i2,j,jz;
	float_tt totOcc,lastOcc, choice;
	long idum = -(long)time(NULL);
	int printLevel = 1;
	

	// printf("Time: %d\n",time(NULL));
	natomsFinal = natoms;
	std::sort(atoms.begin(), atoms.end(), atomCompareZYX());
	//qsort((void *)atoms,natoms,sizeof(atom),atomCompareZYX);

	for(jz = 0,i=0;i<natoms;i++) {
		if (atoms[i].Znum > 0) {
			totOcc = atoms[i].occ;
			for (j=i+1;j<natoms;j++) {
				// if there is anothe ratom that comes close to within 0.1*sqrt(3) A we will increase 
				// the total occupany and the counter j.
				if ((fabs(atoms[i].x-atoms[j].x) < 0.1) && (fabs(atoms[i].y-atoms[j].y) < 0.1) && (fabs(atoms[i].z-atoms[j].z) < 0.1)) {
					totOcc += atoms[j].occ;
				}
				else break;
			} // j-loop
			// if we encountered atoms in the same position, or the occupancy of the current atom is not 1, then
			// do something about it:
			// printf("%d: tocc = %g\n",i,totOcc);
			if ((totOcc < 1) || (j > i+1)) { // found atoms at equal positions!
				// ran1 returns a uniform random deviate between 0.0 and 1.0 exclusive of the endpoint values. 
				// 
				// if the total occupancy is less than 1 -> make sure we keep this
				// if the total occupancy is greater than 1 (unphysical) -> rescale all partial occupancies!
				if (totOcc < 1.0) choice = ran1(&idum);   
				else choice = totOcc*ran1(&idum);
				lastOcc = 0;
				for (i2=i;i2<j;i2++) {
					// if choice does not match the current atom:
					// choice will never be 0 or 1(*totOcc) 
					if ((choice <lastOcc) || (choice >=lastOcc+atoms[i2].occ)) {
						atoms[i2].Znum =  0;  // vacancy
						jz++;					
					}
					lastOcc += atoms[i2].occ;
					// if we keep this atom!
					atoms[i2].occ =  1.0;  // to avoid this atom being reduced in occupancy again when reading the cfg file
				}
			}
			i = j-1;
		}
	}
	if ((jz > 0 ) &&(printLevel)) printf("Removed %d atoms because of occupancies < 1 or multiple atoms in the same place\n",jz);

	// end of managing atoms at equal position
	//////////////////////////////////////////////////////////////////////////////

	natomsFinal = natoms-jz;
	// printf("%d atoms left\n",natomsFinal);
	// We now need to move all the atoms with zero Znum to the very end.
	if (jz > 0) {
		std::sort(atoms.begin(), atoms.end(), atomCompareZnum());
		//qsort((void *)atoms,natoms,sizeof(atom),atomCompareZnum);
		// for (i=0;i<natomsFinal;i++) printf("%3d: %d (%.1f %.1f %.1f): occ=%.1f\n",i,atoms[i].Znum,atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].occ);
	}
	return natomsFinal;
}

void computeCenterofMass() {
	int i;
	float_tt cmx=0.0,cmy=0.0,cmz=0.0;

	for (i=0;i<superCell.natoms;i++) {
		cmx += superCell.atoms[i].x;
		cmy += superCell.atoms[i].y;
		cmz += superCell.atoms[i].z;
	}
	superCell.cmx = cmx/(superCell.natoms*superCell.ax);
	superCell.cmy = cmy/(superCell.natoms*superCell.by);
	superCell.cmz = cmz/(superCell.natoms*superCell.cz);
	printf("Center of mass of frame work: (%g, %g, %g)\n",superCell.cmx,superCell.cmy,superCell.cmz);
}


/* This function returns 1 on success, and 0 on failure
*/
int readParams(char *datFileName) {
	int printFlag = 1;
	char title[128], parStr[128], *str;
	int gCount,i,Nkind;  
	char unitCellFile[64];
	std::vector<atom> tempCell;

	if (!parOpen(datFileName)) {
		printf("Could not open data input file %s\n",datFileName);
		return 0;
	}
	resetParamFile();
	while(readparam("crystal:",parStr,0)) nGrains++;  
	resetParamFile();
	while(readparam("amorph:",parStr,0)) nGrains++;  
	resetParamFile();
	while(readparam("special:",parStr,0)) nGrains++;  
	printf("Found data for %d grain(s) (crystalline frame work and amorphous)\n",nGrains);
	if (nGrains == 0) return 0;

	grains = (grainBox *)malloc(nGrains*sizeof(grainBox));
	/* Now we will loop through all the grains and lok for the necessary 
	* data for each grain 
	*/
	if (readparam("box:",parStr,1)) {
		sscanf(parStr,"%lf %lf %lf",&(superCell.ax),&(superCell.by),&(superCell.cz));
	}
	else {
		printf("Size of super cell box not defined - exit!\n");
		exit(0);
	} 

	/* reset the input file and advance to the next crystal row */
	resetParamFile();
	/*
	readparam("crystal:",parStr,0);
	grains[0].name = (char *)malloc(NAME_BUF_LEN);
	strcpy(grains[0].name,parStr);
	grains[0].nplanes = 0;
	grains[0].planes = NULL;
	*/
	gCount = -1;
	/* We will look for the following tokens:
	* tilt: tiltx,tilty,tiltz;
	* translation: shiftx shifty shiftz
	* plane: vectX vectY vectZ point[0] point[1] point[2]
	*/
	while (readNextParam(title,parStr)) {
		// printf("%s\n",parStr);
		/* if we found a new crystal ... */
		if (strncmp(title,"crystal:",8) == 0) {
			gCount++;
			grains[gCount].name = (char *)malloc(NAME_BUF_LEN);
			grains[gCount].amorphFlag = 0;
			grains[gCount].density = 0;
			grains[gCount].rmin = 0;
			grains[gCount].rFactor = 1.0;
			grains[gCount].sphereRadius = 0;
			/* first we want to extract crystal name and file name for
			* unit cell data file from this one line
			*/
			grains[gCount].name=parStr;
			str = strnext(parStr," \t");
			if (str != NULL) {
				grains[gCount].name[str-parStr-1]='\0';
				// printf("%s, name: %s\n",parStr,grains[gCount].name);
				//str = strnext(str," \t");
				//if (str != NULL) {
				strcpy(unitCellFile,str);
				if ((str = strchr(unitCellFile,' ')) != NULL)
					*str = '\0';
				if ((str = strchr(unitCellFile,'\t')) != NULL)
					*str = '\0';
			}
			else {
				printf("Error: no unit cell data file specified for crystal %s\n",
					grains[gCount].name);
				return 0;       
			}

			// sscanf(parStr,"%s %s",grains[gCount].name,unitCellFile);
			grains[gCount].nplanes = 0;
			grains[gCount].planes.clear();
			muls->nCellX = 1;
			muls->nCellY = 1;
			muls->nCellZ = 1;
			muls->ctiltx = 0;
			muls->ctilty = 0;
			muls->ctiltz = 0;

			muls->tds = 0;
			tempCell = readUnitCell(grains[gCount].natoms, unitCellFile, *muls, 0);

			/*****************************************************
			* Test code
			*
			printf("ax: %g, by: %g, cz: %g\n",muls->ax,muls->by,muls->c);
			for (i=0;i<grains[gCount].natoms;i++) {
			printf("%d: %d (%g,%g,%g)\n",i,tempCell[i].Znum,tempCell[i].x,tempCell[i].y,tempCell[i].z);
			}
			*
			*****************************************************/

			grains[gCount].unitCell = std::vector<atom>(grains[gCount].natoms);
			grains[gCount].unitCell = tempCell;

			grains[gCount].alpha = muls->cAlpha;
			grains[gCount].beta  = muls->cBeta;
			grains[gCount].gamma = muls->cGamma;
			grains[gCount].ax = muls->ax;
			grains[gCount].by = muls->by;
			grains[gCount].cz = muls->c;
		}
		/***************************************************
		* amorphous stuff
		*/
		else if (strncmp(title,"amorph:",7) == 0) {
			gCount++;
			grains[gCount].name = (char *)malloc(NAME_BUF_LEN);
			grains[gCount].amorphFlag = AMORPHOUS;
			grains[gCount].density = 0;
			grains[gCount].rmin = 1000;
			grains[gCount].rFactor = 1.2f;
			grains[gCount].sphereRadius = 0;
			/* first we want to extract crystal name and file name for
			* unit cell data file from this one line
			*/
			grains[gCount].name=parStr;
			str = strnext(parStr," \t");
			if (str != NULL) {
				grains[gCount].name[str-parStr-1]='\0';
				// printf("%s, name: %s\n",parStr,grains[gCount].name);
				//str = strnext(str," \t");
				//if (str != NULL) {
				strcpy(unitCellFile,str);
				if ((str = strchr(unitCellFile,' ')) != NULL)
					*str = '\0';
				if ((str = strchr(unitCellFile,'\t')) != NULL)
					*str = '\0';
			}
			else {
				printf("Error: no unit cell data file specified for crystal %s\n",
					grains[gCount].name);
				return 0;       
			}

			// sscanf(parStr,"%s %s",grains[gCount].name,unitCellFile);
			grains[gCount].nplanes = 0;
			grains[gCount].planes.clear();
			muls->nCellX = 1;
			muls->nCellY = 1;
			muls->nCellZ = 1;
			muls->ctiltx = 0;
			muls->ctilty = 0;
			tempCell = readUnitCell(grains[gCount].natoms, unitCellFile, *muls, 0);
			grains[gCount].unitCell = tempCell;
			//grains[gCount].unitCell = (atom *)malloc(grains[gCount].natoms * sizeof(atom));
			//memcpy(grains[gCount].unitCell, tempCell, grains[gCount].natoms * sizeof(atom));
			grains[gCount].alpha = 0;
			grains[gCount].beta  = 0;
			grains[gCount].gamma = 0;
			grains[gCount].ax = 0;
			grains[gCount].by = 0;
			grains[gCount].cz = 0;
		}
		/***************************************************
		* code for specially distributed amorphous stuff
		*/
		else if (strncmp(title,"special:",8) == 0) {
			gCount++;
			grains[gCount].name = (char *)malloc(NAME_BUF_LEN);
			grains[gCount].amorphFlag = SPECIAL_GRAIN;
			grains[gCount].density = 0;
			grains[gCount].rmin = 1000;      
			grains[gCount].rFactor = 1.2f;
			grains[gCount].sphereRadius = 0;
			/* first we want to extract crystal name and file name for
			* unit cell data file from this one line
			*/
			grains[gCount].name = parStr;  // name of this grain
			// don't need cfg input file
			grains[gCount].natoms = 0;
			grains[gCount].unitCell.clear();
			// sscanf(parStr,"%s %s",grains[gCount].name,unitCellFile);
			grains[gCount].nplanes = 0;
			grains[gCount].planes.clear();
			grains[gCount].alpha = 0;
			grains[gCount].beta  = 0;
			grains[gCount].gamma = 0;
			grains[gCount].ax = 0;
			grains[gCount].by = 0;
			grains[gCount].cz = 0;
		}  // end of "special"
		/* if we found tilt data */
		else if (gCount >= 0) { 
			if (strncmp(title,"tilt:",5) == 0) {
				sscanf(parStr,"%lf %lf %lf",&(grains[gCount].tiltx),
					&(grains[gCount].tilty),&(grains[gCount].tiltz));
				if (strstr(parStr,"degree") != NULL) {
					grains[gCount].tiltx *= static_cast<float_tt>(PI180);
					grains[gCount].tilty *= static_cast<float_tt>(PI180);
					grains[gCount].tiltz *= static_cast<float_tt>(PI180);
				}
			}
			/* assign density */
			else if (strncmp(title,"density:",8) == 0) {
				sscanf(parStr,"%lf",&(grains[gCount].density));
				grains[gCount].rmin = static_cast<float_tt>(pow(sqrt(15.0/144.0)/grains[gCount].density,1.0/3.0));
			}
			/* assign density factor */
			else if (strncmp(title,"rmin:",5) == 0) {
				sscanf(parStr,"%lf",&(grains[gCount].rmin));
				grains[gCount].density = static_cast<float_tt>(sqrt(15.0/144)*pow(grains[gCount].rmin,3));
			}
			else if (strncmp(title,"r-factor:",9) == 0) {
				sscanf(parStr,"%lf",&(grains[gCount].rFactor));
			}
			/* if we found shift data */
			else if (strncmp(title,"translation:",12) == 0) {
				sscanf(parStr,"%lf %lf %lf",&(grains[gCount].shiftx),
					&(grains[gCount].shifty),&(grains[gCount].shiftz));    
			}
			else if (strncmp(title,"sphere:",7) == 0) {
				sscanf(parStr,"%lf %lf %lf %lf",&(grains[gCount].sphereRadius),
					&(grains[gCount].sphereX),&(grains[gCount].sphereY),&(grains[gCount].sphereZ));    
			}
			/* if we found a new plane for this crystal */
			else if (strncmp(title,"plane:",6) == 0) {
				grains[gCount].nplanes++;
				grains[gCount].planes.resize(grains[gCount].nplanes);
				QSf3Vec point = grains[gCount].planes[grains[gCount].nplanes-1].point;
				QSf3Vec vect1 = grains[gCount].planes[grains[gCount].nplanes-1].vect1;
				QSf3Vec vect2 = grains[gCount].planes[grains[gCount].nplanes-1].vect2;
				
				sscanf(parStr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
					&(point[0]), &(point[1]), &(point[2]),
					&(vect1[0]), &(vect1[1]), &(vect1[2]), 
					&(vect2[0]), &(vect2[1]), &(vect2[2]));
				grains[gCount].planes[grains[gCount].nplanes-1].norm = vect1.cross(vect2);
			} /* end of if plane ... */
			else if (strncmp(title,"atom:",5) == 0) {
				grains[gCount].natoms++;
				grains[gCount].unitCell.resize(grains[gCount].natoms);
				// Now read Znum, r (->z), and count (->y)
				sscanf(parStr,"%d %f %d %f",&(grains[gCount].unitCell[grains[gCount].natoms-1].Znum),
					&(grains[gCount].unitCell[grains[gCount].natoms-1].z),&Nkind,
					&(grains[gCount].unitCell[grains[gCount].natoms-1].y));

				// assign number of atoms directly, if specified in input file
				if (Nkind > 0) grains[gCount].unitCell[grains[gCount].natoms-1].y = (float)Nkind;
				for (i=0;i<muls->atomKinds;i++) 	
					if (muls->Znums[i] == grains[gCount].unitCell[grains[gCount].natoms-1].Znum) break;
				if (i == muls->atomKinds) {
					muls->atomKinds++;
					muls->Znums.resize(muls->atomKinds);
					muls->Znums[i] = grains[gCount].unitCell[grains[gCount].natoms-1].Znum;
				}
			} /* end of if "atom:" */
		} /* end of if gCount >=0 */
	}
	parClose();
	if (printFlag) showData();
	return 1;
}


void showData() {
	int g,p;

	printf("Supercell dimensions: %g x %g x %g A\n",
		superCell.ax,superCell.by,superCell.cz);

	for (g=0;g<nGrains;g++) {
		if (grains[g].amorphFlag == 0) {
			printf("%d: Grain %s (ax=%g, by=%g, cz=%g, alpha=%g, beta=%g, gamma=%g)\n",
				g,grains[g].name,grains[g].ax,grains[g].by,grains[g].cz,
				grains[g].alpha,grains[g].beta,grains[g].gamma);
			printf("%d atoms in unit cell\n",grains[g].natoms);
			printf("tilt=(%g, %g, %g)rad shift=(%g, %g, %g)A\n",
				grains[g].tiltx,grains[g].tilty,grains[g].tiltz,
				grains[g].shiftx,grains[g].shifty,grains[g].shiftz);
		}
		else {
			printf("%d: Grain %s (density=%g, rmin=%g, r-factor=%g\n",
				g,grains[g].name,grains[g].density,grains[g].rmin,grains[g].rFactor);
		}
		printf("planes:\n");
		for (p=0;p<grains[g].nplanes;p++) 
			printf("vect1=(%g %g %g) vect2=(%g %g %g) point=(%g %g %g) normal=(%g %g %g)\n",
			grains[g].planes[p].vect1[0],grains[g].planes[p].vect1[1],
			grains[g].planes[p].vect1[2],
			grains[g].planes[p].vect2[0],grains[g].planes[p].vect2[1],
			grains[g].planes[p].vect2[2],
			grains[g].planes[p].point[0],grains[g].planes[p].point[1],
			grains[g].planes[p].point[2],
			grains[g].planes[p].norm[0],grains[g].planes[p].norm[1],
			grains[g].planes[p].norm[2]);
	} // for g=0 ...
}


/********************************************************************
* This function adds only the crystalline atoms to the superCell
* the amorphous ones are handled by makeAmorphous,  and makeSpecial
********************************************************************/
void makeSuperCell() {
	int g,p,iatom,ix,iy,iz,atomCount = 0;
	// atom *atomPtr;
	// atomPtr = (atom *)malloc(sizeof(atom));
	//static float_tt *axCell,*byCell,*czCell=NULL;
	//QSfVec axCell, byCell, czCell;
	QSf3Mat Mm, Mminv, Mrot, Mr, Mr2;
	QSf3Vec a, b;
	float_tt maxLength,dx,dy,dz,d,dxs,dys,dzs;
	atom newAtom;
	// float_tt xpos,ypos,zpos;
	int nxmin,nxmax,nymin,nymax,nzmin,nzmax;

	/* calculate maximum length in supercell box, which is naturally 
	* the room diagonal:
	*/
	// calculation depended on ax being a pointer?
	//maxLength = vectLength(superCell.ax);
	maxLength = sqrt(superCell.ax*superCell.ax+
		superCell.by*superCell.by+
		superCell.cz*superCell.cz);

	// TODO: these are slices.
	//axCell=Mm[0]; byCell=Mm[1]; czCell=Mm[2];

	atomCount = superCell.natoms;
	for (g=0;g<nGrains;g++) {
		/********************************************************
		* if this grain is a crystalline one ... 
		*/
		if (grains[g].amorphFlag == 0) {
			dx = grains[g].shiftx/superCell.ax; 
			dy = grains[g].shifty/superCell.by; 
			dz = grains[g].shiftz/superCell.cz;
			/* find the rotated unit cell vectors .. */
			makeCellVect(grains[g], Mm);
			// showMatrix(Mm,3,3,"M");
			///////////////////////////////////////////////////////////////////

			//memset(Mrot[0],0,3*3*sizeof(float_tt));
			Mrot.setZero();
			Mrot(0,0) = 1.0; Mrot(1,1) = 1.0; Mrot(2,2) = 1.0; 
			Mr2 = Mrot;
			//memcpy(Mr2[0],Mrot[0],3*3*sizeof(float_tt));
			Mr.setZero();
			//memset(Mr[0],0,3*3*sizeof(float_tt));
			Mr(0,0) = 1.0; Mr(1,1) = cos(grains[g].tiltx); Mr(2,1) = sin(grains[g].tiltx); 
			Mr(1,2) = -sin(grains[g].tiltx); Mr(2,2) = cos(grains[g].tiltx);
			//  20130514 correct eigen format, I think.
			Mr2 = Mr*Mrot;
			//matrixProduct(Mrot,3,3,Mr,3,3,Mr2);
			Mrot = Mr2;
			//memcpy(Mrot[0],Mr2[0],3*3*sizeof(float_tt));
			// showMatrix(Mrot,3,3,"Mrotx");

			//memset(Mr[0],0,3*3*sizeof(float_tt));
			Mr.setZero();
			Mr(1,1) = 1.0; Mr(0,0) = cos(grains[g].tilty); Mr(2,0) = -sin(grains[g].tilty); 
			Mr(0,2) = sin(grains[g].tilty); Mr(2,2) = cos(grains[g].tilty);
			//matrixProduct(Mrot,3,3,Mr,3,3,Mr2);
			//  20130514 correct eigen format, I think.
			Mr2 = Mr*Mrot;
			//memcpy(Mrot[0],Mr2[0],3*3*sizeof(float_tt));
			Mrot = Mr2;
			// showMatrix(Mrot,3,3,"Mrotxy");
			Mr.setZero();
			//memset(Mr[0],0,3*3*sizeof(float_tt));
			Mr(2,2) = 1.0; Mr(0,0) = cos(grains[g].tiltz); Mr(1,0) = sin(grains[g].tiltz); 
			Mr(0,1) = -sin(grains[g].tiltz); Mr(1,1) = cos(grains[g].tiltz);
			//  20130514 correct eigen format, I think.
			Mr2 = Mr*Mrot;
			//matrixProduct(Mrot,3,3,Mr,3,3,Mr2);
			Mrot = Mr2;
			//memcpy(Mrot[0],Mr2[0],3*3*sizeof(float_tt));
			// showMatrix(Mrot,3,3,"Mrotxyz");

			///////////////////////////////////////////////////////////////////
			/*
			rotateVect(axCell,axCell,grains[g].tiltx,grains[g].tilty,grains[g].tiltz);
			rotateVect(byCell,byCell,grains[g].tiltx,grains[g].tilty,grains[g].tiltz);
			rotateVect(czCell,czCell,grains[g].tiltx,grains[g].tilty,grains[g].tiltz);
			*/
			//inverse_3x3(Mminv[0],Mm[0]);
			Mminv = Mm.inverse();
			//matrixProduct(Mm,3,3,Mrot,3,3,Mr2);
			//  20130514 correct eigen format, I think.
			Mr2 = Mrot*Mm;
			//memcpy(Mm[0],Mr2[0],3*3*sizeof(float_tt));
			Mm = Mr2;
			// showMatrix(Mm,3,3,"M");
			//inverse_3x3(Mr2[0],Mm[0]);
			Mr2 = Mm.inverse();

			/* find out how far we will have to go in units of unit cell vectors.
			* when creating the supercell by checking the number of unit cell vectors 
			* necessary to reach every corner of the supercell box.
			*/
			a.setZero();
			//matrixProduct(a,1,3,Mr2,3,3,b);
			//  20130514 correct eigen format, I think.
			b=Mr2*a;
			// showMatrix(Mm,3,3,"M");
			// showMatrix(Mminv,3,3,"M");
			nxmin = nxmax = (int)floor(b(0,0)-dx); 
			nymin = nymax = (int)floor(b(1,0)-dy); 
			nzmin = nzmax = (int)floor(b(2,0)-dz);
			for (ix=0;ix<=1;ix++) for (iy=0;iy<=1;iy++) for (iz=0;iz<=1;iz++) {
				a(0,0)=ix*superCell.ax; a(1,0)=iy*superCell.by; a(2,0)=iz*superCell.cz;

				//matrixProduct(a,1,3,Mr2,3,3,b);
				//  20130514 correct eigen format, I think.
				b=Mr2*a;
				// showMatrix(b,1,3,"b");
				if (nxmin > (int)floor(b(0,0)-dx)) nxmin=(int)floor(b(0,0)-dx);
				if (nxmax < (int)ceil( b(0,0)-dx)) nxmax=(int)ceil( b(0,0)-dx);
				if (nymin > (int)floor(b(1,0)-dy)) nymin=(int)floor(b(1,0)-dy);
				if (nymax < (int)ceil( b(1,0)-dy)) nymax=(int)ceil( b(1,0)-dy);
				if (nzmin > (int)floor(b(2,0)-dz)) nzmin=(int)floor(b(2,0)-dz);
				if (nzmax < (int)ceil( b(2,0)-dz)) nzmax=(int)ceil( b(2,0)-dz);	  
			}
			// nxmin--;nxmax++;nymin--;nymax++;nzmin--;nzmax++;
			superCell.atoms.resize((superCell.natoms+1+
									(nxmax-nxmin)*(nymax-nymin)*
									(nzmax-nzmin)*grains[g].natoms));
			// showMatrix(Mm,3,3,"Mm");
			// showMatrix(Mminv,3,3,"Mminv");
			printf("Grain %d: range: (%d..%d, %d..%d, %d..%d)\n",
				g,nxmin,nxmax,nymin,nymax,nzmin,nzmax);

			dx = grains[g].shiftx; 
			dy = grains[g].shifty; 
			dz = grains[g].shiftz;
			for (iatom=0;iatom<grains[g].natoms;iatom++) {
				memcpy(&newAtom,&(grains[g].unitCell[iatom]),sizeof(atom));
				// We need to convert the cartesian coordinates of this atom
				// to fractional ones:
				b(0,0) = newAtom.x;
				b(1,0) = newAtom.y;
				b(2,0) = newAtom.z;
				//  20130514 correct eigen format, I think.
				a = Mminv*b;
				//matrixProduct(b,1,3,Mminv,3,3,a);
				newAtom.x = a(0,0); newAtom.y = a(1,0); newAtom.z = a(2,0);
				//printf("%2d: %d (%g,%g,%g) (%g,%g,%g)\n",iatom,newAtom.Znum,
				//		 b(0,0),b(1,0),b(2,0),a(0,0),a(1,0),a(2,0));
				for (ix=nxmin;ix<=nxmax;ix++) {
					for (iy=nymin;iy<=nymax;iy++) {
						for (iz=nzmin;iz<=nzmax;iz++) {
							/* atom position in reduced coordinates: */
							// a(0,0) = ix+newAtom.x; a(1,0) = iy+newAtom.y; a(2,0) = iz+newAtom.z;		 
							a(0,0) = newAtom.x+ix; a(1,0) = newAtom.y+iy; a(2,0) = newAtom.z+iz;
							//  20130514 correct eigen format, I think.
							b=Mm*a;
							//matrixProduct(a,1,3,Mm,3,3,b);

							/*	  
							b(0,0) = a(0,0)*Mm(0,0)+a(1,0)*Mm(1,0)+a(2,0)*Mm(2,0);
							b(1,0) = a(0,0)*Mm(0,1)+a(1,0)*Mm(1,1)+a(2,0)*Mm(2,1);
							b(2,0) = a(0,0)*Mm(0,2)+a(1,0)*Mm(1,2)+a(2,0)*Mm(2,2);
							*/
							/* // same as matrixProduct: 
							b(0,0) = a(0,0)*Mm(0,0)+a(1,0)*Mm(0,1)+a(2,0)*Mm(0,2);
							b(1,0) = a(0,0)*Mm(1,0)+a(1,0)*Mm(1,1)+a(2,0)*Mm(1,2);
							b(2,0) = a(0,0)*Mm(2,0)+a(1,0)*Mm(2,1)+a(2,0)*Mm(2,2);
							*/
							superCell.atoms[atomCount].x  = b(0,0)+dx; 
							superCell.atoms[atomCount].y  = b(1,0)+dy; 
							superCell.atoms[atomCount].z  = b(2,0)+dz; 
							if ((superCell.atoms[atomCount].x >= 0) && 
								(superCell.atoms[atomCount].x < superCell.ax) &&
								(superCell.atoms[atomCount].y >= 0) && 
								(superCell.atoms[atomCount].y < superCell.by) &&
								(superCell.atoms[atomCount].z >= 0) && 
								(superCell.atoms[atomCount].z < superCell.cz)) {
									// If this is a sphere:
									if (grains[g].sphereRadius > 0) {
										dxs = superCell.atoms[atomCount].x - grains[g].sphereX;
										dys = superCell.atoms[atomCount].y - grains[g].sphereY;
										dzs = superCell.atoms[atomCount].z - grains[g].sphereZ;
										if (dxs*dxs+dys*dys+dzs*dzs < grains[g].sphereRadius*grains[g].sphereRadius) { 
											superCell.atoms[atomCount].dw = newAtom.dw;
											superCell.atoms[atomCount].occ = newAtom.occ;
											superCell.atoms[atomCount].q = newAtom.q;
											superCell.atoms[atomCount].Znum = newAtom.Znum;
											atomCount++;
										}
									}
									// If this is a straight-edged grain
									else {
										for (p=0;p<grains[g].nplanes;p++) {
											/*
											printf("hello %d (%g %g %g)\n",g,
											superCell.atoms[atomCount].x,superCell.atoms[atomCount].y,
											superCell.atoms[atomCount].z);
											*/
											d = findLambda(grains[g].planes[p],&(superCell.atoms[atomCount].z),-1);

											/*
											printf("%3d lambda: %g (%g %g %g), (%g %g %g), %d\n",atomCount,d,
											newAtom.x,newAtom.y,newAtom.z,
											superCell.atoms[atomCount].x,superCell.atoms[atomCount].y,
											superCell.atoms[atomCount].z,grains[g].nplanes);
											*/		
											if (d < 0)
												break;
										}
										/* if all the previous test have been successful, this atom is IN,
										* which means that we also need to copy the other data elements
										* for this atom. 
										*/
										if (p == grains[g].nplanes) {
											superCell.atoms[atomCount].q = newAtom.q;
											superCell.atoms[atomCount].dw = newAtom.dw;
											superCell.atoms[atomCount].occ = newAtom.occ;
											superCell.atoms[atomCount].Znum = newAtom.Znum;
											atomCount++;
										}
									} // if this is a sphere or not ...
							}
						} /* iz ... */
					} /* iy ... */
				} /* ix ... */
			} /* iatom ... */
			superCell.natoms = atomCount;
		} /* end of if !amorph,i.e. crystalline */
	} /* g=0..nGrains .. */
	/*
	atomPtr->x = 0.5; atomPtr->y = 0.2; atomPtr->z = 0.7;
	findLambda(grains[0].planes,&(atomPtr->z),1);
	*/ 

}

/*****************************************************************
* Create the amorphous part of the model
****************************************************************/
void makeAmorphous() {
	int g,p,iatom,ix,iy,iz,ic,atomCount = 0,amorphSites,amorphAtoms;
	//static float_tt *axCell,*byCell,*czCell=NULL;

	QSf3Mat Mm;

	float_tt rCellx,rCelly,rCellz;
	float_tt d,r;
	atom *amorphCell;
	// atom newAtom;
	// float_tt xpos,ypos,zpos;
	int nx,ny,nz;
	int *randArray,randCount;

	for (g=0;g<nGrains;g++) {
		/********************************************************
		* if this grain is an amorphous one ... 
		*/
		if (grains[g].amorphFlag == AMORPHOUS) {
			r = grains[g].rmin/grains[g].rFactor;
			/* create an hexagonally closed packed unit cell for initial amorphous structure 
			* The length of each vector is now 1 
			*/
			// TODO: verify that this matrix is not transposed!
			Mm(0,0) = r;
			Mm(1,0) = 0;
			Mm(2,0) = 0;
			Mm(0,1) = 0.5f*r; 
			Mm(1,1) = 0.5f*sqrt(3.0f)*r;
			Mm(2,1) = 0;
			Mm(0,2) = 0.5f*r; 
			Mm(1,2) = 0.5f/sqrt(3.0f)*r; 
			Mm(2,2) = sqrt(5.0f)/6*r;
			/* size of rectangular cell containing 4 atoms: */
			rCellx = r; rCelly = static_cast<float_tt>(sqrt(3.0)*r); rCellz = static_cast<float_tt>((sqrt(5.0)*r)/3.0);

			/* determine number of unit cells in super cell */
			nx = (int)(superCell.ax/rCellx);
			ny = (int)(superCell.by/rCelly);
			nz = (int)(superCell.cz/rCellz);
			amorphSites = 4*nx*ny*nz;
			amorphCell = (atom *)malloc((amorphSites+1)*sizeof(atom));

			atomCount = 0;
			for (ix=0;ix<=nx;ix++) {
				for (iy=0;iy<=ny;iy++) {
					for (iz=0;iz<=nz;iz++) {
						for (ic=0;ic<4;ic++) {
							/* check if this atom and any of the 4 atoms per rect. unit cell lie within the super cell 
							*/
							amorphCell[atomCount].x  = ix*rCellx-(ic==3)*Mm(0,0)+(ic % 2)*Mm(0,1)+(ic>1)*Mm(0,2); 
							amorphCell[atomCount].y  = iy*rCelly-(ic==3)*Mm(1,0)+(ic % 2)*Mm(1,1)+(ic>1)*Mm(1,2); 
							amorphCell[atomCount].z  = iz*rCellz-(ic==3)*Mm(2,0)+(ic % 2)*Mm(2,1)+(ic>1)*Mm(2,2); 
							if ((amorphCell[atomCount].x >= 0) && 
								(amorphCell[atomCount].x < superCell.ax) &&
								(amorphCell[atomCount].y >= 0) && 
								(amorphCell[atomCount].y < superCell.by) &&
								(amorphCell[atomCount].z >= 0) && 
								(amorphCell[atomCount].z < superCell.cz)) {
									for (p=0;p<grains[g].nplanes;p++) {
										d = findLambda(grains[g].planes[p],&(amorphCell[atomCount].z),-1);
										if (d < 0)
											break;
									}
									/* if all the previous test have been successful, this atom is IN */
									if (p == grains[g].nplanes) atomCount++;
							}
						} /* ic ... */
					} /* iz ... */
				} /* iy ... */
			} /* ix ... */
			amorphSites = atomCount;
			/****************************************************************************************
			* Now we have all the sites within the bounding planes on which we can put atoms
			*/
			/* the true number of amorphous atoms is # of sites / rFactor^3 */
			amorphAtoms = (int)floor((float_tt)amorphSites/pow(grains[g].rFactor,3.0f));
			if (amorphAtoms > amorphSites) amorphAtoms = amorphSites;
			randCount = amorphSites;
			atomCount = superCell.natoms;

			superCell.atoms.resize((atomCount+amorphAtoms+1));

			randArray = (int *)malloc(amorphSites*sizeof(int));
			for (ix=0;ix<amorphSites;ix++) randArray[ix] = ix;
			/*
			printf("memory allocation: sC.atoms: %d .. %d, rArray: %d .. %d\n",
			(int)superCell.atoms,(int)superCell.atoms+(atomCount+amorphAtoms+1)*sizeof(atom),
			(int)randArray,(int)randArray+amorphSites*sizeof(int));
			*/

			for (ix=amorphAtoms;ix>0;ix--) {
				do {
					iy = (int)((float_tt)rand()*(float_tt)(randCount-1)/(float_tt)(RAND_MAX));
					if (iy >= randCount) iy = randCount-1;
					if (iy < 0) iy = 0;
					//	printf("%5d, iy: %d, sites: %d, atoms: %d  ",ix,iy,randCount,amorphAtoms);
					iz = randArray[iy];
					if (iz > amorphSites) {
						printf("%5d, iy: %d, sites: %d, atoms: %d  ",ix,iy,randCount,amorphAtoms);
						printf("iz: %d (%.2f)\n",iz,(superCell.atoms[atomCount].z));
						printf("makeAmorphous: Error because of overlapping memory areas!\n");
						for (iz=0;iz<=amorphAtoms;iz++)
							printf("iz=%d: %d\n",iz,randArray[iz]);
						exit(0);
					}
				} while (iz > amorphSites);

				/* replace already choosen sites with unused ones, so that we don't occupy
				* any site twice 
				*/
				if (iy == randCount-1)  randCount--;	
				else 
					randArray[iy] = randArray[--randCount];


				iatom = ix % grains[g].natoms;
				superCell.atoms[atomCount].q = grains[g].unitCell[iatom].q;
				superCell.atoms[atomCount].dw = grains[g].unitCell[iatom].dw;
				superCell.atoms[atomCount].occ = grains[g].unitCell[iatom].occ;
				superCell.atoms[atomCount].Znum = grains[g].unitCell[iatom].Znum;
				superCell.atoms[atomCount].x = amorphCell[iz].x;
				superCell.atoms[atomCount].y = amorphCell[iz].y;
				superCell.atoms[atomCount].z = amorphCell[iz].z;
				atomCount++;
			} 
			superCell.natoms = atomCount;
			free(randArray);
			free(amorphCell);
		} /* end of if amorph,i.e. crystalline */
	} /* g=0..nGrains .. */

}


/*****************************************************************
* Create the amorphous part with special atomic distributions 
* and add itto the modelements.
****************************************************************/
#define SQR(x) ((x)*(x))
#define TRIAL_COUNT 1000000
void makeSpecial(int distPlotFlag) {
	int p,g,iatom,atomCount = 0,amorphAtoms;
	float_tt d,r,x,y,z,dist,volume;
	int i,j,Znum,count;
	long seed;
	QSf3Vec pos,center;
	QSfVec grainBound(6);
	int trials = 0,type;
	//char *ptr;

	seed = -(long)(time(NULL));  // initialize random number generator.

	for (g=0;g<nGrains;g++) {
		/********************************************************
		* if this grain is a special one ... 
		*/
		if (grains[g].amorphFlag == SPECIAL_GRAIN) {
			//  type = atoi(grains[g].name);
			//ptr = grains[g].name;
			// finds first character after a space in string grains[g]name.
			//while (*ptr != ' ') ptr++;
			size_t type_pos = grains[g].name.find(" ")+1;
			type = atoi(&(grains[g].name[type_pos]));
			//*ptr = 0;

			// printf("Distribution type: %d \n",type);
			printf("%s: distribution type: %d (%s)\n",grains[g].name,type,
				(type == 2) ? "double gaussian" : (type == 1 ? "single gaussian" : "random"));
			/**************************************************************
			* We would like to calculate the Volume of this grain.  
			* We do this by brute force, by simply checking, if randomly
			* distributed points are within the grain, or not.
			*************************************************************/
			grainBound[0] = superCell.ax;  grainBound[1] = 0;
			grainBound[2] = superCell.by;  grainBound[3] = 0;
			grainBound[4] = superCell.cz;  grainBound[5] = 0;
			for (count=0, i=0; i<TRIAL_COUNT;i++) {
				// remember that we have to create a vector with
				// z,y,x, because that is the way the atom struct is 
				pos[2] = superCell.ax*ran(&seed);
				pos[1] = superCell.by*ran(&seed);
				pos[0] = superCell.cz*ran(&seed);
				for (p=0;p<grains[g].nplanes;p++) {
					d = findLambda(grains[g].planes[p],pos,-1);
					if (d < 0) break;
				}
				// if all the previous tests have been successful, this atom is IN 
				if (p == grains[g].nplanes) {
					count++;
					// center of this grain
					center[0] += pos[2]; center[1] += pos[1]; center[2] += pos[0]; 
					// boundaries in X-direction of this grain
					if (grainBound[0] > pos[2]) grainBound[0] = pos[2]; // xmin
					if (grainBound[1] < pos[2]) grainBound[1] = pos[2]; // xmax
					if (grainBound[2] > pos[1]) grainBound[2] = pos[1]; // ymin
					if (grainBound[3] < pos[1]) grainBound[3] = pos[1]; // ymax
					if (grainBound[4] > pos[0]) grainBound[4] = pos[0]; // zmin
					if (grainBound[5] < pos[0]) grainBound[5] = pos[0]; // zmax
				}
			}
			center[0] /= (float_tt)count;
			center[1] /= (float_tt)count;
			center[2] /= (float_tt)count;
			volume = superCell.ax*superCell.by*superCell.cz*(float_tt)count/(float_tt)TRIAL_COUNT;
			printf("Volume: %gA^3, %g %%\n",volume,(float_tt)(100*count)/(float_tt)TRIAL_COUNT);
			printf("boundaries: x: %g..%g, y: %g..%g, z: %g..%g\n",
				grainBound[0],grainBound[1],grainBound[2], 
				grainBound[3],grainBound[4],grainBound[5]); 

			// First we need to find out how much memory we need to reserve for this grain
			amorphAtoms = 0;
			for (iatom=0;iatom<grains[g].natoms;iatom++) {
				if (grains[g].unitCell[iatom].y < 1.0) {  // if this is a concentration, and no count
					grains[g].unitCell[iatom].y *= volume;  // then convert it to number of atoms
				}
				amorphAtoms += (int)(grains[g].unitCell[iatom].y);
			}

			superCell.atoms.resize((amorphAtoms+superCell.natoms+1));
			atomCount = superCell.natoms;  // start adding amorphous atoms, where we left off.

			// Now we can loop through and add these atoms randomly to the grain
			for (iatom=0;iatom<grains[g].natoms;iatom++) {
				r = grains[g].unitCell[iatom].z;             // radius of this atom
				count = (int)(grains[g].unitCell[iatom].y);  // number of atoms of this kind
				Znum = grains[g].unitCell[iatom].Znum;
				covRad[Znum-1] = r;                          // set radius of other atoms also
				for (j=0;j<count;j++) {
					do { // make it lie within the grain bounding planes
						do { // make the atoms not touch eachother
							// z = superCell.cz*ran(&seed);
							// y = superCell.by*ran(&seed);
							z = grainBound[4]+ran(&seed)*(grainBound[5]-grainBound[4]);	    
							y = grainBound[2]+ran(&seed)*(grainBound[3]-grainBound[2]);	    
							if (fabs(superCell.cz-z) < 2e-5) z = 0.0;
							if (fabs(superCell.by-y) < 2e-5) y = 0.0;
							if (iatom > 0) {
								x = grainBound[0]+ran(&seed)*(grainBound[1]-grainBound[0]);	    
							}
							else {
								switch (type) {
		case 0:
			x = grainBound[0]+ran(&seed)*(grainBound[1]-grainBound[0]);	    
			break;
		case 1:
			x = xDistrFun1(center[0],0.08f*(grainBound[1]-grainBound[0]));	    
			break;
		case 2:
			x = xDistrFun2(center[0],0.80f*(grainBound[1]-grainBound[0]),
				0.08f*(grainBound[1]-grainBound[0]));
			break;
		default:
			x = grainBound[0]+ran(&seed)*(grainBound[1]-grainBound[0]);	    
								}
							}
							if (fabs(superCell.ax-x) < 2e-5) x = 0.0;
							// Now we must check, whether atoms overlap
							for (i=0;i<atomCount;i++) {
								for (p=-1;p<=1;p++) {
									dist = sqrt(SQR(x-superCell.atoms[i].x)+
										SQR(y-superCell.atoms[i].y+p*superCell.by)+
										SQR(z-superCell.atoms[i].z));		  
									if (dist < r+covRad[superCell.atoms[i].Znum-1]) break;		
								}	    
								if (p < 2) break;
							}
							trials++;
							if (trials % amorphAtoms == 0)
								printf("Average trials per atom: %d times, success: %g %%\n",
								trials/amorphAtoms,100.0*(atomCount-superCell.natoms)/
								(float_tt)amorphAtoms);
						} while (i < atomCount);  
						// try until we find one that does not touch any other

						// superCell.atoms[atomCount].dw = 0.0;
						superCell.atoms[atomCount].dw = 0.45f*28.0f/(2.0f*Znum);

						superCell.atoms[atomCount].occ  = 1.0f;
						superCell.atoms[atomCount].q    = 0;
						superCell.atoms[atomCount].Znum = Znum;
						superCell.atoms[atomCount].x    = x;
						superCell.atoms[atomCount].y    = y;
						superCell.atoms[atomCount].z    = z;

						for (p=0;p<grains[g].nplanes;p++) {
							d = findLambda(grains[g].planes[p],&(superCell.atoms[atomCount].z),-1);
							if (d < 0) break;
						}
						// if all the previous tests have been successful, this atom is IN 
					} while(p < grains[g].nplanes);
					atomCount++;
				} // for j=0..count
				printf("%d (%d): %d \n",iatom,Znum,count);
			} // for iatom = 0..natoms
			printf("\n%d amorphous atoms, volume: %gA^3 (%g%%), center: %g, width: %g\n",
				atomCount-superCell.natoms,
				volume,100.0*volume/(superCell.ax*superCell.by*superCell.cz),center[0],
				grainBound[1]-grainBound[0]);
			switch (type) {
	  case 2:
		  xDistrFun2(0.0,0.0,0.0);
		  break;
	  case 1:
		  xDistrFun1(0.0,0.0);
		  break;
			}
			superCell.natoms = atomCount;
		} // if special_grain 
	} // g=0..nGrains ..
	/*******************************************************************
	* Now we must produce distribution plots of the different atom kinds
	*/
	if (distPlotFlag) makeDistrPlot(superCell.atoms,superCell.ax);
}

/* This function creates the single gaussian distruted x-values for the dopand atoms
*/
float_tt xDistrFun1(float_tt xcenter,float_tt width) {
	static long idum = 0;
	static int count = 0;
	static float_tt x2=0,xavg=0;
	float_tt x;

	if (idum == 0) idum  = -(long)(time(NULL)); 

	if (width >0) {
		count ++;
		x = xcenter+width*gasdev(&idum);
		xavg += x;
		x2 += SQR(x-xcenter);
		return x;
	}
	else {
		printf("Statistics (%d): xavg=%g, x2dev=%g\n",count,xavg/(float_tt)count,sqrt(x2/(float_tt)count));
		xavg = 0;
		x2 = 0;
		return 0;
	}
}

/* This function creates the double gaussian distruted x-values for the dopand atoms
* width1 i the distance between the 2 peaks
* width 2 is the width of a single peak.
*/
float_tt xDistrFun2(float_tt xcenter,float_tt width1,float_tt width2) {
	static long idum = 0;
	static long seed = 0;
	static int count = 0;
	static float_tt x2=0,xavg=0;
	float_tt x,dx;

	if (idum == 0) idum  = -(long)(time(NULL)); 
	if (seed == 0) seed  = -(long)(time(NULL)); 

	if (width2 >0) {
		count ++;
		dx = width2*gasdev(&idum);
		x = xcenter+width1*((ran(&seed) > 0.5)-0.5)+dx;
		xavg += x;
		x2 += SQR(dx);
		return x;
	}
	else {
		printf("Statistics (%d): xavg=%g, x2dev=%g\n",count,xavg/(float_tt)count,sqrt(x2/(float_tt)count));
		xavg = 0;
		x2 = 0;
		return 0;
	}
}


// one can run "xmgr -nxy disList.dat &" to view the data produced by this function
#define DR 1.1
void makeDistrPlot(std::vector<atom>atoms,float_tt ax) {
	int j,i,count,ind;
	QSiMat list;
	FILE *fp;

	printf("Atom kinds: %d: ",muls->atomKinds);
	for (i=0;i<muls->atomKinds;i++) printf(" %3d ",muls->Znums[i]);
	printf("\n");

	count = (int)(ax/DR+1);
	list = QSiMat(muls->atomKinds,count);
	list.setZero();
	for (j=0;j<atoms.size();j++) {
		ind = (int)(atoms[j].x/DR);
		if (ind < 0) ind = 0;
		if (ind >= count) ind = count;
		for (i=0;i<muls->atomKinds;i++) if (muls->Znums[i] == atoms[j].Znum) break;
		if (i==muls->atomKinds) {
			// printf("Error: wrong Z (%d)\n",atoms[j].Znum);
		}
		else list(ind,i)++;
	}
	fp = fopen("disList.dat","w");
	for (j=0;j<count;j++) {
		fprintf(fp,"%.3f ",j*DR);
		for (i=0;i<muls->atomKinds;i++) fprintf(fp,"%d ",list(j,i));
		fprintf(fp,"\n");
	}
}
