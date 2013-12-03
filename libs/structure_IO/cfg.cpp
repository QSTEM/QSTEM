/*
  QSTEM - image simulation for TEM/STEM/CBED
  Copyright (C) 2000-2010  Christoph Koch
  Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "cfg.hpp"

int CStructureCfg::Write(atom *atoms,int natoms,char *fileName,MULS *muls) {
	FILE *fp;
	int j;
	char elem[16];
	double ax,by,cz;

	if (natoms < 1) {
		printf("Atom array empty - no file written\n");
		return 1;
	}

	fp = fopen(fileName, "w" );
	if( fp == NULL ) {
		printf("Cannot open file %s\n",fileName);
		return 0;
	}

	// ax = muls->ax > MIN_EDGE_LENGTH ? muls->ax : MIN_EDGE_LENGTH;
	// by = muls->by > MIN_EDGE_LENGTH ? muls->by : MIN_EDGE_LENGTH;
	// cz = muls->c  > MIN_EDGE_LENGTH ? muls->c  : MIN_EDGE_LENGTH;
	ax = muls->ax;
	by = muls->by;
	cz = muls->c;

	fprintf(fp,"Number of particles = %d\n",natoms);
	fprintf(fp,"A = 1.0 Angstrom (basic length-scale)\n");
	fprintf(fp,"H0(1,1) = %g A\nH0(1,2) = 0 A\nH0(1,3) = 0 A\n",ax);
	fprintf(fp,"H0(2,1) = 0 A\nH0(2,2) = %g A\nH0(2,3) = 0 A\n",by);
	fprintf(fp,"H0(3,1) = 0 A\nH0(3,2) = 0 A\nH0(3,3) = %g A\n",cz);
	fprintf(fp,".NO_VELOCITY.\nentry_count = 6\n");
	printf("ax: %g, by: %g, cz: %g n: %d\n",muls->ax,muls->by,muls->c,natoms);


	elem[2] = '\0';
	elem[0] = elTable[2*atoms[0].Znum-2];
	elem[1] = elTable[2*atoms[0].Znum-1];
	// printf("ax: %g, by: %g, cz: %g n: %d\n",muls->ax,muls->by,muls->c,natoms);
	if (elem[1] == ' ') elem[1] = '\0';
	fprintf(fp,"%g\n%s\n",2.0*atoms[0].Znum,elem);
	fprintf(fp,"%g %g %g %g %g %g\n",atoms[0].x/ax,atoms[0].y/by,atoms[0].z/cz,
		atoms[0].dw,atoms[0].occ,atoms[0].q);



	for (j=1;j<natoms;j++) {
		if (atoms[j].Znum != atoms[j-1].Znum) {
			elem[0] = elTable[2*atoms[j].Znum-2];
			elem[1] = elTable[2*atoms[j].Znum-1];
			if (elem[1] == ' ') elem[1] = '\0';
			fprintf(fp,"%g\n%s\n",2.0*atoms[j].Znum,elem);
			// printf("%d: %g\n%s\n",j,2.0*atoms[j].Znum,elem);
		}
		fprintf(fp,"%g %g %g %g %g %g\n",atoms[j].x/ax,atoms[j].y/by,atoms[j].z/cz,
			atoms[j].dw,atoms[j].occ,atoms[j].q);
		// if (atoms[j].occ != 1) printf("Atom %d: occ = %g\n",j,atoms[j].occ);
	} 
	fclose(fp);

	return 1;
}


#define MIN_EDGE_LENGTH 5.18 /* minimal allowed edge length in A
* going below this limit will crash 
* AtomEye.
*/

// write CFG file using atomic positions stored in pos, Z's in Znum and DW-factors in dw
// the unit cell is assumed to be cubic
int CStructureCfg::WriteFractCubic(double *pos,int *Znum,double *dw,int natoms,char *fileName,
                                   double a,double b,double c) {
  FILE *fp;
  int j;
  char elem[16];
  double ax,by,cz;
  
  if (natoms < 1) {
    printf("Atom array empty - no file written\n");
    return 1;
  }
  
  fp = fopen(fileName, "w" );
  if( fp == NULL ) {
    printf("Cannot open file %s\n",fileName);
    return 0;
  }
  
  ax = a > MIN_EDGE_LENGTH ? a : MIN_EDGE_LENGTH;
  by = b > MIN_EDGE_LENGTH ? b : MIN_EDGE_LENGTH;
  cz = c > MIN_EDGE_LENGTH ? c : MIN_EDGE_LENGTH;
  
  fprintf(fp,"Number of particles = %d\n",natoms);
  fprintf(fp,"A = 1.0 Angstrom (basic length-scale)\n");
  fprintf(fp,"H0(1,1) = %g A\nH0(1,2) = 0 A\nH0(1,3) = 0 A\n",ax);
  fprintf(fp,"H0(2,1) = 0 A\nH0(2,2) = %g A\nH0(2,3) = 0 A\n",by);
  fprintf(fp,"H0(3,1) = 0 A\nH0(3,2) = 0 A\nH0(3,3) = %g A\n",cz);
  fprintf(fp,".NO_VELOCITY.\nentry_count = 4\n");
  printf("ax: %g, by: %g, cz: %g n: %d\n",ax,by,c,natoms);

  
  elem[2] = '\0';
  elem[0] = elTable[2*Znum[0]-2];
  elem[1] = elTable[2*Znum[0]-1];
  // printf("ax: %g, by: %g, cz: %g n: %d\n",muls->ax,muls->by,muls->c,natoms);
  if (elem[1] == ' ') elem[1] = '\0';
  fprintf(fp,"%g\n%s\n",2.0*Znum[0],elem);
  fprintf(fp,"%g %g %g %g\n",pos[0]*a/ax,pos[1]*b/by,pos[2]*c/cz,dw[0]);

  for (j=1;j<natoms;j++) {
    if (Znum[j] != Znum[j-1]) {
      elem[0] = elTable[2*Znum[j]-2];
      elem[1] = elTable[2*Znum[j]-1];
      if (elem[1] == ' ') elem[1] = '\0';
      fprintf(fp,"%g\n%s\n",2.0*Znum[j],elem);
      // printf("%d: %g\n%s\n",j,2.0*atoms[j].Znum,elem);
    }
    fprintf(fp,"%g %g %g %g\n",pos[3*j+0]*a/ax,pos[3*j+1]*b/by,pos[3*j+2]*c/cz,dw[j]);
  } 
  fclose(fp);

  return 1;
}

/***********************************************************************
* The following function returns the number of atoms in the specified
* CFG file and updates the cell parameters in the muls struct
***********************************************************************/
int CStructureCfg::ReadCellParams(MULS *muls, double **Mm, char *fileName) {
	int ncoord,i;
	char buf[256];
	double lengthScale;

	// printf("Paramter File pointer (1): %d\n",(int)getFp());
	parFpPush(); /* push the current parameter file file pointer 
				 * on the stack to make this function totaly 
				 * transparent 
				 */
	if (!parOpen(fileName)) {
		printf("Could not open CFG input file %s\n",fileName);
		parFpPull();  /* restore old parameter file pointer */
		return 0;
	}
	// printf("Paramter File pointer (2): %d\n",(int)getFp());

	resetParamFile();  
	setComment('#');  
	if (readparam("Number of particles =",buf,1)) sscanf(buf,"%d",&ncoord);
	if (readparam("A =",buf,1)) sscanf(buf,"%lf",&lengthScale);

	if (readparam("H0(1,1) =",buf,1)) sscanf(buf,"%lf",Mm[0]+0);
	if (readparam("H0(1,2) =",buf,1)) sscanf(buf,"%lf",Mm[0]+1);
	if (readparam("H0(1,3) =",buf,1)) sscanf(buf,"%lf",Mm[0]+2);

	if (readparam("H0(2,1) =",buf,1)) sscanf(buf,"%lf",Mm[0]+3);
	if (readparam("H0(2,2) =",buf,1)) sscanf(buf,"%lf",Mm[0]+4);
	if (readparam("H0(2,3) =",buf,1)) sscanf(buf,"%lf",Mm[0]+5);

	if (readparam("H0(3,1) =",buf,1)) sscanf(buf,"%lf",Mm[0]+6);
	if (readparam("H0(3,2) =",buf,1)) sscanf(buf,"%lf",Mm[0]+7);
	if (readparam("H0(3,3) =",buf,1)) sscanf(buf,"%lf",Mm[0]+8);
	/*
	if (readparam(".NO_VELOCITY.",buf,1)) noVelocityFlag = 1; 
	if (readparam("entry_count =",buf,1)) sscanf(buf,"%lf",&entryCount);
	if (!noVelocityFlag) entryCount+=3;
	*/
	setComment('%');
	parClose();   
	parFpPull();  /* restore old parameter file pointer */

	for (i=0;i<9;i++) Mm[0][i] *= lengthScale;

	muls->ax = sqrt(Mm[0][0]*Mm[0][0]+Mm[0][1]*Mm[0][1]+Mm[0][2]*Mm[0][2]);
	muls->by = sqrt(Mm[1][0]*Mm[1][0]+Mm[1][1]*Mm[1][1]+Mm[1][2]*Mm[1][2]);
	muls->c  = sqrt(Mm[2][0]*Mm[2][0]+Mm[2][1]*Mm[2][1]+Mm[2][2]*Mm[2][2]);
	muls->cGamma = atan2(Mm[1][1],Mm[1][0]);
	muls->cBeta = acos(Mm[2][0]/muls->c);
	muls->cAlpha = acos(Mm[2][1]*sin(muls->cGamma)/muls->c+cos(muls->cBeta)*cos(muls->cGamma));
	muls->cGamma /= (float)PI180;
	muls->cBeta  /= (float)PI180;
	muls->cAlpha /= (float)PI180;
	if (ncoord < 1) {
		printf("Number of atoms in CFG file not specified!\n");
		ncoord = 0;
	}
	return ncoord;
}


/*******************************************************************************
* This function reads the atomic position and element data for a single atom
* from a .cfg file.  The atomic positions are given in reduced coordinates.
* Calling this function several times will read one atom after another from the 
* file.
* This function will return -1, if the end of file is reached prematurely, 0 otherwise.
*******************************************************************************/

int CStructureCfg::ReadNextAtom(atom *newAtom, int flag, char *fileName) {
	static FILE *fp=NULL;
	static int noVelocityFlag = 1,entryCount = 3,element = 1;
	static char buf[NCMAX];
	static double *atomData = NULL;
	static double mass = 28;
	char *str = NULL;
	int j;


	if (flag < 0) {
		if (fp != NULL) {
			parClose();   
			parFpPull();  /* restore old parameter file pointer */      
			fp = NULL;
			setComment('%');
		}
		return -1;
	}

	if (fp == NULL) {
		parFpPush();  /* save old parameter file pointer */      
		if (!parOpen(fileName)) {
			printf("Could not open CFG input file %s\n",fileName);
			parFpPull();  /* restore old parameter file pointer */
			return -1;
		}
		resetParamFile();  
		if (readparam(".NO_VELOCITY.",buf,1)) noVelocityFlag = 1; 
		else noVelocityFlag = 0;
		if (readparam("entry_count =",buf,1)) sscanf(buf,"%d",&entryCount);
		if (!noVelocityFlag) entryCount+=3;
		fp = getFp();  /* get the file pointer from the parameter file routines */
		atomData = (double *)malloc((entryCount+1)*sizeof(double));
	}

	if (fgets(buf,NCMAX,fp) == NULL) return -1;
	/* check, if this is a new mass number */
	str = strnext(buf," \t");
	if ((atof(buf) >= 1.0) && ((str==NULL) || (*str == '#'))) {
		mass = atof(buf);
		// printf("nV: %d, eC: %d (%g)\n",noVelocityFlag, entryCount,atof(buf));
		if (fgets(buf,NCMAX,fp) == NULL) return -1;    
		element = getZNumber(buf); 
		// printf("*** found element %d (%s %d) ***\n",element,buf,strlen(buf));
		if (fgets(buf,NCMAX,fp) == NULL) return -1;
	}
	str = buf;
	// skip leading spaces:
	while (strchr(" \t",*str) != NULL) str++; 
	for (j=0;j<entryCount;j++) {
		if (str==NULL) {
			printf("readNextCFGatom: Error: incomplete data line: >%s<\n",buf);
			return -1;
		}
		atomData[j] = atof(str); str=strnext(str," \t");
	}


	newAtom->Znum = element;
	newAtom->x    = atomData[0];
	newAtom->y    = atomData[1];
	newAtom->z    = atomData[2];
	// newAtom->dw   = 0.45*28.0/((double)(2*element));	
	// printf("Element: %d, mass=%g\n",element,mass);
	newAtom->dw   = 0.45*28.0/mass;	
	newAtom->occ  = 1.0;
	newAtom->q    = 0.0;
	// read the DW-factor
	if (entryCount > 3+3*(1-noVelocityFlag)) 
		newAtom->dw = atomData[3+3*(1-noVelocityFlag)];
	// read the atom's occupancy:
	if (entryCount > 4+3*(1-noVelocityFlag)) 
		newAtom->occ = atomData[4+3*(1-noVelocityFlag)];
	// read the atom's charge:
	if (entryCount > 5+3*(1-noVelocityFlag)) 
		newAtom->q = atomData[5+3*(1-noVelocityFlag)];
	// printf("Atom: %d (%g %g %g), occ=%g, q=%g\n",newAtom->Znum,newAtom->x,newAtom->y,newAtom->z,newAtom->occ,newAtom->q);	


	return 0;
}
