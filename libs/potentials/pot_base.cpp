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

#include "pot_base.hpp"

CPotential::CPotential(ConfigReaderPtr &configReader)
{
  configReader->ReadProbeArraySize(unsigned int &nx, unsigned int &ny)
  configReader->ReadResolution(m_dx, m_dy);
  configReader->ReadVoltage(m_v0);
  configReader->ReadPotentialOutputParameters(m_savePotential, m_saveProjectedPotential, m_plotPotential);
  
}

void CPotential::Initialize()
{
	m_ddx = m_dx/(double)OVERSAMPLING;
	m_ddy = m_dy/(double)OVERSAMPLING;
	m_ddz = m_dz/(double)OVERSAMPLINGZ;

	/* For now we don't care, if the box has only small 
		* prime factors, because we will not fourier transform it
		* especially not very often.
		*/
	m_boxNx = (int)(m_radius/m_ddx+2.0);  
	m_boxNy = (int)(m_radius/m_ddy+2.0);  
}

/****************************************************************************
* function: atomBoxLookUp - looks up potential at position x, y, z, relative to atom center
*
* Znum = element
* x,y,z = real space position (in A)
* B = Debye-Waller factor, B=8 pi^2 <u^2>
***************************************************************************/
void CPotential::atomBoxLookUp(complex_tt &val, int Znum, float_tt x, float_tt y, float_tt z, float_tt B) 
{
	int boxNx,boxNy,boxNz;
	float_tt dx,dy,dz,ddx,ddy,ddz;
	int ix,iy,iz;
	float_tt maxRadius2;
	char fileName[256],systStr[256];
	int tZ, tnx, tny, tnz, tzOversample;  
	float_tt tdx, tdy, tdz, tv0, tB;
	FILE *fp;
	int numRead = 0,dummy;

	/* initialize all the atoms to non-used */
	if (!m_atomBoxes.find(Znum)) {
		m_atomBoxes[Znum]=atomBoxPtr(new atomBox());
		m_atomBoxes[Znum]->potential = NULL;
		m_atomBoxes[Znum]->rpotential = NULL;
		m_atomBoxes[Znum]->B = -1.0;

		m_radius2 = m_radius*m_radius;

		
		if (m_printLevel > 2)
			printf("Atombox has real space resolution of %g x %g x %gA (%d x %d x %d pixels)\n",
			ddx,ddy,ddz,boxNx,boxNy,boxNz);
	}
	// printf("Debugging: %d %g %g: %g\n",Znum,m_atomBoxes[Znum]->B,B,fabs(m_atomBoxes[Znum]->B - B));

	/* Creating/Reading a atombox for every new kind of atom, but only as needed */
	if (fabs(m_atomBoxes[Znum]->B - B) > 1e-6) {
		//  printf("Debugging 1 (%d: %.7g-%.7g= %.7g), %d\n",
		//	   Znum,m_atomBoxes[Znum]->B,B,fabs(m_atomBoxes[Znum]->B - B),fabs(m_atomBoxes[Znum]->B - B) > 1e-6);
		m_atomBoxes[Znum]->B = B;
		/* Open the file with the projected potential for this particular element
		*/
		sprintf(fileName,"potential_%d_B%d.prj",Znum,(int)(100.0*B));
		if ((fp=fopen(fileName,"r")) == NULL) {
			sprintf(systStr,"scatpot %s %d %g %d %d %d %g %g %g %d %g",
				fileName,Znum,B,boxNx,boxNy,boxNz,ddx,ddy,ddz,OVERSAMPLINGZ,m_v0);
			if (m_printLevel > 2) {
				printf("Could not find precalculated potential for Z=%d,"
					" will calculate now.\n",Znum);
				printf("Calling: %s\n",systStr);
			}
			system(systStr);
			for (dummy=0;dummy < 10000;dummy++);
			if ((fp=fopen(fileName,"r")) == NULL) {
				printf("cannot calculate projected potential using scatpot - exit!\n");	
				exit(0);
			}  

		}
		fgets(systStr,250,fp);
		sscanf(systStr,"%d %le %d %d %d %le %le %le %d %le\n",
			&tZ, &tB, &tnx, &tny, &tnz, &tdx, &tdy, &tdz, &tzOversample, &tv0);
		/* If the parameters in the file don't match the current ones,
		* we need to create a new potential file
		*/
		if ((tZ != Znum) || (fabs(tB-B)>1e-6) || (tnx != boxNx) || (tny != boxNy) || (tnz != boxNz) ||
			(fabs(tdx-ddx) > 1e-5) || (fabs(tdy-ddy) > 1e-5) || (fabs(tdz-ddz) > 1e-5) || 
			(tzOversample != OVERSAMPLINGZ) || (tv0 != m_v0)) {
				if (m_printLevel > 2) {
					printf("Potential input file %s has the wrong parameters\n",fileName);
					printf("Parameters:\n"
						"file:    Z=%d, B=%.3f A^2 (%d, %d, %d) (%.7f, %.7f %.7f) nsz=%d V=%g\n"
						"program: Z=%d, B=%.3f A^2 (%d, %d, %d) (%.7f, %.7f %.7f) nsz=%d V=%g\n"
						"will create new potential file, please wait ...\n",
						tZ,tB,tnx,tny,tnz,tdx,tdy,tdz,tzOversample,tv0,
						Znum,B,boxNx,boxNy,boxNz,ddx,ddy,ddz,OVERSAMPLINGZ,m_v0);
				}
				/* Close the old file, Create a new potential file now 
				*/
				fclose(fp);
				sprintf(systStr,"scatpot %s %d %g %d %d %d %g %g %g %d %g",
					fileName,Znum,B,boxNx,boxNy,boxNz,ddx,ddy,ddz,OVERSAMPLINGZ,m_v0);
				system(systStr);
				if ((fp=fopen(fileName,"r")) == NULL) {
					printf("cannot calculate projected potential using scatpot - exit!\n");
					exit(0);
				}  
				fgets(systStr,250,fp);
		}

		/* Finally we can read in the projected potential
		*/
		if (B == 0) {
			m_atomBoxes[Znum]->rpotential = float3D(boxNz,boxNx,boxNy,"atomBox");
			numRead = fread(m_atomBoxes[Znum]->rpotential[0][0],sizeof(float_tt),
				(size_t)(boxNx*boxNy*boxNz),fp);	
		}
		else {
			m_atomBoxes[Znum]->potential = complex3D(boxNz,boxNx,boxNy,"atomBox");
			numRead = fread(m_atomBoxes[Znum]->potential[0][0],sizeof(complex_tt),
				(size_t)(boxNx*boxNy*boxNz),fp);
		}

		/* writeImage_old(m_atomBoxes[Znum]->potential[0],boxNx,boxNy, 0.0,"potential.img");
		system("showimage potential.img");
		*/
		fclose(fp);

		if (numRead == boxNx*boxNy*boxNz) {
			if (m_printLevel > 1)
				printf("Sucessfully read in the projected potential\n");
		}
		else {
			printf("error while reading potential file %s: read %d of %d values\n",
				fileName,numRead,boxNx*boxNy*boxNz);
			exit(0);
		}
	}	
}



/*****************************************************
* void make3DSlices()
*
* This function will create a 3D potential from whole
* unit cell, slice it, and make transr/i and propr/i
* Call this function with center = NULL, if you don't
* want the array to be shifted.
****************************************************/
void CPotential::make3DSlices(int nlayer,char *fileName,atom *center) 
{
	// FILE *fpu2;
	char filename[512];
	int natom,iatom,iz;  /* number of atoms */
	atom *atoms;
	float_tt c,atomX,atomY,atomZ;
	int i=0,j,ix,iy,iax,iay,iaz,sliceStep;
	int iAtomX,iAtomY,iAtomZ,iRadX,iRadY,iRadZ,iRad2;
	int iax0,iax1,iay0,iay1,iaz0,iaz1,nyAtBox,nyAtBox2,nxyAtBox,nxyAtBox2,iOffsX,iOffsY,iOffsZ;
	int nzSub,Nr,ir,Nz_lut;
	int iOffsLimHi,iOffsLimLo,iOffsStep;

	float_tt *slicePos;
	float_tt z,x,y,r,ddr,dr,r2sqr,x2,y2,potVal,dOffsZ;
	// char *sliceFile = "slices.dat";
	char buf[BUF_LEN];
	FILE *sliceFp;
	float_tt minX,maxX,minY,maxY,minZ,maxZ;
	float_tt atomRadius2;
	time_t time0,time1;
	float s11,s12,s21,s22;
	fftwf_complex	*atPotPtr;
	float *potPtr=NULL, *ptr;
	static int divCount = 0;
	static float_tt **tempPot = NULL;
#if FLOAT_PRECISION == 1
	static fftwf_complex ***oldTrans = NULL;
	static fftwf_complex ***oldTrans0 = NULL;
#else
	static complex_tt ***oldTrans = NULL;
	static complex_tt ***oldTrans0 = NULL;
#endif
	complex_tt dPot;
#if Z_INTERPOLATION
	float_tt ddz;
#endif
#if USE_Q_POT_OFFSETS
	fftwf_complex	*atPotOffsPtr;
#endif

        // parameters to be saved to output files (e.g. thickness)
        std::map<std::string, double> params;

	if (m_trans == NULL) {
		printf("Severe error: trans-array not allocated - exit!\n");
		exit(0);
	}

	if (oldTrans0 == NULL) {
#if FLOAT_PRECISION == 1
		oldTrans0 = (fftwf_complex ***)fftw_malloc(nlayer * sizeof(fftwf_complex**));
#else
		oldTrans0 = (complex_tt ***)fftw_malloc(nlayer * sizeof(complex_tt**));
#endif
		for (i=0;i<nlayer;i++) {
			// printf("%d %d\n",i,(int)(m_trans));
			oldTrans0[i] = m_trans[i]; 
		}
		oldTrans = m_trans;
	}
	if (oldTrans != m_trans)
		printf("Warning: Transmission function pointer has changed!\n");

	/* return, if there is nothing to do */
	if (nlayer <1)
		return;

	/* we need to keep track of which subdivision of the unit cell we are in
	* If the cell is not subdivided, then muls.cellDiv-1 = 0.
	*/
	if ((divCount == 0) || (muls->equalDivs))
		divCount = muls->cellDiv;
	divCount--;

	/* we only want to reread and shake the atoms, if we have finished the 
	* current unit cell.
	*/
	if (divCount == muls->cellDiv-1) {
		if (muls->avgCount == 0) {
			// if this is the first run, the atoms have already been
			// read during initialization
			natom = (*muls).natom;
			atoms = (*muls).atoms;
		}
		else {
			/* 
			the following function makes an array of natom atoms from
			the input file (with x,y,z,dw,occ);
			*/
			// the last parameter is handleVacancies.  If it is set to 1 vacancies 

			// and multiple occupancies will be handled. 
			atoms = readUnitCell(&natom,fileName,muls,1);
			if (m_printLevel>=3)
				printf("Read %d atoms from %s, tds: %d\n",natom,fileName,muls->tds);
			muls->natom = natom;
			muls->atoms = atoms;
		}
		minX = maxX = atoms[0].x;
		minY = maxY = atoms[0].y;
		minZ = maxZ = atoms[0].z;

		for (i=0;i<natom;i++) {
			if (atoms[i].x < minX) minX = atoms[i].x;
			if (atoms[i].x > maxX) maxX = atoms[i].x;
			if (atoms[i].y < minY) minY = atoms[i].y;
			if (atoms[i].y > maxY) maxY = atoms[i].y;
			if (atoms[i].z < minZ) minZ = atoms[i].z;
			if (atoms[i].z > maxZ) maxZ = atoms[i].z;
		}
		/*
		printf("Root of mean square TDS displacement: %f A (wobble=%g at %gK) %g %g %g\n",
		sqrt(u2/natom),wobble,(*muls).tds_temp,ux,uy,uz);
		*/
		if (m_printLevel >= 2) {
			printf("range of thermally displaced atoms (%d atoms): \n",natom);
			printf("X: %g .. %g\n",minX,maxX);
			printf("Y: %g .. %g\n",minY,maxY);
			printf("Z: %g .. %g\n",minZ,maxZ);
		}
		/* 
		define the center of our unit cell by moving the atom specified
		by "center" at position (0.5,0.5,0.0) 
		*/

		if (center != NULL) {
			dx = (*muls).ax/2.0f - (*center).x;	
			dy = (*muls).by/2.0f - (*center).y;	
			dz = -(*center).z;
			for (i=0;i<natom;i++) {
				atoms[i].x += dx;
				if (atoms[i].x < 0.0f) atoms[i].x += (*muls).ax;
				else if (atoms[i].x > (*muls).ax) atoms[i].x -= (*muls).ax;
				atoms[i].y += dy;
				if (atoms[i].y < 0.0f) atoms[i].y += (*muls).by;
				else if (atoms[i].y > (*muls).by) atoms[i].y -= (*muls).by;
				atoms[i].z += dz;
				if (atoms[i].z < 0.0f) atoms[i].z += (*muls).c;
				else if (atoms[i].z > (*muls).c) atoms[i].z -= (*muls).c;
			}
		}

		/**********************************************************
		* Sort the atoms in z.
		*********************************************************/
		qsort(atoms,natom,sizeof(atom),atomCompare);
		if ((*muls).cfgFile != NULL) {
			sprintf(buf,"%s/%s",muls->folder,muls->cfgFile);
			// append the TDS run number
			if (strcmp(buf+strlen(buf)-4,".cfg") == 0) *(buf+strlen(buf)-4) = '\0';
			if (muls->tds) sprintf(buf+strlen(buf),"_%d.cfg",muls->avgCount);
			else sprintf(buf+strlen(buf),".cfg");
		
			// printf("Will write CFG file <%s> (%d)\n",buf,muls->tds)
			writeCFG(atoms,natom,buf,muls);

			if (muls->readPotential) {
				sprintf(buf,"nanopot %s/%s %d %d %d %s",muls->folder,muls->cfgFile,
					ny,nx,muls->slices*muls->cellDiv,muls->folder);
				system(buf);
			}
		}
	} /* end of if divCount==cellDiv-1 ... */
	else {
		natom = muls->natom;
		atoms = muls->atoms;
	}

	/************************************************************** 
	*	setup the slices with their start and end positions
	*	then loop through all the atoms and add their potential to
	*	the slice that their potential reaches into (up to RMAX)
	*************************************************************/
	// c = (*muls).c/(float_tt)((*muls).cellDiv);
	c = muls->sliceThickness * muls->slices;
	dx = resolutionX;
	dy = resolutionY;
	dr   = dx/OVERSAMP_X;  // define step width in which radial V(r,z) is defined 
	iRadX = (int)ceil(m_radius/dx);
	iRadY = (int)ceil(m_radius/dy);
	iRadZ = (int)ceil(m_radius/muls->sliceThickness);
	iRad2 = iRadX*iRadX+iRadY*iRadY;

	if (m_printLevel >= 3) {
		printf("Slab thickness: %gA z-offset: %gA (cellDiv=%d)\n",
			c,c*(float_tt)(muls->cellDiv-divCount-1),divCount);
	}	 
	/*******************************************************
	* initializing slicPos, cz, and transr
	*************************************************************/
	if ((*muls).cz == NULL) {
		(*muls).cz = float1D(nlayer,"cz");
	}
	// sliceFp = fopen(sliceFile,"r");
	sliceFp = NULL;
	slicePos = float1D(nlayer,"slicePos");


	if (muls->sliceThickness == 0)
		(*muls).cz[0] = c/(float_tt)nlayer;
	else
		(*muls).cz[0] = muls->sliceThickness;
	slicePos[0] = (*muls).czOffset;  
	/*
	************************************************************/


	for (i=1;i<nlayer;i++) {
		if (sliceFp == NULL) (*muls).cz[i] = (*muls).cz[0];  
		/* don't need to all be the same, yes they do for fast 3D-FFT method! */
		else {
			fgets(buf,BUF_LEN,sliceFp);
			(*muls).cz[i] = atof(buf);
		}
		slicePos[i] = slicePos[i-1]+(*muls).cz[i-1]/2.0+(*muls).cz[i]/2.0;
	}

	memset(m_trans[0][0],0,nlayer*nx*ny*sizeof(fftwf_complex));
	/* check whether we have constant slice thickness */

	if (muls->fftpotential) {
		for (i = 0;i<nlayer;i++)  if ((*muls).cz[0] != (*muls).cz[i]) break;
		if (i<nlayer) printf("Warning: slice thickness not constant, will give wrong results (iz=%d)!\n",i);

	}

	/*************************************************************************
	* read the potential that has been created externally!
	*/
	if (muls->readPotential) {
          std::vector<unsigned> slice(1);
          for (i=(divCount+1)*muls->slices-1,j=0;i>=(divCount)*muls->slices;i--,j++) {
            slice[0]=i;
            sprintf(buf,"%s/potential",muls->folder);
            imageIO->ReadImage((void **)tempPot,buf, slice);
            for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
                m_trans[j][ix][iy][0] = tempPot[ix][iy];
                m_trans[j][ix][iy][1] = 0.0;
              }
          }
          return;
	}

	// reset the potential to zero:  
	memset((void *)&(m_trans[0][0][0][0]),0,
		muls->slices*muls->potNx*muls->potNy*sizeof(complex_tt));
	nyAtBox   = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionY);
	nxyAtBox  = nyAtBox*(2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionX));
	nyAtBox2  = 2*nyAtBox;
	nxyAtBox2 = 2*nxyAtBox;
	sliceStep = 2*muls->potNx*muls->potNy;

	/*
	for (i=0;i<nlayer;i++)
	printf("slice center: %g width: %g\n",slicePos[i],(*muls).cz[i]);
	*/ 

	/****************************************************************
	* Loop through all the atoms and add their potential 
	* to the slices:									 
	***************************************************************/

	time(&time0);
	for (iatom = 0;iatom<natom;iatom++) {
		// make sure we skip vacancies:
		while (atoms[iatom].Znum == 0) iatom++;
		if (iatom >=natom) break;

		if ((m_printLevel >= 4) && (muls->displayPotCalcInterval > 0)) {
			if (((iatom+1) % (muls->displayPotCalcInterval)) == 0) {
				printf("Adding potential for atom %d (Z=%d, pos=[%.1f, %.1f, %.1f])\n",iatom+1,atoms[iatom].Znum,atoms[iatom].x,atoms[iatom].y,atoms[iatom].z);
			}
		}
		// printf("c=%g, slice thickness=%g, slices=%d, %d\n",c,muls->sliceThickness,muls->slices,muls->displayPotCalcInterval);
		/*
		* c = the thickness of the current slab.
		*
		* if the z-position of this atom is outside the potential slab
		* we won't consider it and skip to the next
		*/
		/* cellDiv = number of times that the big super cell is divided into
		* less big ones, yet often still bigger than a single unit cell
		* (for saving memory)
		* divCount = counter of how many such semi-super cells we already 
		* passed through.
		* c = height in A of one such semi-super cell.
		* Since cellDiv is always >=1, and divCount starts at 0, the actual position
		* of this atom with the super cell is given by:
		*/
		/* c*(float_tt)((*muls).cellDiv-divCount-1) will pick the right super-cell
		* division in the big super-cell
		* The z-offset 0.5*cz[0] will position atoms at z=0 into the middle of the first 
		* slice.
		*/
		atomZ = atoms[iatom].z-c*(float_tt)(muls->cellDiv-divCount-1) + muls->czOffset 
			-(0.5*muls->sliceThickness*(1-muls->centerSlices));
		// make sure that slices are centered for 2D and 3D differently:
		if (muls->potential3D==0)	atomZ += 0.5*muls->sliceThickness;
		else atomZ -= muls->sliceThickness;

		/* Now we need to find the first atom that contributes to this slice */
		/* Make use of the fact that we sorted the atoms in z */
		if ((*muls).nonPeriodZ) {
			if (((*muls).potential3D) && (atomZ -m_radius > c)) break;	 
			if (((*muls).potential3D==0) && (atomZ > c)) break;		
			do {
				// printf("z: %g c: %g\n",atomZ,c);
				if (((*muls).potential3D) && (atomZ+m_radius+muls->sliceThickness >=0)) break;
				if (((*muls).potential3D==0) && (atomZ >=0)) break;			  
				// atomZ = atoms[++iatom].z-c*(float_tt)((*muls).cellDiv-divCount-1)+ (0.5*(*muls).cz[0]*muls->centerSlices);	
				atomZ = atoms[++iatom].z-c*(float_tt)(muls->cellDiv-divCount-1) + muls->czOffset 
					-(0.5*muls->sliceThickness*(1-muls->centerSlices));
				if (muls->potential3D==0)	atomZ += 0.5*muls->sliceThickness;
				else atomZ -= muls->sliceThickness;
			}
			while (iatom < natom-1);
		}
		/* atom coordinates in cartesian coords
		* The x- and y-position will be offset by the starting point
		* of the actually needed array of projected potential
		*/
		atomX = atoms[iatom].x -(*muls).potOffsetX;
		atomY = atoms[iatom].y -(*muls).potOffsetY;

		/* so far we need periodicity in z-direction.
		* This requirement can later be removed, if we 
		* use some sort of residue slice which will contain the 
		* proj. potential that we need to add to the first slice of
		* the next stack of slices
		*	
		* 
		*/

		/*************************************************************
		* real space potential lookup table summation
		************************************************************/
		if (!muls->fftpotential) {
			/* Warning: will assume constant slice thickness ! */
			/* do not round here: atomX=0..dx -> iAtomX=0 */
			/*
			iAtomX = (int)(atomX/dx);	
			if (atomX/dx < (float)iAtomX) iAtomX--; // in case iAtomX is negative
			iAtomY = (int)(atomY/dy);
			if (atomY/dy < (float)iAtomY) iAtomY--;
			iAtomZ = (int)(atomZ/(*muls).cz[0]);
			if (atomZ/(*muls).cz[0] < (float)iAtomZ) iAtomZ--;
			*/
			iAtomX = (int)floor(atomX/dx);  
			iAtomY = (int)floor(atomY/dy);
			iAtomZ = (int)floor(atomZ/muls->cz[0]);

			// printf("atomZ(%d)=%g(%d)\t",iatom,atomZ,iAtomZ);

			if (muls->displayPotCalcInterval > 0) {
				if ((m_printLevel>=3) && ((iatom+1) % muls->displayPotCalcInterval == 0)) {
					printf("adding atom %d [%.3f %.3f %.3f (%.3f)], Z=%d\n",
						iatom+1,atomX+(*muls).potOffsetX,atomY+(*muls).potOffsetY,
						atoms[iatom].z,atomZ,atoms[iatom].Znum);
					/*    (*muls).ax,(*muls).by,(*muls).potOffsetX,(*muls).potOffsetY); */
				}
			}

			for (iax = -iRadX;iax<=iRadX;iax++) {
				if ((*muls).nonPeriod) {
					if (iax+iAtomX < 0) {
						iax = -iAtomX;
						if (abs(iax)>iRadX) break;
					} 
					if (iax+iAtomX >= nx)	break;
				}
				x = (float_tt)(iAtomX+iax)*dx-atomX;
				ix = (iax+iAtomX+16*nx) % nx;	/* shift into the positive range */
				for (iay=-iRadY;iay<=iRadY;iay++) {
					if ((*muls).nonPeriod) {
						if (iay+iAtomY < 0) {
							iay = -iAtomY;
							if (abs(iay)>iRadY) break;
						} 
						if (iay+iAtomY >= ny)	break;
					}
					y = (float_tt)(iAtomY+iay)*dy-atomY;
					iy = (iay+iAtomY+16*ny) % ny;	  /* shift into the positive range */
					r2sqr = x*x + y*y;

					if (r2sqr <= atomRadius2) {

						if (muls->potential3D) {
							/* calculate the range which we have left to cover with z-variation */
							/* iRadZ is the number of slices (rounded up) that this atom
							* will contribute to, given its current x,y-radius
							*/ 
							iRadZ = (int)(sqrt(atomRadius2-r2sqr)/(*muls).cz[0]+1.0);
							/* loop through the slices that this atoms contributes to */
							for (iaz=-iRadZ;iaz <=iRadZ;iaz++) {
								if ((*muls).nonPeriodZ) {
									if (iaz+iAtomZ < 0)  {
										if (-iAtomZ <= iRadZ) iaz = -iAtomZ;
										else break;
										if (abs(iaz)>nlayer) break;
									} 
									if (iaz+iAtomZ >= nlayer)	break;
								}
								z = (float_tt)(iAtomZ+iaz+0.5)*(*muls).cz[0]-atomZ;
								/* shift into the positive range */
								iz = (iaz+iAtomZ+32*nlayer) % nlayer;	  
								/* x,y,z is the true vector from the atom center
								* We can look up the proj potential at that spot
								* using trilinear extrapolation.
								*/
								atomBoxLookUp(&dPot,muls,atoms[iatom].Znum,x,y,z,
									muls->tds ? 0 : atoms[iatom].dw);
								//    printf("access: %d %d %d\n",iz,ix,iy);
								m_trans[iz][ix][iy][0] += dPot[0];
								m_trans[iz][ix][iy][1] += dPot[1]; 	
							} /* end of for iaz=-iRadZ .. iRadZ */
						} /* end of if potential3D */

						/********************************************************************/ 

						else { /* if 2D potential */
							if ((*muls).nonPeriodZ) {
								if (iAtomZ < 0)  break;			
								if (iAtomZ >= nlayer)	break;	
							}		 
							iz = (iAtomZ+32*nlayer) % nlayer;	  /* shift into the positive range */
							atomBoxLookUp(&dPot,muls,atoms[iatom].Znum,x,y,0,
								muls->tds ? 0 : atoms[iatom].dw);
							z = (float_tt)(iAtomZ+1)*(*muls).cz[0]-atomZ;

							/* 
							printf("iz=%d,%d, z=%g, atomZ=%g(%g), iAtomZ=%d, c=%g (%d, %d))\n",
							iz,iatom,z,atomZ,atoms[iatom].z,iAtomZ,c, divCount,(*muls).cellDiv);
							*/
							/* split the atom if it is close to the top edge of the slice */
							//    printf("access: %d %d %d\n",iz,ix,iy);
							if ((z<0.15*(*muls).cz[0]) && (iz >0)) {
								m_trans[iz][ix][iy][0] += 0.5*dPot[0];
								m_trans[iz][ix][iy][1] += 0.5*dPot[1]; 
								m_trans[iz-1][ix][iy][0] += 0.5*dPot[0];
								m_trans[iz-1][ix][iy][1] += 0.5*dPot[1];			
							}
							/* split the atom if it is close to the bottom edge of the slice */
							else {
								if ((z>0.85*(*muls).cz[0]) && (iz < nlayer-1)) {
									m_trans[iz][ix][iy][0] += 0.5*dPot[0];
									m_trans[iz][ix][iy][1] += 0.5*dPot[1];	
									m_trans[iz+1][ix][iy][0] += 0.5*dPot[0];
									m_trans[iz+1][ix][iy][1] += 0.5*dPot[1]; 		
								}
								else {
									m_trans[iz][ix][iy][0] += dPot[0];
									m_trans[iz][ix][iy][1] += dPot[1];	
								}
							}
							/*
							if ((iatom % 100 == 0) && (fabs(x) <0.04) && (fabs(y) < 0.04))
							printf("i: %d iz: %d atomZ: %g (%g)\n",iatom,iz,atomZ,z);
							*/ 	  
						}
					}
				}
			}
		}


		/**************************************************************************
		* Newer, even faster method based on FFT of tabulated scattering factors
		**************************************************************************/ 
		else {	/* fftpotential */
			iAtomX = (int)floor(atomX/dx);	
			iAtomY = (int)floor(atomY/dy);
			if (muls->potential3D) atomZ+=muls->sliceThickness;  // why ??? !!!!!
			// printf("%d: pos=[%d, %d, %.1f]\n",iatom,iAtomX,iAtomY,atomZ);
			// atomZ is z-distance with respect to the start of the current stack of slices.
			// ddz = atomZ-dz*iAtomZ;


			/////////////////////////////////////////////////////////////////////////////
			// if we need to cut away at the edges, i.e. non-periodic potential arrays:
			if (muls->nonPeriod) {
				if (muls->potential3D) {
					iAtomZ = (int)floor(atomZ/muls->sliceThickness+0.5);
					// printf("iAtomZ: %d\n",iAtomZ);

					iax0 = iAtomX-iRadX <  0 ? 0 : iAtomX-iRadX;
					iax1 = iAtomX+iRadX >= muls->potNx ? muls->potNx-1 : iAtomX+iRadX;
					iay0 = iAtomY-iRadY <  0 ? 0 : iAtomY-iRadY;
					iay1 = iAtomY+iRadY >= muls->potNy ? muls->potNy-1 : iAtomY+iRadY;
					// if within the potential map range:
					if ((iax0 <  muls->potNx) && (iax1 >= 0) && (iay0 <  muls->potNy) && (iay1 >= 0)) {


						// define range of sampling from atomZ-/+atomRadius
						iaz0 = iAtomZ-iRadZ <  0 ? -iAtomZ : -iRadZ;
						iaz1 = iAtomZ+iRadZ >= muls->slices ?  muls->slices-iAtomZ-1 : iRadZ;
						// iaz0 = 0;  iaz1 = 0;
						// printf("iatomZ: %d, %d..%d cz=%g, %g (%d), dOffsZ=%g (%d)\n",iAtomZ,iaz0,iaz1,muls->sliceThickness,atomZ,(int)atomZ,(iAtomZ+iaz0-atomZ/muls->sliceThickness)*nzSub,(int)(iAtomZ+iaz0-atomZ/muls->sliceThickness)*nzSub+0.5);
						if ((iAtomZ+iaz0 <	muls->slices) && (iAtomZ+iaz1 >= 0)) {
							// retrieve the pointer for this atom
							atPotPtr     = getAtomPotential3D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw,&nzSub,&Nr,&Nz_lut);
#if USE_Q_POT_OFFSETS
							// retrieve the pointer to the array of charge-dependent potential offset
							// This function will return NULL; if the charge of this atom is zero:
							atPotOffsPtr = getAtomPotentialOffset3D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw,&nzSub,&Nr,&Nz_lut,atoms[iatom].q);
#endif // USE_Q_POT_OFFSETS
							iOffsLimHi   =  Nr*(Nz_lut-1);
							iOffsLimLo   = -Nr*(Nz_lut-1);
							iOffsStep    = nzSub*Nr;

							// Slices around the slice that this atom is located in must be affected by this atom:
							// iaz must be relative to the first slice of the atom potential box.
							for (iax=iax0; iax <= iax1; iax++) {
								potPtr = &(m_trans[iAtomZ+iaz0][iax][iay0][0]);
								// potPtr = &(m_trans[iAtomZ-iaz0+iaz][iax][iay0][0]);
								// printf("access: %d %d %d (%d)\n",iAtomZ+iaz0,iax,iay0,(int)potPtr);							

								//////////////////////////////////////////////////////////////////////
								// Computation of Radius must be made faster by using pre-calculated ddx
								// and LUT for sqrt:
								// use of sqrt slows down from 120sec to 180 sec.
								x2 = iax*dx - atomX;  x2 *= x2;
								for (iay=iay0; iay <= iay1; iay++) {
									// printf("iax=%d, iay=%d\n",iax,iay); 
									y2 = iay*dy - atomY;  y2 *= y2;
									r = sqrt(x2+y2);
									// r = (x2+y2);
									ddr = r/dr;
									ir	= (int)floor(ddr);
									// add in different slices, once r has been defined
									if (ir < Nr-1) {
										ddr = ddr-(float_tt)ir;
										ptr = potPtr;

										dOffsZ = (iAtomZ+iaz0-atomZ/muls->sliceThickness)*nzSub;
#if Z_INTERPOLATION
										iOffsZ = (int)dOffsZ;
										ddz    = fabs(dOffsZ - (float_tt)iOffsZ);
#else
										iOffsZ = (int)(dOffsZ+0.5);
#endif
										iOffsZ *= Nr;

										for (iaz=iaz0; iaz <= iaz1; iaz++) {
											potVal = 0;
											// iOffsZ = (int)(fabs(iAtomZ+iaz-atomZ/muls->sliceThickness)*nzSub+0.5);
											if (iOffsZ < 0) {
												if (iOffsZ > iOffsLimLo) {
													// do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
													potVal = (1-ddz)*((1-ddr)*atPotPtr[ir-iOffsZ+Nr][0]+ddr*atPotPtr[ir+1-iOffsZ+Nr][0])+
														ddz *((1-ddr)*atPotPtr[ir-iOffsZ   ][0]+ddr*atPotPtr[ir+1-iOffsZ   ][0]);
#if USE_Q_POT_OFFSETS
													// add the charge-dependent potential offset
													if (atPotOffsPtr != NULL) {
														potVal += atoms[iatom].q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir-iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1-iOffsZ+Nr][0])+
															ddz *((1-ddr)*atPotOffsPtr[ir-iOffsZ   ][0]+ddr*atPotOffsPtr[ir+1-iOffsZ   ][0]));		
													}	
#endif  // USE_Q_POT_OFFSETS
#else   // Z_INTERPOLATION
													potVal = (1-ddr)*atPotPtr[ir-iOffsZ+Nr][0]+ddr*atPotPtr[ir+1-iOffsZ+Nr][0];
#if USE_Q_POT_OFFSETS
													// add the charge-dependent potential offset
													if (atPotOffsPtr != NULL) {
														potVal += atoms[iatom].q*((1-ddr)*atPotOffsPtr[ir-iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1-iOffsZ+Nr][0]);
													}	
#endif  // USE_Q_POT_OFFSETS
#endif  // Z_INTERPOLATION
													// if (r < dr) printf("%d: iaz=%d,pot=%g,ddr=%g, iOffsZ=%d\n",iatom,iaz,potVal,ddr,iOffsZ);												
												}
											} // if iOffZ < 0
											else {
												// select the pointer to the right layer in the lookup box
												// printf("%4d: iOffsZ: %d, iaz: %d (slice: %d, pos: %g [%d .. %d])\n",iatom,iOffsZ,iaz,iAtomZ+iaz,atomZ,iaz0,iaz1);
												if (iOffsZ < iOffsLimHi) {
													// do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
													potVal = (1-ddz)*((1-ddr)*atPotPtr[ir+iOffsZ][0]+ddr*atPotPtr[ir+1+iOffsZ][0])+
														ddz *((1-ddr)*atPotPtr[ir+iOffsZ+Nr][0]+ddr*atPotPtr[ir+1+iOffsZ+Nr][0]);
#if USE_Q_POT_OFFSETS
													// add the charge-dependent potential offset
													if (atPotOffsPtr != NULL) {
														potVal += atoms[iatom].q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0])+
															ddz *((1-ddr)*atPotOffsPtr[ir+iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1+iOffsZ+Nr][0]));
													}
#endif  // USE_Q_POT_OFFSETS
#else   // Z_INTERPOLATION
													potVal = (1-ddr)*atPotPtr[ir+iOffsZ][0]+ddr*atPotPtr[ir+1+iOffsZ][0];
#if USE_Q_POT_OFFSETS
													// add the charge-dependent potential offset
													if (atPotOffsPtr != NULL) {
														potVal += atoms[iatom].q*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0]);												
													}
#endif  // USE_Q_POT_OFFSETS
#endif  // Z_INTERPOLATION
													// if (r < dr) printf("%d: iaz=%d,pot=%g,ddr=%g, iOffsZ=%d\n",iatom,iaz,potVal,ddr,iOffsZ);

												}
											} // if iOffsZ >=0
											*ptr += potVal;  // ptr = potPtr = m_trans[...]

											ptr  += sliceStep;	// advance to the next slice
											// add the remaining potential to the next slice:
											// if (iaz < iaz1)	*ptr += (1-ddz)*potVal;
											iOffsZ += iOffsStep;
										} // end of iaz-loop
									} // if ir < Nr
									// advance pointer to next complex potential point:
									potPtr+=2;
									// make imaginary part zero for now
									// *potPtr = 0;
									// potPtr++;

									// wrap around when end of y-line is reached:
									// if (++iay % muls->potNy == 0) potPtr -= 2*muls->potNy;		
								} // iay=iay0 .. iay1	  
							} // iax=iax0 .. iax1
						} // iaz0+iAtomZ < muls->slices
						// dOffsZ = (iAtomZ-atomZ/muls->sliceThickness)*nzSub;
						// printf("%5d (%2d): iAtomZ=%2d, offsZ=%g, diff=%g, (%g)\n",
						//	  iatom,atoms[iatom].Znum,iAtomZ,dOffsZ,dOffsZ - (int)(dOffsZ+0.5),iAtomZ-atomZ/muls->sliceThickness);
					} // if within bounds	
				}  // muls->potential3D and non-periodic in x-y
				////////////////////////////////////////////////////////////////////
				// 2D potential calculation already seems to work!
				// However, the potential is calculated wrongly, it is not integrated 
				// in z-direction.	This must be done in the potential slice initialization
				// procedure.
				else {
					iAtomZ = (int)floor(atomZ/muls->sliceThickness);
					iax0 = iAtomX-iRadX <  0 ? 0 : iAtomX-iRadX;
					iax1 = iAtomX+iRadX >= muls->potNx ? muls->potNx-1 : iAtomX+iRadX;
					iay0 = iAtomY-iRadY <  0 ? 0 : iAtomY-iRadY;
					iay1 = iAtomY+iRadY >= muls->potNy ? muls->potNy-1 : iAtomY+iRadY;
					// if within the potential map range:
					if ((iax0 <  muls->potNx) && (iax1 >= 0) && (iay0 <  muls->potNy) && (iay1 >= 0)) {
						ddx = (-(float_tt)iax0+(atomX/dx-(float_tt)iRadX))*(float_tt)OVERSAMP_X;
						ddy = (-(float_tt)iay0+(atomY/dy-(float_tt)iRadY))*(float_tt)OVERSAMP_X;
						iOffsX = (int)floor(ddx);
						iOffsY = (int)floor(ddy);
						ddx -= (float_tt)iOffsX;
						ddy -= (float_tt)iOffsY;
						s11 = (1-ddx)*(1-ddy);
						s12 = (1-ddx)*ddy;
						s21 = ddx*(1-ddy);
						s22 = ddx*ddy;
						atPotPtr = getAtomPotential2D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw);

						for (iax=iax0; iax < iax1; iax++) {
							// printf("(%d, %d): %d,%d\n",iax,nyAtBox,(iOffsX+OVERSAMP_X*(iax-iax0)),iOffsY+iay1-iay0);
							// potPtr and ptr are of type (float *)
							potPtr = &(m_trans[iAtomZ][iax][iay0][0]);
							ptr = &(atPotPtr[(iOffsX+OVERSAMP_X*(iax-iax0))*nyAtBox+iOffsY][0]);
							for (iay=iay0; iay < iay1; iay++) {
								*potPtr += s11*(*ptr)+s12*(*(ptr+2))+s21*(*(ptr+nyAtBox2))+s22*(*(ptr+nyAtBox2+2));

								potPtr++;
								// *potPtr = 0;
								potPtr++;
								ptr += 2*OVERSAMP_X;
							}
						}

					}  // not muls->potential3D
				}
			} // end of if not periodic

			/////////////////////////////////////////////////////////////////////////////
			// if the potential array is periodic:
			else {


				// printf("Z=%d (z=%d), iOffs: %d, %d (%g, %g) %g, %g, %g %g, atom: %g, %g\n",
				//		atoms[iatom].Znum,iAtomZ,iOffsX,iOffsY,ddx,ddy,s11,s12,s21,s22,atomX,atomY);
				////////////////////////////////////////////////////////////////////
				// add code here!
				if (muls->potential3D) {
					iAtomZ = (int)floor(atomZ/muls->sliceThickness+0.5);
					iax0 = iAtomX-iRadX;
					iax1 = iAtomX+iRadX; 
					iay0 = iAtomY-iRadY;
					iay1 = iAtomY+iRadY;

					// define range of sampling from atomZ-/+atomRadius
					iaz0 = iAtomZ-iRadZ <  0 ? -iAtomZ : -iRadZ;
					iaz1 = iAtomZ+iRadZ >= muls->slices ?  muls->slices-iAtomZ-1 : iRadZ;
					// if (iatom < 2) printf("iatomZ: %d, %d cz=%g, %g: %d, %d\n",iAtomZ,iaz0,muls->sliceThickness,atomZ,(int)(-2.5-atomZ),(int)(atomZ+2.5));
					// printf("%d: iatomZ: %d, %d cz=%g, %g\n",iatom,iAtomZ,iaz0,muls->sliceThickness,atomZ);
					if ((iAtomZ+iaz0 <  muls->slices) && (iAtomZ+iaz1 >= 0)) {
						// retrieve the pointer for this atom
						atPotPtr = getAtomPotential3D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw,&nzSub,&Nr,&Nz_lut);
#if USE_Q_POT_OFFSETS
						// retrieve the pointer to the array of charge-dependent potential offset
						// This function will return NULL; if the charge of this atom is zero:
						atPotOffsPtr = getAtomPotentialOffset3D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw,&nzSub,&Nr,&Nz_lut,atoms[iatom].q);
#endif // USE_Q_POT_OFFSETS
						iOffsLimHi =	Nr*(Nz_lut-1);
						iOffsLimLo = -Nr*(Nz_lut-1);
						iOffsStep  = nzSub*Nr;

						// Slices around the slice that this atom is located in must be affected by this atom:
						// iaz must be relative to the first slice of the atom potential box.
						for (iax=iax0; iax < iax1; iax++) {
							potPtr = &(m_trans[iAtomZ+iaz0][(iax+2*muls->potNx) % muls->potNx][(iay0+2*muls->potNy) % muls->potNy][0]);
							// potPtr = &(m_trans[iAtomZ-iaz0+iaz][iax][iay0][0]);
							x2 = iax*dx - atomX;	x2 *= x2;
							for (iay=iay0; iay < iay1; ) {
								// printf("iax=%d, iay=%d\n",iax,iay); 
								y2 = iay*dy - atomY;	y2 *= y2;
								r = sqrt(x2+y2);
								ddr = r/dr;
								ir  = (int)floor(ddr);
								// add in different slices, once r has been defined
								if (ir < Nr-1) {
									ddr = ddr-(float_tt)ir;
									ptr = potPtr;
									// Include interpolation in z-direction as well (may do it in a very smart way here !):

									dOffsZ = (iAtomZ+iaz0-atomZ/muls->sliceThickness)*nzSub;
#if Z_INTERPOLATION
									iOffsZ = (int)dOffsZ;
									ddz	 = fabs(dOffsZ - (float_tt)iOffsZ);
#else  // Z_INTERPOLATION
									iOffsZ = (int)(dOffsZ+0.5);
#endif  // Z_INTERPOLATION
									iOffsZ *= Nr;

									for (iaz=iaz0; iaz <= iaz1; iaz++) {
										potVal = 0;
										// iOffsZ = (int)(fabs(iAtomZ+iaz-atomZ/muls->sliceThickness)*nzSub+0.5);
										if (iOffsZ < 0) {
											if (iOffsZ > iOffsLimLo) {
												// do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
												potVal = (1-ddz)*((1-ddr)*atPotPtr[ir-iOffsZ+Nr][0]+ddr*atPotPtr[ir+1-iOffsZ+Nr][0])+
													ddz *((1-ddr)*atPotPtr[ir-iOffsZ	 ][0]+ddr*atPotPtr[ir+1-iOffsZ	 ][0]);
#if USE_Q_POT_OFFSETS
												// add the charge-dependent potential offset
												if (atPotOffsPtr != NULL) {
													potVal += atoms[iatom].q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir-iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1-iOffsZ+Nr][0])+
														ddz *((1-ddr)*atPotOffsPtr[ir-iOffsZ   ][0]+ddr*atPotOffsPtr[ir+1-iOffsZ   ][0]));		
												}	
#endif  // USE_Q_POT_OFFSETS
#else  // Z_INTERPOLATION
												potVal = (1-ddr)*atPotPtr[ir-iOffsZ+Nr][0]+ddr*atPotPtr[ir+1-iOffsZ+Nr][0];
#if USE_Q_POT_OFFSETS
												// add the charge-dependent potential offset
												if (atPotOffsPtr != NULL) {
													potVal += atoms[iatom].q*((1-ddr)*atPotOffsPtr[ir-iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1-iOffsZ+Nr][0]);
												}	
#endif  // USE_Q_POT_OFFSETS
#endif  // Z_INTERPOLATION
											}
										}
										else {
											// select the pointer to the right layer in the lookup box
											// printf("%4d: iOffsZ: %d, iaz: %d (slice: %d, pos: %g [%d .. %d])\n",iatom,iOffsZ,iaz,iAtomZ+iaz,atomZ,iaz0,iaz1);
											if (iOffsZ < iOffsLimHi) {
												// do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
												potVal = (1-ddz)*((1-ddr)*atPotPtr[ir+iOffsZ][0]+ddr*atPotPtr[ir+1+iOffsZ][0])+
													ddz *((1-ddr)*atPotPtr[ir+iOffsZ+Nr][0]+ddr*atPotPtr[ir+1+iOffsZ+Nr][0]);
#if USE_Q_POT_OFFSETS
												// add the charge-dependent potential offset
												if (atPotOffsPtr != NULL) {
													potVal += atoms[iatom].q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0])+
														ddz *((1-ddr)*atPotOffsPtr[ir+iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1+iOffsZ+Nr][0]));
												}
#endif  // USE_Q_POT_OFFSETS
#else  // Z_INTERPOLATION
												potVal = (1-ddr)*atPotPtr[ir+iOffsZ][0]+ddr*atPotPtr[ir+1+iOffsZ][0];
#if USE_Q_POT_OFFSETS
												// add the charge-dependent potential offset
												if (atPotOffsPtr != NULL) {
													potVal += atoms[iatom].q*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0]);												
												}
#endif  // USE_Q_POT_OFFSETS
#endif  // Z_INTERPOLATION
											}
										}
										*ptr += potVal;

										ptr  += sliceStep;  // advance to the next slice
										// add the remaining potential to the next slice:
										// if (iaz < iaz1)  *ptr += (1-ddz)*potVal;
										iOffsZ += iOffsStep;
									}
								} // if ir < Nr-1
								// advance pointer to next complex potential point:
								potPtr+=2;
								// make imaginary part zero for now
								// *potPtr = 0;
								// potPtr++;

								// wrap around when end of y-line is reached:
								if (++iay % muls->potNy == 0) potPtr -= 2*muls->potNy;		
							}   
						}
					} // iaz0+iAtomZ < muls->slices
				}  // muls->potential3D	
				////////////////////////////////////////////////////////////////////
				// 2D potential (periodic) calculation already seems to work!
				else {
					iAtomZ = (int)floor(atomZ/muls->sliceThickness);
					iax0 = iAtomX-iRadX+2*muls->potNx;
					iax1 = iAtomX+iRadX+2*muls->potNx; 
					iay0 = iAtomY-iRadY+2*muls->potNy;
					iay1 = iAtomY+iRadY+2*muls->potNy;

					ddx = (-(float_tt)iax0+(atomX/dx-(float_tt)(iRadX-2*muls->potNx)))*(float_tt)OVERSAMP_X;
					ddy = (-(float_tt)iay0+(atomY/dy-(float_tt)(iRadY-2*muls->potNy)))*(float_tt)OVERSAMP_X;
					iOffsX = (int)floor(ddx);
					iOffsY = (int)floor(ddy);
					ddx -= (float_tt)iOffsX;
					ddy -= (float_tt)iOffsY;
					/*
					s11 = (1-ddx)*(1-ddy);
					s12 = (1-ddx)*ddy;
					s21 = ddx*(1-ddy);
					s22 = ddx*ddy;
					*/									
				    s22 = (1-ddx)*(1-ddy);
					s21 = (1-ddx)*ddy;
					s12 = ddx*(1-ddy);
					s11 = ddx*ddy;

					atPotPtr = getAtomPotential2D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw);

					// if (iatom < 3) printf("atom #%d: ddx=%g, ddy=%g iatomZ=%d, atomZ=%g, %g\n",iatom,ddx,ddy,iAtomZ,atomZ,atoms[iatom].z);
					for (iax=iax0; iax < iax1; iax++) {  // TODO: should use ix += OVERSAMP_X
						// printf("(%d, %d): %d,%d\n",iax,nyAtBox,(iOffsX+OVERSAMP_X*(iax-iax0)),iOffsY+iay1-iay0);
						// potPtr and ptr are of type (float *)
						//////////////////
						// Only the exact slice that this atom is located in is affected by this atom:
						int atPosX = (OVERSAMP_X*(iax-iax0)-iOffsX);
						if ((atPosX >= 0) && (atPosX < nyAtBox-1)) {
							ptr = &(atPotPtr[atPosX*nyAtBox-iOffsY][0]);
						for (iay=iay0; iay < iay1; iay++) {
							// wrap around when end of y-line is reached:
								int atPosY = (iay-iay0)*OVERSAMP_X-iOffsY; 
								if ((atPosY < nyAtBox-1) && (atPosY >=0)) {
							// do the real part
									m_trans[iAtomZ][iax % muls->potNx][iay % muls->potNy][0] +=
									     s11*(*ptr)+s12*(*(ptr+2))+s21*(*(ptr+nyAtBox2))+s22*(*(ptr+nyAtBox2+2));
								}
							// make imaginary part zero for now
								// *potPtr = 0;  potPtr++;
							ptr += 2*OVERSAMP_X;
						}
						} // if atPosX within limits
					}   
				}  // muls->potential3D	== 0
			}
			////////////////////////////////////////////////////////////////////
		} /* end of if (fftpotential) */
	} /* for iatom =0 ... */
	time(&time1);
	if (iatom > 0)
	if (m_printLevel) printf("%g sec used for real space potential calculation (%g sec per atom)\n",difftime(time1,time0),difftime(time1,time0)/iatom);
	else
	if (m_printLevel) printf("%g sec used for real space potential calculation\n",difftime(time1,time0));


	/*************************************************/
	/* Save the potential slices					   */

	if (muls->savePotential) {
		for (iz = 0;iz<nlayer;iz++){
			/*
			muls->thickness = iz;
			showCrossSection(muls,m_transr[iz],nx,1,0);
			*/	
			// find the maximum value of each layer:
			potVal = m_trans[iz][0][0][0];
			for (ddx=potVal,ddy = potVal,ix=0;ix<muls->potNy*muls->potNx;potVal = m_trans[iz][0][++ix][0]) {
				if (ddy<potVal) ddy = potVal; 
				if (ddx>potVal) ddx = potVal; 
			}

#ifndef WIN32
			sprintf(filename,"%s/%s%d.img",muls->folder,muls->fileBase,iz);
#else
			sprintf(filename,"%s\\%s%d.img",muls->folder,muls->fileBase,iz);
#endif
			if (m_printLevel >= 3)
				printf("Saving (complex) potential layer %d to file %s (r: %g..%g)\n",iz,filename,ddx,ddy); 

                        params["Thickness"]=muls->sliceThickness;
			sprintf(buf,"Projected Potential (slice %d)",iz);		 
			imageIO->WriteComplexImage((void **)m_trans[iz],filename, params, std::string(buf));
		} // loop through all slices
	} /* end of if savePotential ... */
	if (muls->saveTotalPotential) {
		if (tempPot == NULL) tempPot = float2D(muls->potNx,muls->potNy,"total projected potential");

		for (ix=0;ix<muls->potNx;ix++) for (iy=0;iy<muls->potNy;iy++) {
			tempPot[ix][iy] = 0;
			for (iz=0;iz<nlayer;iz++) tempPot[ix][iy] += m_trans[iz][ix][iy][0];
		}

		for (ddx=tempPot[0][0],ddy = potVal,ix=0;ix<muls->potNy*muls->potNx;potVal = tempPot[0][++ix]) {
			if (ddy<potVal) ddy = potVal; 
			if (ddx>potVal) ddx = potVal; 
		}
#ifndef WIN32
		sprintf(filename,"%s/%sProj",muls->folder,muls->fileBase);	
#else
		sprintf(filename,"%s\\%sProj",muls->folder,muls->fileBase);	
#endif
		if (m_printLevel >= 2)
			printf("Saving total projected potential to file %s (r: %g..%g)\n",filename,ddx,ddy); 
                params["Thickness"]=muls->sliceThickness;
		sprintf(buf,"Projected Potential (sum of %d slices)",muls->slices);
		imageIO->WriteRealImage((void **)tempPot, fileName, params, std::string(buf));
	}

} // end of make3DSlices



/********************************************************************************
* Create Lookup table for 3D potential due to neutral atoms
********************************************************************************/
#define PHI_SCALE 47.87658
complex_tt *getAtomPotential3D(int Znum, MULS *muls,float_tt B,int *nzSub,int *Nr,int *Nz_lut) {
	int ix,iy,iz,iiz,ind3d,iKind,izOffset;
	float_tt zScale,kzmax,zPos,xPos;
#if FLOAT_PRECISION == 1
	fftwf_plan plan;
#else
	fftw_plan plan;
#endif
	static float_tt f,phase,s2,s3,kmax2,smax2,kx,kz,dkx,dky,dkz; // ,dx2,dy2,dz2;
	static int nx,ny,nz,nzPerSlice;
	static complex_tt **atPot = NULL;
	static complex_tt *temp = NULL;
#if SHOW_SINGLE_POTENTIAL == 1
	ImageIOPtr imageio = ImageIOPtr();
	static complex_tt *ptr = NULL;
	static char fileName[256];
#endif 
	static double *splinb=NULL;
	static double *splinc=NULL;
	static double *splind=NULL;


	// scattering factors in:
	// float scatPar[4][30]
	if (atPot == NULL) {
		splinb = double1D(N_SF, "splinb" );
		splinc = double1D(N_SF, "splinc" );
		splind = double1D(N_SF, "splind" );


		nx = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionX);
		ny = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionY);
		// The FFT-resolution in the z-direction must be high enough to avoid 
		// artifacts due to premature cutoff of the rec. space scattering factor 
		// we will therefore make it roughly the same as the x-resolution
		// However, we will make sure that a single slice contains an integer number 
		// of sampling points.
		nzPerSlice = (int)floor(OVERSAMP_X*muls->sliceThickness/muls->resolutionX);
		// make nzPerSlice odd:
		if (2.0*(nzPerSlice >> 1) == nzPerSlice) nzPerSlice += 1;
		// Total number of z-positions should be twice that of atomRadius/sliceThickness 
		nz = (2*(int)ceil(muls->atomRadius/muls->sliceThickness))*nzPerSlice;
		if (m_printLevel > 1) printf("Will use %d sampling points per slice, total nz=%d (%d)\n",nzPerSlice,nz,nzPerSlice >> 1);

		dkx = 0.5*OVERSAMP_X/(nx*muls->resolutionX);  // nx*muls->resolutionX is roughly 2*muls->atomRadius
		dky = dkx;                                    
		dkz = nzPerSlice/(double)(nz*muls->sliceThickness);
		kmax2 = 0.5*nx*dkx/(double)OVERSAMP_X;  // largest k that we'll admit
		smax2 = kmax2;

		printf("dkx = %g, nx = %d, kmax2 = %g\n",dkx,nx,kmax2);
		if (m_printLevel > 1) printf("Cutoff scattering angle: kmax=%g (1/A), dk=(%g,%g %g)\n",kmax2,dkx,dky,dkz);
		scatPar[0][N_SF-1] = 1.2*smax2;
		scatPar[0][N_SF-2] = 1.1*smax2;
		scatPar[0][N_SF-3] = smax2;
		// adjust the resolution of the lookup table if necessary
		if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3]) {

			if (1) {
				// set additional scattering parameters to zero:
				for (ix = 0;ix < N_SF-10;ix++) {
					if (scatPar[0][N_SF-4-ix] < scatPar[0][N_SF-3]-0.001*(ix+1)) break;
					scatPar[0][N_SF-4-ix] = scatPar[0][N_SF-3]-0.001*(ix+1);
					for (iy=1; iy<N_ELEM;iy++) scatPar[iy][N_SF-4-ix] = 0; 
				}
			}
			else {
				for (ix = 0;ix < 20;ix++) {
					if (scatPar[0][N_SF-4-ix] < scatPar[0][N_SF-3]) break;
					scatPar[0][N_SF-4-ix] = scatPar[0][N_SF-3];
					for (iy=1; iy<N_ELEM;iy++) scatPar[iy][N_SF-4-ix] = scatPar[iy][N_SF-3];	
				}
			}		
			if (m_printLevel > 1) printf("getAtomPotential3D: set resolution of scattering factor to %g/A!\n",
				scatPar[0][N_SF-4-ix]);
		}	// end of if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3])
		smax2 *= smax2;
		kmax2 *= kmax2;
		// allocate a list of pointers for the element-specific potential lookup table
		atPot = (fftwf_complex **)malloc((NZMAX+1)*sizeof(fftwf_complex *));
		for (ix=0;ix<=NZMAX;ix++) atPot[ix] = NULL;
		temp  = (fftwf_complex*) fftwf_malloc(nx*nz*sizeof(fftwf_complex));
	}
	// initialize this atom, if it has not been done yet:
	if (atPot[Znum] == NULL) {
		iKind = Znum;

		// setup cubic spline interpolation:
		splinh(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF);

		// allocate a 3D array:
		atPot[Znum] = (complex_tt*) fftw_malloc(nx*nz/4*sizeof(complex_tt));
		memset(temp,0,nx*nz*sizeof(complex_tt));
		kzmax	  = dkz*nz/2.0; 
		// define x-and z-position of atom center:
		// The atom should sit in the top-left corner, 
		// however (nzPerSlice+1)/2 above zero in z-direction
		xPos = -2.0*PI*0.0;  // or muls->resolutionX*nx/(OVERSAMP_X), if in center
		izOffset = (nzPerSlice-1)/2;
		zPos = -2.0*PI*(muls->sliceThickness/nzPerSlice*(izOffset));

		// What this look-up procedure will do is to supply V(r,z) computed from fe(q).
		// Since V(r,z) is rotationally symmetric we might as well compute 
		// V(x,y,z) at y=0, i.e. V(x,z).  
		// In order to do the proper 3D inverse FT without having to do a complete 3D FFT
		// we will pre-compute the qy-integral for y=0.

		// kzborder = dkz*(nz/(2*OVERSAMP_Z) -1); 
		for (iz=0;iz<nz;iz++) {
			kz = dkz*(iz<nz/2 ? iz : iz-nz);	
			// We also need to taper off the potential in z-direction
			// in order to avoid cutoff artifacts.
			// zScale = fabs(kz) <= kzborder ? 1.0 : 0.5+0.5*cos(M_PI*(fabs(kz)-kzborder)/(kzmax-kzborder));
			// printf("iz=%d, kz=%g, zScale=%g ",iz,kz,zScale);
			for (ix=0;ix<nx;ix++) {
				kx = dkx*(ix<nx/2 ? ix : ix-nx);	   
				s2 = (kx*kx+kz*kz);
				// if this is within the allowed circle:
				if (s2<smax2) {
					ind3d = ix+iz*nx;
					// f = fe3D(Znum,k2,muls->tds,1.0,muls->scatFactor);
					// multiply scattering factor with Debye-Waller factor:
					// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
					f = seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF,sqrt(s2))*exp(-s2*B*0.25);
					// perform the qy-integration for qy <> 0:
					for (iy=1;iy<nx;iy++) {
						s3 = dkx*iy;
						s3 = s3*s3+s2;
						if (s3<smax2) {
							f += 2*seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF,sqrt(s3))*exp(-s3*B*0.25);
						}
						else break;
					}
					f *= dkx;  
					// note that the factor 2 is missing in the phase (2pi k*r)
					// this places the atoms in the center of the box.
					phase	= kx*xPos + kz*zPos;
					temp[ind3d][0] = f*cos(phase);  // *zScale
					temp[ind3d][1] = f*sin(phase);  // *zScale
					// if ((kx==0) && (ky==0)) printf(" f=%g (%g, [%g, %g])\n",f,f*zScale,atPot[Znum][ind3d][0],atPot[Znum][ind3d][1]);
				}
			}
		} // for iz ...


#if SHOW_SINGLE_POTENTIAL
		// 0 thickness
		imageio = ImageIOPtr(new CImageIO(nz, nx, 0, dkz, dkx));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		sprintf(fileName,"pot_rec_%d.img",Znum);
		imageio->SetThickness(muls->sliceThickness);
		imageio->WriteComplexImage((void**)temp, fileName);
#endif	  
		// This converts the 2D kx-kz  map of the scattering factor to a 2D real space map.
#if FLOAT_PRECISION ==1
		plan = fftwf_plan_dft_2d(nz,nx,temp,temp,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
#else
                plan = fftw_plan_dft_2d(nz,nx,temp,temp,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
#endif
		// dx2 = muls->resolutionX*muls->resolutionX/(OVERSAMP_X*OVERSAMP_X);
		// dy2 = muls->resolutionY*muls->resolutionY/(OVERSAMP_X*OVERSAMP_X);
		// dz2 = muls->sliceThickness*muls->sliceThickness/(OVERSAMP_Z*OVERSAMP_Z);
		// We also make sure that the potential touches zero at least somewhere.  This will avoid 
		// sharp edges that could produce ringing artifacts.  
		// It is certainly debatable whether this is a good apprach, or not. 
		// printf("Setting up %d x %d potential for Znum=%d, (atom kind %d), Omega=%g\n",nx,nz,Znum,iKind,dkx*dky*dkz);
		// min = atPot[Znum][ny/2+nx/2*ny+nz/2*nx*ny][0];
		for (ix=0;ix<nx/2;ix++)  for (iz=0;iz<nz/2;iz++) {
			ind3d = ix+iz*nx/2;
			// Integrate over nzPerSlice neighboring layers here:::::::::
			for (zScale=0,iiz=-izOffset;iiz<=izOffset;iiz++) {
				if (iz+izOffset+iiz < nz/2) zScale += temp[ix+(iz+izOffset+iiz)*nx][0];
			}
			if (zScale < 0) zScale = 0;
			// assign the iz-th slice the sum of the 3 other slices:
			// and divide by unit cell volume (since this is in 3D):
			// Where does the '1/2' come from???  OVERSAMP_X*OVERSAMP_Y/8 = 1/2
			// if nothing has changed, then OVERSAMP_X=2 OVERSAMP_Z=18.
			// remember, that s=0.5*k; 	
			// This potential will later again be scaled by lambda*gamma (=0.025*1.39139)
			atPot[Znum][ind3d][0] = 47.8658*dkx*dkz/(nz)*zScale; 

			// *8*14.4*0.529=4*a0*e (s. Kirkland's book, p. 207)
			// 2*pi*14.4*0.529 = 7.6176;
			// if (atPot[Znum][ind3d][0] < min) min = atPot[Znum][ind3d][0];	  
			atPot[Znum][ind3d][1]= 0;
		}
		// make sure we don't produce negative potential:
		// if (min < 0) for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) atPot[Znum][iy+ix*ny][0] -= min;
#if SHOW_SINGLE_POTENTIAL
		imageio = ImageIOPtr(new CImageIO(nz/2, nx/2, 0, muls->sliceThickness/nzPerSlice, 
			muls->resolutionX/OVERSAMP_X));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		imageio->SetThickness(nz*muls->sliceThickness/nzPerSlice);
		sprintf(fileName,"potential_rz_%d.img",Znum);
		ptr = atPot[Znum];
		imageio->WriteComplexImage((void**)ptr, fileName);
#endif	  
		if (m_printLevel > 1) printf("Created 3D (r-z) %d x %d potential array for Z=%d (%d, B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)\n",
			nx/2,nz/2,Znum,iKind,B,dkx,dky,dkz,izOffset);
	}
	*Nz_lut = nz/2;
	*nzSub  = nzPerSlice;
	*Nr	  = nx/2;
	return atPot[Znum];
}



/********************************************************************************
* Lookup function for 3D potential offset due to charged atoms (ions)
********************************************************************************/
fftwf_complex *getAtomPotentialOffset3D(int Znum, MULS *muls,double B,int *nzSub,int *Nr,int *Nz_lut,float q) {
	int ix,iy,iz,iiz,ind3d,iKind,izOffset;
	double zScale,kzmax,zPos,xPos;
	fftwf_plan plan;
	static double f,phase,s2,s3,kmax2,kx,kz,dkx,dky,dkz; // ,dx2,dy2,dz2;
	static int nx,ny,nz,nzPerSlice;
	static fftwf_complex **atPot = NULL;
	static fftwf_complex *temp = NULL;
#if SHOW_SINGLE_POTENTIAL == 1
	ImageIOPtr imageio = ImageIOPtr();
	static fftwf_complex *ptr = NULL;
	static char fileName[256];
#endif 
	static double *splinb=NULL;
	static double *splinc=NULL;
	static double *splind=NULL;


	// if there is no charge to this atom, return NULL:
	if (q == 0) return NULL;


	// scattering factors in:
	// float scatPar[4][30]
	if (atPot == NULL) {
		splinb = double1D(N_SF, "splinb" );
		splinc = double1D(N_SF, "splinc" );
		splind = double1D(N_SF, "splind" );


		nx = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionX);
		ny = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionY);
		// The FFT-resolution in the z-direction must be high enough to avoid 
		// artifacts due to premature cutoff of the rec. space scattering factor 
		// we will therefore make it roughly the same as the x-resolution
		// However, we will make sure that a single slice contains an integer number 
		// of sampling points.
		nzPerSlice = (int)floor(OVERSAMP_X*muls->sliceThickness/muls->resolutionX);
		// make nzPerSlice odd:
		if (2.0*floor((double)(nzPerSlice >> 1)) == nzPerSlice) nzPerSlice += 1;
		// Total number of z-positions should be twice that of atomRadius/sliceThickness 
		nz = (2*(int)ceil(muls->atomRadius/muls->sliceThickness))*nzPerSlice;
		if (m_printLevel > 1) printf("Potential offset: will use %d sampling points per slice, total nz=%d (%d)\n",nzPerSlice,nz,nzPerSlice >> 1);

		dkx = 0.5*OVERSAMP_X/(nx*muls->resolutionX);  
		dky = 0.5*OVERSAMP_X/(ny*muls->resolutionY);
		dkz = nzPerSlice/(double)(nz*muls->sliceThickness);
		kmax2 = 0.5*nx*dkx/(double)OVERSAMP_X;	// largest k that we'll admit

		// printf("Cutoff scattering angle:kmax=%g, smax=%g (1/A), dk=(%g,%g %g)\n",kmax2,S_SCALE*kmax2,dkx,dky,dkz);
		scatParOffs[0][N_SF-1] = 1.2*kmax2;
		scatParOffs[0][N_SF-2] = 1.1*kmax2;
		scatParOffs[0][N_SF-3] = kmax2;
		if (scatParOffs[0][N_SF-4] > scatParOffs[0][N_SF-3]) {

			if (1) {
				// set additional scattering parameters to zero:
				for (ix = 0;ix < N_SF-10;ix++) {
					if (scatParOffs[0][N_SF-4-ix] < scatParOffs[0][N_SF-3]-0.001*(ix+1)) break;
					scatParOffs[0][N_SF-4-ix] = scatParOffs[0][N_SF-3]-0.001*(ix+1);
					for (iy=1; iy<N_ELEM;iy++) scatParOffs[iy][N_SF-4-ix] = 0;	
				}
			}
			else {
				for (ix = 0;ix < 20;ix++) {
					if (scatParOffs[0][N_SF-4-ix] < scatParOffs[0][N_SF-3]) break;
					scatParOffs[0][N_SF-4-ix] = scatParOffs[0][N_SF-3];
					for (iy=1; iy<N_ELEM;iy++) scatParOffs[iy][N_SF-4-ix] = scatParOffs[iy][N_SF-3];	
				}
			}	   
			if (m_printLevel > 1) printf("getAtomPotentialOffset3D: reduced angular range of scattering factor to %g/A!\n",
				scatParOffs[0][N_SF-4-ix]);
		}  // end of if (scatParOffs[0][N_SF-4] > scatParOffs[0][N_SF-3])
		kmax2 *= kmax2;

		atPot = (fftwf_complex **)malloc((NZMAX+1)*sizeof(fftwf_complex *));
		for (ix=0;ix<=NZMAX;ix++) atPot[ix] = NULL;
		temp  = (fftwf_complex*)fftwf_malloc(nx*nz*sizeof(fftwf_complex));
	}
	// initialize this atom, if it has not been done yet:
	if (atPot[Znum] == NULL) {
#if USE_REZ_SFACTS
		iKind = Znum;
#else
		printf("Using charged atoms only works with scattering factors by Rez et al!\n",Znum);
		exit(0);
#endif


		// setup cubic spline interpolation:
		splinh(scatParOffs[0],scatParOffs[iKind],splinb,splinc,splind,N_SF);

		atPot[Znum] = (fftwf_complex*)fftwf_malloc(nx*nz/4*sizeof(fftwf_complex));
		memset(temp,0,nx*nz*sizeof(fftwf_complex));
		kzmax	 = dkz*nz/2.0; 
		// define x-and z-position of atom center:
		// The atom should sit in the top-left corner, 
		// however (nzPerSlice+1)/2 above zero in z-direction
		xPos = -2.0*PI*0.0;  // or muls->resolutionX*nx/(OVERSAMP_X), if in center
		izOffset = (nzPerSlice-1)/2;
		zPos = -2.0*PI*(muls->sliceThickness/nzPerSlice*(izOffset));

		// kzborder = dkz*(nz/(2*OVERSAMP_Z) -1); 
		for (iz=0;iz<nz;iz++) {
			kz = dkz*(iz<nz/2 ? iz : iz-nz);   
			// We also need to taper off the potential in z-direction
			// in order to avoid cutoff artifacts.
			// zScale = fabs(kz) <= kzborder ? 1.0 : 0.5+0.5*cos(M_PI*(fabs(kz)-kzborder)/(kzmax-kzborder));
			// printf("iz=%d, kz=%g, zScale=%g ",iz,kz,zScale);
			for (ix=0;ix<nx;ix++) {
				kx = dkx*(ix<nx/2 ? ix : ix-nx);	  
				s2 = (kx*kx+kz*kz);
				// if this is within the allowed circle:
				if (s2<kmax2) {
					ind3d = ix+iz*nx;
					// f = fe3D(Znum,k2,muls->tds,1.0,muls->scatFactor);
					// multiply scattering factor with Debye-Waller factor:
					// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
					f = seval(scatParOffs[0],scatParOffs[iKind],splinb,splinc,splind,N_SF,sqrt(s2))*exp(-s2*B*0.25);
					// perform the qy-integration for qy <> 0:
					for (iy=1;iy<nx;iy++) {
						s3 = dky*iy;
						s3 = s3*s3+s2;
						if (s3<kmax2) {
							f += 2*seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF,sqrt(s3))*exp(-s3*B*0.25);
						}
						else break;
					}
					f *= dkx;  
					// note that the factor 2 is missing in the phase (2pi k*r)
					// this places the atoms in the center of the box.
					phase  = kx*xPos + kz*zPos;
					temp[ind3d][0] = f*cos(phase);	// *zScale
					temp[ind3d][1] = f*sin(phase);	// *zScale
					// if ((kx==0) && (ky==0)) printf(" f=%g (%g, [%g, %g])\n",f,f*zScale,atPot[Znum][ind3d][0],atPot[Znum][ind3d][1]);
				}
			}
		} // for iz ...




#if SHOW_SINGLE_POTENTIAL
		imageio = ImageIOPtr(new CImageIO(nz, nx, 0, dkx, dkz, std::vector<double>(), 
			"rec. space potential"));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		imageio->SetThickness(muls->sliceThickness);
		imageio->WriteComplexImage((void**)temp, fileName);
#endif	  

		plan = fftwf_plan_dft_2d(nz,nx,temp,temp,FFTW_BACKWARD,FFTW_ESTIMATE);
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
		// dx2 = muls->resolutionX*muls->resolutionX/(OVERSAMP_X*OVERSAMP_X);
		// dy2 = muls->resolutionY*muls->resolutionY/(OVERSAMP_X*OVERSAMP_X);
		// dz2 = muls->sliceThickness*muls->sliceThickness/(OVERSAMP_Z*OVERSAMP_Z);
		// We also make sure that the potential touches zero at least somewhere.  This will avoid 
		// sharp edges that could produce ringing artifacts.  
		// It is certainly debatable whether this is a good apprach, or not. 
		// printf("Setting up %d x %d potential for Znum=%d, (atom kind %d), Omega=%g\n",nx,nz,Znum,iKind,dkx*dky*dkz);
		// min = atPot[Znum][ny/2+nx/2*ny+nz/2*nx*ny][0];
		for (ix=0;ix<nx/2;ix++)  for (iz=0;iz<nz/2;iz++) {
			ind3d = ix+iz*nx/2;
			// Integrate over nzPerSlice neighboring layers here:::::::::
			for (zScale=0,iiz=-izOffset;iiz<=izOffset;iiz++) {
				if (iz+izOffset+iiz < nz/2) zScale += temp[ix+(iz+izOffset+iiz)*nx][0];
			}
			if (zScale < 0) zScale = 0;
			// assign the iz-th slice the sum of the 3 other slices:
			// and divide by unit cell volume (since this is in 3D):
			// Where does the '1/2' come from???  OVERSAMP_X*OVERSAMP_Y/8 = 1/2
			// if nothing has changed, then OVERSAMP_X=2 OVERSAMP_Z=18.
			// remember, that s=0.5*k;	   
			// This potential will later again be scaled by lambda*gamma (=0.025*1.39139)
			atPot[Znum][ind3d][0] = 47.8658*dkx*dkz/(nz)*zScale; 
			atPot[Znum][ind3d][1] = 0;
		}
#if SHOW_SINGLE_POTENTIAL
		imageio = ImageIOPtr(new CImageIO(nz/2, nx/2, 0, muls->resolutionX/OVERSAMP_X, 
		muls->sliceThickness/nzPerSlice));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		imageio->SetThickness(nz*muls->sliceThickness/nzPerSlice);
		sprintf(fileName,"potentialOffs_rz_%d.img",Znum);
		ptr = atPot[Znum];
		imageio->WriteComplexImage((void**)ptr, fileName);
#endif	  
		if (m_printLevel > 1) printf("Created 3D (r-z) %d x %d potential offset array for Z=%d (%d, B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)\n",
			nx/2,nz/2,Znum,iKind,B,dkx,dky,dkz,izOffset);
	}
	*Nz_lut = nz/2;
	*nzSub = nzPerSlice;
	*Nr	 = nx/2;
	return atPot[Znum];
}
// #undef SHOW_SINGLE_POTENTIAL
// end of fftwf_complex *getAtomPotential3D(...)

#define PHI_SCALE 47.87658
// #define SHOW_SINGLE_POTENTIAL 0
////////////////////////////////////////////////////////////////////////////
// This function should be used yet, because it computes the projected
// potential wrongly, since it doe not yet perform the projection!!!
fftwf_complex *getAtomPotential2D(int Znum, MULS *muls,double B) {
	int ix,iy,iz,ind,iKind;
	double min;
	fftwf_plan plan;
	static double f,phase,s2,s3,kmax2,kx,ky,dkx,dky;
	static int nx,ny;
	static fftwf_complex **atPot = NULL;
#if SHOW_SINGLE_POTENTIAL == 1
	ImageIOPtr imageio = ImageIOPtr();
	static char fileName[256];
#endif 
	static double *splinb=NULL;
	static double *splinc=NULL;
	static double *splind=NULL;


	// scattering factors in:
	// float scatPar[4][30]
	if (atPot == NULL) {
		splinb = double1D(N_SF, "splinb" );
		splinc = double1D(N_SF, "splinc" );
		splind = double1D(N_SF, "splind" );

		nx = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionX);
		ny = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionY);
		dkx = 0.5*OVERSAMP_X/((nx)*muls->resolutionX);  
		dky = 0.5*OVERSAMP_X/((ny)*muls->resolutionY);
		kmax2 = 0.5*nx*dkx/(double)OVERSAMP_X;  // largest k that we'll admit

		printf("Cutoff scattering angle:kmax=%g (1/A)\n",kmax2);
		scatPar[0][N_SF-1] = 1.2*kmax2;
		scatPar[0][N_SF-2] = 1.1*kmax2;
		scatPar[0][N_SF-3] = kmax2;
		if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3]) {

			if (1) {
				// set additional scattering parameters to zero:
				for (ix = 0;ix < N_SF-10;ix++) {
					if (scatPar[0][N_SF-4-ix] < scatPar[0][N_SF-3]-0.001*(ix+1)) break;
					scatPar[0][N_SF-4-ix] = scatPar[0][N_SF-3]-0.001*(ix+1);
					for (iy=1; iy<N_ELEM;iy++) scatPar[iy][N_SF-4-ix] = 0;	
				}
			}

			if (m_printLevel > 1) printf("getAtomPotential2D: reduced angular range of scattering factor to %g/A!\n",
				scatPar[0][N_SF-4-ix]);
		}  // end of if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3])
		kmax2 *= kmax2;



		atPot = (fftwf_complex **)malloc((NZMAX+1)*sizeof(fftwf_complex *));
		for (ix=0;ix<=NZMAX;ix++) atPot[ix] = NULL;
	}
	// initialize this atom, if it has not been done yet:
	if (atPot[Znum] == NULL) {
		iKind = Znum;
		// setup cubic spline interpolation:
		splinh(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF);

		atPot[Znum] = (fftwf_complex*) fftwf_malloc(nx*ny*sizeof(fftwf_complex));
		// memset(temp,0,nx*nz*sizeof(fftwf_complex));
		memset(atPot[Znum],0,nx*ny*sizeof(fftwf_complex));
		for (ix=0;ix<nx;ix++) {
			kx = dkx*(ix<nx/2 ? ix : nx-ix);      
			for (iy=0;iy<ny;iy++) {
				ky = dky*(iy<ny/2 ? iy : ny-iy);      
				s2 = (kx*kx+ky*ky);
				// if this is within the allowed circle:
				if (s2<kmax2) {
					ind = iy+ix*ny;
					// f = fe3D(Znum,k2,muls->tds,1.0,muls->scatFactor);
					// multiply scattering factor with Debye-Waller factor:
					// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
					f = seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF,sqrt(s2))*exp(-s2*B*0.25);
					phase = PI*(kx*muls->resolutionX*nx+ky*muls->resolutionY*ny);
					atPot[Znum][ind][0] = f*cos(phase);
					atPot[Znum][ind][1] = f*sin(phase);
				}
			}
		}
#if SHOW_SINGLE_POTENTIAL == 1
		imageio = ImageIOPtr(new CImageIO(ny, nx, 0, dkx, dky, std::vector<double>(), 
		"potential"));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		imageio->SetThickness(muls->sliceThickness);
		imageio->WriteComplexImage((void**)atPot[Znum], fileName);
#endif    
		plan = fftwf_plan_dft_2d(nx,ny,atPot[Znum],atPot[Znum],FFTW_BACKWARD,FFTW_ESTIMATE);
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
		for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
				atPot[Znum][iy+ix*ny][0] *= dkx*dky*(OVERSAMP_X*OVERSAMP_X);  
		}
		// make sure we don't produce negative potential:
		// if (min < 0) for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) atPot[Znum][iy+ix*ny][0] -= min;
#if SHOW_SINGLE_POTENTIAL == 1
		imageio = ImageIOPtr(new CImageIO(nx, ny, 0, muls->resolutionX/OVERSAMP_X, 
			muls->resolutionY/OVERSAMP_X, std::vector<double>(), "potential"));
		// This scattering factor agrees with Kirkland's scattering factor fe(q)
		//imageio->SetThickness(nz*muls->sliceThickness/nzPerSlice);
		sprintf(fileName,"potential_%d.img",Znum);
		imageio->WriteComplexImage((void**)atPot[Znum], fileName);
#endif    
		printf("Created 2D %d x %d potential array for Z=%d (%d, B=%g A^2)\n",nx,ny,Znum,iKind,B);
	}
	return atPot[Znum];
}
#undef PHI_SCALE
#undef SHOW_SINGLE_POTENTIAL

/**************************************************************
* The imaginary part of the trans arrays is already allocated
* The projected potential is already located in trans[][][][0]
*
**************************************************************/
#define PHI_SCALE 47.87658
void initSTEMSlices(MULS *muls, int nlayer) {
	int printFlag = 0;
	int ilayer;
	int nbeams;
	double scale,vz,vzscale,mm0,wavlen;
	int nx,ny,ix,iy; // iz;
	float_tt temp,k2max,k2,kx,ky;
	static float_tt *kx2= NULL,*ky2 = NULL; /* *kx= NULL,*ky= NULL, */
	float_tt pi;
	double fftScale;
	double timer1,timer2,time2=0,time1=0;
	// char filename[32];

	pi = (float)PI;

	nx = (*muls).potNx;
	ny = (*muls).potNy;

	/**************************************************/
	/* Setup all the reciprocal lattice vector arrays */
	if ((kx2 == NULL)||(ky2 == NULL)) {
		kx2    = float1D(nx, "kx2" );
		ky2    = float1D(ny, "ky2" );
		/*
		kx     = float1D(nx, "kx" );
		ky     = float1D(ny, "ky" );
		*/
		for(ix=0; ix<nx; ix++) {
			kx = (ix>nx/2) ? (float_tt)(ix-nx)/(*muls).potSizeX : 
				(float_tt)ix/(*muls).potSizeX;
		kx2[ix] = kx*kx;
		/*
		kx[ix] = (ix>nx/2) ? (float_tt)(ix-nx)/(*muls).potSizeX : 
		(float_tt)ix/(*muls).potSizeX;
		kx2[ix] = kx[ix]*kx[ix];
		*/
		}
		for( iy=0; iy<ny; iy++) {
			ky = (iy>ny/2) ? 
				(float_tt)(iy-ny)/(*muls).by : (float_tt)iy/(*muls).potSizeY;
			ky2[iy] = ky*ky;
		}    
		if (m_printLevel > 2) printf("Reciprocal lattice vector arrays initialized ... \n");
	}
	/**************************************************/
	/* setup of all the necessary parameters: */
	wavlen = (float_tt) wavelength((*muls).v0);

	nbeams = 0;

	/* setup of k-vector arrays, etc... */
	k2max = nx/(2.0F*(*muls).potSizeX);
	temp = ny/(2.0F*(*muls).potSizeY);
	if( temp < k2max ) k2max = temp;
	k2max = (float_tt)(BW * k2max);
	if (printFlag) {
		printf("Bandwidth limited to a real space resolution of %f Angstroms\n",
			1.0F/k2max);
		printf("   (= %.2f mrad)  for symmetrical anti-aliasing.",
			wavlen*k2max*1000.0F);
		printf(" (a: %g, b: %g)\n",(*muls).potSizeX,(*muls).potSizeY);
	}
	k2max = k2max*k2max;
	(*muls).k2max = k2max;

	if(m_trans==NULL) {
		printf("Memory for trans has not been allocated\n");
		exit(0);
	}


	/**************************************************************
	* read in the potential layers and exponentiate them
	* to form the transmission function exp(-i*phi(x)) = transr/i
	* 
	* This is independent of kx, ky - so we can load these arrays once,
	* and keep them for all inc. beam angles
	*
	**************************************************************/


	mm0 = 1.0F + muls->v0/511.0F;  // relativistic corr. factor gamma
	// scale = mm0*sigma(muls->v0) / 1000.0;  /* in 1/(volt-Angstroms) */
	// we will need to increase lambda and decrease gamm by the mean inner crystal potential effect.
	scale = mm0*wavelength(muls->v0);
	if (m_printLevel > 1) printf("Making phase gratings for %d layers (scale=%g rad/VA, gamma=%g, sigma=%g) ... \n",nlayer,scale,mm0,sigma(muls->v0));


	fftScale = 1.0/(nx*ny);
	vzscale= 1.0;
	timer1 = cputim();    
	for( ilayer=0;  ilayer<nlayer; ilayer++ ) {     
		timer2 = cputim();    
		for( iy=0; iy<ny; iy++) for( ix=0; ix<nx; ix++) {
			vz= m_trans[ilayer][ix][iy][0]*scale;  // scale = lambda*gamma
			// include absorption:
			// vzscale= exp(-m_trans[ilayer][ix][iy][1]*scale);
			/* printf("vz(%d %d) = %g\n",ix,iy,vz); */
			m_trans[ilayer][ix][iy][0] =  cos(vz);
			m_trans[ilayer][ix][iy][1] =  sin(vz);
		}
	}

	/*******************************************************************
	* FFT/IFFT the transmit functions in order to bandwidth limit them
	*******************************************************************/ 
	if (muls->bandlimittrans) {
		timer2 = cputim();    
#if FLOAT_PRECISION == 1
		fftwf_execute(muls->fftPlanPotForw);
#else
		fftw_execute(muls->fftPlanPotForw);
#endif
		time2 = cputim()-timer2;
		//     printf("%g sec used for 1st set of FFTs\n",time2);  
		for( ilayer=0;  ilayer<nlayer; ilayer++ ) {     
			for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
				k2= ky2[iy] + kx2[ix];
				if (k2 < k2max) {
					nbeams++;
					m_trans[ilayer][ix][iy][0] *= fftScale;
					m_trans[ilayer][ix][iy][1] *= fftScale;
				}
				else {
					m_trans[ilayer][ix][iy][0] = 0.0F;
					m_trans[ilayer][ix][iy][1] = 0.0F;
				}	
			}
		}  /* end for(ilayer=... */
		timer2 = cputim();    
		// old code: fftwnd_one((*muls).fftPlanPotInv, m_trans[ilayer][0], NULL);
#if FLOAT_PRECISION == 1
		fftwf_execute(muls->fftPlanPotInv);
#else
		fftw_execute(muls->fftPlanPotInv);
#endif
		time2 += cputim()-timer2;
	}  /* end of ... if bandlimittrans */
	time1 = cputim()-timer1;

	if (m_printLevel > 1) {
		if ((*muls).bandlimittrans) {
			printf("%g sec used for making phase grating, %g sec for %d FFTs.\n",
				time1,time2,2*nlayer);
		}
		else
			printf("%g sec used for making phase grating, no bandwidth limiting\n",
			time1);

		if (printFlag)
			printf("Number of symmetrical non-aliasing beams = %d\n", nbeams);
	}  

	/*
	printf("Size in pixels Nx x Ny= %d x %d = %d beams\n",
	nx,ny, nx*ny);
	printf("Lattice constant a = %.4f, b = %.4f\n", (*muls).ax,(*muls).by);
	*/
}  // initSTEMSlices

#undef PHI_SCALE
