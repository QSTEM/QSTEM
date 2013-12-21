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

/* file: fileio_fftw3.c */

#include <stdio.h>	/*  ANSI-C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

// #include "../lib/floatdef.hpp"
#include "stemtypes_fftw3.hpp"
#include "memory_fftw3.hpp"	/* memory allocation routines */
// #include "tiffsubs.hpp"
#include "matrixlib.hpp"
#include "readparams.hpp"
#include "fileio_fftw3.hpp"
// #include "stemlib.hpp"

#include "readparams.hpp"

#include "elTable.hpp"

#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#ifdef WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif

#define NCINMAX  500	/* max number of characers in stacking spec */
#define NRMAX	50	/* number of values in look-up-table in vzatomLUT */
#define RMIN	0.01	/* min r (in Ang) range of LUT for vzatomLUT() */
#define RMAX	5

#define NPDTMAX 8       /* number of parameters for doyle turner sfacts */
#define NPMAX	12	/* number of parameters for each Z */
#define NZMIN	1	/* min Z in featom.tab */
//#define NZMAX	98      /* max Z (in featom.dat ZMAX=103) */
#define NPARAM	64	/* number of parameters in tiff files */

#define N_A 6.022e+23
#define K_B 1.38062e-23      /* Boltzman constant */
#define PID 3.14159265358979 /* pi */
#define SQRT_PI 1.77245385090552 /* sqrt(pi) */
#define HPLANCK 6.6262e-34   /* planck's constant */
#define AMU 1.66053e-27      /* atomic mass unit */
#define THZ_AMU_PI_HBAR 0.05012012918415 /*  A°^2*THz*amu/(pi*hbar) */
#define THZ_AMU_HBAR 0.15745702964189    /*   A°^2*THz*amu/(hbar)   */
// 4.46677327584453 /* 1e10/sqrt(THz*amu/(pi*hbar)) */ 
#define THZ_HBAR_2KB  3.81927135604119     /* THz*hbar/(2*kB) */
#define THZ_HBAR_KB   1.90963567802059     /* THz*hbar/kB */
#define AMU_THZ2_A2_KB   1.20274224623720     /* AMU*THz^2*A^2/kB */

#define MAX_MASS_INDEX 59
int idArraySize=0;
int *idArray = NULL;
int idArrayPtr = 0;
double massArray[MAX_MASS_INDEX]={1.008,                                         4.0026,
6.939,9.012,      10.81,12.01,14.01,16.00,19.00,20.18,
23.00,24.31,      26.98,28.09,30.97,32.06,35.45,39.95,
39.10,40.08,
44.96,47.90,50.94,52.00,54.94,55.85,58.93,58.71,63.55,65.38,
69.72,72.59,74.92,78.96,79.90,83.80,
85.47,87.62,88.91,91.22,92.91,95.94,98.91,10.07,102.9,106.4,107.9,112.4,
114.8,118.7,121.8,127.6,126.9,131.3,
132.9054,137.33,138.9055,178.49,180.9479};
/* so far this list goes up to Xe (from Gerthsen) .. Ta (webelememts) */
double chargeTable[MAX_MASS_INDEX];



// #define NCINMAX 500
// #define NPARAM	64    /* number of parameters */

//MCS - why do we have this function and initMuls in stem3.cpp?
MULS initMu() {
	MULS muls;
	int sCount,i,slices = 2;
	char waveFile[32];
	char *waveFileBase = "w";

	muls.slices = slices;

	/* general setup: */
	muls.lpartl = 0;

	// muls.wave = NULL;
	muls.atomRadius = 5.0;  /* radius in A for making the potential boxes */

	for (sCount =0;sCount<slices;sCount++)
		muls.cin2[sCount] = 'a'+sCount;
	for (sCount = slices;sCount < NCINMAX;sCount++)
		muls.cin2[sCount] = 0;
	muls.nlayer = slices;
	muls.saveFlag = 0;

	muls.sigmaf = 0;
	muls.dfdelt = 0;
	muls.acmax = 0;
	muls.acmin = 0;
	muls.aobj = 0;
	muls.aAIS = 0;
	// muls.areaAIS = 1.0;

	// Tomography parameters:
	muls.tomoTilt = 0;
	muls.tomoStart = 0;
	muls.tomoStep = 0;
	muls.tomoCount = 0;  // indicate: NO Tomography simulation.


	/* make multislice read the inout files and assign transr and transi: */
	muls.trans = NULL;
	muls.cz = NULL;  // (float_t *)malloc(muls.slices*sizeof(float_t));

	muls.onlyFresnel = 0;
	muls.showPhaseplate = 0;
	muls.czOffset = 0;  /* defines the offset for the first slice in 
						fractional coordinates        */
	muls.normHolog = 0;
	muls.gaussianProp = 0;

	muls.sparam = (float *)malloc(NPARAM*sizeof(float));
	for (i=0;i<NPARAM;i++)
		muls.sparam[i] = 0.0;

	/****************************************************/
	/* copied from slicecell.c                          */
	muls.pendelloesung = NULL;
	return muls;
}
// #undef NCINMAX 500
// #undef NPARAM	64    /* number of parameters */


#define CHARGE 0.0

void writeFrameWork(FILE *fp,superCellBox superCell) {
	int i,id = 0,newId = 1;
	double charge=0.0;

	if (fp != NULL) fprintf(fp,"frame1 1 framework\n");
	for (i=0;i<superCell.natoms;i++) {
		if (i==0) {
			if (idArray == NULL) {
				idArraySize = 2;
				idArray=(int *)malloc(idArraySize*sizeof(int));
				id = 0;
				idArrayPtr = 0;
			}
			else {
				/* check, if this Znum is already somewhere in the list, and if not,
				* add it to the list, extending its size, if necessary
				*/
				for (id=0;id<=idArrayPtr;id++) if (superCell.atoms[i].Znum == idArray[id]) break;
				if (id>idArrayPtr) {
					if (idArrayPtr >= idArraySize-1) {
						idArraySize *=2;
						idArray = (int *)realloc(idArray,idArraySize*sizeof(int));
					}
					idArray[++idArrayPtr] = superCell.atoms[i].Znum;
					newId = 1;
				}           
			}
			idArray[id] = superCell.atoms[i].Znum;
		}
		else {
			newId = 0;
			if (superCell.atoms[i].Znum != superCell.atoms[i-1].Znum) {
				/* check, if this Znum is already somewhere in the list, and if not,
				* add it to the list, extending its size, if necessary
				*/
				for (id=0;id<=idArrayPtr;id++) if (superCell.atoms[i].Znum == idArray[id]) break;
				if (id>idArrayPtr) {
					if (idArrayPtr >= idArraySize-1) {
						idArraySize *=2;
						idArray = (int *)realloc(idArray,idArraySize*sizeof(int));
					}
					idArray[++idArrayPtr] = superCell.atoms[i].Znum;
					newId = 1;
				}
			}
		}
		if (newId) {
			switch (superCell.atoms[i].Znum) {
	  case 14: charge = 4.0; break;
	  case 7:  charge = -3.0; break;
	  case 8:  charge = -2.0; break;
	  case 39: charge = 3.0; break;
	  default: charge = 0.0; 
			}

			if (fp != NULL) fprintf(fp,"%d %7.3f %7.3f %7.3f %6.3f %6.3f %c%c\n",id+1,superCell.atoms[i].x,
				superCell.atoms[i].y,superCell.atoms[i].z,
				massArray[superCell.atoms[i].Znum-1],charge,
				elTable[2*superCell.atoms[i].Znum-2],elTable[2*superCell.atoms[i].Znum-1]);
		}
		else
			if (fp != NULL) fprintf(fp,"%d %7.3f %7.3f %7.3f\n",id+1,superCell.atoms[i].x,
			superCell.atoms[i].y,superCell.atoms[i].z);
	}
}


#define A 0.0     // 975.857
#define B 43684.0   // 425379 
#define C 3.4483  // 5.04748

/* This function will write the amorphous data to the MD input file starting at 
* atoms nstart and ending just before atom nstop.   
*/
void writeAmorphous(FILE *fp,superCellBox superCell,int nstart,int nstop) {
	int i,j,id;
	int *idCountArray = NULL;
	double charge,b,x,y,z;
	// char elem[8];

	memset(chargeTable,0,MAX_MASS_INDEX*sizeof(double));

	printf("amorph: %d ..%d-1\n",nstart,nstop);

	if (idArray == NULL) {
		idArraySize = 2;
		idArray=(int *)malloc(idArraySize*sizeof(int));
		id = 0;
		idArrayPtr = 0;
		idArray[0] = superCell.atoms[nstart].Znum;
	}
	else {
		/* check, if this Znum is already somewhere in the list, and if not,
		* add it to the list, extending its size, if necessary
		*/
		for (id=0;id<=idArrayPtr;id++) if (superCell.atoms[nstart].Znum == idArray[id]) break;
		if (id>idArrayPtr) {
			if (idArrayPtr >= idArraySize-1) {
				idArraySize *=2;
				idArray = (int *)realloc(idArray,idArraySize*sizeof(int));
			}
			idArray[++idArrayPtr] = superCell.atoms[nstart].Znum;
		}           
	}
	idCountArray = (int *)malloc(idArraySize*sizeof(int));
	memset(idCountArray,0,idArraySize*sizeof(int));
	idCountArray[id] = 1;

	for (i=nstart+1;i<nstop;i++) {
		if (superCell.atoms[i].Znum != superCell.atoms[i-1].Znum) {
			/* check, if this Znum is already somewhere in the list, and if not,
			* add it to the list, extending its size, if necessary
			*/

			for (id=0;id<=idArrayPtr;id++) {
				if (superCell.atoms[i].Znum == idArray[id]) break;
			}
			if (id>idArrayPtr) {
				if (idArrayPtr >= idArraySize-1) {
					idArraySize *=2;
					idArray = (int *)realloc(idArray,idArraySize*sizeof(int));
					idCountArray = (int *)realloc(idCountArray,idArraySize*sizeof(int));
				}
				idArray[++idArrayPtr] = superCell.atoms[i].Znum;
			}
		}
		idCountArray[id] += 1;
	}

	/**************************************************************************************
	* Now that we found all the species present in the amorphous phase we can go ahead and 
	* make a list of them, together with the number of atoms present of each kind
	*/
	for (id=0;id<=idArrayPtr;id++) {
		if (idCountArray[id] > 0) { 
			if (idArray[id] >= MAX_MASS_INDEX) {
				printf("mass exceeds array!!! - extend array!\n");
				exit(0);
			}
			switch (idArray[id]) {
	  case 14: charge = 4.0; break;
	  case 7:  charge = -3.0; break;
	  case 8:  charge = -2.0; break;
	  case 39: charge = 3.0; break;
	  default: charge = 0.0; 
			}
			if (fp != NULL) {
				fprintf(fp,"%c%c %d\n",elTable[2*idArray[id]-2],elTable[2*idArray[id]-1],idCountArray[id]);
				fprintf(fp,"%d %7.3f %7.3f %7.3f %6.3f %6.3f %c%c\n",id+1,0.0,0.0,0.0,
					massArray[idArray[id]-1],charge,
					elTable[2*idArray[id]-2],elTable[2*idArray[id]-1]);    
			}
		}
	}
	if (fp != NULL) fprintf(fp,"end\n");
	/***********************************************************************
	* The list of the potential parameters is next
	*/
	if (fp != NULL) {
	fprintf(fp,"buckingham\n");
	for (id=0;id<=idArrayPtr;id++)
		for (j=id;j<=idArrayPtr;j++) {
			// The following parameters are from S. Garofalini, J. Am. Cer. Soc. 67, 133 (1987)
			switch (10*(id+1) + (j+1)) {
	  case 11: b = 0.0; break;  // N -N
	  case 12: b = 4.5; break;  // N -Si
	  case 13: b = 4.5; break;  // N -Y
	  case 14: b = 0.0; break;  // N -O
	  case 22: b = 1.8770; break;  // Si-Si
	  case 23: b = 4.4200; break;  // Si-Y
	  case 24: b = 2.9620; break;  // Si-O
	  case 33: b = 9.9706; break;  // Y -Y
	  case 34: b = 6.0802; break;  // Y -O
	  case 44: b = 0.7254; break;  // O -O		
	  default : b=0.0;
			}
			b = b*6.022e4;  // *1e-16*6.022e23*1e-3
			fprintf(fp,"%d %d %g %g %g\n",id+1,j+1,A,b,C);

		}
		fprintf(fp,"end\n");

		/*************************************************************************
		* We can now write down the size of the MD cell:
		*/
		fprintf(fp,"%g %g %g 90 90 90 1 1 1\n",superCell.ax,superCell.by,superCell.cz);

		/****************************************************************************
		* now list the position of each one of the atoms in the amorphous phase
		* in FRACTIONAL COORDINATES !!!
		*/
		for(i=nstart;i<nstop;i++) {   
			x = superCell.atoms[i].x/superCell.ax;
			y = superCell.atoms[i].y/superCell.by;
			z = superCell.atoms[i].z/superCell.cz;
			if (fabs(x-1.0) < 1e-5) x = 0.0;
			if (fabs(y-1.0) < 1e-5) y = 0.0;
			if (fabs(z-1.0) < 1e-5) z = 0.0;
			fprintf(fp,"%c%c %7.5f %7.5f %7.5f\n",elTable[2*superCell.atoms[i].Znum-2],
				elTable[2*superCell.atoms[i].Znum-1],x,y,z);
		}
		if (nstart > 0)
			fprintf(fp,"frame1 %7.5f %7.5f %7.5f 1 0 0 0\n",superCell.cmx,superCell.cmy,superCell.cmz);
		fprintf(fp,"end\n");
	}  // end of if fp != NULL
}
