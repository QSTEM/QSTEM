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

#include <stdio.h>	/*  ANSI-C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "stemlib.hpp"
#include "memory_fftw3.hpp"	/* memory allocation routines */
#include "stemutil.hpp"
// #include "tiffsubs.hpp"
#include "imagelib_fftw3.hpp"
#include "fileio_fftw3.hpp"
// #include "floatdef.hpp"
// #include "imagelib.hpp"


#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#ifdef WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif

#define OVERSAMP_X 2
#define OVERSAMP_Z 18

#define NSMAX 1000	/* max number of slices */
#define NLMAX	52	/* maximum number of layers */
#define NCINMAX  500	/* max number of characers in stacking spec */
#define NCMAX 256	/* max characters in file names */
#define NPARAM	64	/* number of parameters */
#define NZMIN	1	/* min Z in featom.tab */
#define NZMAX	103	/* max Z in featom.tab */
#define EXTRA_LAYERS 3  /* number of extra layers for potential overlap */

#define SUB_SLICES  5    /* number of sub slices per slice (for integration) */
#define POTENTIAL_3D
#define INTEGRAL_TOL 1e-5
#define MAX_INTEGRAL_STEPS 15
#define MIN_INTEGRAL_STEPS 2
/*#define USE_VATOM_LUT */ /* set if you want to use vzatomLUT/v3DatomLUT */
/*#define USE_VZATOM_IN_CENTER */
/////////////////////////////////////////////////
// for debugging:
#define SHOW_SINGLE_POTENTIAL 0
/////////////////////////////////////////////////



#define BUF_LEN 256
#define PI 3.14159265358979
#define USE_REZ_SFACTS    1  // used in getAtomPotential3D and getAtomPotentialOffset3D 
#define Z_INTERPOLATION   0  // used in make3DSlices (central function for producing atom potential slices)
#define USE_Q_POT_OFFSETS 1  // used in make3DSlices (central function for producing atom potential slices)


const char cname[] = "abcdefghijklmnopqrstuvwxyz"
"ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const double pid=PI;
const double pid2=2*PI;

// define all the scat facts (or just those for Sr, Ti, O, In, P, He,Cl, Si,Ca,Ba,Fe) from Rez et al / Doyle and Turner
// the final 3 values must be zero to achieve a nice tapering off
#if USE_REZ_SFACTS
// provide array variable names 
// - scatPar[N_ELEM][N_SF] and
// - scatParOffs[N_ELEM][N_SF]
// and also define N_SF and N_ELEM:
#include "scatfactsRez.hpp"
#else

#define N_SF 30
#define N_ELEM 14
// tabulated scattering factors from Doyle and Turner (copied by hand!)
double scatPar[N_ELEM][N_SF] = {{0.0000,0.0500,0.1000,0.1500,0.2000,0.2500,0.3000,0.3500,0.4000,0.4500,0.5000,
0.6000,0.7000,0.8000,0.9000,1.0000,1.2000,1.4000,1.6000,1.8000,2.0000,2.5000,
3.0000,3.5000,4.0000,5.0000,6.0000,80,90,100},
{13.1090,11.4760,8.4780,6.2000,4.7940,3.8820,3.2240,2.7180,2.3150,1.9930,1.7330,
1.3500,1.0890,0.9020,0.7620,0.6510,0.4880,0.3750,0.2970,0.2390,0.1940,0.1290,
0.0920,0.0690,0.0540,0.0350,0.0240,0,0,0},
{8.7760,7.9370,6.1990,4.6430,3.5640,2.8440,2.3410,1.9640,1.6680,1.4280,1.2300,
0.9300,0.7210,0.5730,0.4670,0.3890,0.2850,0.2190,0.1750,0.1430,0.1170,0.0780,
0.0550,0.0410,0.0310,0.0200,0.0140,0,0,0},
{1.9830,1.9370,1.8080,1.6250,1.4220,1.2220,1.0400,0.8810,0.7470,0.6350,0.5420,
0.4030,0.3070,0.2410,0.1930,0.1590,0.1130,0.0850,0.0660,0.0530,0.0440,0.0290,
0.0200,0.0150,0.0120,0.0080,0.0050,0,0,0},
{10.434,9.7680,8.2970,6.8050,5.6010,4.6820,3.9730,3.4120,2.9550,2.5760,2.2570,
1.7580,1.3970,1.1320,0.9360,0.7890,0.5870,0.4580,0.3680,0.3020,0.2500,0.1660,
0.1180,0.0880,0.0680,0.0345,0.0310,0,0,0},
{5.4880,5.1920,4.4570,3.5860,2.7960,2.1690,1.7020,1.3620,1.1150,0.9330,0.7970,
0.6100,0.4870,0.4010,0.3350,0.2840,0.2100,0.1600,0.1250,0.1000,0.0820,0.0530,
0.0370,0.0280,0.0210,0.0140,0.0100,0,0,0},
{0.4180,0.4100,0.3900,0.3590,0.3230,0.2860,0.2500,0.2170,0.1890,0.1640,0.1430,
0.1100,0.0860,0.0680,0.0550,0.0460,0.0320,0.0240,0.0190,0.0150,0.0120,0.0080,
0.0050,0.0040,0.0030,0.0020,0.0010,0,0,0},
{4.8570,4.6850,4.2270,3.6200,2.9970,2.4380,1.9740,1.6060,1.3190,1.0980,0.8280,
0.6920,0.5410,0.4400,0.3660,0.3110,0.2320,0.1780,0.1410,0.1130,0.0930,0.0600,
0.0420,0.0310,0.0240,0.0160,0.0110,0,0,0},
{5.8280,5.4210,4.4670,3.4370,2.5890,1.9690,1.5340,1.2310,1.0170,0.8610,0.7430,
0.5780,0.4650,0.3830,0.3200,0.2700,0.1980,0.1500,0.1170,0.0930,0.0760,0.0500,
0.0350,0.0260,0.0200,0.0130,0.0090,0,0,0},
{9.9130,8.7030,6.3880,4.5500,3.4080,2.6950,2.2060,1.8380,1.5480,1.3140,1.1230,
0.8380,0.6470,0.5150,0.4220,0.3540,0.2620,0.2020,0.1620,0.1320,0.1070,0.0710,
0.0500,0.0370,0.0280,0.0180,0.0130,0,0,0},
{18.267,15.854,11.675,8.6820,6.8290,5.5700,4.6280,3.8950,3.3180,2.8610,2.4940,
1.9510,1.5700,1.2880,1.0730,0.9040,0.6660,0.5110,0.4110,0.3370,0.2770,0.1890,
0.1340,0.1000,0.0780,0.0510,0.0360,0,0,0},
{7.1650,6.6690,5.5580,4.4360,3.5620,2.9280,2.4610,2.1040,1.8180,1.5840,1.3880,
1.0800,0.8540,0.6860,0.5610,0.4660,0.3360,0.2550,0.2020,0.1650,0.1360,0.0910,
0.0650,0.0480,0.0370,0.0240,0.0170,0,0,0},
// Al:
{5.8990,5.3710,4.2370,3.1280,2.2990,1.7370,1.3630,1.1110,0.9320,0.8010,0.7000,
0.5510,0.4450,0.3660,0.3040,0.2550,0.1850,0.1390,0.1090,0.0870,0.0700,0.0460,
0.0320,0.0240,0.0190,0.0120,0.0090,0,0,0},
// Y:
{12.307,10.968,8.3983,6.3131,4.9361,4.0087,3.3336,2.8186,2.4107,2.0840,1.8121,
1.4173,1.1500,0.9536,0.8049,0.6849,0.5114,0.3917,0.3098,0.2480,0.2031,0.1363,
0.0957,0.0727,0.0569,0.0369,0.0258,0,0,0}};
#endif  // USE_REZ_SFACTS

void writePix(char *outFile,complex_tt **pict,MULS *muls,int iz) {
	float_tt *sparam;
	float_tt rmin,rmax;
	int i,j, result;

	rmin  = pict[0][0][0];
	rmax  = rmin;

	sparam = float1D( NPARAM, "sparam" );    

	for( i=0; i<(*muls).nx; i++)
		for(j=0; j<(*muls).ny; j++) 
		{
			if(pict[i][j][0] < rmin ) rmin = pict[i][j][0];
			if(pict[i][j][0] > rmax ) rmax = pict[i][j][0];
		}
		printf("min: %g  max: %g\n",rmin,rmax);

		if (rmin==rmax)
			rmax = rmin+1.0;
		sparam[pRMAX]  = rmax;
		sparam[pIMAX]  = rmax;
		sparam[pRMIN]  = rmin;
		sparam[pIMIN]  = rmin;
		sparam[pXCTILT]  = 0.0f;
		sparam[pYCTILT] = 0.0f;
		sparam[pENERGY] = (*muls).v0;
		sparam[pDX] = (*muls).ax/(float_tt)(*muls).nx;
		sparam[pDY] = (*muls).by/(float_tt)(*muls).ny;
		sparam[pWAVEL] = (float_tt)wavelength((*muls).v0);
		sparam[pNSLICES] = 0.0F;  /* ??? */
		sparam[pDEFOCUS] = 0.0;
		sparam[pOAPERT] = 0.0;
		sparam[pCS] = 0.0;
		sparam[pCAPERT] = 0.0;
		sparam[pDDF] = 0.0;
		sparam[pC]= (*muls).cz[iz];
		result = 1;
		/*
		result = tcreateFloatPixFile(outFile,pict,(long)(*muls).nx,
		(long)(*muls).ny,1,sparam);
		*/

		if (result != 1)
			printf("\ncould not write output file %s\n",outFile);
}




/******************************************************************
* runMulsSTEM() - do the multislice propagation in STEM/CBED mode
* 
*    Each probe position is running this function.  Each CPU is thus
*      running a separate instance of the function.  It is nested in
*      the main OpenMP parallel region - specifying critical, single, and
*      barrier OpenMP pragmas should be OK.
*
* waver, wavei are expected to contain incident wave function 
* they will be updated at return
*****************************************************************/
int runMulsSTEM(MULS *muls, WavePtr wave, PotPtr pot) {
	int printFlag = 0; 
	int showEverySlice=1;
	int islice,i,ix,iy,mRepeat;
	float_tt cztot=0.0;
	float_tt wavlen,scale,sum=0.0; //,zsum=0.0
	// static int *layer=NULL;
	float_tt x,y;
	int absolute_slice;

	char outStr[64];
	double fftScale;

	printFlag = (muls->printLevel > 3);
	fftScale = 1.0/(muls->nx*wave->m_ny);

	wavlen = (float_tt)wavelength((*muls).v0);

	/*  calculate the total specimen thickness and echo */
	cztot=0.0;
	for( islice=0; islice<(*muls).slices; islice++) {
		cztot += (*muls).cz[islice];
	}
	if (printFlag)
		printf("Specimen thickness: %g Angstroms\n", cztot);

	scale = 1.0F / (((float_tt)wave->m_nx) * ((float_tt)wave->m_ny));

	for (mRepeat = 0; mRepeat < muls->mulsRepeat1; mRepeat++) 
	{
		for( islice=0; islice < muls->slices; islice++ ) 
		{
			absolute_slice = (muls->totalSliceCount+islice);

			/***********************************************************************
			* Transmit is a simple multiplication of wave with trans in real space
			**********************************************************************/
			wave->Transmit(pot, islice);   
			/***************************************************** 
			* remember: prop must be here to anti-alias
			* propagate is a simple multiplication of wave with prop
			* but it also takes care of the bandwidth limiting
			*******************************************************/
#if FLOAT_PRECISION == 1
			fftwf_execute(wave->m_fftPlanWaveForw);
#else
			fftw_execute(wave->m_fftPlanWaveForw);
#endif
			wave->Propagate();
			//propagate_slow(wave, muls->nx, muls->ny, muls);

            muls->detectors->CollectIntensity(wave, muls->totalSliceCount+islice*(1+mRepeat));
			//collectIntensity(muls, wave, muls->totalSliceCount+islice*(1+mRepeat));

			if (muls->mode != STEM) {
				/* write pendelloesung plots, if this is not STEM */
				writeBeams(muls,wave,islice, absolute_slice);
			}

			// go back to real space:
#if FLOAT_PRECISION == 1
			fftwf_execute(wave->fftPlanWaveInv);
#else
			fftw_execute(wave->fftPlanWaveInv);
#endif
			// old code: fftwnd_one((*muls).fftPlanInv,(complex_tt *)wave[0][0], NULL);
			fft_normalize((void **)wave->wave,muls->nx,muls->ny);

			// write the intermediate TEM wave function:

			/********************************************************************
			* show progress:
			********************************************************************/
			wave->thickness = (absolute_slice+1)*muls->sliceThickness;
			if ((printFlag)) {
				sum = 0.0;
				for( ix=0; ix<(*muls).nx; ix++)  for( iy=0; iy<(*muls).ny; iy++) {
					sum +=  wave->wave[ix][iy][0]* wave->wave[ix][iy][0] +
						wave->wave[ix][iy][1]* wave->wave[ix][iy][1];
				}
				sum *= scale;

				sprintf(outStr,"position (%3d, %3d), slice %4d (%.2f), int. = %f", 
					wave->detPosX, wave->detPosY,
					muls->totalSliceCount+islice,wave->thickness,sum );
				if (showEverySlice)
					printf("%s\n",outStr);
				else {
					printf("%s",outStr);
					for (i=0;i<(int)strlen(outStr);i++) printf("\b");
				}
			}
			
			if ((muls->mode == TEM) || ((muls->mode == CBED)&&(muls->saveLevel > 1))) 
			{
                          // TODO (MCS 2013/04): this restructure probably broke this file saving - 
                          //   need to rewrite a function to save things for TEM/CBED?
                          // This used to call interimWave(muls,wave,muls->totalSliceCount+islice*(1+mRepeat));
                          interimWave(muls,wave,absolute_slice*(1+mRepeat)); 
                          muls->detectors->CollectIntensity(wave, absolute_slice*(1+mRepeat));

                          //collectIntensity(muls,wave,absolute_slice*(1+mRepeat));
			}
		} /* end for(islice...) */
		// collect intensity at the final slice
		//collectIntensity(muls, wave, muls->totalSliceCount+muls->slices*(1+mRepeat));
	} /* end of mRepeat = 0 ... */
	if (printFlag) printf("\n***************************************\n");

	/****************************************************
	****************************************************
	**           this -WAS- the big loop              **
	****************************************************
	***************************************************/

	// TODO: modifying shared value from multiple threads?
	//#pragma omp single
	muls->rmin  = wave->wave[0][0][0];
	//#pragma omp single
	muls->rmax  = (*muls).rmin;
	//#pragma omp single
	muls->aimin = wave->wave[0][0][1];
	//#pragma omp single
	muls->aimax = (*muls).aimin;

	sum = 0.0;
	for( ix=0; ix<muls->nx; ix++)  for( iy=0; iy<muls->ny; iy++) {
		x =  wave->wave[ix][iy][0];
		y =  wave->wave[ix][iy][1];
		if( x < (*muls).rmin ) (*muls).rmin = x;
		if( x > (*muls).rmax ) (*muls).rmax = x;
		if( y < (*muls).aimin ) (*muls).aimin = y;
		if( y > (*muls).aimax ) (*muls).aimax = y;
		sum += x*x+y*y;
	}
	// TODO: modifying shared value from multiple threads?
	//  Is this sum supposed to be across multiple pixels?
	//#pragma omp critical
	wave->intIntensity = sum*scale;

	if (printFlag) {
		printf( "pix range %g to %g real,\n"
			"          %g to %g imag\n",  
			(*muls).rmin,(*muls).rmax,(*muls).aimin,(*muls).aimax);

	}
	if (muls->saveFlag) {
		if ((muls->saveLevel > 1) || (muls->cellDiv > 1)) {
			wave->WriteWave();
			if (printFlag)
                          printf("Created complex image file %s\n",(*wave).fileout.c_str());
		}
	}
	return 0;
}  // end of runMulsSTEM


////////////////////////////////////////////////////////////////
// save the current wave function at this intermediate thickness:
void interimWave(MULS *muls,WavePtr wave,int slice) {
	int t;
	char fileName[256]; 
        std::map<std::string, double> params;

	if ((slice < muls->slices*muls->cellDiv-1) && ((slice+1) % muls->outputInterval != 0)) return;

	t = (int)((slice)/muls->outputInterval);
	// write the high tension, too:
        
        params["HT"] = muls->v0;
        params["Cs"] = muls->Cs;
        params["Defocus"] = muls->df0;
        params["Astigmatism Magnitude"] = muls->astigMag;
        params["Astigmatism Angle"] = muls->astigAngle;
        params["Focal Spread"] = muls->Cc * sqrt(muls->dE_E*muls->dE_E+muls->dV_V*muls->dV_V+muls->dI_I*muls->dI_I);
        params["Convergence Angle"] = muls->alpha;
        params["Beam Tilt X"] = muls->btiltx;
        params["Beam Tilt Y"] = muls->btilty;
	
	// produce the following filename:
	// wave_avgCount_thicknessIndex.img or
	// wave_thicknessIndex.img if tds is turned off
	if (muls->tds) wave->WriteWave(muls->avgCount, t, "Wave Function", params);
        else wave->WriteWave(t, "Wave Function", params);
}

/*****  saveSTEMImages *******/
// Saves all detector images (STEM images) that are defined in muls.
//   When saving intermediate STEM images is enabled, this also saves
//   the intermediate STEM images for each detector.
void saveSTEMImages(MULS *muls)
{
  std::map<std::string, double> params;
  params["Runs Averaged"]=(double)muls->avgCount+1;
  muls->detectors->SaveDetectors(params);  
}


void fft_normalize(void **array,int nx, int ny) {
	int ix,iy;
	double fftScale;

        complex_tt **carray;
	carray = (complex_tt **)array;

	fftScale = 1.0/(double)(nx*ny);
	for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
		carray[ix][iy][0] *= fftScale;
		carray[ix][iy][1] *= fftScale;
	}
}

void showPotential(fftw_complex ***pot,int nz,int nx,int ny,double dx,double dy,double dz) {
	char *fileName = "potential.dat";
	char systStr[256];
	FILE *fp;
	int ix,iz;
	static fftw_complex *data = NULL;
	int length;
	float r;


	/*
	make data array:
	*/
	length = (nx < ny) ? nx : ny ;
	length +=2;
	if (data == NULL)
		data = (fftw_complex *)malloc(length * sizeof(fftw_complex));

	/* 
	copy data to array
	*/

	for (ix=0;ix<nx;ix++) {
		data[ix][0] = pot[0][ix][ix][0];
		data[ix][1] = pot[0][ix][ix][1];
		for (iz=1;iz<nz;iz++) {
			data[ix][0] += 2*pot[iz][ix][ix][0];
			data[ix][1] += 2*pot[iz][ix][ix][1];
		}
		/*    printf("ix: %d, pot: %g\n",ix,data[ix]); */
	}

	if ((fp = fopen(fileName,"w")) == NULL) {
		printf("Could not open %s for writing!\n",fileName);
		return;
	}
	for (ix=0;ix<nx;ix++) {
		r = sqrt(ix*ix*(dx*dx+dy*dy));
		fprintf(fp,"%g",r);
		fprintf(fp,"\t%g\t%g",data[ix][0],data[ix][1]);
		/*    for (iz = 0;iz < ((nz>10) ? 10 : nz);iz++) {
		r = sqrt(ix*ix*(dx*dx+dy*dy)+iz*iz*dz*dz);      
		fprintf(fp,"\t%g",pot[iz][ix][ix][0]*r);      
		}
		*/
		fprintf(fp,"\n");
	}
	fclose(fp);

	sprintf(systStr,"xmgr -nxy %s &",fileName);
	system(systStr);
}


/*****************************************************************
* This function will write a data file with the pendeloesungPlot 
* for selected beams
****************************************************************/
void writeBeams(MULS *muls, WavePtr wave, int ilayer, int absolute_slice) {
	static char fileAmpl[32];
	static char filePhase[32];
	static char fileBeam[32];
	static FILE *fp1 = NULL,*fpAmpl = NULL,*fpPhase=NULL;
	int ib;
	static std::vector<int> hbeam,kbeam;
	static float_tt zsum = 0.0f,scale;
	float_tt rPart,iPart,ampl,phase;
	static char systStr[64];
	// static int counter=0;

	if (!muls->lbeams)
		return;

	if ((muls->mode != REFINE) && ((*muls).mode != CBED)) {
		if (ilayer < 0) {
			if (fp1 != NULL) fclose(fp1);
			if (fpAmpl != NULL) fclose(fpAmpl);
			if (fpPhase != NULL) fclose(fpPhase);
			fp1 = fpAmpl = fpPhase = NULL;
			sprintf(systStr,"xmgr -nxy %s &",fileAmpl);
			system(systStr);
			return;
		}

		if ((fp1 == NULL) || (fpAmpl == NULL) || (fpPhase == NULL)) {
			scale = 1.0F / ( ((float_tt)muls->nx) * ((float_tt)muls->ny) );
			hbeam = (*muls).hbeams;
			kbeam = (*muls).kbeams;
			if ((hbeam.empty()) || (kbeam.empty())) {
				printf("ERROR: hbeam or kbeam == NULL!\n");
				exit(0);
			}

			sprintf(fileAmpl,"%s/beams_amp.dat",(*muls).folder.c_str());
			sprintf(filePhase,"%s/beams_phase.dat",(*muls).folder.c_str());
			sprintf(fileBeam,"%s/beams_all.dat",(*muls).folder.c_str());
			fp1 = fopen(fileBeam, "w" );
			fpAmpl = fopen( fileAmpl, "w" );
			fpPhase = fopen( filePhase, "w" );
			if(fp1==NULL) {
				printf("can't open file %s\n", fileBeam);
				exit(0);
			}
			if(fpAmpl==NULL) {
				printf("can't open amplitude file %s\n",fileAmpl);
				exit(0);
			}
			if(fpPhase==NULL) {
				printf("can't open phase file %s\n", filePhase);
				exit(0);
			}
			fprintf(fp1, " (h,k) = ");
			for(ib=0; ib<(*muls).nbout; ib++) {
				fprintf(fp1," (%d,%d)", muls->hbeams[ib],  muls->kbeams[ib]);
			}
			fprintf( fp1, "\n" );
			fprintf( fp1, "nslice, (real,imag) (real,imag) ...\n\n");
			for( ib=0; ib<muls->nbout; ib++)
			{
				// printf("beam: %d [%d,%d]",ib,hbeam[ib],kbeam[ib]);			
				if(hbeam[ib] < 0 ) hbeam[ib] = muls->nx + hbeam[ib];
				if(kbeam[ib] < 0 ) kbeam[ib] = muls->ny + kbeam[ib];
				if(hbeam[ib] < 0 ) hbeam[ib] = 0;
				if(kbeam[ib] < 0 ) kbeam[ib] = 0;
				if(hbeam[ib] > muls->nx-1 ) hbeam[ib] = muls->nx-1;
				if(kbeam[ib] > muls->ny-1 ) kbeam[ib] = muls->ny-1;
				// printf(" => [%d,%d] %d %d\n",hbeam[ib],kbeam[ib],muls->nx,muls->ny);			
			}
			/****************************************************/
			/* setup of beam files, include the t=0 information */
			fprintf( fpAmpl, "%g",0.0);
			fprintf( fpPhase, "%g",0.0);
			for( ib=0; ib<muls->nbout; ib++) {
				ampl = 0.0;
				if ((hbeam[ib] == 0) && (kbeam[ib]==0))
					ampl = 1.0;
				fprintf(fpAmpl,"\t%g",ampl);
				fprintf(fpPhase,"\t%g",0.0);
			}
			fprintf( fpAmpl, "\n");
			fprintf( fpPhase, "\n");
		} /* end of if fp1 == NULL ... i.e. setup */


		zsum += (*muls).cz[ilayer];

		fprintf( fp1, "%g", zsum);
		fprintf( fpAmpl, "%g",zsum);
		fprintf( fpPhase, "%g",zsum);
		for( ib=0; ib<(*muls).nbout; ib++) {
			fprintf(fp1, "\t%g\t%g",
				rPart = scale*(*wave).wave[hbeam[ib]][kbeam[ib]][0],
				iPart = scale*(*wave).wave[hbeam[ib]][kbeam[ib]][1]);
			ampl = (float_tt)sqrt(rPart*rPart+iPart*iPart);
			phase = (float_tt)atan2(iPart,rPart);	
			fprintf(fpAmpl,"\t%g",ampl);
			fprintf(fpPhase,"\t%g",phase);
		}
		fprintf( fp1, "\n");
		fprintf( fpAmpl, "\n");
		fprintf( fpPhase, "\n");
	} /* end of if muls.mode != REFINE */
	
	if (muls->mode == TEM) {
		if (muls->pendelloesung == NULL) {
			muls->pendelloesung = 
				float2D((*muls).nbout,
				(*muls).slices*(*muls).mulsRepeat1*(*muls).mulsRepeat2*(*muls).cellDiv,
				"pendelloesung");
			scale = 1.0/(muls->nx*muls->ny); 
			printf("Allocated memory for pendelloesung plot (%d x %d)\n",
				(*muls).nbout,(*muls).slices*(*muls).mulsRepeat1*(*muls).mulsRepeat2);
		}
		for( ib=0; ib<muls->nbout; ib++) {
			rPart = (*wave).wave[muls->hbeams[ib]][muls->kbeams[ib]][0];
			iPart = (*wave).wave[muls->hbeams[ib]][muls->kbeams[ib]][1];
			muls->pendelloesung[ib][absolute_slice] = scale*(float_tt)(rPart*rPart+iPart*iPart);
			// printf("slice: %d beam: %d [%d,%d], intensity: %g\n",muls->nslic0,ib,muls->hbeam[ib],muls->kbeam[ib],muls->pendelloesung[ib][muls->nslic0]);			
		} // end of ib=0 ... 
	}
	
}



//////////////////////////////////////////////////////////////////////////////////
// DO NOT use this function, gives wrong results.
// 
// This method of sampling z at the slice interval gives very wrong HOLZ results,
// since this means a linear interpolation of a highly non-linear function.  Atoms 
// are quasi split up into two layers.  This causes a cutting in half of the z-periodicity
// for some unit cells (e.g. STO sampled 2 slices per unit cell)
//
/*
#define S_SCALE 0.5
#define PHI_SCALE 47.87658
#define SHOW_SINGLE_POTENTIAL 1
fftwf_complex *getAtomPotential3D_3DFFT(int Znum, MULS *muls,double B) {
	int ix,iy,iz,iiz,ind3d,ind3dd,iKind;
	double zScale,kzmax,kzborder;
	fftwf_plan plan;
	static double f,phase,s2,kmax2,kx,ky,kz,dkx,dky,dkz,dx2,dy2,dz2;
	static int nx,ny,nz;
	static fftwf_complex **atPot = NULL;
#if SHOW_SINGLE_POTENTIAL == 1
	//static imageStruct *header = NULL;
	char fileName[256];
	ImageIOPtr imageIO = ImageIOPtr(new CImageIO(muls->potNx,muls->potNy,
				muls->sliceThickness,muls->resolutionX/OVERSAMP_X,
				muls->resolutionY/OVERSAMP_X));;
#endif 
	static double *splinb=NULL;
	static double *splinc=NULL;
	static double *splind=NULL;


	// scattering factors in:
	// float scatPar[4][30]
	if (atPot == NULL) {
		splinb = double1D(30, "splinb" );
		splinc = double1D(30, "splinc" );
		splind = double1D(30, "splind" );


		// Why do I use nx+2 and make dkx=1/nx? ... don't know anymore. :(    
		nx = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionX)+2;
		ny = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionY)+2;
		// make the atom sphere have an even number of slices in z-direction
		// the 2 center slices will be equal.
		// The FFT-resolution in the z-direction must be high enough to avoid 
		// artifacts due to premature cutoff of the rec. space scattering factor 
		nz = 2*OVERSAMP_Z*(int)ceil(muls->atomRadius/muls->sliceThickness);  
		dkx = OVERSAMP_X/((nx-2)*muls->resolutionX);  // 0.5-factor s->k
		dky = OVERSAMP_X/((ny-2)*muls->resolutionY);
		dkz = OVERSAMP_Z/(double)(nz*muls->sliceThickness);
		// if (nz<=2) dkz=0;
		// else dkz = 1/((nz-2)*muls->sliceThickness);
		kmax2 = 0.5*nx*dkx/(double)OVERSAMP_X;  // largest k that we'll admit

		printf("Cutoff scattering angle:kmax=%g, smax=%g (1/A)\n",kmax2,S_SCALE*kmax2);
		scatPar[0][29] = 1.2*S_SCALE*kmax2;
		scatPar[0][28] = 1.1*S_SCALE*kmax2;
		scatPar[0][27] = S_SCALE*kmax2;
		if (scatPar[0][26] > scatPar[0][27]) {
			// set additional scattering parameters to zero:
			for (ix = 0;ix < 20;ix++) {
				if (scatPar[0][26-ix] < scatPar[0][27]-0.1*(ix+1)) break;
				scatPar[0][26-ix] = scatPar[0][27]-0.1*(ix+1);
				for (iy=1; iy<=8;iy++) scatPar[iy][26-ix] = 0;	
			}
			if (muls->printLevel > 1)
				printf("getAtomPotential3D: reduced angular range of scattering factor to %g/A!\n",scatPar[0][26-ix]);
		} 
		kmax2 *= kmax2;

		atPot = (complex_tt **)malloc((NZMAX+1)*sizeof(complex_tt *));
		for (ix=0;ix<=NZMAX;ix++) atPot[ix] = NULL;
	}
	// initialize this atom, if it has not been done yet:
	if (atPot[Znum] == NULL) {
		switch (Znum) {
	case 38: iKind = 1; break;  // Sr
	case 22: iKind = 2; break;  // Ti
	case  8: iKind = 3; break;  // O
	case 49: iKind = 4; break;  // In
	case 15: iKind = 5; break;  // P
	case  2: iKind = 6; break;  // He
	case 17: iKind = 7; break;  // Cl
	case 14: iKind = 8; break;  // Si
	case 20: iKind = 9; break;  // Ca
	case 56: iKind = 10; break;  // Ba
	case 26: iKind = 11; break;  // Fe
	default: 
		printf("This atom kind (%d) is not supported yet - sorry!\n",Znum);
		exit(0);
		}


		// setup cubic spline interpolation:
		splinh(scatPar[0],scatPar[iKind],splinb,splinc,splind,30);
                // TODO: handle multiple precision possibility (need fftwmalloc?)
		atPot[Znum] = (complex_tt*)fftwf_malloc(nx*ny*nz*sizeof(complex_tt));
		memset(atPot[Znum],0,nx*ny*nz*sizeof(complex_tt));
		kzmax    = dkz*nz/2.0; 
		kzborder = dkz*(nz/(2*OVERSAMP_Z) -1); 
		for (iz=0;iz<nz;iz++) {
			kz = dkz*(iz<nz/2 ? iz : iz-nz);    
			// We also need to taper off the potential in z-direction
			// in order to avoid cutoff artifacts.
			zScale = fabs(kz) <= kzborder ? 1.0 : 
				0.5+0.5*cos(PI*(fabs(kz)-kzborder)/(kzmax-kzborder));
			// printf("iz=%d, kz=%g, zScale=%g ",iz,kz,zScale);
			for (ix=0;ix<nx;ix++) {
				kx = dkx*(ix<nx/2 ? ix : ix-nx);      
				for (iy=0;iy<ny;iy++) {
					ky = dky*(iy<ny/2 ? iy : iy-ny);      
					s2 = S_SCALE*S_SCALE*(kx*kx+ky*ky+kz*kz);
					// if this is within the allowed circle:
					if (s2<S_SCALE*S_SCALE*kmax2) {
						ind3d = iy+ix*ny+iz*nx*ny;
						// f = fe3D(Znum,k2,muls->tds,1.0,muls->scatFactor);
						// multiply scattering factor with Debye-Waller factor:
						// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
						f = seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,30,sqrt(s2))*exp(-s2*B);
						// note that the factor 2 is missing in the phase (2pi k*r)
						// this places the atoms in the center of the box.
						phase = PI*(kx*muls->resolutionX*nx/(OVERSAMP_X)+ky*muls->resolutionY*ny/(OVERSAMP_X));
						phase += PI*kz/OVERSAMP_Z*(muls->sliceThickness*(nz+1));
						atPot[Znum][ind3d][0] = zScale*f*cos(phase);
						atPot[Znum][ind3d][1] = zScale*f*sin(phase);
						// if ((kx==0) && (ky==0)) printf(" f=%g (%g, [%g, %g])\n",f,f*zScale,atPot[Znum][ind3d][0],atPot[Znum][ind3d][1]);
					}
				}
			}
		}

		plan = fftwf_plan_dft_3d(nz,nx,ny,atPot[Znum],atPot[Znum],FFTW_BACKWARD,FFTW_ESTIMATE);
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
		dx2 = muls->resolutionX*muls->resolutionX/(OVERSAMP_X*OVERSAMP_X);
		dy2 = muls->resolutionY*muls->resolutionY/(OVERSAMP_X*OVERSAMP_X);
		dz2 = muls->sliceThickness*muls->sliceThickness/(OVERSAMP_Z*OVERSAMP_Z);
		// Here we make sure that our atom is not bigger than the desired radius, i.e. that 
		// we really have round blobs.
		// We also make sure that the potential touches zero at least somewhere.  This will avoid 
		// sharp edges that could produce ringing artifacts.  
		// It is certainly debatable whether this is a good apprach, or not. 
		// printf("Setting up %d x %d x %d potential for Znum=%d, (atom kind %d), Omega=%g\n",nx,ny,nz,Znum,iKind,dkx*dky*dkz);
		// min = atPot[Znum][ny/2+nx/2*ny+nz/2*nx*ny][0];
		for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) for (iz=0;iz<nz/OVERSAMP_Z;iz++) {
			ind3d = iy+ix*ny+iz*nx*ny;
			// Integrate over 3 neighboring layers here:::::::::
			for (zScale=0,iiz=0;iiz<OVERSAMP_Z;iiz++) {
				ind3dd = iy+ix*ny+(OVERSAMP_Z*iz+iiz)*nx*ny;
				// if (sqrt((ix-nx/2)*(ix-nx/2)*dx2+(iy-ny/2)*(iy-ny/2)*dy2+(iz-0.5-nz/2)*(iz-0.5-nz/2)*dz2) > muls->atomRadius) {
				if ((sqrt((ix-nx/2)*(ix-nx/2)*dx2+(iy-ny/2)*(iy-ny/2)*dy2+
					(OVERSAMP_Z*iz+iiz+0.5-nz/2)*(OVERSAMP_Z*iz+iiz+0.5-nz/2)*dz2)) > muls->atomRadius) {
						// atPot[Znum][ind3d][0] *= dkx*dky*dkz;
						atPot[Znum][ind3dd][0] = 0;
				}
				zScale += atPot[Znum][ind3dd][0];
			}
			// assign the iz-th slice the sum of the 3 other slices:
			// and divide by unit cell area (volume, if in 3D):
			atPot[Znum][ind3d][0] = zScale*dkx*dky/nz;	
			if (atPot[Znum][ind3d][0] < 0) atPot[Znum][ind3d][0] = 0;
			// if (atPot[Znum][ind3d][0] < min) min = atPot[Znum][ind3d][0];

			atPot[Znum][ind3d][1]= 0;
		}
		// printf("Found minimum potential value of %g ... subtracting it from 3D potential.\n",min);
		for (zScale = 0,ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++)  for (iz=0;iz<nz/OVERSAMP_Z;iz++){
			if ((sqrt((ix-nx/2)*(ix-nx/2)*dx2+(iy-ny/2)*(iy-ny/2)*dy2+
				(iz+0.5-nz/(2*OVERSAMP_Z))*(iz+0.5-nz/(2*OVERSAMP_Z))*dz2)) <= muls->atomRadius) {
					//      if (sqrt((ix-nx/2)*(ix-nx/2)*dx2+(iy-ny/2)*(iy-ny/2)*dy2+(iz-nz/2)*(iz-nz/2)*dz2) <= muls->atomRadius) {	
					ind3d = iy+ix*ny+iz*nx*ny;
					if ((ix == nx/2) && (iy == ny/2)) zScale += atPot[Znum][ind3d][0];
					// atPot[Znum][ind3d][0] -= min;
			}
		}
		// make sure we don't produce negative potential:
		// if (min < 0) for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) atPot[Znum][iy+ix*ny][0] -= min;
#if SHOW_SINGLE_POTENTIAL == 1
		for (iz=0;iz<nz/OVERSAMP_Z;iz++) {
			sprintf(fileName,"potential_%d_%d.img",Znum,iz);
			imageIO->SetThickness(iz);
			ptr = &(atPot[Znum][iz*nx*ny]);
			imageIO->WriteRealImage((void **)ptr,fileName);
		}
#endif    
		if (muls->printLevel > 0)
			printf("Created 3D %d x %d x %d potential array for Z=%d (%d, B=%g, sum=%g)\n",nx,ny,nz/OVERSAMP_Z,Znum,iKind,B,zScale);
	}
	return atPot[Znum];
}
#undef SHOW_SINGLE_POTENTIAL

*/

