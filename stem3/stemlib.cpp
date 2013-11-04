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
#define OVERSAMPLING 3
#define OVERSAMPLINGZ (3*OVERSAMPLING)
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


/**********************************************
* This function creates a incident STEM probe 
* at position (dx,dy)
* with parameters given in muls
*
* The following Abberation functions are being used:
* 1) ddf = Cc*dE/E + Cc2*(dE/E)^2,    
*    Cc, Cc2 = chrom. Abber. (1st, 2nd order) [1]
* 2) chi(qx,qy) = (2*pi/lambda)*{0.5*C1*(qx^2+qy^2)+
*                 0.5*C12a*(qx^2-qy^2)+
*                 C12b*qx*qy+
*                 C21a/3*qx*(qx^2+qy^2)+
*                 ... 
*                 +0.5*C3*(qx^2+qy^2)^2
*                 +0.125*C5*(qx^2+qy^2)^3
*                 ... (need to finish)
*
*
*    qx = acos(kx/K), qy = acos(ky/K) 
*
* References:
* [1] J. Zach, M. Haider, 
*    "Correction of spherical and Chromatic Abberation 
*     in a low Voltage SEM", Optik 98 (3), 112-118 (1995)
* [2] O.L. Krivanek, N. Delby, A.R. Lupini,
*    "Towards sub-Angstroem Electron Beams", 
*    Ultramicroscopy 78, 1-11 (1999)
*
*********************************************/

#define SMOOTH_EDGE 5 // make a smooth edge on AIS aperture over +/-SMOOTH_EDGE pixels
void probe(MULS *muls, WavePtr wave, double dx, double dy)
{
	// static char *plotFile = "probePlot.dat",systStr[32];
	int ix, iy, nx, ny, ixmid, iymid;
	int CsDefAstOnly = 0;
	float rmin, rmax, aimin, aimax;
	// float **pixr, **pixi;
	double  kx, ky, ky2,k2, ktheta2, ktheta, k2max, v0, wavlen,ax,by,x,y,
		rx2, ry2,rx,ry, pi, scale, pixel,alpha,
		df, df_eff, chi1, chi2,chi3, sum, chi, time,r,phi;
	double gaussScale = 0.05;
	double envelope,delta,avgRes,edge;

	// FILE *fp=NULL;

	/* temporary fix, necessary, because fftw has rec. space zero 
	in center of image:
	*/
	nx = (*muls).nx;
	ny = (*muls).ny;
	ax = nx*(*muls).resolutionX; 
	by = ny*(*muls).resolutionY; 
	dx = ax-dx;
	dy = by-dy;
	gaussScale = (*muls).gaussScale;
	// average resolution:
	avgRes = sqrt(0.5*(muls->resolutionX*muls->resolutionX+muls->resolutionY*muls->resolutionY));
	edge = SMOOTH_EDGE*avgRes;

	/********************************************************
	* formulas from:
	* http://cimesg1.epfl.ch/CIOL/asu94/ICT_8.html
	*
	* dE_E = dE/E = energy spread of emitted electrons
	* dV_V = dV/V = acc. voltage fluctuations
	* dI_I = dI/I = lens current fluctuations
	* delta defocus in Angstroem (Cc in A)
	*******************************************************/
	delta = muls->Cc*muls->dE_E;
	if (muls->printLevel > 2) printf("defocus offset: %g nm (Cc = %g)\n",delta,muls->Cc);

	if (wave->wave == NULL) {
		printf("Error in probe(): Wave not allocated!\n");
		exit(0);
	}

	/**********************************************************
	*  Calculate misc constants  
	*********************************************************/  
	time = cputim( );
	pi = 4.0 * atan( 1.0 );

	rx = 1.0/ax;
	rx2 = rx * rx;
	ry = 1.0/by;
	ry2 = ry * ry;

	ixmid = nx/2;
	iymid = ny/2;

	// df = muls->df0;
	v0 = muls->v0;
	wavlen = 12.26/ sqrt( v0*1.e3 + v0*v0*0.9788 );

	/*  printf("Wavelength: %g A\n",wavlen);
	*/


	// chi2 = (*muls).Cs*0.5*wavlen*wavlen;
	// chi3 = (*muls).C5*0.25*wavlen*wavlen*wavlen*wavlen;
	/* delta *= 0.5*delta*pi*pi*wavlen*wavlen; */

	/* convert convergence angle from mrad to rad */
	alpha = 0.001*muls->alpha;
	k2max = sin(alpha)/wavlen;  /* = K0*sin(alpha) */
	k2max = k2max * k2max;

	/*   Calculate MTF 
	NOTE zero freg is in the bottom left corner and
	expandes into all other corners - not in the center
	this is required for FFT

	PIXEL = diagonal width of pixel squared
	if a pixel is on the apertur boundary give it a weight
	of 1/2 otherwise 1 or 0
	*/
	pixel = ( rx2 + ry2 );
	scale = 1.0/sqrt((double)nx*(double)ny);

	/*
	if ((muls.a33 == 0) && (muls.a31 == 0) && (muls.a44 == 0) && (muls.a42 == 0) &&
	(muls.a55 == 0) && (muls.a53 == 0) && (muls.a51 == 0) && 
	(muls.a66 == 0) && (muls.a64 == 0) && (muls.a62 == 0) && (muls.C5 == 0)) {
	CsDefAstOnly = 1;
	}
	*/

	for( iy=0; iy<ny; iy++) {
		ky = (double) iy;
		if( iy > iymid ) ky = (double) (iy-ny);
		ky2 = ky*ky*ry2;
		for( ix=0; ix<nx; ix++) {
			kx = (double) ix;
			if( ix > ixmid ) kx = (double) (ix-nx);
			k2 = kx*kx*rx2 + ky2;
			ktheta2 = k2*(wavlen*wavlen);
			ktheta = sqrt(ktheta2);
			phi = atan2(ry*ky,rx*kx);
			// compute the effective defocus from the actual defocus and the astigmatism: 
			// df_eff = df + muls->astigMag*cos(muls->astigAngle+phi);

			// chi = chi1*k2*(df_eff +chi2*k2)-2.0*pi*( (dx*kx/ax) + (dy*ky/by) );
			// defocus, astigmatism, and shift:
			chi = ktheta2*(muls->df0+delta + muls->astigMag*cos(2.0*(phi-muls->astigAngle)))/2.0;
			ktheta2 *= ktheta;  // ktheta^3 
			if ((muls->a33 > 0) || (muls->a31 > 0)) {
				chi += ktheta2*(muls->a33*cos(3.0*(phi-muls->phi33))+muls->a31*cos(phi-muls->phi31))/3.0;
			}	
			ktheta2 *= ktheta;   // ktheta^4
			if ((muls->a44 > 0) || (muls->a42 > 0) || (muls->Cs != 0)) {
				// chi += ktheta2*(muls->a33*cos(3*(phi-muls->phi33))+muls->a31*cos(phi-muls->phi31))/3.0;
				chi += ktheta2*(muls->a44*cos(4.0*(phi-muls->phi44))+muls->a42*cos(2.0*(phi-muls->phi42))+muls->Cs)/4.0;  
				//                     1/4*(a(4,4).*cos(4*(kphi-phi(4,4)))+a(4,2).*cos(2*(kphi-phi(4,2)))+c(4)).*ktheta.^4+...
			}
			ktheta2 *= ktheta;    // ktheta^5
			if ((muls->a55 > 0) || (muls->a53 > 0) || (muls->a51 > 0)) {
				chi += ktheta2*(muls->a55*cos(5.0*(phi-muls->phi55))+muls->a53*cos(3.0*(phi-muls->phi53))+muls->a51*cos(phi-muls->phi51))/5.0;
				//                     1/5*(a(5,5).*cos(5*(kphi-phi(5,5)))+a(5,3).*cos(3*(kphi-phi(5,3)))+a(5,1).*cos(1*(kphi-phi(5,1)))).*ktheta.^5+...
			}
			ktheta2 *= ktheta;    // ktheta^6
			if ((muls->a66 > 0) || (muls->a64 > 0) || (muls->a62 = 0) || (muls->C5 != 0)) {
				chi += ktheta2*(muls->a66*cos(6.0*(phi-muls->phi66))+muls->a64*cos(4.0*(phi-muls->phi64))+muls->a62*cos(2.0*(phi-muls->phi62))+muls->C5)/6.0;
				//                     1/6*(a(6,6).*cos(6*(kphi-phi(6,6)))+a(6,4).*cos(4*(kphi-phi(6,4)))+a(6,2).*cos(2*(kphi-phi(6,2)))+c(6)).*ktheta.^6);
			}

			chi *= 2*pi/wavlen;
			chi -= 2.0*pi*( (dx*kx/ax) + (dy*ky/by) );
			// include higher order aberrations


			if ( ( (*muls).ismoth != 0) && 
				( fabs(k2-k2max) <= pixel)) {
					wave->wave[ix][iy][0]= (float) ( 0.5*scale * cos(chi));
					wave->wave[ix][iy][1]= (float) (-0.5*scale* sin(chi));
			} 
			else if ( k2 <= k2max ) {
				wave->wave[ix][iy][0]= (float)  scale * cos(chi);
				wave->wave[ix][iy][1]= (float) -scale * sin(chi);
			} 
			else {
				wave->wave[ix][iy][0] = wave->wave[ix][iy][1] = 0.0f;
			}
		}
	}
	/* Fourier transform into real space */
	// fftwnd_one(muls->fftPlanInv, &(muls->wave[0][0]), NULL);
#if FLOAT_PRECISION == 1
	fftwf_execute(wave->fftPlanWaveInv);
#else
	fftw_execute(wave->fftPlanWaveInv);
#endif
	/**********************************************************
	* display cross section of probe intensity
	*/

	/* multiply with gaussian in Real Space in order to avoid artifacts */
	if (muls->gaussFlag) {
		for( ix=0; ix<nx; ix++) {
			for( iy=0; iy<ny; iy++) {
				r = exp(-((ix-nx/2)*(ix-nx/2)+(iy-ny/2)*(iy-ny/2))/(nx*nx*gaussScale));
				wave->wave[ix][iy][0] *= (float)r;
				wave->wave[ix][iy][1] *= (float)r;
			}
		}  
	}

	/* Apply AIS aperture in Real Space */
	// printf("center: %g,%g\n",dx,dy);
	if (muls->aAIS > 0) {
		for( ix=0; ix<nx; ix++) {
			for( iy=0; iy<ny; iy++) {
				x = ix*muls->resolutionX-dx;
				y = iy*muls->resolutionY-dy;
				r = sqrt(x*x+y*y);
				delta = r-0.5*muls->aAIS+edge;
				if (delta > 0) {
					wave->wave[ix][iy][0] = 0;
					wave->wave[ix][iy][1] = 0;
				}
				else if (delta >= -edge) {
					scale = 0.5*(1-cos(pi*delta/edge));
					wave->wave[ix][iy][0] = scale*wave->wave[ix][iy][0];
					wave->wave[ix][iy][1] = scale*wave->wave[ix][iy][1];
				}
			}
		}
	}

	/*  Normalize probe intensity to unity  */

	sum = 0.0;
	for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) 
		sum +=  wave->wave[ix][iy][0]*wave->wave[ix][iy][0]
	+ wave->wave[ix][iy][1]*wave->wave[ix][iy][1];

	scale = 1.0 / sum;
	scale = scale * ((double)nx) * ((double)ny);
	scale = (double) sqrt( scale );

	for( ix=0; ix<nx; ix++) 
		for( iy=0; iy<ny; iy++) {
			wave->wave[ix][iy][0] *= (float) scale;
			wave->wave[ix][iy][1] *= (float) scale;
		}

		/*  Output results and find min and max to echo
		remember that complex pix are stored in the file in FORTRAN
		order for compatability
		*/

		rmin = wave->wave[0][0][0];
		rmax = rmin;
		aimin = wave->wave[0][0][1];
		aimax = aimin;
		for( iy=0; iy<ny; iy++) {
			for( ix=0; ix<nx; ix++) {
				if( wave->wave[ix][iy][0] < rmin ) rmin = wave->wave[ix][iy][0];
				if( wave->wave[ix][iy][0] > rmax ) rmax = wave->wave[ix][iy][0];
				if( wave->wave[ix][iy][1] < aimin ) aimin = wave->wave[ix][iy][1];
				if( wave->wave[ix][iy][1] > aimax ) aimax = wave->wave[ix][iy][1];
			}
		}
		(*muls).rmin = rmin;
		(*muls).rmax = rmax;
		(*muls).aimin = aimin;
		(*muls).aimax = aimax;

		/**********************************************************/

}  /* end probe() */


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
int runMulsSTEM(MULS *muls, WavePtr wave) {
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
	fftScale = 1.0/(muls->nx*muls->ny);

	wavlen = (float_tt)wavelength((*muls).v0);

	/*  calculate the total specimen thickness and echo */
	cztot=0.0;
	for( islice=0; islice<(*muls).slices; islice++) {
		cztot += (*muls).cz[islice];
	}
	if (printFlag)
		printf("Specimen thickness: %g Angstroms\n", cztot);

	scale = 1.0F / (((float_tt)muls->nx) * ((float_tt)muls->ny));

	for (mRepeat = 0; mRepeat < muls->mulsRepeat1; mRepeat++) 
	{
		for( islice=0; islice < muls->slices; islice++ ) 
		{
			absolute_slice = (muls->totalSliceCount+islice);

			// if ((muls->cubez > 0) && (muls->thickness >= muls->cubez)) break;
			//  else if ((muls->cubez == 0) && (muls->thickness >= muls->c)) break;

			/***********************************************************************
			* Transmit is a simple multiplication of wave with trans in real space
			**********************************************************************/
			transmit((void **)wave->wave, (void **)(muls->trans[islice]), muls->nx,muls->ny, wave->iPosX, wave->iPosY);
			//    writeImage_old(wave,(*muls).nx,(*muls).ny,(*muls).thickness,"wavet.img");      
			/***************************************************** 
			* remember: prop must be here to anti-alias
			* propagate is a simple multiplication of wave with prop
			* but it also takes care of the bandwidth limiting
			*******************************************************/
#if FLOAT_PRECISION == 1
			fftwf_execute(wave->fftPlanWaveForw);
#else
			fftw_execute(wave->fftPlanWaveForw);
#endif
			propagate_slow((void **)wave->wave, muls->nx, muls->ny, muls);

			collectIntensity(muls, wave, muls->totalSliceCount+islice*(1+mRepeat));

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

			/*
			sprintf(outStr,"wave%d.img",islice);
			writeImage_old(wave,(*muls).nx,(*muls).ny,(*muls).thickness,"wavep.img");
			*/    

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
				collectIntensity(muls,wave,absolute_slice*(1+mRepeat));
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
			wave->WriteWave(wave->fileout);
			if (printFlag)
				printf("Created complex image file %s\n",(*wave).fileout);    
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
	if (muls->tds) sprintf(fileName,"%s/wave_%d_%d.img",muls->folder,muls->avgCount,t);
	else sprintf(fileName,"%s/wave_%d.img",muls->folder,t);
	wave->WriteWave(fileName, "Wave Function", params);
}

/********************************************************************
* collectIntensity(muls, wave, slice)
* collect the STEM signal on the annular detector(s) defined in muls
* and write the appropriate pixel in the image for each detector and thickness
* The number of images is determined by the following formula:
* muls->slices*muls->cellDiv/muls->outputInterval 
* There are muls->detectorNum different detectors
*******************************************************************/
void collectIntensity(MULS *muls, WavePtr wave, int slice) 
{
	int i,ix,iy,ixs,iys,t;
	float_tt k2;
	double intensity,scale,scaleCBED,scaleDiff,intensity_save;
	char fileName[256],avgName[256]; 
	float_tt **diffpatAvg = NULL;
	int tCount = 0;

	std::vector<std::vector<DetectorPtr> > detectors;

	scale = muls->electronScale/((double)(muls->nx*muls->ny)*(muls->nx*muls->ny));
	// scaleCBED = 1.0/(scale*sqrt((double)(muls->nx*muls->ny)));
	scaleDiff = 1.0/sqrt((double)(muls->nx*muls->ny));

	tCount = (int)(ceil((double)((muls->slices * muls->cellDiv) / muls->outputInterval)));

	// we write directly to the shared muls object.  This is safe only because 
	//    each thread is accessing different pixels in the output images.
	detectors = muls->detectors;

	if (muls->outputInterval == 0) t = 0;
	else if (slice < ((muls->slices*muls->cellDiv)-1))
	{
		t = (int)((slice) / muls->outputInterval);
		//if (t > tCount)
		//{
			// printf("t = %d, which is greater than tCount (%d)\n",t,tCount);
			//t = tCount-1;
		//}
	}
	else
	{
		t = tCount;
	}

	int position_offset = wave->detPosY * muls->scanXN + wave->detPosX;

	// Multiply each image by its number of averages and divide by it later again:
	for (i=0;i<muls->detectorNum;i++) 
	{
		detectors[t][i]->image[wave->detPosX][wave->detPosY]  *= detectors[t][i]->Navg;	
		detectors[t][i]->image2[wave->detPosX][wave->detPosY] *= detectors[t][i]->Navg;	
		detectors[t][i]->error = 0;
	}
	/* add the intensities in the already 
	fourier transformed wave function */
	for (ix = 0; ix < muls->nx; ix++) 
	{
		for (iy = 0; iy < muls->ny; iy++) 
		{
			k2 = muls->kx2[ix]+muls->ky2[iy];
			intensity = (wave->wave[ix][iy][0]*wave->wave[ix][iy][0]+
				wave->wave[ix][iy][1]*wave->wave[ix][iy][1]);
			wave->diffpat[(ix+muls->nx/2)%muls->nx][(iy+muls->ny/2)%muls->ny] = intensity*scaleDiff;
			intensity *= scale;
			for (i=0;i<muls->detectorNum;i++) {
				if ((k2 >= detectors[t][i]->k2Inside) && (k2 <= detectors[t][i]->k2Outside)) 
				{
					// detector in center of diffraction pattern:
					if ((detectors[t][i]->shiftX == 0) && (detectors[t][i]->shiftY == 0)) 
					{
						detectors[t][i]->image[wave->detPosX][wave->detPosY] += intensity;
						// misuse the error number for collecting this pixels raw intensity
						detectors[t][i]->error += intensity;
					}
					/* special case for shifted detectors: */		
					else 
					{
						intensity_save = intensity;
						ixs = (ix+(int)detectors[t][i]->shiftX+muls->nx) % muls->nx;
						iys = (iy+(int)detectors[t][i]->shiftY+muls->ny) % muls->ny;	    
						intensity = scale * (wave->wave[ixs][iys][0]*wave->wave[ixs][iys][0]+
							wave->wave[ixs][iys][1]*wave->wave[ixs][iys][1]);
						detectors[t][i]->image[wave->detPosX][wave->detPosY] += intensity;
						// repurpose the error number for collecting this pixels raw intensity
						detectors[t][i]->error += intensity;
						/* restore intensity, so that it will not be shifted for the other detectors */
						intensity = intensity_save;
					}
				} /* end of if k2 ... */
			} /* end of for i=0 ... detectorNum */
		} /* end of for iy=0... */
	} /* end of for ix = ... */

	////////////////////////////////////////////////////////////////////////////
	// write the diffraction pattern to disc in case we are working in CBED mode
	if ((muls->mode == CBED) && (muls->saveLevel > 0)) {
		sprintf(avgName,"%s/diff_%d.img",muls->folder,t);
		// for (ix=0;ix<muls->nx*muls->ny;ix++) wave->diffpat[0][ix] *= scaleCBED;
		if (muls->avgCount == 0) {
                  wave->WriteDiffPat(avgName);
		}
		else {
			wave->ReadAvgArray(avgName);
			for (ix=0;ix<muls->nx*muls->ny;ix++) {
				wave->avgArray[0][ix] = (muls->avgCount*wave->avgArray[0][ix]+wave->diffpat[0][ix])/(muls->avgCount+1);
			}
			wave->WriteAvgArray(avgName);
		}
	}

	// Divide each image by its number of averages again:
	for (i=0;i<muls->detectorNum;i++) {
		// add intensity squared to image2 for this detector and pixel, then rescale:
		detectors[t][i]->image2[wave->detPosX][wave->detPosY] += detectors[t][i]->error*detectors[t][i]->error;
		detectors[t][i]->image2[wave->detPosX][wave->detPosY] /= detectors[t][i]->Navg+1;	

		// do the rescaling for the average image:
		detectors[t][i]->image[wave->detPosX][wave->detPosY] /= detectors[t][i]->Navg+1;	
	}
}

/*****  saveSTEMImages *******/
// Saves all detector images (STEM images) that are defined in muls.
//   When saving intermediate STEM images is enabled, this also saves
//   the intermediate STEM images for each detector.
void saveSTEMImages(MULS *muls)
{
	int i, ix, islice;
	double intensity;
	static char fileName[256]; 
	//imageStruct *header = NULL;
	std::vector<DetectorPtr> detectors;
	float t;
    std::map<std::string, double> params;
    std::vector<unsigned> position(1,0);

	int tCount = (int)(ceil((double)((muls->slices * muls->cellDiv) / muls->outputInterval)));

	// Loop over slices (intermediates)
	for (islice=0; islice <= tCount; islice++)
	{
		if (islice<tCount)
		{
			t = ((islice+1) * muls->outputInterval ) * muls->sliceThickness;
		}
		else
		{
			t = muls->slices*muls->cellDiv*muls->sliceThickness;
		}
		detectors = muls->detectors[islice];
		// write the output STEM images:
		// This is done only after all pixels have completed, so that image is complete.
		for (i=0; i<muls->detectorNum; i++) 
		{
			// calculate the standard error for this image:
			detectors[i]->error = 0;
			intensity             = 0;
			for (ix=0; ix<muls->scanXN * muls->scanYN; ix++) 
			{
				detectors[i]->error += (detectors[i]->image2[0][ix]-
					detectors[i]->image[0][ix] * detectors[i]->image[0][ix]);
				intensity += detectors[i]->image[0][ix] * detectors[i]->image[0][ix];
			}
			detectors[i]->error /= intensity;
            position[0]=islice;
            sprintf(fileName,"%s/%s", muls->folder, detectors[i]->name);
            params["Thickness"]=t;
            params["Runs Averaged"]=(double)muls->avgCount+1;
            params["Error"]=(double)detectors[i]->error;

            // TODO: why is this here?  It's a second image.  Why aren't we just saving another image?
            /*
			for (ix=0; ix<muls->scanXN * muls->scanYN; ix++) 
			{
				detectors[i]->SetParameter(2+ix, (double)detectors[i]->image2[0][ix]);
			}
            */

			// exclude the suffix if this is the last detector (at the final thickness)
			if (islice==tCount)
				detectors[i]->WriteImage(fileName, detectors[i]->name, params);
			else
				detectors[i]->WriteImage(fileName, detectors[i]->name, params, position);
		}
	}
}


/******************************************************************
* propagate_slow() 
* replicates the original way, mulslice did it:
*****************************************************************/
void propagate_slow(void **w,int nx, int ny,MULS *muls)
{
	int ixa, iya;
	float_tt wr, wi, tr, ti,ax,by;
	float_tt scale,t,dz; 
	static float_tt dzs=0;
	static float_tt *propxr=NULL,*propyr=NULL;
	static float_tt *propxi=NULL,*propyi=NULL;
	static float_tt *kx2,*ky2;
	static float_tt *kx,*ky;
	static float_tt k2max=0,wavlen;
	complex_tt **wave;
	wave = (complex_tt **)w;

	ax = (*muls).resolutionX*nx;
	by = (*muls).resolutionY*ny;
	dz = (*muls).cz[0];

	if (dz != dzs) {
		if (propxr == NULL) {
			propxr = float1D(nx, "propxr" );
			propxi = float1D(nx, "propxi" );
			propyr = float1D(ny, "propyr" );
			propyi = float1D(ny, "propyi" );
			kx2    = float1D(nx, "kx2" );
			kx     = float1D(nx, "kx" );
			ky2    = float1D(ny, "ky2" );
			ky     = float1D(ny, "ky" );
		}
		dzs = dz;
		scale = dz*PI;
		wavlen = wavelength((*muls).v0);

		for( ixa=0; ixa<nx; ixa++) {
			kx[ixa] = (ixa>nx/2) ? (float_tt)(ixa-nx)/ax : 
				(float_tt)ixa/ax;
			kx2[ixa] = kx[ixa]*kx[ixa];
			t = scale * (kx2[ixa]*wavlen);
			propxr[ixa] = (float_tt)  cos(t);
			propxi[ixa] = (float_tt) -sin(t);
		}
		for( iya=0; iya<ny; iya++) {
			ky[iya] = (iya>ny/2) ? 
				(float_tt)(iya-ny)/by : 
			(float_tt)iya/by;
			ky2[iya] = ky[iya]*ky[iya];
			t = scale * (ky2[iya]*wavlen);
			propyr[iya] = (float_tt)  cos(t);
			propyi[iya] = (float_tt) -sin(t);
		}
		k2max = nx/(2.0F*ax);
		if (ny/(2.0F*by) < k2max ) k2max = ny/(2.0F*by);
		k2max = 2.0/3.0 * k2max;
		// TODO: modifying shared value from multiple threads?
		k2max = k2max*k2max;
		(*muls).kx2=kx2;
		(*muls).ky2=ky2;
		(*muls).kx=kx;
		(*muls).ky=ky;
	} 
	/* end of: if dz != dzs */
	/*************************************************************/

	/*************************************************************
	* Propagation
	************************************************************/
	for( ixa=0; ixa<nx; ixa++) {
		if( kx2[ixa] < k2max ) {
			for( iya=0; iya<ny; iya++) {
				if( (kx2[ixa] + ky2[iya]) < k2max ) {

					wr = wave[ixa][iya][0];
					wi = wave[ixa][iya][1];
					tr = wr*propyr[iya] - wi*propyi[iya];
					ti = wr*propyi[iya] + wi*propyr[iya];
					wave[ixa][iya][0] = tr*propxr[ixa] - ti*propxi[ixa];
					wave[ixa][iya][1] = tr*propxi[ixa] + ti*propxr[ixa];

				} else
					wave[ixa][iya][0] = wave[ixa][iya][1] = 0.0F;
			} /* end for(iy..) */

		} else for( iya=0; iya<ny; iya++)
			wave[ixa][iya][0] = wave[ixa][iya][1] = 0.0F;
	} /* end for(ix..) */
} /* end propagate_slow() */


/*------------------------ transmit() ------------------------*/
/*
transmit the wavefunction thru one layer 
(simply multiply wave by transmission function)

waver,i[ix][iy]  = real and imaginary parts of wavefunction
transr,i[ix][iy] = real and imag parts of transmission functions

nx, ny = size of array

on entrance waver,i and transr,i are in real space

only waver,i will be changed by this routine
*/
void transmit(void **wave, void **trans,int nx, int ny,int posx,int posy) {
	int ix, iy;
	double wr, wi, tr, ti;

	complex_tt **w,**t;
	w = (complex_tt **)wave;
	t = (complex_tt **)trans;

	/*  trans += posx; */
	for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
		wr = w[ix][iy][0];
		wi = w[ix][iy][1];
		tr = t[ix+posx][iy+posy][0];
		ti = t[ix+posx][iy+posy][1];
		w[ix][iy][0] = wr*tr - wi*ti;
		w[ix][iy][1] = wr*ti + wi*tr;
	} /* end for(iy.. ix .) */
} /* end transmit() */

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
	static int *hbeam=NULL,*kbeam=NULL;
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
			hbeam = (*muls).hbeam;
			kbeam = (*muls).kbeam;
			if ((hbeam == NULL) || (kbeam == NULL)) {
				printf("ERROR: hbeam or kbeam == NULL!\n");
				exit(0);
			}

			sprintf(fileAmpl,"%s/beams_amp.dat",(*muls).folder);
			sprintf(filePhase,"%s/beams_phase.dat",(*muls).folder);
			sprintf(fileBeam,"%s/beams_all.dat",(*muls).folder);
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
				fprintf(fp1," (%d,%d)", muls->hbeam[ib],  muls->kbeam[ib]);
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
			rPart = (*wave).wave[muls->hbeam[ib]][muls->kbeam[ib]][0];
			iPart = (*wave).wave[muls->hbeam[ib]][muls->kbeam[ib]][1];
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
