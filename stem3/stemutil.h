#ifndef STEMUTIL_H
#define STEMUTIL_H

#include "data_containers.h"
#include "stemtypes_fftw3.h"

#define pNPIX		0   /* number of pix 1 for real and 2 for complex */
#define	pRMAX		1	/* maximum value of the real part of the image */
#define	pIMAX		2	/* maximum value of the imaginary part of the image */
#define pRMIN		3	/* minimum value of the real part of the image */
#define	pIMIN		4	/* minimum value of the imaginary part of the image */
#define pXBTILT		5	/* x beam tilt in rad */
#define pYBTILT		6	/* y beam tilt in rad */
#define pC		7	/* c unit cell dimension in Angstroms */
#define pRES		8	/* real space resolution in atompot */
#define pXCTILT		9	/* x crystal tilt in rad */
#define pYCTILT		10	/* y crystal tilt in rad */
#define pDEFOCUS	11	/* defocus in Angstroms */
#define pASTIG		12	/* astigmatism in Angstroms */
#define pTHETA		13	/* angle of astigmatism in radians */
#define pDX		14	/* dimension of pixel in x direction in Angstroms */
#define pDY		15	/* dimension of pixel in y direction in Angstroms */
#define pENERGY		16	/* beam energy in keV */
#define pOAPERT		17	/* objective aperture semi-angle in radians */
#define pCS		18	/* spherical aberration in Angstroms */
#define pWAVEL		19	/* electron wavelength in Angstroms */
#define pCAPERT		21	/* condenser (CTEM) illumination angle in radians */
#define pDDF		22	/* defocus spread in Angstroms */
#define pNSLICES	29	/* number of slices */
#define pMINDET		31	/* minimum detector angle (STEM) in radians */
#define pMAXDET		32	/* maximum detector angle (STEM) in radians */

typedef struct timeval timev;
typedef struct timezone timez;

int phononDisplacement(float_tt *u,MULS *muls,int id,int icx,int icy,
		       int icz,int atomCount,float_tt dw,int maxAtom, int Znum);

void *memcopy(void *dest, const void *src, size_t n);
// void saveSTEMimages(MULS *muls);
// atom *readUnitCell(int &natom,char *fileName,MULS *muls,int handleVacancies);
atom *readCFGUnitCell(int &natom,char *fileName,MULS *muls);

float_tt wavelength( float_tt kev );

float_tt gasdev(long *idum);
float_tt rangauss( unsigned long *iseed );
float_tt sigma( float_tt kev );
//void splinh( QSfVec &x, QSfVec &y, QSfMat &coeffs, int n);
//void splinh( float_tt x[], float_tt y[],
	     //float_tt b[], float_tt c[], float_tt d[], int n);
//float_tt seval( float_tt *x, float_tt *y, float_tt *b, float_tt *c,
	     //float_tt *d, int n, float_tt x0 );
float_tt ranflat( unsigned long *iseed );
int parlay( const char c[], int islice[], int nsmax, int lmax,
	    int *nslice, int fperr );
/* long powerof2( long n ); */
float_tt sfLUT(float_tt s,int atKind, MULS *muls);
float_tt bicubic(QSfMat ff,int Nz, int Nx,float_tt z,float_tt x);


float_tt getTime();
float_tt cputim();



#endif


