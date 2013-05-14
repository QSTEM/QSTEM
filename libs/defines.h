////////////////////////////////////////////////////////////////////////
// define whether to use single or double precision
///////////////////////////////////////////////////////////////////////
#define FLOAT_PRECISION 1


#define BW (2.0F/3.0F)	/* bandwidth limit */
#define DOYLE_TURNER 0
#define WEICK_KOHL 1
#define CUSTOM 2
#define STEM    1
#define CBED    2
#define TEM     3
#define REFINE  4
#define MSCBED  5
#define TOMO    6

////////////////////////////////////////////////////////////////////////
// Define physical constants
////////////////////////////////////////////////////////////////////////
#define ELECTRON_CHARGE (1.6021773e-19)
#define PICO_AMPERE (1e-12/ELECTRON_CHARGE)
#define MILLISEC_PICOAMP (1e-3*PICO_AMPERE)

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

#define NCINMAX  500	/* max number of characers in stacking spec */
#define NRMAX	50	/* number of values in look-up-table in vzatomLUT */
#define RMIN	0.01	/* min r (in Ang) range of LUT for vzatomLUT() */
#define RMAX	5

#define NPDTMAX 8       /* number of parameters for doyle turner sfacts */
#define NPMAX	12	/* number of parameters for each Z */
//#define NZMIN	1	/* min Z in featom.tab */
//#define NZMAX	98      /* max Z (in featom.dat ZMAX=103) */
//#define	NCMAX	132	/* characters per line to read */
//#define NPARAM	64	/* number of parameters in tiff files */

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
#define PI180 1.7453292519943e-2
#define USE_REZ_SFACTS    1  // used in getAtomPotential3D and getAtomPotentialOffset3D 
#define Z_INTERPOLATION   0  // used in make3DSlices (central function for producing atom potential slices)
#define USE_Q_POT_OFFSETS 1  // used in make3DSlices (central function for producing atom potential slices)

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
