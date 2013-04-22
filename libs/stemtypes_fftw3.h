#ifndef STEMTYPES_H
#define STEMTYPES_H

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

// #include "floatdef.h"
#include "fftw3.h"




////////////////////////////////////////////////////////////////
#if FLOAT_PRECISION == 1
#define fftw_real float
#ifndef float_tt
#define float_tt  float
#endif
#define real      float
#else  // FLOAT_PRECISION
#define fftw_real double
#ifndef float_tt
#define float_tt  float
#endif
#define real      float
#endif  // FLOAT_PRECISION
////////////////////////////////////////////////////////////////

typedef struct atomStruct {
  float z,y,x;
  // float dx,dy,dz;  // thermal displacements
  float dw;      // Debye-Waller factor
  float occ;     // occupancy
  float q;       // charge 
  int Znum;
} atom;

/* Planes will be defined by the standard equation for a plane, i.e.
 * a point (point) and 2 vectors (vect1, vect2)
 */
typedef struct planeStruct {
  double normX,normY,normZ;
  double vect1X,vect1Y,vect1Z;
  double vect2X,vect2Y,vect2Z;
  double pointX,pointY,pointZ;
} plane;

typedef struct grainBoxStruct {
  int amorphFlag;
  double density,rmin, rFactor;  /* density, atomic distance, reduced atomic distance
				  * for amorphous material.  Red. r is for making a hex.
				  * closed packed structure, which will fill all space, 
				  * but will only be sparsely filled, and later relaxed.
				  */ 
  char *name;
  atom *unitCell; /* definition of unit cell */
  int natoms;     /* number of atoms in unit cell */
  double ax,by,cz; /* unit cell parameters */
  double alpha, beta, gamma; /* unit cell parameters */
  double tiltx,tilty,tiltz;
  double shiftx,shifty,shiftz;
  plane *planes;   /* pointer to array of bounding planes */
  double sphereRadius, sphereX,sphereY,sphereZ; /* defines a sphere instead of a grain with straight edges */
  int nplanes; /* number of planes in array planes */
} grainBox;

typedef struct superCellBoxStruct {
  double cmx,cmy,cmz;  /* fractional center of mass coordinates */
  double ax,by,cz;
  int natoms;
  atom *atoms; /* contains all the atoms within the super cell */
} superCellBox;

typedef struct atomBoxStruct {
  int used;   /* indicate here whether this atom is used in the
		 particular problem */
  int nx,ny,nz;
  float_tt dx,dy,dz;
  double B;
#if FLOAT_PRECISION == 1
  fftwf_complex ***potential;   /* 3D array containg 1st quadrant of real space potential */
  float_tt ***rpotential;   /* 3D array containg 1st quadrant of real space potential */
#else
  fftw_complex ***potential;   /* 3D array containg 1st quadrant of real space potential */
  float_tt ***rpotential;   /* 3D array containg 1st quadrant of real space potential */
#endif
} atomBox;

typedef struct detectorStruct {
  float_tt rInside,rOutside;
  float_tt k2Inside,k2Outside;
  char name[32];
  float_tt **image;        // place for storing avg image = sum(data)/Navg
  float_tt **image2;        // we will store sum(data.^2)/Navg 
  float_tt error;  
  float_tt shiftX,shiftY;
  int Navg;
} DETECTOR;

// a structure for a probe/parallel beam wavefunction.
// Separate from mulsliceStruct for parallelization.
typedef struct waveStruct {
	int iPosX,iPosY;      /* integer position of probe position array */
	int detPosX,detPosY;
	char *fileStart;
	char *fileout;
	real **diffpat;
	real **avgArray;
	char *avgName;
#if FLOAT_PRECISION == 1
  fftwf_plan fftPlanWaveForw,fftPlanWaveInv;
  fftwf_complex  **wave; /* complex wave function */
#else
  fftw_plan fftPlanWaveForw,fftPlanWaveInv;
  fftw_complex  **wave; /* complex wave function */
#endif
  
} WAVEFUNC;

typedef struct mulsliceStruct {
  int mode;                             /* determine the mode that this program runs in
					 * can be STEM, TEM, CBED ... */
  int printLevel;                       /* Flag indicating how much output should appear
					 * in the window. */
  int saveLevel;
  int complete_pixels;  //the number of pixels completed so far

#if FLOAT_PRECISION == 1
  fftwf_plan fftPlanPotInv,fftPlanPotForw;
  // wave moved to probeStruct
  //fftwf_complex  **wave; /* complex wave function */
  fftwf_complex ***trans;
#else
  fftw_plan fftPlanPotInv,fftPlanPotForw;
  // wave moved to probeStruct
  //fftw_complex  **wave; /* complex wave function */
  fftw_complex ***trans;
#endif

  real **diffpat;
  real czOffset;
  real xOffset;
  real yOffset;
  
  char *cin2;				/* stacking sequence */
  char *fileBase;
  char **filein;			/* array of input potential files */
  int lpartl, lstartl;	                /* flags indicating partial 
					   coherence */
  char *atomPosFile;
                                        /* and start wavefunction */	
  float_tt v0;				/* inc. beam energy */
  float_tt resolutionX;                  /* real space pixelsize for wave function and potential */
  float_tt resolutionY;                  /* real space pixelsize for wave function and potential */
  float_tt ctiltx,ctilty,ctiltz;	        /* crystal tilt in mrad */
  char *cfgFile;                        /* file name for writing tilted atomic configuration */
  float_tt cubex,cubey,cubez;            /* dimension of crystal cube, if zero, then nx,ny,nz *
					 * will be used */
  int adjustCubeSize;
  float_tt btiltx,btilty;   	        /* beam tilt in mrad*/
	int tiltBack;               /* tilt back the wave below the specimen */
  int *hbeam,*kbeam;		        /* arrays to hold recorded 
					   beam indicies */
  int lbeams;				/* flag indicating, whether 
					   to record beams */	
  char *filebeam;		 	/* file, that beams get recorded in */
  int nbout;				/* number of recorded beams */

  int nslic0;				/* slice counter */
  int mulsRepeat1;                      /* # of times to repeat structure */
  int mulsRepeat2;                      /* for REFINE mode # of mulsRun repeats */
  int slices;                           /* number of different slices */
  int centerSlices;                     /* flag indicating how to cut the sample */
  float_tt **pendelloesung;              /* pendelloesung plot for REFINE mode */
  float_tt ax,by,c;	                /* lattice parameters */
  float_tt cAlpha,cBeta,cGamma;
  double **Mm;                          /* metric matrix Mm(ax,by,cz,alpha,beta,gamma) */
  int nCellX,nCellY,nCellZ;             /* number of unit cells in x-y-z dir*/
  int natom;				/* number of atoms in "atoms" */
  atom *atoms;				/* 3D atoms array */	
  float_tt atomRadius;                   /* for atom potential boxes */
  float_tt potOffsetX,potOffsetY;        /* offset of potential array from zero */
  float_tt potSizeX,potSizeY;            /* real space dimensions of potential array in A */
  int potNx,potNy;                      /* size of projected potential in pixels */
  int nx,ny;				/* size of wave function arrays */
  int avgCount;
  float_tt thickness;

  float_tt C5;
  float_tt dE_E;
  float_tt dV_V;
  float_tt dI_I;
  float_tt alpha;
  float_tt sourceRadius;
  float_tt Cc;
  float_tt df0;				/* defocus */
  float_tt astigMag;				/* astigmatism*/
  float_tt astigAngle;				/* angle of astigmatism */

  int ismoth;                          /* smoothen the probe wave function */
  int gaussFlag;
  float_tt gaussScale;
  int showProbe;
  int displayProgInterval;             /* show progress every .. beam positions */
  int displayPotCalcInterval;             /* show progress every .. beam positions */

  float_tt beamCurrent;  // pico Ampere
  float_tt dwellTime;    // msec
  float_tt electronScale;  // number of electrons

  
  int totalSliceCount;
  int outputInterval;    // output results every n slices

  float_tt aobj;				/* obj aperture */
  float_tt aAIS;                         /* condensor aperture in A (projected size, */
                                        /* for Koehler illumination)                */
  // float_tt areaAIS;                      /* fractional area illuminated by AIS (def=1) */
  float_tt Cs;			      	/* spher. aberration */
  /////////////////////////////////////////////
  // more aberrations:
  float_tt a33;
  float_tt a31;
  float_tt a44;
  float_tt a42;
  float_tt a55;
  float_tt a53;
  float_tt a51;
  float_tt a66;
  float_tt a64;
  float_tt a62;
  float_tt phi33;
  float_tt phi31;
  float_tt phi44;
  float_tt phi42;
  float_tt phi55;
  float_tt phi53;
  float_tt phi51;
  float_tt phi66;
  float_tt phi64;
  float_tt phi62;


  float_tt acmax,acmin;
  float_tt sigmaf;
  float_tt dfdelt;
  float_tt dfa2,dfa3;
  float_tt dfa2phi,dfa3phi;
  float_tt chi,phi;
  float_tt *sparam;

  int saveFlag;			/* flag indicating, whether to save the result */
  float_tt rmin,rmax;		/* min and max of real part */
  float_tt aimin,aimax;		/* min and max of imag part */
  float_tt *kx2,*ky2,k2max,*kx,*ky;

  int nlayer;
  float_tt *cz;
  float_tt sliceThickness;
  int onlyFresnel;
  int startSpherical;
  float_tt startDistance;
  float_tt maxAngle;
  float_tt gaussWidth;
  float_tt defInfinity;
  int accumulateIntensity;
  int deconvolute;
  int showPhaseplate;
  int normHolog;
  int gaussianProp;    /* convolute fresnel propagator with gaussian or not */
  int nonPeriodZ;      /* for slicecell (make non periodic in Z */
  int nonPeriod;       /* for slicecell (make non periodic in x,y */
  int bandlimittrans;  /* flag for bandwidth limiting transmission function */
  int fftpotential;    /* flag indicating that we should use FFT for V_proj calculation */
  int plotPotential;
  int storeSeries;
  int tds;
  int Einstein;        /* if set (default=set), the Einstein model will be used */
  char *phononFile;    /* file name for detailed phonon modes */
  int atomKinds;
  int *Znums;
  double **rPotential;   /* array containing real space potential LUT for each atom kind present */
  double *sfkArray;
  double **sfTable;
  int sfNk;              /* number of k-points in sfTable and sfkArray */
  double *u2,*u2avg;     /* (current/averaged) rms displacement of atoms */
  float_tt tds_temp;
  int savePotential;
  int saveTotalPotential;
  int readPotential;
  float_tt scanXStart,scanXStop,scanYStart,scanYStop;
  int scanXN,scanYN;
  float_tt intIntensity;
  double imageGamma;
  char *folder;
  int avgRuns;
  int potential3D;
  int scatFactor;
  int Scherzer;
  double *chisq;
  int webUpdate;
  int cellDiv;
  int equalDivs;           // this flag indicates whether we can reuse already pre-calculated potential data

  /* Parameters for STEM-detectors */
  int detectorNum;
  DETECTOR *detectors;  /* we will alow as many detector 
			   definitions as the user wants */
  int save_output_flag;
  
  double *dE_EArray;

  // Tomography parameters:
  double tomoTilt;  // current tilt in tomography series
  double tomoStart; // in rad
  double tomoStep;  // in rad
  int    tomoCount;  // number of diffraction patterns.
  double zoomFactor; // increases the size of the super-box in x,y, in order to
                     // make full use of atoms present, creates vacuum edge around sample.

} MULS; 

#endif // STEMTYPES_H
