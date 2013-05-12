#ifndef DATA_CONTAINERS_H
#define DATA_CONTAINERS_H

#include <vector>
// #include "stemtypes_fftw3.h"
#include "memory_fftw3.h"

// a structure for a probe/parallel beam wavefunction.
// Separate from mulsliceStruct for parallelization.
class WAVEFUNC 
{
public:
	int iPosX,iPosY;      /* integer position of probe position array */
	int nx, ny;			/* size of diffpat arrays */
	int detPosX,detPosY;
	char fileStart[512];
	char fileout[512];
	float2DArray diffpat;
	//real **diffpat;
	float2DArray avgArray;
	//real **avgArray;
	char avgName[512];
	float_tt thickness;
	float_tt intIntensity;

	float_tt rmin,rmax;		/* min and max of real part */
	float_tt aimin,aimax;		/* min and max of imag part */
	float1DArray kx2, ky2, kx, ky;

	fftwf_plan fftPlanWaveForw,fftPlanWaveInv;
	complex2DArray wave; /* complex wave function */

public:
	// initializing constructor:
	WAVEFUNC(int nx, int ny);
	// define a copy constructor to create new arrays
	WAVEFUNC( WAVEFUNC& other );

	void ZeroWave(void);
};

class MULS {
public:
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
  complex3DArray trans;
#else
  fftw_plan fftPlanPotInv,fftPlanPotForw;
  // wave moved to probeStruct
  //fftw_complex  **wave; /* complex wave function */
  fftw_complex ***trans;
#endif

  float2DArray diffpat;
  real czOffset;
  real xOffset;
  real yOffset;
  
  char cin2[1024];				/* stacking sequence */
  char fileBase[512];
  char **filein;			/* array of input potential files */
  int lpartl, lstartl;	                /* flags indicating partial 
					   coherence */
  char atomPosFile[512];
                                        /* and start wavefunction */	
  float_tt v0;				/* inc. beam energy */
  float_tt resolutionX;                  /* real space pixelsize for wave function and potential */
  float_tt resolutionY;                  /* real space pixelsize for wave function and potential */
  float_tt ctiltx,ctilty,ctiltz;	        /* crystal tilt in mrad */
  char cfgFile[512];                        /* file name for writing tilted atomic configuration */
  float_tt cubex,cubey,cubez;            /* dimension of crystal cube, if zero, then nx,ny,nz *
					 * will be used */
  int adjustCubeSize;
  float_tt btiltx,btilty;   	        /* beam tilt in mrad*/
  int tiltBack;               /* tilt back the wave below the specimen */
  int *hbeam,*kbeam;		        /* arrays to hold recorded 
					   beam indicies */
  int lbeams;				/* flag indicating, whether 
					   to record beams */	
  char filebeam[512];		 	/* file, that beams get recorded in */
  int nbout;				/* number of recorded beams */

  //int nslic0;				/* slice counter */
  int mulsRepeat1;                      /* # of times to repeat structure */
  int mulsRepeat2;                      /* for REFINE mode # of mulsRun repeats */
  int slices;                           /* number of different slices */
  int centerSlices;                     /* flag indicating how to cut the sample */
  float2DArray pendelloesung;              /* pendelloesung plot for REFINE mode */
  float_tt ax,by,c;	                /* lattice parameters */
  float_tt cAlpha,cBeta,cGamma;
  double2DArray Mm;                          /* metric matrix Mm(ax,by,cz,alpha,beta,gamma) */
  int nCellX,nCellY,nCellZ;             /* number of unit cells in x-y-z dir*/
  int natom;				/* number of atoms in "atoms" */
  atom *atoms;				/* 3D atoms array */	
  float_tt atomRadius;                   /* for atom potential boxes */
  float_tt potOffsetX,potOffsetY;        /* offset of potential array from zero */
  float_tt potSizeX,potSizeY;            /* real space dimensions of potential array in A */
  int potNx,potNy;                      /* size of projected potential in pixels */
  int nx,ny;				/* size of wave function arrays */
  int avgCount;
  //float_tt thickness;

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
  float1DArray sparam;

  int saveFlag;			/* flag indicating, whether to save the result */
  float_tt k2max;

  // propagators used in propagate_slow
  float1DArray propxr, propyr;
  float1DArray propxi, propyi;

  int nlayer;
  float1DArray cz;
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
  char phononFile[512];    /* file name for detailed phonon modes */
  int atomKinds;
  int *Znums;
  //boost::shared_ptr<boost::shared_ptr<double>> rPotential;   /* array containing real space potential LUT for each atom kind present */
  double2DArray sfkArray;
  double2DArray sfTable;
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
  char folder[1024];
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
  /* we will alow as many detector 
			   definitions as the user wants */
  std::vector<std::vector<DETECTOR>> detectors;
  //DETECTOR *detectors;
  int save_output_flag;
  
  double1DArray dE_EArray;

  // Tomography parameters:
  double tomoTilt;  // current tilt in tomography series
  double tomoStart; // in rad
  double tomoStep;  // in rad
  int    tomoCount;  // number of diffraction patterns.
  double zoomFactor; // increases the size of the super-box in x,y, in order to
                     // make full use of atoms present, creates vacuum edge around sample.

}; 

#endif