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

#ifndef DATA_CONTAINERS_H
#define DATA_CONTAINERS_H

#include <vector>
#include "stemtypes_fftw3.hpp"
#include "imagelib_fftw3.hpp"

#include "detectors.hpp"

class MULS {
public:
  MULS();
  
  int mode;                             /* determine the mode that this program runs in
					 * can be STEM, TEM, CBED ... */
  int printLevel;                       /* Flag indicating how much output should appear
					 * in the window. */
  int saveLevel;
  int complete_pixels;  //the number of pixels completed so far

  complex_tt ***trans;

#if FLOAT_PRECISION == 1
  fftwf_plan fftPlanPotInv,fftPlanPotForw;
#else
  fftw_plan fftPlanPotInv,fftPlanPotForw;
#endif

  float_tt **diffpat;
  float_tt czOffset;
  float_tt xOffset;
  float_tt yOffset;
  
  std::string input_ext, output_ext; /* file extensions for data input/output */

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
  bool adjustCubeSize;
  float_tt btiltx,btilty;   	        /* beam tilt in mrad*/
  bool tiltBack;               /* tilt back the wave below the specimen */
  std::vector<int> hbeams,kbeams;	/* arrays to hold recorded 
					   beam indicies */
  bool lbeams;				/* flag indicating, whether 
					   to record beams */	
  char filebeam[512];		 	/* file, that beams get recorded in */
  unsigned nbout;				/* number of recorded beams */

  //int nslic0;				/* slice counter */
  int mulsRepeat1;                      /* # of times to repeat structure */
  int mulsRepeat2;                      /* for REFINE mode # of mulsRun repeats */
  unsigned slices;                           /* number of different slices */
  bool centerSlices;                     /* flag indicating how to cut the sample */
  float_tt **pendelloesung;              /* pendelloesung plot for REFINE mode */
  float_tt ax,by,c;	                /* lattice parameters */
  float_tt cAlpha,cBeta,cGamma;
  double **Mm;                          /* metric matrix Mm(ax,by,cz,alpha,beta,gamma) */
  unsigned nCellX,nCellY,nCellZ;             /* number of unit cells in x-y-z dir*/
  int natom;				/* number of atoms in "atoms" */
  atom *atoms;				/* 3D atoms array */	
  float_tt atomRadius;                   /* for atom potential boxes */
  float_tt potOffsetX,potOffsetY;        /* offset of potential array from zero */
  float_tt potSizeX,potSizeY;            /* real space dimensions of potential array in A */
  unsigned potNx,potNy;                      /* size of projected potential in pixels */
  unsigned nx,ny;				/* size of wave function arrays */
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

  bool ismoth;                          /* smoothen the probe wave function */
  bool gaussFlag;
  float_tt gaussScale;
  int showProbe;
  unsigned displayProgInterval;             /* show progress every .. beam positions */
  unsigned displayPotCalcInterval;             /* show progress every .. beam positions */

  float_tt beamCurrent;  // pico Ampere
  float_tt dwellTime;    // msec
  float_tt electronScale;  // number of electrons

  
  int totalSliceCount;
  unsigned outputInterval;    // output results every n slices

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
  bool gaussianProp;    /* convolute fresnel propagator with gaussian or not */
  bool periodicZ;      /* for slicecell (make non periodic in Z */
  bool periodicXY;       /* for slicecell (make non periodic in x,y */
  bool bandlimittrans;  /* flag for bandwidth limiting transmission function */
  bool fftpotential;    /* flag indicating that we should use FFT for V_proj calculation */
  bool plotPotential;
  bool storeSeries;
  bool tds;
  bool Einstein;        /* if set (default=set), the Einstein model will be used */
  std::string phononFile;    /* file name for detailed phonon modes */
  int atomKinds;
  int *Znums;
  double **rPotential;   /* array containing real space potential LUT for each atom kind present */
  double *sfkArray;
  double **sfTable;
  int sfNk;              /* number of k-points in sfTable and sfkArray */
  double *u2,*u2avg;     /* (current/averaged) rms displacement of atoms */
  float_tt tds_temp;
  bool savePotential;
  bool saveTotalPotential;
  bool readPotential;
  float_tt scanXStart,scanXStop,scanYStart,scanYStop;
  unsigned scanXN,scanYN;
  float_tt intIntensity;
  std::string folder;
  unsigned avgRuns;
  bool potential3D;
  int scatFactor;
  int Scherzer;
  std::vector<double> chisq;
  int cellDiv;
  bool equalDivs;           // this flag indicates whether we can reuse already pre-calculated potential data

  /* Parameters for STEM-detectors */
  int detectorNum;
  /* we will alow as many detector 
			   definitions as the user wants */
  DetectorMgrPtr detectors;
  // std::vector<std::vector<DetectorPtr> > detectors;
  //DETECTOR *detectors;
  bool save_output_flag;
  
  double *dE_EArray;

  // Tomography parameters:
  float_tt tomoTilt;  // current tilt in tomography series
  float_tt tomoStart; // in rad
  float_tt tomoStep;  // in rad
  int    tomoCount;  // number of diffraction patterns.
  float_tt zoomFactor; // increases the size of the super-box in x,y, in order to
                     // make full use of atoms present, creates vacuum edge around sample.

}; 

#endif
