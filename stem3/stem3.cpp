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

#define VERSION 2.30
#define VIB_IMAGE_TEST

#ifndef _WIN32
#define UNIX 
#endif
/* #define USE_FFT_POT */
// for memory leak checking in windows.  Should not affect speed of release builds.
#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#ifdef _WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif
#include <string.h>
#ifndef _WIN32
#ifdef __cplusplus
#include <cmath>
#else
#include <math.h>
#endif
#else
#include <math.h>
#endif

#include <time.h>
#include <ctype.h>
#include <sys/stat.h>
// #include <stat.h>

#include <omp.h>

#include "stem3.hpp"

#include "memory_fftw3.hpp"	/* memory allocation routines */
#include "readparams.hpp"
#include "imagelib_fftw3.hpp"
#include "fileio_fftw3.hpp"
#include "matrixlib.hpp"
#include "stemlib.hpp"
#include "stemutil.hpp"
// #include "weblib.hpp"
// 20131222 - customslice doesn't seem to actually be used.
//#include "customslice.hpp"
#include "data_containers.hpp"

//#include "readparams.hpp"

#define NCINMAX 1024
#define NPARAM	64    /* number of parameters */
#define MAX_SCANS 1   /* maximum number of linescans per graph window */
#define PHASE_GRATING 0
#define BUF_LEN 256

#define DELTA_T 1     /* number of unit cells between pictures */
#define PICTS 5      /* number of different thicknesses */
#define NBITS 8	       /* number of bits for writeIntPix */

/* global variable: */
MULS muls;

void usage() {
	printf("usage: stem [input file='stem.dat']\n\n");
}


/***************************************************************
***************************************************************
* MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN *
***************************************************************
**************************************************************/


int main(int argc, char *argv[]) {
  int i; 
  double timerTot;
  char cinTemp[BUF_LEN];

  timerTot = cputim();
  for (i=0;i<BUF_LEN;i++)
    cinTemp[i] = 0;

#ifdef UNIX
  system("date");
#endif

  std::string fileName;
  /*************************************************************
   * read in the parameters
   ************************************************************/  
  if (argc < 2)
    fileName = "stem.dat";
  else
    fileName=argv[1];
  // Initialize the config file reader
  ConfigReaderPtr configReader = GetConfigReader(fileName);
  
  if (!configReader->IsValid())
    {
      printf("could not open input file %s!\n",fileName.c_str());
      exit(0);
      usage();
    }
        
  // Read potential parameters and initialize a pot object
  WavePtr initialWave = WavePtr(new WAVEFUNC(configReader));
  PotPtr potential = GetPotential(configReader);
  readFile(configReader);

  displayParams(initialWave, potential);
#ifdef _OPENMP
  omp_set_dynamic(1);
#endif

  int mode;
  configReader->ReadMode(mode);
  switch (mode) {
  case CBED:   doCBED(initialWave, potential);   break;	
  case STEM:   doSTEM(initialWave, potential);   break;
  case TEM:    doTEM(initialWave, potential);    break;
    //case MSCBED: doMSCBED(initialWave, potential); break;
    //case TOMO:   doTOMO(initialWave, potential);   break;
    // case REFINE: doREFINE(); break;
  default:
    printf("Mode not supported\n");
  }

#if _DEBUG
  _CrtDumpMemoryLeaks();
#endif

  return 0;
}

/***************************************************************
***************************************************************
** End of MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN  ** 
***************************************************************
**************************************************************/

void initMuls() {
	int sCount,i,slices;

	slices = muls.slices;

       	for (sCount =0;sCount<slices;sCount++)
		muls.cin2[sCount] = 'a'+sCount;
	for (sCount = slices;sCount < NCINMAX;sCount++)
		muls.cin2[sCount] = 0;
	muls.nlayer = slices;
	
       	muls.sparam = (float *)malloc(NPARAM*sizeof(float));
	for (i=0;i<NPARAM;i++)
		muls.sparam[i] = 0.0;
}

/************************************************************************
*
***********************************************************************/

void displayParams(WavePtr &wave, PotPtr &pot) {
  FILE *fpDir;
  char systStr[64];
  double k2max,temp;
  int i,j;
  static char Date[16],Time[16];
  time_t caltime;
  struct tm *mytime;
  const double pi=3.1415926535897;

  /*
  if (wave->printLevel < 1) {
    if ((fpDir = fopen(muls.folder.c_str(),"r"))) {
      fclose(fpDir);
      // printf(" (already exists)\n");
    }
    else {
      sprintf(systStr,"mkdir %s",muls.folder.c_str());
      system(systStr);
      // printf(" (created)\n");
    }	  
    return;
  }
  */
  caltime = time( NULL );
  mytime = localtime( &caltime );
  strftime( Date, 12, "%Y:%m:%d", mytime );
  strftime( Time, 9, "%H:%M:%S", mytime );
  
  printf("\n*****************************************************\n");
  printf("* Running program STEM3 (version %.2f) in %s mode\n",VERSION,
         (muls.mode == STEM) ? "STEM" : (muls.mode==TEM) ? "TEM" : 
         (muls.mode == CBED) ? "CBED" : (muls.mode==TOMO)? "TOMO" : 
         "???"); 
  printf("* Date: %s, Time: %s\n",Date,Time);
	
  printf("* Input file:           %s\n",muls.atomPosFile);
  
  /* create the data folder ... */
  printf("* Data folder:          ./%s/ ",muls.folder.c_str()); 
	
  printf("* Super cell divisions: %d (in z direction) %s\n",muls.cellDiv,muls.equalDivs ? "equal" : "non-equal");
  printf("* Slices per division:  %d (%gA thick slices [%scentered])\n",
         muls.slices,muls.sliceThickness,(muls.centerSlices) ? "" : "not ");
  printf("* Output every:         %d slices\n",muls.outputInterval);
  

  /* 
     if (muls.ismoth) printf("Type 1 (=smooth aperture), ");
     if (muls.gaussFlag) printf("will apply gaussian smoothing"); 
     printf("\n");
  */

  /***************************************************/
  /*  printf("Optimizing fftw plans according to probe array (%d x %dpixels = %g x %gA) ...\n",
      muls.nx,muls.ny,muls.nx*muls.resolutionX,muls.ny*muls.resolutionY);
  */
  
  printf("* TDS:                  %d runs)\n",muls.avgRuns);

  printf("*\n*****************************************************\n");
}


void readArray(FILE *fp, char *title,double *array,int N) {
	int i;
	char buf[512],*str;

	if (!readparam(fp, title,buf,1)) printf("%s array not found - exit\n",title), exit(0);
	i=0;
	str = buf;
	if (strchr(" \t\n",*str) != NULL) str = strnext(str," \t");
	while (i<N) {
		array[i++] = atof(str);
		str = strnext(str," \t\n");
		while (str == NULL) {
                  if (!readNextLine(fp, buf,511)) 
                    printf("Incomplete reading of %s array - exit\n",title), exit(0);
			str = buf;
			if (strchr(" \t\n",*str) != NULL) str = strnext(str," \t\n");
		}  
	}
}

/***********************************************************************
* readSFactLUT() reads the scattering factor lookup table from the 
* input file
**********************************************************************/
/*
// TODO: 20131222 - does not seem to be used anywhere
void readSFactLUT() {
	int Nk,i,j;
	double **sfTable=NULL;
	double *kArray = NULL;
	char buf[256], elem[8];

	if (readparam("Nk:",buf,1))
		Nk = atoi(buf);
	else {
		printf("Could not find number of k-points for custom scattering factors (Nk)\n");
		exit(0);
	}

	// allocate memory for sfTable and kArray:
	sfTable = double2D(muls.atomKinds,Nk+1,"sfTable");
	kArray  = double1D(Nk+1,"kArray");

	// read the k-values:
	readArray("k:",kArray,Nk);
	kArray[Nk] = 2.0*kArray[Nk-1];

	for (j=0;j<muls.atomKinds;j++) {
		elem[3] = '\0';
		elem[2] = ':';
		elem[0] = elTable[2*muls.Znums[j]-2];
		elem[1] = elTable[2*muls.Znums[j]-1];
		if (elem[1] == ' ') {
			elem[1] = ':';
			elem[2] = '\0';
		}
		// printf("%s\n",elem);
		readArray(elem,sfTable[j],Nk);
		sfTable[j][Nk] = 0.0;
	}

	if (0) {
		printf("k: ");
		for (i=0;i<=Nk;i++) printf("%.3f ",kArray[i]);
		for (j=0;j<muls.atomKinds;j++) {
			printf("\n%2d: ",muls.Znums[j]);
			for (i=0;i<=Nk;i++) printf("%.3f ",sfTable[j][i]);
		}
		printf("\n");
	}
	muls.sfTable = sfTable;
	muls.sfkArray = kArray;
	muls.sfNk = Nk+1;
}
*/

/************************************************************************
* readFile() 
*
* reads the parameters from the input file and does some 
* further setup accordingly
*
***********************************************************************/
/*
void readFile(ConfigReaderPtr &configReader) {
	char answer[256],*strptr;
	FILE *fpTemp;
	float ax,by,c;
	char buf[BUF_LEN],*strPtr;
	int i,ix;
	int potDimensions[2];
	long ltime;
	unsigned iseed;
	double dE_E0,x,y,dx,dy;
	const double pi=3.1415926535897;


	ltime = (long) time(NULL);
	iseed = (unsigned) ltime;
	muls.cubex = 0.0;
	muls.cubey = 0.0;
	muls.cubez = 0.0;

        configReader->ReadMode(muls.mode);
        
        if (muls.mode == STEM)
          configReader->ReadOutputLevel(muls.printLevel, muls.saveLevel, muls.displayPotCalcInterval,
                                        muls.displayProgInterval);
        else
          configReader->ReadOutputLevel(muls.printLevel, muls.saveLevel, muls.displayPotCalcInterval);            

        configReader->ReadNCells(muls.nCellX, muls.nCellY, muls.nCellZ, muls.cellDiv);

	// Read the beam tilt parameters	
        configReader->ReadBeamTilt(muls.btiltx, muls.btilty, muls.tiltBack);

	//Read the crystal tilt parameters
        configReader->ReadCrystalCubeAndTilt(muls.ctiltx, muls.ctilty, muls.ctiltz,
                                             muls.cubex, muls.cubey, muls.cubez,
                                             muls.adjustCubeSize);

	// temperature related data must be read before reading the atomic positions:
	//
        configReader->ReadTemperatureData(muls.tds, muls.tds_temp, muls.phononFile, muls.Einstein);

	// Read the atomic model positions !!!
	sprintf(muls.atomPosFile,muls.fileBase);
	// remove directory in front of file base: 
	while ((strptr = strchr(muls.fileBase,'\\')) != NULL) strcpy(muls.fileBase,strptr+1);

	// add a '_' to fileBase, if not existent 
	if (strrchr(muls.fileBase,'_') != muls.fileBase+strlen(muls.fileBase)-1) {
		if ((strPtr = strchr(muls.fileBase,'.')) != NULL) sprintf(strPtr,"_");
		else strcat(muls.fileBase,"_");
	}
	if (strchr(muls.atomPosFile,'.') == NULL) {
          // take atomPosFile as is, or add an ending to it, if it has none yet
          if (strrchr(muls.atomPosFile,'.') == NULL) {
            sprintf(buf,"%s.cssr",muls.atomPosFile);
            if ((fpTemp=fopen(buf,"r")) == NULL) {
              sprintf(buf,"%s.cfg",muls.atomPosFile);
              if ((fpTemp=fopen(buf,"r")) == NULL) {
                printf("Could not find input file %s.cssr or %s.cfg\n",
                       muls.atomPosFile,muls.atomPosFile);
                exit(0);
              }
              strcat(muls.atomPosFile,".cfg");
              fclose(fpTemp);
            }
            else {
              strcat(muls.atomPosFile,".cssr");
              fclose(fpTemp);
            }
          }
	}
        configReader->ReadSliceOffset(muls.xOffset, muls.yOffset);

	// the last parameter is handleVacancies.  If it is set to 1 vacancies 
	// and multiple occupancies will be handled. 
	// _CrtSetDbgFlag  _CRTDBG_CHECK_ALWAYS_DF();
	// printf("memory check: %d, ptr= %d\n",_CrtCheckMemory(),(int)malloc(32*sizeof(char)));

	muls.atoms = readUnitCell(&(muls.natom),muls.atomPosFile,&muls,1);

	// printf("memory check: %d, ptr= %d\n",_CrtCheckMemory(),(int)malloc(32*sizeof(char)));

	if (muls.atoms == NULL) {
		printf("Error reading atomic positions!\n");
		exit(0);
	}
	if (muls.natom == 0) {
		printf("No atom within simulation boundaries!\n");
		exit(0);
	}
	// printf("hello!\n");
	ax = muls.ax/muls.nCellX;
	by = muls.by/muls.nCellY;;
	c =  muls.c/muls.nCellZ;


	// Done reading atomic positions 

        configReader->ReadProbeArraySize(muls.nx, muls.ny);
        configReader->ReadResolution(muls.resolutionX, muls.resolutionY);
        configReader->ReadVoltage(muls.v0);
        configReader->ReadSliceParameters(muls.centerSlices, muls.sliceThickness,
                                          muls.slices, muls.outputInterval, muls.czOffset);

        // TODO: this section is horribly broken
	if (readparam("slice-thickness:",buf,1)) {
		sscanf(buf,"%g",&(muls.sliceThickness));
		if (readparam("slices:",buf,1)) {
			sscanf(buf,"%d",&(muls.slices));
		}
		else {
			if (muls.cubez >0)
				muls.slices = (int)(muls.cubez/(muls.cellDiv*muls.sliceThickness)+0.99);
			else
				muls.slices = (int)(muls.c/(muls.cellDiv*muls.sliceThickness)+0.99);
		}
		muls.slices += muls.centerSlices;
	}
	else {
		muls.slices = 0; 
		if (readparam("slices:",buf,1)) {
			sscanf(buf,"%d",&(muls.slices));
			if (muls.sliceThickness == 0.0) {
				if ((muls.slices == 1) && (muls.cellDiv == 1)) {
					if (muls.cubez >0)
						muls.sliceThickness = (muls.centerSlices) ? 2.0*muls.cubez/muls.cellDiv : muls.cubez/muls.cellDiv;
					else
						muls.sliceThickness = (muls.centerSlices) ? 1.0*muls.c/(muls.cellDiv) : muls.c/(muls.cellDiv);
				}
				else {
					if (muls.cubez >0) {
						muls.sliceThickness = muls.cubez/(muls.cellDiv*muls.slices-muls.centerSlices);
					}
					else {
						muls.sliceThickness = muls.c/(muls.cellDiv*muls.slices);
					}
				}
			}
			else {
				muls.cellDiv = (muls.cubez >0) ? (int)ceil(muls.cubez/(muls.slices*muls.sliceThickness)) :
					(int)ceil(muls.c/(muls.slices*muls.sliceThickness));
			if (muls.cellDiv < 1) muls.cellDiv = 1;
			}
		}
	}
	if (muls.slices == 0) {
		if (muls.printLevel > 0) printf("Error: Number of slices = 0\n");
		exit(0);
	}
	// Find out whether we need to recalculate the potential every time, or not

	muls.equalDivs = ((!muls.tds)  && (muls.nCellZ % muls.cellDiv == 0) && 
		(fabs(muls.slices*muls.sliceThickness-muls.c/muls.cellDiv) < 1e-5));

	initMuls();  

	// Fit the resolution to the wave function array, if not specified different
	if (muls.resolutionX == 0.0)
		muls.resolutionX = muls.ax / (double)muls.nx;
	if (muls.resolutionY == 0.0)
		muls.resolutionY = muls.by / (double)muls.ny;



	// Optional parameters:
	// determine whether potential periodic or not, etc.:
        configReader->ReadPeriodicParameters(muls.periodicXY, muls.periodicZ);
	if ((muls.periodicZ) && (muls.cellDiv > 1)) {
		printf("****************************************************************\n"
			"* Warning: cannot use cell divisions >1 and Z-periodic potential\n"
			"* periodicZ = NO\n"
			"****************************************************************\n");
		muls.periodicZ = false;
	}

        configReader->ReadBandLimitTrans(muls.bandlimittrans);
        configReader->ReadLoadPotential(muls.readPotential);
        configReader->ReadPotentialOutputParameters(muls.savePotential, muls.saveTotalPotential,
                                                    muls.plotPotential);
        configReader->ReadPotentialCalculationParameters(muls.fftpotential, muls.potential3D);

        configReader->ReadAverageParameters(muls.avgRuns, muls.storeSeries);

	if (!muls.tds) muls.avgRuns = 1;

	muls.scanXStart = muls.ax/2.0;
	muls.scanYStart = muls.by/2.0;
	muls.scanXN = 1;
	muls.scanYN = 1;
	muls.scanXStop = muls.scanXStart;
	muls.scanYStop = muls.scanYStart;


	switch (muls.mode) {
		/////////////////////////////////////////////////////////
		// read the position for doing CBED: 
	case CBED:
          configReader->ReadScanParameters(muls.scanXStart, muls.scanYStart);
          muls.scanXStop = muls.scanXStart;
          muls.scanYStop = muls.scanYStart;
          break;

		/////////////////////////////////////////////////////////
		// Read STEM scanning parameters 

	case STEM:
          configReader->ReadScanParameters(muls.scanXStart, muls.scanXStop, muls.scanXN,
                                           muls.scanYStart, muls.scanYStop, muls.scanYN);

		if (muls.scanXN < 1) muls.scanXN = 1;
		if (muls.scanYN < 1) muls.scanYN = 1;
		// if (muls.scanXStart > muls.scanXStop) muls.scanXN = 1;

		muls.displayProgInterval = muls.scanYN*muls.scanYN;
	}
	// printf("Potential progress interval: %d\n",muls.displayPotCalcInterval);

	// Read STEM/CBED probe parameters 
	muls.Cc = 0.0;
        configReader->ReadDoseParameters(muls.beamCurrent, muls.dwellTime);
        configReader->ReadProbeParameters(muls.dE_E, muls.dI_I, muls.dV_V, muls.alpha, muls.aAIS,
                                          muls.sourceRadius, muls.ismoth, muls.gaussScale, muls.gaussFlag);

	muls.electronScale = muls.beamCurrent*muls.dwellTime*MILLISEC_PICOAMP;
	//////////////////////////////////////////////////////////////////////

	// memorize dE_E0, and fill the array of well defined energy deviations
	dE_E0 = sqrt(muls.dE_E*muls.dE_E+
		muls.dI_I*muls.dI_I+
		muls.dV_V*muls.dV_V);
	muls.dE_EArray = (double *)malloc((muls.avgRuns+1)*sizeof(double));
	muls.dE_EArray[0] = 0.0;

	//**********************************************************
	// quick little fix to calculate gaussian energy distribution
	// without using statistics (better for only few runs)
	
	if (muls.printLevel > 0) printf("avgRuns: %d\n",muls.avgRuns);
	// serious bug in Visual C - dy comes out enormous.
	//dy = sqrt((double)pi)/((double)2.0*(double)(muls.avgRuns));
	// using precalculated sqrt(pi):
	dy = 1.772453850905/((double)2.0*(double)(muls.avgRuns));
	dx = pi/((double)(muls.avgRuns+1)*20);
	for (ix=1,x=0,y=0;ix<muls.avgRuns;x+=dx) {
		y += exp(-x*x)*dx;
		if (y>=ix*dy) {
			muls.dE_EArray[ix++] = x*2*dE_E0/pi;
			if (muls.printLevel > 2) printf("dE[%d]: %g eV\n",ix,muls.dE_EArray[ix-1]*muls.v0*1e3);
			if (ix < muls.avgRuns) {
				muls.dE_EArray[ix] = -muls.dE_EArray[ix-1];
				ix ++;
				if (muls.printLevel > 2) printf("dE[%d]: %g eV\n",ix,muls.dE_EArray[ix-1]*muls.v0*1e3);
			}
		}
	}

        configReader->ReadAberrationAmplitudes(muls.Cs, muls.C5, muls.Cc,
                                               muls.df0, muls.Scherzer, muls.astigMag,
                                               muls.a33, muls.a31,
                                               muls.a44, muls.a42,
                                               muls.a55, muls.a53, muls.a51,
                                               muls.a66, muls.a64, muls.a62);
        configReader->ReadAberrationAngles(muls.astigAngle,
                                           muls.phi33, muls.phi31,
                                           muls.phi44, muls.phi42,
                                           muls.phi55, muls.phi53, muls.phi51,
                                           muls.phi66, muls.phi64, muls.phi62);

        // Convert magnitudes from mm into A
        muls.Cs*=1e7;
        muls.C5*=1e7;
        muls.Cc*=1e7;
        
        switch(muls.Scherzer)
          {
          case 1:
            muls.df0 = -(float)sqrt(1.5*muls.Cs*(wavelength(muls.v0)));
            break;
          case 2:
            muls.df0 = -(float)sqrt(muls.Cs*(wavelength(muls.v0)));
            break;
          default:
            // convert defocus provided in nm into Angstrom
            muls.df0*=10;
            break;
          }

	// convert to A from nm:
	muls.astigMag *= 10.0;
	// convert astigAngle from deg to rad:
	muls.astigAngle *= pi/180.0;

	muls.phi33 /= (float)RAD2DEG;
	muls.phi31 /= (float)RAD2DEG;
	muls.phi44 /= (float)RAD2DEG;
	muls.phi42 /= (float)RAD2DEG;
	muls.phi55 /= (float)RAD2DEG;
	muls.phi53 /= (float)RAD2DEG;
	muls.phi51 /= (float)RAD2DEG;
	muls.phi66 /= (float)RAD2DEG;
	muls.phi64 /= (float)RAD2DEG;
	muls.phi62 /= (float)RAD2DEG;



	//**********************************************************************
	// Parameters for image display and directories, etc.

        configReader->ReadOutputName(muls.folder);

	//************************************************************************
        // read the different detector configurations                           
	resetParamFile();
	muls.detectorNum = 0;

	if (muls.mode == STEM) 
	{
          // instantiating the detector manager creates the detectors
          muls.detectors = DetectorMgrPtr(new DetectorManager(configReader));
        }
        
        // ************************************************************************
        // Tomography Parameters:
	if (muls.mode == TOMO) {     
          configReader->ReadTomoParameters(muls.tomoTilt, muls.tomoStart, muls.tomoStep, muls.tomoCount, 
                                           muls.zoomFactor);
          if ((muls.tomoStep == 0) && (muls.tomoStep > 1))
            muls.tomoStep = -2.0*muls.tomoStart/(double)(muls.tomoCount - 1);
	}


	// *******************************************************************
	// Read in parameters related to the calculation of the projected
	// Potential
        configReader->ReadAtomRadius(muls.atomRadius);
        configReader->ReadStructureFactorType(muls.scatFactor);


	// ***************************************************************
	// We now need to determine the size of the potential array, 
	// and the offset from its edges in A.  We only need to calculate
	// as much potential as we'll be illuminating later with the 
	// electron beam.

	if ((muls.mode == STEM) || (muls.mode == CBED)) {
		// we are assuming that there is enough atomic position data: 
		muls.potOffsetX = muls.scanXStart - 0.5*muls.nx*muls.resolutionX;
		muls.potOffsetY = muls.scanYStart - 0.5*muls.ny*muls.resolutionY;
                // adjust scanStop so that it coincides with a full pixel: 
		muls.potNx = (int)((muls.scanXStop-muls.scanXStart)/muls.resolutionX);
		muls.potNy = (int)((muls.scanYStop-muls.scanYStart)/muls.resolutionY);
		muls.scanXStop = muls.scanXStart+muls.resolutionX*muls.potNx;
		muls.scanYStop = muls.scanYStart+muls.resolutionY*muls.potNy;
		muls.potNx+=muls.nx;
		muls.potNy+=muls.ny;
		muls.potSizeX = muls.potNx*muls.resolutionX;
		muls.potSizeY = muls.potNy*muls.resolutionY;
	}
	else {
		muls.potNx = muls.nx;
		muls.potNy = muls.ny;
		muls.potSizeX = muls.potNx*muls.resolutionX;
		muls.potSizeY = muls.potNy*muls.resolutionY;
		muls.potOffsetX = muls.scanXStart - 0.5*muls.potSizeX;
		muls.potOffsetY = muls.scanYStart - 0.5*muls.potSizeY;
	}  
	// **************************************************************
	// Check to see if the given scan parameters really fit in cell 
	// dimensions:
	if ((muls.scanXN <=0) ||(muls.scanYN <=0)) {
		printf("The number of scan pixels must be >=1\n");
		exit(0);
	}
	if ((muls.scanXStart<0) || (muls.scanYStart<0) ||
		(muls.scanXStop<0) || (muls.scanYStop<0) ||
		(muls.scanXStart>muls.ax) || (muls.scanYStart>muls.by) ||
		(muls.scanXStop>muls.ax) || (muls.scanYStop>muls.by)) {
			printf("Scanning window is outside model dimensions (%g,%g .. %g,%g) [ax = %g, by = %g]!\n",muls.scanXStart,muls.scanYStart,muls.scanXStop,muls.scanYStop,muls.ax,muls.by);
			exit(0);
	}
	// *************************************************************
	// read in the beams we want to plot in the pendeloesung plot
	// Only possible if not in STEM or CBED mode 
	muls.lbeams = 0;   // flag for beam output 
        muls.nbout = 0;    // number of beams 
	resetParamFile();
	if ((muls.mode != STEM) && (muls.mode != CBED)) {
          configReader->ReadPendelloesungParameters(muls.hbeams, muls.kbeams, muls.lbeams, muls.nbout);
          muls.hbeams[i] *= muls.nCellX;
          muls.kbeams[i] *= muls.nCellY;

          muls.hbeams[i] = (muls.hbeams[i]+muls.nx) % muls.nx;
          muls.kbeams[i] = (muls.kbeams[i]+muls.ny) % muls.ny;
          for (size_t i=0; i<muls.hbeams.size(); i++)
            {
              printf("beam %d [%d %d]\n",i,muls.hbeams[i],muls.kbeams[i]);
            }
	}

	// TODO: possible breakage here - MCS 2013/04 - made muls.cfgFile be allocated on the struct
	//       at runtim - thus this null check doesn't make sense anymore.  Change cfgFile set
	//   Old comment:
	//	if cfgFile != NULL, the program will later write a the atomic config to this file 
	//muls.cfgFile = NULL;
	if (readparam("CFG-file:",buf,1)) 
	{
		sscanf(buf,"%s",muls.cfgFile);
	}

	// allocate memory for wave function 

	potDimensions[0] = muls.potNx;
	potDimensions[1] = muls.potNy;
	muls.trans = complex3D(muls.slices,muls.potNx,muls.potNy,"trans");
#if FLOAT_PRECISION == 1
	// printf("allocated trans %d %d %d\n",muls.slices,muls.potNx,muls.potNy);
	muls.fftPlanPotForw = fftwf_plan_many_dft(2,potDimensions, muls.slices,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy, FFTW_FORWARD, fftMeasureFlag);
	muls.fftPlanPotInv = fftwf_plan_many_dft(2,potDimensions, muls.slices,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy, FFTW_BACKWARD, fftMeasureFlag);
#else
	muls.fftPlanPotForw = fftw_plan_many_dft(2,potDimensions, muls.slices,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy, FFTW_FORWARD, fftMeasureFlag);
	muls.fftPlanPotInv = fftw_plan_many_dft(2,potDimensions, muls.slices,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy,muls.trans[0][0], NULL,
		1, muls.potNx*muls.potNy, FFTW_BACKWARD, fftMeasureFlag);
#endif

	////////////////////////////////////
	if (muls.printLevel >= 4) 
		printf("Memory for transmission function (%d x %d x %d) allocated and plans initiated\n",muls.slices,muls.potNx,muls.potNy);


	// printf("%d %d %d %d\n",muls.nx,muls.ny,sizeof(complex_tt),(int)(&muls.wave[2][2])-(int)(&muls.wave[2][1]));


} // end of readFile()
*/


