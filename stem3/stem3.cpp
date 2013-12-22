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

#include "memory_fftw3.hpp"	/* memory allocation routines */
#include "readparams.hpp"
#include "imagelib_fftw3.hpp"
#include "fileio_fftw3.hpp"
#include "matrixlib.hpp"
#include "stemlib.hpp"
#include "stemutil.hpp"
// #include "weblib.hpp"
#include "customslice.hpp"
#include "data_containers.hpp"

#include "readparams.hpp"

#include "config_readers.hpp"

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

void makeAnotation(float_tt **pict,int nx,int ny,char *text);
void initMuls();
void writeIntPix(char *outFile,float_tt **pict,int nx,int ny);
void runMuls(int lstart);
void saveLineScan(int run);
void readBeams(FILE *fp);
void doCBED(WavePtr initialWave, PotPtr pot);
void doSTEM(WavePtr initialWave, PotPtr pot);
void doTEM(WavePtr initialWave, PotPtr pot);
void doMSCBED(WavePtr initialWave, PotPtr pot);
void doTOMO(WavePtr initialWave, PotPtr pot);
void readFile(ConfigReaderPtr &configReader);
void displayParams();

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
  muls.nCellX = 1; muls.nCellY = 1; muls.nCellZ = 1; 

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

  displayParams();
#ifdef _OPENMP
  omp_set_dynamic(1);
#endif

        int mode;
        configReader->ReadMode(mode);
	switch (mode) {
        case CBED:   doCBED(initialWave, potential);   break;	
        case STEM:   doSTEM(initialWave, potential);   break;
        case TEM:    doTEM(initialWave, potential);    break;
        case MSCBED: doMSCBED(initialWave, potential); break;
        case TOMO:   doTOMO(initialWave, potential);   break;
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
void displayProgress(int flag) {
	// static double timer;
	static double timeAvg = 0;
	static double intensityAvg = 0;
	static time_t time0,time1;
	double curTime;
	int jz;

	if (flag < 0) {
		time(&time0);
		// timer = cputim();
		return;
	}
	time(&time1);  
	curTime = difftime(time1,time0);
	/*   curTime = cputim()-timer;
	if (curTime < 0) {
	printf("timer: %g, curr. time: %g, diff: %g\n",timer,cputim(),curTime);
	}
	*/
	if (muls.printLevel > 0) {

		if (muls.tds) {
			timeAvg = ((muls.avgCount)*timeAvg+curTime)/(muls.avgCount+1);
			intensityAvg = ((muls.avgCount)*intensityAvg+muls.intIntensity)/(muls.avgCount+1);
			printf("\n********************** run %3d ************************\n",muls.avgCount+1);
			// if (muls.avgCount < 1) {
			printf("* <u>: %3d |",muls.Znums[0]);
			for (jz=1;jz<muls.atomKinds;jz++) printf(" %8d |",muls.Znums[jz]);  
			printf(" intensity | time(sec) |    chi^2  |\n");
			// }
			/*
			printf("* %9g | %9g | %9g \n",muls.u2,muls.intIntensity,curTime);  
			}
			else {
			*/
			printf("*");
			for (jz=0;jz<muls.atomKinds;jz++) printf(" %8f |",(float)(muls.u2[jz]));  
			printf(" %9f | %9f | %9f |\n",muls.intIntensity,curTime,muls.avgCount > 0 ? muls.chisq[muls.avgCount-1] : 0);
			printf("*");
			for (jz=0;jz<muls.atomKinds;jz++) printf(" %8f |",(float)(muls.u2avg[jz]));  
			printf(" %9f | %9f \n",intensityAvg,timeAvg);
		}
		else {
			printf("\n**************** finished after %.1f sec ******************\n",curTime);
		}
	}  // end of printLevel check.

	time(&time0);
	//  timer = cputim();

}

void displayParams() {
	FILE *fpDir;
	char systStr[64];
	double k2max,temp;
	int i,j;
	static char Date[16],Time[16];
	time_t caltime;
	struct tm *mytime;
	const double pi=3.1415926535897;

	if (muls.printLevel < 1) {
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

	
	printf("* Beam tilt:            x=%g deg, y=%g deg (tilt back == %s)\n",muls.btiltx*RAD2DEG,muls.btilty*RAD2DEG,
		(muls.tiltBack == 1 ? "on" : "off"));
        printf("* Super cell divisions: %d (in z direction) %s\n",muls.cellDiv,muls.equalDivs ? "equal" : "non-equal");
	printf("* Slices per division:  %d (%gA thick slices [%scentered])\n",
		muls.slices,muls.sliceThickness,(muls.centerSlices) ? "" : "not ");
	printf("* Output every:         %d slices\n",muls.outputInterval);

	
	printf("* Beams:                %d x %d \n",muls.nx,muls.ny);  

	printf("* Aperture half angle:  %g mrad\n",muls.alpha);
	printf("* AIS aperture:         ");
	if (muls.aAIS > 0) printf("%g A\n",muls.aAIS);
	else printf("none\n");
	printf("* beam current:         %g pA\n",muls.beamCurrent);
	printf("* dwell time:           %g msec (%g electrons)\n",
		muls.dwellTime,muls.electronScale);

	printf("* Damping dE/E: %g / %g \n",sqrt(muls.dE_E*muls.dE_E+muls.dV_V*muls.dV_V+muls.dI_I*muls.dI_I)*muls.v0*1e3,muls.v0*1e3);

	/* 
	if (muls.ismoth) printf("Type 1 (=smooth aperture), ");
	if (muls.gaussFlag) printf("will apply gaussian smoothing"); 
	printf("\n");
	*/

	/***************************************************/
	/*  printf("Optimizing fftw plans according to probe array (%d x %dpixels = %g x %gA) ...\n",
	muls.nx,muls.ny,muls.nx*muls.resolutionX,muls.ny*muls.resolutionY);
	*/
	k2max = muls.nx/(2.0*muls.potSizeX);
	temp = muls.ny/(2.0*muls.potSizeY);
	if( temp < k2max ) k2max = temp;
	k2max = (BW * k2max);

	printf("* Real space res.:      %gA (=%gmrad)\n",
		1.0/k2max,wavelength(muls.v0)*k2max*1000.0);
	printf("* Reciprocal space res: dkx=%g, dky=%g\n",
		1.0/(muls.nx*muls.resolutionX),1.0/(muls.ny*muls.resolutionY));
	if (muls.mode == STEM) {
		printf("*\n"
			"* STEM parameters:\n");
		printf("* Maximum scattering angle:  %.0f mrad\n",
			0.5*2.0/3.0*wavelength(muls.v0)/muls.resolutionX*1000);    
                muls.detectors->PrintDetectors();
		
		printf("* Scan window:          (%g,%g) to (%g,%g)A, %d x %d = %d pixels\n",
			muls.scanXStart,muls.scanYStart,muls.scanXStop,muls.scanYStop,
			muls.scanXN,muls.scanYN,muls.scanXN*muls.scanYN);
        } /* end of if mode == STEM */

	/***********************************************************************
	* TOMOGRAPHY Mode
	**********************************************************************/
        if (muls.mode == TOMO) {
		printf("*\n"
			"* TOMO parameters:\n");
		printf("* Starting angle:       %g mrad (%g deg)\n",
			muls.tomoStart,muls.tomoStart*0.18/pi);
		printf("* Angular increase:     %g mrad (%g deg)\n",
			muls.tomoStep,muls.tomoStep*0.180/pi);
		printf("* Number of dp's:       %d\n",muls.tomoCount);
		printf("* Zoom factor:          %g\n",muls.zoomFactor);
	}

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


/************************************************************************
* doTOMO performs a Diffraction Tomography simulation
*
* This routine creates a script file which can then be run by a separate 
* command.
* For now this routine will only allow tomography about the y-axis.
* To do anything else one can start with a previously rotated super-cell.
*
* Important parameters: tomoStart, tomoStep, tomoCount, zoomFactor
***********************************************************************/
void doTOMO(WavePtr initialWave, PotPtr pot) {
	double boxXmin=0,boxXmax=0,boxYmin=0,boxYmax=0,boxZmin=0,boxZmax=0;
	double mAx,mBy,mCz;
	int ix,iy,iz,iTheta,i;
	double u[3],**Mm = NULL;
	double theta = 0;
	atom *atoms = NULL;
	char cfgFile[64],stemFile[128],scriptFile[64],diffAnimFile[64];
	FILE *fpScript,*fpDiffAnim;


	Mm = muls.Mm;
	atoms = (atom *)malloc(muls.natom*sizeof(atom));

	boxXmin = boxXmax = muls.ax/2.0;
	boxYmin = boxYmax = muls.by/2.0;
	boxZmin = boxZmax = muls.c/2.0;

	// For all tomography tilt angles:
	for (iTheta = 0;iTheta < muls.tomoCount;iTheta++) {
		theta = muls.tomoStart + iTheta*muls.tomoStep; // theta in mrad
		// Try different corners of the box, and see, how far they poke out.
		for (ix=-1;ix<=1;ix++) for (iy=-1;iy<=1;iy++) for (iz=-1;iz<=1;iz++) {
			// Make center of unit cell rotation center
			u[0]=ix*muls.ax/2; u[1]=iy*muls.by/2.0; u[2]=iz*muls.c/2.0;

			// rotate about y-axis
			rotateVect(u,u,0,theta*1e-3,0);

			// shift origin back to old (0,0,0):
			u[0]+=muls.ax/2; u[1]+=muls.by/2.0; u[2]+=muls.c/2.0;

			boxXmin = boxXmin>u[0] ? u[0] : boxXmin; boxXmax = boxXmax<u[0] ? u[0] : boxXmax; 
			boxYmin = boxYmin>u[1] ? u[1] : boxYmin; boxYmax = boxYmax<u[1] ? u[1] : boxYmax; 
			boxZmin = boxZmin>u[2] ? u[2] : boxZmin; boxZmax = boxZmax<u[2] ? u[2] : boxZmax; 

		}
	} /* for iTheta ... */

	// find max. box size:
	boxXmax -= boxXmin;
	boxYmax -= boxYmin;
	boxZmax -= boxZmin;
	printf("Minimum box size for tomography tilt series: %g x %g x %gA, zoom Factor: %g\n",
		boxXmax,boxYmax,boxZmax,muls.zoomFactor);
	boxXmax /= muls.zoomFactor;
	boxYmax = boxXmax*muls.by/muls.ax;

	// boxMin will now be boxCenter:
	boxXmin = 0.5*boxXmax;
	boxYmin = 0.5*boxYmax;
	boxZmin = 0.5*boxZmax;

	// We have to save the original unit cell dimensions
	mAx = muls.ax; mBy = muls.by; mCz = muls.c;
	muls.ax=boxXmax; muls.by=boxYmax; muls.c=boxZmax;

	printf("Will use box sizes: %g x %g x %gA (kept original aspect ratio). \n"
		"Writing structure files now, please wait ...\n",
		boxXmax,boxYmax,boxZmax);

	// open the script file and write the stem instructions in there
	sprintf(scriptFile,"%s/run_tomo",muls.folder.c_str());
	if ((fpScript=fopen(scriptFile,"w")) == NULL) {
		printf("doTOMO: unable to open scriptFile %s for writing\n",scriptFile);
		exit(0);
	}
	fprintf(fpScript,"#!/bin/bash\n\n");

	// open the diffraction animation file 
	sprintf(diffAnimFile,"%s/diff_anim",muls.folder.c_str());
	if ((fpDiffAnim=fopen(diffAnimFile,"w")) == NULL) {
		printf("doTOMO: unable to open diffraction animation %s for writing\n",diffAnimFile);
		exit(0);
	}
	fprintf(fpDiffAnim,"#!/bin/bash\n\n");


	// We will now rotate all the unit cell and write the answer to
	// separate structure output files.
	for (iTheta = 0;iTheta < muls.tomoCount;iTheta++) {
		theta = muls.tomoStart + iTheta*muls.tomoStep; // theta in mrad
		muls.tomoTilt = theta;

		// rotate the structure and write result to local atom array
		for(i=0;i<(muls.natom);i++) {	
			u[0] = muls.atoms[i].x - mAx/2.0; 
			u[1] = muls.atoms[i].y - mBy/2.0; 
			u[2] = muls.atoms[i].z - mCz/2.0; 
			rotateVect(u,u,0,theta*1e-3,0);
			atoms[i].x = u[0]+boxXmin;
			atoms[i].y = u[1]+boxYmin; 
			atoms[i].z = u[2]+boxZmin; 
			atoms[i].Znum = muls.atoms[i].Znum;
			atoms[i].occ = muls.atoms[i].occ;
			atoms[i].dw = muls.atoms[i].dw;
		}

		sprintf(cfgFile,"%s/tomo_%dmrad.cfg",muls.folder.c_str(),(int)theta);
		sprintf(stemFile,"%s/tomo_%dmrad.dat",muls.folder.c_str(),(int)theta);
		printf("Writing file %s | ",cfgFile);
		writeCFG(atoms,muls.natom,cfgFile,&muls);
		sprintf(cfgFile,"tomo_%dmrad.cfg",(int)theta);
		writeSTEMinput(stemFile,cfgFile,&muls);

		// add to script files:
		fprintf(fpScript,"stem tomo_%dmrad.dat\n",(int)theta);
	}
	muls.ax = mAx; muls.by = mBy; muls.c = mCz;
	sprintf(stemFile,"copy fparams.dat %s/",muls.folder.c_str());
	system(stemFile);

	// close script files again
	fprintf(fpDiffAnim,"convert -delay 20 diff*.jpg diff.gif\n");  
	fclose(fpScript);
	fclose(fpDiffAnim);
	sprintf(stemFile,"chmod +x %s",scriptFile);
	system(stemFile);
	sprintf(stemFile,"chmod +x %s",diffAnimFile);
	system(stemFile);

	exit(0);
}

/************************************************************************
* doMSCBED performs a super-CBED calculation based on the MS algorithm
* including phonons.
*
***********************************************************************/
void doMSCBED(WavePtr initialWave, PotPtr pot) {


}

/************************************************************************
* doCBED performs a CBED calculation
*
***********************************************************************/

void doCBED(WavePtr initialWave, PotPtr pot) {
  int ix,iy,i,pCount,result;
  FILE *avgFp,*fp,*fpPos=0;
  double timer,timerTot;
  double probeCenterX,probeCenterY,probeOffsetX,probeOffsetY;
  char buf[BUF_LEN];
  float_tt t=0;
  float_tt **avgPendelloesung = NULL;
  int oldMulsRepeat1 = 1;
  int oldMulsRepeat2 = 1;
  long iseed=0;
  WavePtr wave = WavePtr(new WAVEFUNC(muls.nx,muls.ny, muls.resolutionX, muls.resolutionY, muls.input_ext, muls.output_ext));
  ImageIOPtr imageIO = ImageIOPtr(new CImageIO(muls.nx, muls.ny, muls.input_ext, muls.output_ext));
  std::map<std::string, double> params;
  params["dx"]=muls.resolutionX;
  params["dy"]=muls.resolutionY;
  std::string comment;

  std::vector<unsigned> position(1);         // Used to indicate the number of averages

  muls.chisq = std::vector<double>(muls.avgRuns);

  if (iseed == 0) iseed = -(long) time( NULL );

  if (muls.lbeams) {
    muls.pendelloesung = NULL;
    if (avgPendelloesung == NULL) {
      avgPendelloesung = float2D(muls.nbout,
                                 muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv,
                                 "pendelloesung");
    }    
  }
  probeCenterX = muls.scanXStart;
  probeCenterY = muls.scanYStart;

  timerTot = 0; /* cputim();*/
  displayProgress(-1);

  for (muls.avgCount = 0;muls.avgCount < muls.avgRuns;muls.avgCount++) {
    muls.totalSliceCount = 0;
    pCount = 0;
    /* make sure we start at the beginning of the file 
       so we won't miss any line that contains a sequence,
       because we will not do any EOF wrapping
    */
    resetParamFile();
    
    /* probe(&muls,xpos,ypos); */
    /* make incident probe wave function with probe exactly in the center */
    /* if the potential array is not big enough, the probe can 
     * then also be adjusted, so that it is off-center
     */

    probeOffsetX = muls.sourceRadius*gasdev(&iseed)*SQRT_2;
    probeOffsetY = muls.sourceRadius*gasdev(&iseed)*SQRT_2;
    muls.scanXStart = probeCenterX+probeOffsetX;
    muls.scanYStart = probeCenterY+probeOffsetY;
    probe(&muls, wave,muls.scanXStart-muls.potOffsetX,muls.scanYStart-muls.potOffsetY);
    if (muls.saveLevel > 2) {
      wave->WriteProbe();
    } 	
    // printf("Probe: (%g, %g)\n",muls.scanXStart,muls.scanYStart);
    /*****************************************************************
     * For debugging only!!!
     *
		muls->WriteWave("probe.img")
    *****************************************************************/


    if (muls.sourceRadius > 0) {
      if (muls.avgCount == 0) fpPos = fopen("probepos.dat","w");
      else fpPos = fopen("probepos.dat","a");
      if (fpPos == NULL) {
        printf("Was unable to open file probepos.dat for writing\n");
      }
      else {
        fprintf(fpPos,"%g %g\n",muls.scanXStart,muls.scanYStart);
        fclose(fpPos);
      }
    }

    if ((muls.showProbe) && (muls.avgCount == 0)) {
#ifndef WIN32
      //probePlot(&muls);
      sprintf(buf,"ee %s/probePlot_0.jpg &",muls.folder.c_str());
      system(buf);
#endif
    }
    //muls.nslic0 = 0;

    result = readparam("sequence: ",buf,0);
    while (result) {
      if (((buf[0] < 'a') || (buf[0] > 'z')) && 
          ((buf[0] < '1') || (buf[0] > '9')) &&
          ((buf[0] < 'A') || (buf[0] > 'Z'))) {
        // printf("Stacking sequence: %s\n",buf);
        printf("Can only work with old stacking sequence\n");
        break;
      }
      muls.mulsRepeat1 = 1;
      muls.mulsRepeat2 = 1;
      sscanf(buf,"%d %d",&muls.mulsRepeat1,&muls.mulsRepeat2);
      for (i=0;i<(int)strlen(buf);i++) buf[i] = 0;
      if (muls.mulsRepeat2 < 1) muls.mulsRepeat2 = 1;
      sprintf(muls.cin2,"%d",muls.mulsRepeat1);

      
      /***********************************************************
       * make sure we have enough memory for the pendelloesung plot
       */
      if ((muls.lbeams) && 
          ((oldMulsRepeat1 !=muls.mulsRepeat1) ||
           (oldMulsRepeat2 !=muls.mulsRepeat2))) {
        oldMulsRepeat1 = muls.mulsRepeat1;
        oldMulsRepeat2 = muls.mulsRepeat2;
        if (muls.pendelloesung != NULL)
          free(muls.pendelloesung[0]);
        free(avgPendelloesung[0]);
        muls.pendelloesung = NULL;
        avgPendelloesung = float2D(muls.nbout,
                                   muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv,
                                   "pendelloesung");
      }
      /*********************************************************/
      
      // printf("Stacking sequence: %s\n",buf);

      muls.saveFlag = 0;
      /****************************************
       * do the (small) loop
       *****************************************/
      for (pCount = 0;pCount<muls.mulsRepeat2*muls.cellDiv;pCount++) {
        
        pot->Refresh();
        
        timer = cputim();
        // what probe should runMulsSTEM use here?
        runMulsSTEM(&muls, wave, pot); 
        
        printf("Thickness: %gA, int.=%g, time: %gsec\n",
               wave->thickness,wave->intIntensity,cputim()-timer);

        /***************** Only if Save level > 2: ****************/
        if ((muls.avgCount == 0) && (muls.saveLevel > 2)) {
          wave->WriteWave();
        } 	
#ifdef VIB_IMAGE_TEST_CBED
        wave->WriteWave()
#endif 
          muls.totalSliceCount += muls.slices;
        
      } // end of for pCount = 0... 
      result = readparam("sequence: ",buf,0);
    }
    /*    printf("Total CPU time = %f sec.\n", cputim()-timerTot ); */
    
    // TODO: Why are we reading in a DP at this point?  Do we have one yet?  
    //     What happens if it isn't there?
    wave->ReadDiffPat();
    
    if (muls.avgCount == 0) {
      memcpy((void *)wave->avgArray[0],(void *)wave->diffpat[0],
             (size_t)(muls.nx*muls.ny*sizeof(float_tt)));
      /* move the averaged (raw data) file to the target directory as well */
      // TODO: make sure that DP average gets created properly
      //sprintf(avgName,"%s/diffAvg_%d.img",muls.folder.c_str(),muls.avgCount+1);
      //sprintf(systStr,"mv %s/diff.img %s",muls.folder.c_str(),avgName);
      //system(systStr);
      if (muls.lbeams) {
        for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
          for (ix=0;ix<muls.nbout;ix++) {
            avgPendelloesung[ix][iy] = muls.pendelloesung[ix][iy];
          }
        }
      }
    } // of if muls.avgCount == 0 ...
    else {
      muls.chisq[muls.avgCount-1] = 0.0;
      for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
          t = ((float_tt)muls.avgCount*wave->avgArray[ix][iy]+
               wave->diffpat[ix][iy])/((float_tt)(muls.avgCount+1));
          muls.chisq[muls.avgCount-1] += (wave->avgArray[ix][iy]-t)*(wave->avgArray[ix][iy]-t);
          wave->avgArray[ix][iy] = t;

        }
      muls.chisq[muls.avgCount-1] = muls.chisq[muls.avgCount-1]/(double)(muls.nx*muls.ny);
      params["Tilt"] = muls.tomoTilt;
      params["1/Wavelength"] = 1.0/wavelength(muls.v0);
      wave->WriteDiffPat("Averaged Diffraction pattern, unit: 1/A", params);
                        
      muls.storeSeries = 1;
      if (muls.saveLevel == 0)	muls.storeSeries = 0;
      else if (muls.avgCount % muls.saveLevel != 0) muls.storeSeries = 0;

      if (muls.storeSeries) 
        wave->WriteAvgArray(muls.avgCount+1, "Averaged Diffraction pattern, unit: 1/A", params);


      /* write the data to a file */
      // TODO: vestigial code?  does anyone use this?
      if (muls.saveFlag >-1) {
        char systStr[255];
        sprintf(systStr,"%s/avgresults.dat",muls.folder.c_str());
        if ((avgFp = fopen(systStr,"w")) == NULL )
          printf("Sorry, could not open data file for averaging\n");
        else {
          for (ix =0;ix<muls.avgCount;ix++) {
            fprintf(avgFp,"%d %g\n",ix+1,muls.chisq[ix]);
          }
          fclose(avgFp);
        }
      }
      /*************************************************************/

      /***********************************************************
       * Average over the pendelloesung plot as well
       */
      if (muls.lbeams) {
        for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
          for (ix=0;ix<muls.nbout;ix++) {
            avgPendelloesung[ix][iy] = 
              ((float_tt)muls.avgCount*avgPendelloesung[ix][iy]+
               muls.pendelloesung[ix][iy])/(float_tt)(muls.avgCount+1);
          }
        }
      }
    } /* else ... if avgCount was greater than 0 */
    
    if (muls.lbeams) {
      /**************************************************************
       * The diffraction spot intensities of the selected 
       * diffraction spots are now stored in the 2 dimensional array
       * muls.pendelloesung[beam][slice].
       * We can write the array to a file and display it, just for 
       * demonstration purposes
       *************************************************************/
      char systStr[255];
      sprintf(systStr,"%s/pendelloesung.dat",muls.folder.c_str());
      if ((fp=fopen(systStr,"w")) !=NULL) {
        printf("Writing Pendelloesung data\n");
        for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
          /* write the thicknes in the first column of the file */
          fprintf(fp,"%g",iy*muls.c/((float)(muls.slices*muls.cellDiv)));
          /* write the beam intensities in the following columns */
          for (ix=0;ix<muls.nbout;ix++) {
            fprintf(fp,"\t%g",avgPendelloesung[ix][iy]);
          }
          /* close the line, and start a new one for the next set of
           * intensities
           */
          fprintf(fp,"\n");
        }
        fclose(fp);
      }
      else {
        printf("Could not open file for pendelloesung plot\n");
      }  
    } /* end of if lbeams ... */
    displayProgress(1);
  } /* end of for muls.avgCount=0.. */
  //delete(wave);
}
/************************************************************************
* End of doCBED()
***********************************************************************/

/************************************************************************
* doTEM performs a TEM calculation (with through focus reconstruction)
*
***********************************************************************/

void doTEM(WavePtr initialWave, PotPtr pot) {
	const double pi=3.1415926535897;
	int ix,iy,i,pCount,result;
	FILE *avgFp,*fp; // *fpPos=0;
	double timer,timerTot;
	double x,y,ktx,kty;
	char buf[BUF_LEN];//,avgName[256],systStr[512];
	std::string comment;
	float_tt t;
	float_tt **avgPendelloesung = NULL;
	int oldMulsRepeat1 = 1;
	int oldMulsRepeat2 = 1;
	long iseed=0;
	std::map<std::string, double> params;
	WavePtr wave = WavePtr(new WAVEFUNC(muls.nx,muls.ny,muls.resolutionX,muls.resolutionY, muls.input_ext, muls.output_ext));
	fftwf_complex **imageWave = NULL;

	if (iseed == 0) iseed = -(long) time( NULL );

	muls.chisq=std::vector<double>(muls.avgRuns);

	if (muls.lbeams) {
          muls.pendelloesung = NULL;
          if (avgPendelloesung == NULL) {
            avgPendelloesung = float2D(muls.nbout,
                                       muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv,
                                       "pendelloesung");
          }	  
	}

	timerTot = 0; /* cputim();*/
	displayProgress(-1);
	for (muls.avgCount = 0;muls.avgCount < muls.avgRuns;muls.avgCount++) {
          muls.totalSliceCount = 0;
          
          pCount = 0;

          /* make incident probe wave function with probe exactly in the center */
          /* if the potential array is not big enough, the probe can 
           * then also be adjusted, so that it is off-center
           */

          // muls.scanXStart = muls.nx/2*muls.resolutionX+muls.sourceRadius*gasdev(&iseed)*sqrt(2);
          // muls.scanYStart = muls.ny/2*muls.resolutionY+muls.sourceRadius*gasdev(&iseed)*sqrt(2);
          // probe(&muls,muls.scanXStart,muls.scanYStart);

          //muls.nslic0 = 0;
          // produce an incident plane wave:
          if ((muls.btiltx == 0) && (muls.btilty == 0)) {
            for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
                wave->wave[ix][iy][0] = 1;	wave->wave[ix][iy][1] = 0;
              }
          }
          else {
            // produce a tilted wave function (btiltx,btilty):
            ktx = 2.0*pi*sin(muls.btiltx)/wavelength(muls.v0);
            kty = 2.0*pi*sin(muls.btilty)/wavelength(muls.v0);
            for (ix=0;ix<muls.nx;ix++) {
              x = muls.resolutionX*(ix-muls.nx/2);
              for (iy=0;iy<muls.ny;iy++) {
                y = muls.resolutionY*(ix-muls.nx/2);
                wave->wave[ix][iy][0] = (float)cos(ktx*x+kty*y);	
                wave->wave[ix][iy][1] = (float)sin(ktx*x+kty*y);
              }
            }
          }

          result = readparam("sequence: ",buf,0);
          while (result) {
            if (((buf[0] < 'a') || (buf[0] > 'z')) && 
                ((buf[0] < '1') || (buf[0] > '9')) &&
                ((buf[0] < 'A') || (buf[0] > 'Z'))) {
              // printf("Stacking sequence: %s\n",buf);
              printf("Can only work with old stacking sequence\n");
              break;
            }

            muls.mulsRepeat1 = 1;
            muls.mulsRepeat2 = 1;
            sscanf(buf,"%d %d",&muls.mulsRepeat1,&muls.mulsRepeat2);
            for (i=0;i<(int)strlen(buf);i++)
              buf[i] = 0;
            if (muls.mulsRepeat2 < 1) muls.mulsRepeat2 = 1;
            sprintf(muls.cin2,"%d",muls.mulsRepeat1);
            /***********************************************************
             * make sure we have enough memory for the pendelloesung plot
             */
            if ((muls.lbeams) && 
                ((oldMulsRepeat1 !=muls.mulsRepeat1) ||
                 (oldMulsRepeat2 !=muls.mulsRepeat2))) {
              oldMulsRepeat1 = muls.mulsRepeat1;
              oldMulsRepeat2 = muls.mulsRepeat2;
              if (muls.pendelloesung != NULL)  free(muls.pendelloesung[0]);
              free(avgPendelloesung[0]);
              muls.pendelloesung = NULL;
              avgPendelloesung = float2D(muls.nbout,
                                         muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv,
                                         "pendelloesung");
            }
            /*********************************************************/
            
            // printf("Stacking sequence: %s\n",buf);

            if (muls.equalDivs) {
              if (muls.printLevel > 1) printf("found equal unit cell divisions\n");
              pot->Refresh();
            }

            muls.saveFlag = 0;
            /****************************************
             * do the (small) loop through the slabs
             *****************************************/
            for (pCount = 0;pCount<muls.mulsRepeat2*muls.cellDiv;pCount++) {

              /*******************************************************
               * build the potential slices from atomic configuration
               ******************************************************/
              // if ((muls.tds) || (muls.nCellZ % muls.cellDiv != 0)) {
              if (!muls.equalDivs) {
                pot->Refresh();
              }

              timer = cputim();
              runMulsSTEM(&muls,wave, pot); 
              muls.totalSliceCount += muls.slices;

              if (muls.printLevel > 0) {
                printf("t=%gA, int.=%g time: %gsec (avgCount=%d)\n",
                       wave->thickness,wave->intIntensity,cputim()-timer,muls.avgCount);
              }

              /***************** FOR DEBUGGING ****************/		
              if ((muls.avgCount == 0) && (muls.saveLevel >=0) && (pCount+1==muls.mulsRepeat2*muls.cellDiv)) {
                if (muls.tds) comment = "Test wave function for run 0";
                else comment = "Exit face wave function for no TDS";
                if ((muls.tiltBack) && ((muls.btiltx != 0) || (muls.btilty != 0))) {
                  ktx = -2.0*pi*sin(muls.btiltx)/wavelength(muls.v0);
                  kty = -2.0*pi*sin(muls.btilty)/wavelength(muls.v0);
                  for (ix=0;ix<muls.nx;ix++) {
                    x = muls.resolutionX*(ix-muls.nx/2);
                    for (iy=0;iy<muls.ny;iy++) {
                      y = muls.resolutionY*(ix-muls.nx/2);
                      wave->wave[ix][iy][0] *= cos(ktx*x+kty*y);	
                      wave->wave[ix][iy][1] *= sin(ktx*x+kty*y);
                    }
                  }
                  if (muls.printLevel > 1) printf("** Applied beam tilt compensation **\n");
                }
                
                wave->WriteWave();
              }	
#ifdef VIB_IMAGE_TEST  // doTEM
              if ((muls.tds) && (muls.saveLevel > 2)) {
                params["HT"] = muls.v0;
                params["Cs"] = muls.Cs;
                params["Defocus"] = muls.df0;
                params["Astigmatism Magnitude"] = muls.astigMag;
                params["Astigmatism Angle"] = muls.astigAngle;
                params["Focal Spread"] = muls.Cc * sqrt(muls.dE_E*muls.dE_E+muls.dV_V*muls.dV_V+muls.dI_I*muls.dI_I);
                params["Convergence Angle"] = muls.alpha;
                params["Beam Tilt X"] = muls.btiltx;
                params["Beam Tilt Y"] = muls.btilty;
                comment = "complex exit face Wave function";

                wave->WriteWave(muls.avgCount, comment, params);
              }
#endif 

            } 
            result = readparam("sequence: ",buf,0);
          } 
          /////////////////////////////////////////////////////////////////////////////
          // finished propagating through whole sample, we're at the exit surface now.
          // This means the wave function is used for nothing else than producing image(s)
          // and diffraction patterns.
          //////////////////////////////////////////////////////////////////////////////

          wave->ReadDiffPat();

          if (muls.avgCount == 0) {
            /***********************************************************
             * Save the diffraction pattern
             **********************************************************/	
            memcpy((void *)wave->avgArray[0],(void *)wave->diffpat[0],(size_t)(muls.nx*muls.ny*sizeof(float_tt)));
            /* move the averaged (raw data) file to the target directory as well */
            wave->WriteAvgArray(muls.avgCount+1);

            /***********************************************************
             * Save the Pendelloesung Plot
             **********************************************************/	
            if (muls.lbeams) {
              for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
                for (ix=0;ix<muls.nbout;ix++) {
                  avgPendelloesung[ix][iy] = muls.pendelloesung[ix][iy];
                }
              }
            }
            /***********************************************************
             * Save the defocused image, we can do with the wave what 
             * we want, since it is not used after this anymore. 
             * We will therefore multiply with the transfer function for
             * all the different defoci, inverse FFT and save each image.
             * diffArray will be overwritten with the image.
             **********************************************************/ 
            if (imageWave == NULL) imageWave = complex2D(muls.nx,muls.ny,"imageWave");
            // multiply wave (in rec. space) with transfer function and write result to imagewave
            fftwf_execute(wave->fftPlanWaveForw);
            for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
                // here, we apply the CTF:
                imageWave[ix][iy][0] = wave->wave[ix][iy][0];
                imageWave[ix][iy][1] = wave->wave[ix][iy][1];
              }
            fftwf_execute_dft(wave->fftPlanWaveInv,imageWave[0],imageWave[0]);
            // get the amplitude squared:
            for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
                wave->diffpat[ix][iy] = imageWave[ix][iy][0]*imageWave[ix][iy][0]+imageWave[ix][iy][1]*imageWave[ix][iy][1];
              }
            wave->WriteWaveIntensity();
            // End of Image writing (if avgCount = 0)
            //////////////////////////////////////////////////////////////////////
            
          } // of if muls.avgCount == 0 ...
          else {
            /* 	 readRealImage_old(avgArray,muls.nx,muls.ny,&t,"diffAvg.img"); */
            muls.chisq[muls.avgCount-1] = 0.0;
            for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
                t = ((float_tt)muls.avgCount*wave->avgArray[ix][iy]+
                     wave->diffpat[ix][iy])/((float_tt)(muls.avgCount+1));
                muls.chisq[muls.avgCount-1] += (wave->avgArray[ix][iy]-t)*(wave->avgArray[ix][iy]-t);
                wave->avgArray[ix][iy] = t;
              }
            muls.chisq[muls.avgCount-1] = muls.chisq[muls.avgCount-1]/(double)(muls.nx*muls.ny);
            wave->WriteAvgArray(muls.avgCount+1);

            /* write the data to a file */
            if ((avgFp = fopen("avgresults.dat","w")) == NULL )
              printf("Sorry, could not open data file for averaging\n");
            else {
              for (ix =0;ix<muls.avgCount;ix++) {
                fprintf(avgFp,"%d %g\n",ix+1,muls.chisq[ix]);
              }
              fclose(avgFp);
            }
            /*************************************************************/
            
            /***********************************************************
             * Average over the pendelloesung plot as well
             */
            if (muls.lbeams) {
              for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
                for (ix=0;ix<muls.nbout;ix++) {
                  avgPendelloesung[ix][iy] = 
                    ((float_tt)muls.avgCount*avgPendelloesung[ix][iy]+
                     muls.pendelloesung[ix][iy])/(float_tt)(muls.avgCount+1);
                }
              }
            }
            /***********************************************************
             * Save the defocused image, we can do with the wave what 
             * we want, since it is not used after this anymore. 
             * We will therefore multiply with the transfer function for
             * all the different defoci, inverse FFT and save each image.
             * diffArray will be overwritten with the image.
             **********************************************************/ 
            if (imageWave == NULL) imageWave = complex2D(muls.nx,muls.ny,"imageWave");
            // multiply wave (in rec. space) with transfer function and write result to imagewave
#if FLOAT_PRECISION == 1
            fftwf_execute(wave->fftPlanWaveForw);
#elif FLOAT_PRECISION == 2
            fftw_execute(wave->fftPlanWaveForw);
#endif
            for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
                imageWave[ix][iy][0] = wave->wave[ix][iy][0];
                imageWave[ix][iy][1] = wave->wave[ix][iy][1];
              }
#if FLOAT_PRECISION == 1
            fftwf_execute_dft(wave->fftPlanWaveInv,imageWave[0],imageWave[0]);
#elif FLOAT_PRECISION == 2
            fftw_execute_dft(wave->fftPlanWaveInv,imageWave[0],imageWave[0]);
#endif

            // save the amplitude squared:
            wave->ReadImage();
            for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
                t = ((float_tt)muls.avgCount*wave->diffpat[ix][iy]+
                     imageWave[ix][iy][0]*imageWave[ix][iy][0]+imageWave[ix][iy][1]*imageWave[ix][iy][1])/(float_tt)(muls.avgCount+1);
                wave->diffpat[ix][iy] = t;
              }
            wave->WriteImage();
            // End of Image writing (if avgCount > 0)
            //////////////////////////////////////////////////////////////////////

          } /* else ... if avgCount was greater than 0 */


          /////////////////////////////////////////////////////
          // Save the Pendelloesung plot:
          if (muls.lbeams) {
            /**************************************************************
             * The diffraction spot intensities of the selected 
             * diffraction spots are now stored in the 2 dimensional array
             * muls.pendelloesung[beam][slice].
             * We can write the array to a file and display it, just for 
             * demonstration purposes
             *************************************************************/
            char avgName[255];
            sprintf(avgName,"%s/pendelloesung.dat",muls.folder.c_str());
            if ((fp=fopen(avgName,"w")) !=NULL) {
              printf("Writing Pendelloesung data\n");
              for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
                /* write the thicknes in the first column of the file */
                fprintf(fp,"%g",iy*muls.c/((float)(muls.slices*muls.cellDiv)));
                /* write the beam intensities in the following columns */
                for (ix=0;ix<muls.nbout;ix++) {
                  // store the AMPLITUDE:
                  fprintf(fp,"\t%g",sqrt(avgPendelloesung[ix][iy]/(muls.nx*muls.ny)));
                }
                /* close the line, and start a new one for the next set of
                 * intensities
                 */
                fprintf(fp,"\n");
              }
              fclose(fp);
            }
            else {
              printf("Could not open file for pendelloesung plot\n");
            }	
          } /* end of if lbeams ... */		 
          displayProgress(1);
	} /* end of for muls.avgCount=0.. */  
}
/************************************************************************
* end of doTEM
***********************************************************************/



/************************************************************************
* doSTEM performs a STEM calculation
*
***********************************************************************/

void doSTEM(WavePtr initialWave, PotPtr pot) {
	int ix=0,iy=0,i,pCount,picts,ixa,iya,totalRuns;
	double timer, total_time=0;
	char buf[BUF_LEN];
	float_tt t;
	static float_tt **avgArray=NULL;
	double collectedIntensity;

	std::vector<WavePtr> waves;
	WavePtr wave;

	//pre-allocate several waves (enough for one row of the scan.  
	for (int th=0; th<omp_get_max_threads(); th++)
	{
		waves.push_back(WavePtr(new WAVEFUNC(muls.nx, muls.ny, muls.resolutionX, muls.resolutionY, muls.input_ext, muls.output_ext)));
	}

	muls.chisq = std::vector<double>(muls.avgRuns);
	totalRuns = muls.avgRuns;
	timer = cputim();

	/* average over several runs of for TDS */
	displayProgress(-1);

	for (muls.avgCount = 0;muls.avgCount < totalRuns; muls.avgCount++) {
		total_time = 0;
		collectedIntensity = 0;
		muls.totalSliceCount = 0;
		muls.dE_E = muls.dE_EArray[muls.avgCount];


		/****************************************
		* do the (big) loop
		*****************************************/
		pCount = 0;
		/* make sure we start at the beginning of the file 
		so we won't miss any line that contains a sequence,
		because we will not do any EOF wrapping
		*/
		resetParamFile();
		while (readparam("sequence: ",buf,0)) {
			if (((buf[0] < 'a') || (buf[0] > 'z')) && 
				((buf[0] < '1') || (buf[0] > '9')) &&
				((buf[0] < 'A') || (buf[0] > 'Z'))) {
					printf("Stacking sequence: %s\n",buf);
					printf("Can only work with old stacking sequence\n");
					break;
			}

			// printf("Stacking sequence: %s\n",buf);

			picts = 0;
			/* for the dislocation models picts will be 1, because the atomcoordinates
			* are expressed explicitly for every atom in the whole specimen
			* For perfect Si samples mulsRepeat1 will remain 1, but picts will give
			* the number of unit cells in Z-direction, where the unit cell is defined by 
			* the unit cell in the cssr file multiplied by NCELLZ.  
			* cellDiv will usually be 1 in that case.
			*/
			sscanf(buf,"%d %d",&muls.mulsRepeat1,&picts);
			for (i=0;i<(int)strlen(buf);i++) buf[i] = 0;
			if (picts < 1) picts = 1;
			muls.mulsRepeat2 = picts;
			sprintf(muls.cin2,"%d",muls.mulsRepeat1);
			/* if the unit cell is divided into slabs, we need to multiply
			* picts by that number
			*/
			if ((picts > 1)&& (muls.cubex >0) && (muls.cubey >0) && (muls.cubez>0)) {
				printf("Warning: cube size of height %gA has been defined, ignoring sequence\n",muls.cubez);
				picts = 1;
			}
			picts *= muls.cellDiv;

			if (muls.equalDivs) {
                          pot->Refresh();
                          timer = cputim();
			}

			/****************************************
			* do the (small) loop over slabs
			*****************************************/
			for (pCount=0;pCount<picts;pCount++) {
				/*******************************************************
				* build the potential slices from atomic configuration
				******************************************************/
				if (!muls.equalDivs) {
                                  pot->Refresh();
                                  timer = cputim();
				}

				muls.complete_pixels=0;
				/**************************************************
				* scan through the different probe positions
				*************************************************/
				// default(none) forces us to specify all of the variables that are used in the parallel section.  
				//    Otherwise, they are implicitly shared (and this was cause of several bugs.)
#pragma omp parallel \
	private(ix, iy, ixa, iya, wave, t, timer) \
  shared(pot, pCount, picts, muls, collectedIntensity, total_time, waves) \
	default(none)
#pragma omp for
				for (i=0; i < (muls.scanXN * muls.scanYN); i++)
				{
					timer=cputim();
					ix = i / muls.scanYN;
					iy = i % muls.scanYN;

					wave = waves[omp_get_thread_num()];

					//printf("Scanning: %d %d %d %d\n",ix,iy,pCount,muls.nx);

					/* if this is run=0, create the inc. probe wave function */
					if (pCount == 0) 
					{
						probe(&muls, wave, muls.nx/2*muls.resolutionX, muls.ny/2*muls.resolutionY);

						// TODO: modifying shared value from multiple threads?
						//muls.nslic0 = 0;
						//wave->thickness = 0.0;
					}
                                          
					else 
					{
        					/* load incident wave function and then propagate it */
                                          
                                          wave->ReadWave(ix, iy); /* this also sets the thickness!!! */
						// TODO: modifying shared value from multiple threads?
						//muls.nslic0 = pCount;
					}
					/* run multislice algorithm
					   and save exit wave function for this position 
					   (done by runMulsSTEM), 
					   but we need to define the file name */
					muls.saveFlag = 1;

					wave->iPosX =(int)(ix*(muls.scanXStop-muls.scanXStart)/
									  ((float)muls.scanXN*muls.resolutionX));
					wave->iPosY = (int)(iy*(muls.scanYStop-muls.scanYStart)/
									   ((float)muls.scanYN*muls.resolutionY));
					if (wave->iPosX > muls.potNx-muls.nx)
					{
						wave->iPosX = muls.potNx-muls.nx;  
					}
					if (wave->iPosY > muls.potNy-muls.ny)
					{
						wave->iPosY = muls.potNy-muls.ny;
					}

					// MCS - update the probe wavefunction with its position

					runMulsSTEM(&muls,wave, pot); 


					/***************************************************************
					* In order to save some disk space we will add the diffraction 
					* patterns to their averages now.  The diffraction pattern 
					* should be stored in wave->diffpat (which each thread has independently), 
					* if collectIntensity() has been executed correctly.
					***************************************************************/

					#pragma omp atomic
					collectedIntensity += wave->intIntensity;

					if (pCount == picts-1)  /* if this is the last slice ... */
					{
                                          if (muls.saveLevel > 0) 
						{
							if (muls.avgCount == 0)  
							{
								// initialize the avgArray from the diffpat
								for (ixa=0;ixa<muls.nx;ixa++) 
								{
									for (iya=0;iya<muls.ny;iya++)
									{
										wave->avgArray[ixa][iya]=wave->diffpat[ixa][iya];
									}
								}
							}
							else 
							{
								// printf("Will read image %d %d\n",muls.nx, muls.ny);	
                                                          wave->ReadAvgArray(ix, iy);
								for (ixa=0;ixa<muls.nx;ixa++) for (iya=0;iya<muls.ny;iya++) {
									t = ((float_tt)muls.avgCount * wave->avgArray[ixa][iya] +
										wave->diffpat[ixa][iya]) / ((float_tt)(muls.avgCount + 1));
									if (muls.avgCount>1)
									{
										#pragma omp atomic
										muls.chisq[muls.avgCount-1] += (wave->avgArray[ixa][iya]-t)*
											(wave->avgArray[ixa][iya]-t);
									}
									wave->avgArray[ixa][iya] = t;
								}
							}
							// Write the array to a file, resize and crop it, 
							wave->WriteAvgArray(ix, iy);
							}	
							else {
								if (muls.avgCount > 0)	muls.chisq[muls.avgCount-1] = 0.0;
							}
					} /* end of if pCount == picts, i.e. conditional code, if this
						  * was the last slice
						  */

					#pragma omp atomic
					++muls.complete_pixels;

					if (muls.displayProgInterval > 0) if ((muls.complete_pixels) % muls.displayProgInterval == 0) 
					{
						#pragma omp atomic
						total_time += cputim()-timer;
						printf("Pixels complete: (%d/%d), int.=%.3f, avg time per pixel: %.2fsec\n",
							muls.complete_pixels, muls.scanXN*muls.scanYN, wave->intIntensity,
							(total_time)/muls.complete_pixels);
						timer=cputim();
					}
				} /* end of looping through STEM image pixels */
				/* save STEM images in img files */
				saveSTEMImages(&muls);
				muls.totalSliceCount += muls.slices;
			} /* end of loop through thickness (pCount) */
		} /* end of  while (readparam("sequence: ",buf,0)) */
		// printf("Total CPU time = %f sec.\n", cputim()-timerTot ); 

		/*************************************************************/
		if (muls.avgCount>1)
			muls.chisq[muls.avgCount-1] = muls.chisq[muls.avgCount-1]/(double)(muls.nx*muls.ny);
		muls.intIntensity = collectedIntensity/(muls.scanXN*muls.scanYN);
		displayProgress(1);
	} /* end of loop over muls.avgCount */

}

