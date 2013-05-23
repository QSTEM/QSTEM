/* good settings:
3x3 unit cells
50 kV: gaussian parameter: 8

5 kV: gaussian parameter: 11
good energies: 327, 360,393,520 keV
*/
#define VERSION 3.0
#define VIB_IMAGE_TEST
// #define VIB_IMAGE_TEST_CBED

#ifndef WIN32
#define UNIX 
#endif
/* #define USE_FFT_POT */
// #define WINDOWS
// for memory leak checking
#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#ifdef WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif
//#include <string.h>
#include <string>
#ifndef WIN32
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "memory_fftw3.h"	/* memory allocation routines */
#include "readparams.h"
#include "imagelib_fftw3.h"
#include "fileio_fftw3.h"
#include "matrixlib.h"
#include "stemlib.h"
#include "stemutil.h"
// #include "weblib.h"
#include "customslice.h"
#include "data_containers.h"
#include "stringutils.h"

const char *resultPage = "result.html";
/* global variable: */
MULS muls;
// int fftMeasureFlag = FFTW_MEASURE;
int fftMeasureFlag = FFTW_ESTIMATE;

void makeAnotation(float_tt **pict,int nx,int ny,char *text);
//void initMuls();
void writeIntPix(char *outFile,float_tt **pict,int nx,int ny);
void runMuls(int lstart);
void saveLineScan(int run);
void readBeams(FILE *fp);
void doCBED();
void doSTEM();
void doTEM();
void doMSCBED();
void doTOMO();
void readFile();
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
	float_tt timerTot;
	char fileName[256]; 
	char cinTemp[BUF_LEN];

	timerTot = cputim();
	for (i=0;i<BUF_LEN;i++)
		cinTemp[i] = 0;
	muls.nCellX = 1; muls.nCellY = 1; muls.nCellZ = 1; 

#ifdef UNIX
	system("date");
#endif

	/*************************************************************
	* read in the parameters
	************************************************************/  
	if (argc < 2)
		sprintf(fileName,"stem.dat");
	else
		strcpy(fileName,argv[1]);
	if (parOpen(fileName) == 0) {
		printf("could not open input file %s!\n",fileName);
		usage();
		exit(0);
	}
	readFile();

	displayParams();
#ifdef _OPENMP
	omp_set_dynamic(1);
#endif

	/* report the result on the web page 
	add2Webpage("diff_0_0.img"); 
	*/
	if (muls.mode == STEM) {
		// sprintf(systStr,"mkdir %s",muls.folder);
		// system(systStr);
		for (muls.avgCount=0;muls.avgCount<muls.avgRuns;muls.avgCount++) {
			muls.dE_E = muls.dE_EArray[muls.avgCount];
			// printf("dE/E: %g\n",muls.dE_E);
#ifndef WIN32
			//probePlot(&muls, wave);
#endif
		}    

		// drawStructure();

		muls.showProbe = 0;
		muls.avgCount = -1;
		/*
		wave->thickness = 600;
		makeSTEMPendelloesungPlot();
		*/    
#ifndef WIN32
		//add2STEMWebpage();
#endif
	}

	switch (muls.mode) {
	  case CBED:   doCBED();   break;	
	  case STEM:   doSTEM();   break;
	  case TEM:    doTEM();    break;
	  case MSCBED: doMSCBED(); break;
	  case TOMO:   doTOMO();   break;
	  // case REFINE: doREFINE(); break;
	  default:
		  printf("Mode not supported\n");
	}

	parClose();
#if _DEBUG
	_CrtDumpMemoryLeaks();
#endif

	return 0;
}

int DirExists(char *filename) { 
  struct stat status; 
  status.st_mode = 0;
  
  // Attempt to get the file attributes 
  if (stat(filename,&status) == 0) return 0;
  if (status.st_mode & S_IFDIR )  return 1;
  return 0;
}


/************************************************************************
*
***********************************************************************/
void displayProgress(int flag) {
	// static float_tt timer;
	static float_tt timeAvg = 0;
	static float_tt intensityAvg = 0;
	static time_t time0,time1;
	float_tt curTime;
	int jz;

	if (flag < 0) {
		time(&time0);
		// timer = cputim();
		return;
	}
	time(&time1);  
	curTime = static_cast<float_tt>(difftime(time1,time0));
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
	float_tt k2max,temp;
	int i,j;
	static char Date[16],Time[16];
	time_t caltime;
	struct tm *mytime;

	if (muls.printLevel < 1) {
		if ((fpDir = fopen(muls.folder,"r"))) {
			fclose(fpDir);
			// printf(" (already exists)\n");
		}
		else {
			sprintf(systStr,"mkdir %s",muls.folder);
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
	printf("*****************************************************\n");
	printf("* Print level:          %d\n",muls.printLevel);
	printf("* Save level:           %d\n",muls.saveLevel);
	printf("* Input file:           %s\n",muls.atomPosFile);
	if (muls.savePotential)
		printf("* Potential file name:  %s\n",muls.fileBase);
	/* create the data folder ... */
	printf("* Data folder:          ./%s/ ",muls.folder); 
	if (DirExists(muls.folder)) {
	// if ((fpDir = fopen(muls.folder,"r"))!= NULL) {
	// 	fclose(fpDir);
		printf(" (already exists)\n");
	}
	else {
		sprintf(systStr,"mkdir %s",muls.folder);
		system(systStr);
		printf(" (created)\n");
	}

	if ((muls.cubex == 0) || (muls.cubey == 0) || (muls.cubez == 0))
		printf("* Unit cell:            ax=%g by=%g cz=%g\n",
		muls.cellDims[0],muls.cellDims[1],muls.cellDims[2]);
	else {
		printf("* Size of Cube:         ax=%g by=%g cz=%g\n",
			muls.cubex,muls.cubey,muls.cubez);
		printf("* Cube size adjusted:   %s\n",muls.adjustCubeSize ? "yes" : "no");
	}

	printf("* Super cell:           %d x %d x %d unit cells\n",muls.nCellX,muls.nCellY,muls.nCellZ);
	printf("* Number of atoms:      %d (super cell)\n",muls.natom);
	printf("* Crystal tilt:         x=%g deg, y=%g deg, z=%g deg\n",
		muls.ctiltx*RAD2DEG,muls.ctilty*RAD2DEG,muls.ctiltz*RAD2DEG);
	printf("* Beam tilt:            x=%g deg, y=%g deg (tilt back == %s)\n",muls.btiltx*RAD2DEG,muls.btilty*RAD2DEG,
		(muls.tiltBack == 1 ? "on" : "off"));
	printf("* Model dimensions:     ax=%gA, by=%gA, cz=%gA (after tilt)\n"
		"*                       sampled every %g x %g x %g A\n",
		muls.cellDims[0],muls.cellDims[1],muls.cellDims[2],muls.resolutionX,muls.resolutionY,muls.sliceThickness);
	printf("* Atom species:         %d (Z=%d",muls.atomKinds,muls.Znums[0]);
	for (i=1;i<muls.atomKinds;i++) printf(", %d",muls.Znums[i]); printf(")\n");
	printf("* Super cell divisions: %d (in z direction) %s\n",muls.cellDiv,muls.equalDivs ? "equal" : "non-equal");
	printf("* Slices per division:  %d (%gA thick slices [%scentered])\n",
		muls.slices,muls.sliceThickness,(muls.centerSlices) ? "" : "not ");
	printf("* Output every:         %d slices\n",muls.outputInterval);

	printf("* Potential:            ");
	if (muls.potential3D) printf("3D"); else printf("2D");
	if (muls.fftpotential) printf(" (fast method)\n"); else printf(" (slow method)\n");	
	printf("* Pot. array offset:    (%g,%g,%g)A\n",muls.potOffsetX,muls.potOffsetY,muls.czOffset);
	printf("* Potential periodic:   (x,y): %s, z: %s\n",
		(muls.nonPeriod) ? "no" : "yes",(muls.nonPeriodZ) ? "no" : "yes");

	printf("* Beams:                %d x %d \n",muls.nx,muls.ny);  
	printf("* Acc. voltage:         %g (lambda=%gA)\n",muls.v0,wavelength(muls.v0));
	printf("* C_3 (C_s):            %g mm\n",muls.Cs*1e-7);
	printf("* C_1 (Defocus):        %g nm%s\n",0.1*muls.df0,
		(muls.Scherzer == 1) ? " (Scherzer)" : (muls.Scherzer==2) ? " (opt.)":"");
	printf("* Astigmatism:          %g nm, %g deg\n",0.1*muls.astigMag,RAD2DEG*muls.astigAngle);

	// more aberrations:
	if (muls.a33 > 0)
		printf("* a_3,3:                %g nm, phi=%g deg\n",muls.a33*1e-1,muls.phi33*RAD2DEG);
	if (muls.a31 > 0)
		printf("* a_3,1:                %g nm, phi=%g deg\n",muls.a31*1e-1,muls.phi31*RAD2DEG);

	if (muls.a44 > 0)
		printf("* a_4,4:                %g um, phi=%g deg\n",muls.a44*1e-4,muls.phi44*RAD2DEG);
	if (muls.a42 > 0)
		printf("* a_4,2:                %g um, phi=%g deg\n",muls.a42*1e-4,muls.phi42*RAD2DEG);

	if (muls.a55 > 0)
		printf("* a_5,5:                %g um, phi=%g deg\n",muls.a55*1e-4,muls.phi55*RAD2DEG);
	if (muls.a53 > 0)
		printf("* a_5,3:                %g um, phi=%g deg\n",muls.a53*1e-4,muls.phi53*RAD2DEG);
	if (muls.a51 > 0)
		printf("* a_5,1:                %g um, phi=%g deg\n",muls.a51*1e-4,muls.phi51*RAD2DEG);

	if (muls.a66 > 0)
		printf("* a_6,6:                %g um, phi=%g deg\n",muls.a66*1e-7,muls.phi66*RAD2DEG);
	if (muls.a64 > 0)
		printf("* a_6,4:                %g um, phi=%g deg\n",muls.a64*1e-7,muls.phi64*RAD2DEG);
	if (muls.a62 > 0)
		printf("* a_6,2:                %g um, phi=%g deg\n",muls.a62*1e-7,muls.phi62*RAD2DEG);
	if (muls.C5 != 0)
		printf("* C_5:                  %g mm\n",muls.C5*1e-7);

	printf("* C_c:                  %g mm\n",muls.Cc*1e-7);
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
	printf("* Temperature:          %gK\n",muls.tds_temp);
	if (muls.tds)
		printf("* TDS:                  yes (%d runs)\n",muls.avgRuns);
	else
		printf("* TDS:                  no\n"); 
	if (muls.imageGamma == 0)
		printf("* Gamma for diff. patt: logarithmic\n");
	else
		printf("* Gamma for diff. patt: %g\n",muls.imageGamma);

	/**********************************************************
	* FFTW specific data structures (stores in row major order)
	*/
	if (fftMeasureFlag == FFTW_MEASURE)
		printf("* Probe array:          %d x %d pixels (optimized)\n",muls.nx,muls.ny);
	else
		printf("* Probe array:          %d x %d pixels (estimated)\n",muls.nx,muls.ny);
	printf("*                       %g x %gA\n",
		muls.nx*muls.resolutionX,muls.ny*muls.resolutionY);

	if (fftMeasureFlag == FFTW_MEASURE)
		printf("* Potential array:      %d x %d (optimized)\n",muls.potNx,muls.potNy);
	else
		printf("* Potential array:      %d x %d (estimated)\n",muls.potNx,muls.potNy);
	printf("*                       %g x %gA\n",muls.potSizeX,muls.potSizeY);
	printf("* Scattering factors:   %d\n",muls.scatFactor);
	/***************************************************/
	/*  printf("Optimizing fftw plans according to probe array (%d x %dpixels = %g x %gA) ...\n",
	muls.nx,muls.ny,muls.nx*muls.resolutionX,muls.ny*muls.resolutionY);
	*/
	k2max = muls.nx/(2.0f*muls.potSizeX);
	temp = muls.ny/(2.0f*muls.potSizeY);
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
		printf("* Number of detectors:  %d\n",muls.detectorNum);
		for (i=0;i<muls.detectorNum;i++) {
			printf("* %d (\"%s\"):",i+1,muls.detectors[0][i].name.c_str());
			// TODO: some fixed-length formatting going on here.
			for (j=0;j<14-strlen(muls.detectors[0][i].name.c_str());j++) printf(" ");
			printf(" %g .. %g mrad = (%.2g .. %.2g 1/A)\n",
				muls.detectors[0][i].rInside,
				muls.detectors[0][i].rOutside,
				muls.detectors[0][i].k2Inside,
				muls.detectors[0][i].k2Outside);
			if ((muls.detectors[0][i].shiftX != 0) ||(muls.detectors[0][i].shiftY != 0))
				printf("*   center shifted:     dkx=%g, dky=%g\n",
				muls.detectors[0][i].shiftX,muls.detectors[0][i].shiftY);
		}
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
			muls.tomoStart,muls.tomoStart*0.18/PI);
		printf("* Angular increase:     %g mrad (%g deg)\n",
			muls.tomoStep,muls.tomoStep*0.180/PI);
		printf("* Number of dp's:       %d\n",muls.tomoCount);
		printf("* Zoom factor:          %g\n",muls.zoomFactor);


	}

	printf("*\n*****************************************************\n");

	/* k2max = muls.alpha * 0.001;  k2max = aperture in rad */
}


void readArray(const char *title,float_tt *array,int N) {
	int i;
	char buf[512],*str;

	if (!readparam(title,buf,1)) printf("%s array not found - exit\n",title), exit(0);
	i=0;
	str = buf;
	if (strchr(" \t\n",*str) != NULL) str = strnext(str," \t");
	while (i<N) {
		array[i++] = static_cast<float_tt>(atof(str));
		str = strnext(str," \t\n");
		while (str == NULL) {
			if (!readNextLine(buf,511)) 
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
	QSfMat sfTable;
	QSfVec kArray;


	char buf[256];
	std::string elem;

	if (readparam("Nk:",buf,1))
		Nk = atoi(buf);
	else {
		printf("Could not find number of k-points for custom scattering factors (Nk)\n");
		exit(0);
	}

	// Table has atomKinds columns (x), with sf's for each atom kind stored
	//    as a column.
	sfTable = QSfMat(muls.atomKinds,Nk+1);
	kArray  = QSfVec(Nk+1);

	// read the k-values:
	readArray("k:",kArray.data(),Nk);
	kArray[Nk] = 2.0f*kArray[Nk-1];

	for (j=0;j<muls.atomKinds;j++) {
		elem = ElTable::GetSymbol(muls.Znums[j]).c_str();
		// printf("%s\n",elem);
		readArray(elem.c_str(), sfTable.col(j).data(),Nk);
		sfTable(Nk,j) = 0.0;
	}

	if (0) {
		printf("k: ");
		for (i=0;i<=Nk;i++) printf("%.3f ",kArray[i]);
		for (j=0;j<muls.atomKinds;j++) {
			printf("\n%2d: ",muls.Znums[j]);
			for (i=0;i<=Nk;i++) printf("%.3f ",sfTable(i,j));
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
void readFile() {
	char answer[256],*strptr;
	FILE *fpTemp;
	float ax,by,cz;
	char buf[BUF_LEN],*strPtr;
	int i,ix;
	int potDimensions[2];
	long ltime;
	unsigned long iseed;
	float_tt dE_E0,x,y,dx,dy;

	ltime = (long) time(NULL);
	iseed = (unsigned) ltime;
	muls.cubex = 0.0;
	muls.cubey = 0.0;
	muls.cubez = 0.0;

	muls.mode = STEM;
	if (readparam("mode:",buf,1)) {
		if (strstr(buf,"STEM")) muls.mode = STEM;
		else if (strstr(buf,"TEM")) muls.mode = TEM;
		else if (strstr(buf,"CBED")) muls.mode = CBED;
		else if (strstr(buf,"TOMO")) muls.mode = TOMO;
		else if (strstr(buf,"REFINE")) muls.mode = REFINE;
	}

	muls.printLevel = 2;
	if (readparam("print level:",buf,1)) sscanf(buf,"%d",&(muls.printLevel));
	muls.saveLevel = 0;
	if (readparam("save level:",buf,1)) sscanf(buf,"%d",&(muls.saveLevel));


	/************************************************************************
	* Basic microscope/simulation parameters, 
	*/ 
	//muls.fileBase = (char *)malloc(512);
	//muls.atomPosFile = (char *)malloc(512);
	if (!readparam("filename:",buf,1)) exit(0); sscanf(buf,"%s",muls.fileBase);
	// printf("buf: %s\n",buf);
	// printf("fileBase: %s\n",muls.fileBase);
	/*
	while ((strptr = strchr(muls.fileBase,'"')) != NULL) {
		if (strptr == muls.fileBase) sscanf(strptr+1,"%s",muls.fileBase);	
		else *strptr = '\0';
	}
	*/
	// search for second '"', in case the filename is in quotation marks:
	// remove any quotation marks from the filename:
	std::string fname = std::string(buf);
	fname.erase(std::remove( fname.begin(), fname.end(), '\"' ), fname.end());
	muls.fileBase += fname;
	// append quotation marks to the end, if the muls filebase starts with them.
	if (muls.fileBase[0] == '"') {
		muls.fileBase+='"';
	}
	muls.fileBase = StripSpaces(muls.fileBase);
		//strPtr = strchr(buf,'"');
		//strcpy(muls.fileBase,strPtr+1);
		//strPtr = strchr(muls.fileBase.c_str(),'"');
		//*strPtr = '\0';
	//}


	// printf("fileBase: %s\n",muls.fileBase);

	// printf("File Base: %s\n",muls.fileBase);
	if (readparam("NCELLX:",buf,1)) sscanf(buf,"%d",&(muls.nCellX));
	if (readparam("NCELLY:",buf,1)) sscanf(buf,"%d",&(muls.nCellY));

	muls.cellDiv = 1;
	if (readparam("NCELLZ:",buf,1)) {
		sscanf(buf,"%s",answer);
		if ((strPtr = strchr(answer,'/')) != NULL) {
			strPtr[0] = '\0';
			muls.cellDiv = atoi(strPtr+1);
		}
		muls.nCellZ = atoi(answer);
	}

	/*************************************************
	* Read the beam tilt parameters
	*/
	muls.btiltx = 0.0;
	muls.btilty = 0.0;
	muls.tiltBack = 1;
	answer[0] = '\0';
	if (readparam("Beam tilt X:",buf,1)) { 
		sscanf(buf,"%g %s",&(muls.btiltx),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.btiltx *= static_cast<float_tt>(PI/180.0);
	}
	answer[0] = '\0';
	if (readparam("Beam tilt Y:",buf,1)) { 
		sscanf(buf,"%g %s",&(muls.btilty),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.btilty *= static_cast<float_tt>(PI/180.0);
	}  
	if (readparam("Tilt back:",buf,1)) { 
		sscanf(buf,"%s",answer);
		muls.tiltBack  = (tolower(answer[0]) == (int)'y');
	}


	/*************************************************
	* Read the crystal tilt parameters
	*/
	muls.ctiltx = 0.0;  /* tilt around X-axis in mrad */
	muls.ctilty = 0.0;  /* tilt around y-axis in mrad */	
	muls.ctiltz = 0.0;  /* tilt around z-axis in mrad */	
	answer[0] = '\0';
	if (readparam("Crystal tilt X:",buf,1)) { 
		sscanf(buf,"%g %s",&(muls.ctiltx),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.ctiltx *= static_cast<float_tt>(PI/180.0);
	}
	answer[0] = '\0';
	if (readparam("Crystal tilt Y:",buf,1)) { 
		sscanf(buf,"%g %s",&(muls.ctilty),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.ctilty *= static_cast<float_tt>(PI/180.0);
	}  
	answer[0] = '\0';
	if (readparam("Crystal tilt Z:",buf,1)) { 
		sscanf(buf,"%g %s",&(muls.ctiltz),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.ctiltz *= static_cast<float_tt>(PI/180.0);
	}
	muls.cubex = 0; muls.cubey = 0; muls.cubez = 0;
	if (readparam("Cube:",buf,1)) { 
		sscanf(buf,"%g %g %g",&(muls.cubex),&(muls.cubey),&(muls.cubez)); /* in A */
	}

	muls.adjustCubeSize = 0;
	if (readparam("Adjust cube size with tilt:",buf,1)) { 
		sscanf(buf,"%s",answer);
		muls.adjustCubeSize  = (tolower(answer[0]) == (int)'y');
	}

	/***************************************************************************
	* temperature related data must be read before reading the atomic positions:
	***************************************************************************/
	if (readparam("tds:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.tds = (tolower(answer[0]) == (int)'y');
	}
	else muls.tds = 0;
	if (readparam("temperature:",buf,1)) sscanf(buf,"%g",&(muls.tds_temp));
	else muls.tds_temp = 300.0;
	muls.Einstein = 1;
	//muls.phononFile = NULL;
	if (readparam("phonon-File:",buf,1)) {
		sscanf(buf,"%s",muls.phononFile);
		muls.Einstein = 0;
	}

	/**********************************************************************
	* Read the atomic model positions !!!
	*********************************************************************/
	muls.atomPosFile = muls.fileBase;
	
	if (!muls.fileBase.empty())
	{
		char lastChar = *muls.fileBase.rbegin();
		// if our fileBase is a directory, but no file:
		if (lastChar == '\\' || lastChar == '/')
		{
			muls.fileBase = "_";
		}
		// Otherwise, try to split off the directory part, keeping only the filename.
		muls.fileBase = ExtractFilename(muls.fileBase);
	}
	else
	{
		muls.fileBase = "_";
	}

	/* take atomPosFile as is, or add an ending to it, if it has none yet
	*/
	if (muls.atomPosFile.find('.')==muls.atomPosFile.npos) {
		std::string tmp;
		tmp = muls.atomPosFile += ".cssr";
		if ((fpTemp=fopen(tmp.c_str(),"r")) == NULL) {
			tmp = muls.atomPosFile += ".cfg";
			if ((fpTemp=fopen(tmp.c_str(),"r")) == NULL) {
				printf("Could not find input file %s.cssr or %s.cfg\n",
					muls.atomPosFile, muls.atomPosFile);
				exit(0);
			}
			muls.atomPosFile+=".cfg";
			fclose(fpTemp);
		}
		else {
			muls.atomPosFile+=".cssr";
			fclose(fpTemp);
		}
	}
	// We need to initialize a few variables, before reading the atomic 
	// positions for the first time.
	muls.natom = 0;
	muls.atomKinds = 0;
	//muls.atoms = NULL;
	//muls.Znums = NULL;
	//muls.u2 = NULL;
	//muls.u2avg = NULL;

	muls.xOffset = 0.0; /* slize z-position offset in cartesian coords */
	if (readparam("xOffset:",buf,1)) sscanf(buf,"%g",&(muls.xOffset));
	muls.yOffset = 0.0; /* slize z-position offset in cartesian coords */
	if (readparam("yOffset:",buf,1)) sscanf(buf,"%g",&(muls.yOffset));
	// printf("Reading Offset: %f, %f\n",muls.xOffset,muls.yOffset);

	// the last parameter is handleVacancies.  If it is set to 1 vacancies 
	// and multiple occupancies will be handled. 
	// _CrtSetDbgFlag  _CRTDBG_CHECK_ALWAYS_DF();
	// printf("memory check: %d, ptr= %d\n",_CrtCheckMemory(),(int)malloc(32*sizeof(char)));

	muls.atoms = readUnitCell(muls.natom, muls.atomPosFile, muls, 1);

	// printf("memory check: %d, ptr= %d\n",_CrtCheckMemory(),(int)malloc(32*sizeof(char)));


	if (muls.atoms.empty()) 
	{
		printf("Error reading atomic positions!\n");
		exit(0);
	}
	if (muls.natom == 0) {
		printf("No atom within simulation boundaries!\n");
		exit(0);
	}
	// printf("hello!\n");
	ax = muls.cellDims[0]/muls.nCellX;
	by = muls.cellDims[1]/muls.nCellY;;
	cz =  muls.cellDims[2]/muls.nCellZ;


	/*****************************************************************
	* Done reading atomic positions 
	****************************************************************/

	if (!readparam("nx:",buf,1)) exit(0); sscanf(buf,"%d",&(muls.nx));
	if (readparam("ny:",buf,1)) sscanf(buf,"%d",&(muls.ny));
	else muls.ny = muls.nx;

	muls.resolutionX = 0.0;
	muls.resolutionY = 0.0;
	if (readparam("resolutionX:",buf,1)) sscanf(buf,"%g",&(muls.resolutionX));
	if (readparam("resolutionY:",buf,1)) sscanf(buf,"%g",&(muls.resolutionY));
	if (!readparam("v0:",buf,1)) exit(0); sscanf(buf,"%g",&(muls.v0));

	muls.centerSlices = 0;
	if (readparam("center slices:",buf,1)) {
		// answer[0] =0;
		sscanf(buf,"%s",answer);
		// printf("center: %s (%s)\n",answer,buf);
		muls.centerSlices = (tolower(answer[0]) == (int)'y');
	}
	// just in case the answer was not exactly 1 or 0:
	// muls.centerSlices = (muls.centerSlices) ? 1 : 0;

	muls.sliceThickness = 0.0;
	if (readparam("slice-thickness:",buf,1)) {
		sscanf(buf,"%g",&(muls.sliceThickness));
		if (readparam("slices:",buf,1)) {
			sscanf(buf,"%d",&(muls.slices));
		}
		else {
			if (muls.cubez >0)
				muls.slices = (int)(muls.cubez/(muls.cellDiv*muls.sliceThickness)+0.99);
			else
				muls.slices = (int)(muls.cellDims[2]/(muls.cellDiv*muls.sliceThickness)+0.99);
		}
		muls.slices += muls.centerSlices;
	}
	else {
		muls.slices = 0; 
		if (readparam("slices:",buf,1)) {
			sscanf(buf,"%d",&(muls.slices));
			// muls.slices = (int)(muls.slices*muls.nCellZ/muls.cellDiv);
			if (muls.sliceThickness == 0.0) {
				if ((muls.slices == 1) && (muls.cellDiv == 1)) {
					if (muls.cubez >0)
						muls.sliceThickness = (muls.centerSlices) ? 2.0f*muls.cubez/muls.cellDiv : muls.cubez/muls.cellDiv;
					else
						// muls.sliceThickness = (muls.centerSlices) ? 2.0*muls.c/(muls.cellDiv) : muls.c/(muls.cellDiv);
						muls.sliceThickness = (muls.centerSlices) ? 1.0f*muls.cellDims[2]/(muls.cellDiv) : muls.cellDims[2]/(muls.cellDiv);
				}
				else {
					if (muls.cubez >0) {
						muls.sliceThickness = muls.cubez/(muls.cellDiv*muls.slices-muls.centerSlices);
					}
					else {
						muls.sliceThickness = muls.cellDims[2]/(muls.cellDiv*muls.slices);
					}
				}
			}
			else {
				muls.cellDiv = (muls.cubez >0) ? (int)ceil(muls.cubez/(muls.slices*muls.sliceThickness)) :
					(int)ceil(muls.cellDims[2]/(muls.slices*muls.sliceThickness));
			if (muls.cellDiv < 1) muls.cellDiv = 1;
			}
		}
	}
	if (muls.slices == 0) {
		if (muls.printLevel > 0) printf("Error: Number of slices = 0\n");
		exit(0);
	}
	/* Find out whether we need to recalculate the potential every time, or not
	*/

	muls.equalDivs = ((!muls.tds)  && (muls.nCellZ % muls.cellDiv == 0) && 
		(fabs(muls.slices*muls.sliceThickness-muls.cellDims[2]/muls.cellDiv) < 1e-5));

	// read the output interval:
	muls.outputInterval = muls.slices;
	if (readparam("slices between outputs:",buf,1)) sscanf(buf,"%d",&(muls.outputInterval));
	if (muls.outputInterval < 1) muls.outputInterval= muls.slices;

	// muls is initialized in its constructor.
	//initMuls();  
	muls.czOffset = 0.0; /* slize z-position offset in cartesian coords */
	if (readparam("zOffset:",buf,1)) sscanf(buf,"%g",&(muls.czOffset));



	/***********************************************************************
	* Fit the resolution to the wave function array, if not specified different
	*/
	if (muls.resolutionX == 0.0)
		muls.resolutionX = muls.cellDims[0] / (float_tt)muls.nx;
	if (muls.resolutionY == 0.0)
		muls.resolutionY = muls.cellDims[1] / (float_tt)muls.ny;



	/************************************************************************
	* Optional parameters:
	* determine whether potential periodic or not, etc.:
	*/
	muls.nonPeriodZ = 1;
	muls.nonPeriod = 1;
	if (readparam("periodicXY:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.nonPeriod = (tolower(answer[0]) != (int)'y');
	}
	if (readparam("periodicZ:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.nonPeriodZ = (tolower(answer[0]) != (int)'y'); /* if 'y' -> nonPeriodZ=0 */
	}
	if ((muls.nonPeriodZ == 0) && (muls.cellDiv > 1)) {
		printf("****************************************************************\n"
			"* Warning: cannot use cell divisions >1 and Z-periodic potential\n"
			"* periodicZ = NO\n"
			"****************************************************************\n");
		muls.nonPeriodZ = 1;
	}

	muls.bandlimittrans = 1;
	if (readparam("bandlimit f_trans:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.bandlimittrans = (tolower(answer[0]) == (int)'y');
	}    
	muls.readPotential = 0;
	if (readparam("read potential:",buf,1)) {
		sscanf(buf," %s",answer);
		muls.readPotential = (tolower(answer[0]) == (int)'y');
	}  
	muls.savePotential = 0;
	if (readparam("save potential:",buf,1)) {
		sscanf(buf," %s",answer);
		muls.savePotential = (tolower(answer[0]) == (int)'y');
	}  
	muls.saveTotalPotential = 0;
	if (readparam("save projected potential:",buf,1)) {
		sscanf(buf," %s",answer);
		muls.saveTotalPotential = (tolower(answer[0]) == (int)'y');
	}  
	muls.plotPotential = 0;
	if (readparam("plot V(r)*r:",buf,1)) {
		sscanf(buf," %s",answer);
		muls.plotPotential = (tolower(answer[0]) == (int)'y');
	}  
	muls.fftpotential = 1;
	if (readparam("one time integration:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.fftpotential = (tolower(answer[0]) == (int)'y');
	}
	muls.potential3D = 1;
	if (readparam("potential3D:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.potential3D = (tolower(answer[0]) == (int)'y');
	}
	muls.avgRuns = 10;
	if (readparam("Runs for averaging:",buf,1))
		sscanf(buf,"%d",&(muls.avgRuns));

	muls.storeSeries = 0;
	if (readparam("Store TDS diffr. patt. series:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.storeSeries = (tolower(answer[0]) == (int)'y');
	}  

	if (!muls.tds) muls.avgRuns = 1;

	muls.scanXStart = muls.cellDims[0]/2.0f;
	muls.scanYStart = muls.cellDims[1]/2.0f;
	muls.scanXN = 1;
	muls.scanYN = 1;
	muls.scanXStop = muls.scanXStart;
	muls.scanYStop = muls.scanYStart;


	switch (muls.mode) {
		/////////////////////////////////////////////////////////
		// read the position for doing CBED: 
	case CBED:
		if (readparam("scan_x_start:",buf,1)) sscanf(buf,"%g",&(muls.scanXStart));
		if (readparam("scan_y_start:",buf,1)) sscanf(buf,"%g",&(muls.scanYStart));
		muls.scanXStop = muls.scanXStart;
		muls.scanYStop = muls.scanYStart;
		break;

		/////////////////////////////////////////////////////////
		// Read STEM scanning parameters 

	case STEM:
		/* Read in scan coordinates: */
		if (!readparam("scan_x_start:",buf,1)) exit(0); 
		sscanf(buf,"%g",&(muls.scanXStart));
		if (!readparam("scan_x_stop:",buf,1)) exit(0); 
		sscanf(buf,"%g",&(muls.scanXStop));
		if (!readparam("scan_x_pixels:",buf,1)) exit(0); 
		sscanf(buf,"%d",&(muls.scanXN));
		if (!readparam("scan_y_start:",buf,1)) exit(0); 
		sscanf(buf,"%g",&(muls.scanYStart));
		if (!readparam("scan_y_stop:",buf,1)) exit(0); 
		sscanf(buf,"%g",&(muls.scanYStop));
		if (!readparam("scan_y_pixels:",buf,1)) exit(0); 
		sscanf(buf,"%d",&(muls.scanYN));

		if (muls.scanXN < 1) muls.scanXN = 1;
		if (muls.scanYN < 1) muls.scanYN = 1;
		// if (muls.scanXStart > muls.scanXStop) muls.scanXN = 1;

		muls.displayProgInterval = muls.scanYN*muls.scanYN;
		if (readparam("propagation progress interval:",buf,1)) 
			sscanf(buf,"%d",&(muls.displayProgInterval));
	}
	muls.displayPotCalcInterval = 1000;
	if (readparam("potential progress interval:",buf,1)) 
		sscanf(buf,"%d",&(muls.displayPotCalcInterval));
	// printf("Potential progress interval: %d\n",muls.displayPotCalcInterval);

	/**********************************************************************
	* Read STEM/CBED probe parameters 
	*/
	muls.dE_E = 0.0;
	muls.dI_I = 0.0;
	muls.dV_V = 0.0;
	muls.Cc = 0.0;
	if (readparam("dE/E:",buf,1))
		muls.dE_E = static_cast<float_tt>(atof(buf));
	if (readparam("dI/I:",buf,1))
		muls.dI_I = static_cast<float_tt>(atof(buf));
	if (readparam("dV/V:",buf,1))
		muls.dV_V = static_cast<float_tt>(atof(buf));
	if (readparam("Cc:",buf,1))
		muls.Cc = static_cast<float_tt>(1e7*atof(buf));


	/* memorize dE_E0, and fill the array of well defined energy deviations */
	dE_E0 = sqrt(muls.dE_E*muls.dE_E+
		muls.dI_I*muls.dI_I+
		muls.dV_V*muls.dV_V);
	muls.dE_EArray = QSfVec(muls.avgRuns+1);//(float_tt *)malloc((muls.avgRuns+1)*sizeof(float_tt));
	muls.dE_EArray[0] = 0.0;
	/***********************************************************
	* Statistical gaussian energy spread
	* (takes too long for the statistics to become gaussian)
	*/

	/* for (i = 1;i <= muls.avgRuns*muls.tds; i++) {
	muls.dE_EArray[i] = rangauss(&iseed)*dE_E0;     
	printf("dE/E[%d]: %g\n",i,muls.dE_EArray[i]); 
	}
	*/

	/**********************************************************
	* quick little fix to calculate gaussian energy distribution
	* without using statistics (better for only few runs)
	*/
	if (muls.printLevel > 0) printf("avgRuns: %d\n",muls.avgRuns);
	// serious bug in Visual C - dy comes out enormous.
	//dy = sqrt((float_tt)pi)/((float_tt)2.0*(float_tt)(muls.avgRuns));
	// using precalculated sqrt(pi):
	dy = 1.772453850905f/(2.0f*(muls.avgRuns));
	dx = static_cast<float_tt>(PI/((muls.avgRuns+1)*20.0f));
	for (ix=1,x=0,y=0;ix<muls.avgRuns;x+=dx) {
		y += exp(-x*x)*dx;
		if (y>=ix*dy) {
			muls.dE_EArray[ix++] = static_cast<float_tt>(x*2*dE_E0/PI);
			if (muls.printLevel > 2) printf("dE[%d]: %g eV\n",ix,muls.dE_EArray[ix-1]*muls.v0*1e3f);
			if (ix < muls.avgRuns) {
				muls.dE_EArray[ix] = -muls.dE_EArray[ix-1];
				ix ++;
				if (muls.printLevel > 2) printf("dE[%d]: %g eV\n",ix,muls.dE_EArray[ix-1]*muls.v0*1e3f);
			}
		}
	}


	if (!readparam("Cs:",buf,1))  exit(0); 
	sscanf(buf,"%g",&(muls.Cs)); /* in mm */
	muls.Cs *= 1.0e7; /* convert Cs from mm to Angstroem */

	muls.C5 = 0;
	if (readparam("C5:",buf,1)) { 
		sscanf(buf,"%g",&(muls.C5)); /* in mm */
		muls.C5 *= 1.0e7; /* convert C5 from mm to Angstroem */
	}

	/* assume Scherzer defocus as default */
	muls.df0 = -(float)sqrt(1.5*muls.Cs*(wavelength(muls.v0))); /* in A */
	muls.Scherzer = 1;
	if (readparam("defocus:",buf,1)) { 
		sscanf(buf,"%s",answer);
		/* if Scherzer defocus */
		if (tolower(answer[0]) == 's') {
			muls.df0 = -(float)sqrt(1.5*muls.Cs*(wavelength(muls.v0)));
			muls.Scherzer = 1;
		}
		else if (tolower(answer[0]) == 'o') {
			muls.df0 = -(float)sqrt(muls.Cs*(wavelength(muls.v0)));
			muls.Scherzer = 2;
		}
		else {
			sscanf(buf,"%g",&(muls.df0)); /* in nm */
			muls.df0 = 10.0f*muls.df0;       /* convert defocus to A */
			muls.Scherzer = (-(float)sqrt(1.5*muls.Cs*(wavelength(muls.v0)))==muls.df0);
		}
	}
	// Astigmatism:
	muls.astigMag = 0;
	if (readparam("astigmatism:",buf,1)) sscanf(buf,"%g",&(muls.astigMag)); 
	// convert to A from nm:
	muls.astigMag = 10.0f*muls.astigMag;
	muls.astigAngle = 0;
	if (readparam("astigmatism angle:",buf,1)) sscanf(buf,"%g",&(muls.astigAngle)); 
	// convert astigAngle from deg to rad:
	muls.astigAngle *= static_cast<float_tt>(PI/180.0);

	////////////////////////////////////////////////////////
	// read in more aberrations:
	muls.a33 = 0;
	muls.a31 = 0;
	muls.a44 = 0;
	muls.a42 = 0;
	muls.a55 = 0;
	muls.a53 = 0;
	muls.a51 = 0;
	muls.a66 = 0;
	muls.a64 = 0;
	muls.a62 = 0;

	muls.phi33 = 0;
	muls.phi31 = 0;
	muls.phi44 = 0;
	muls.phi42 = 0;
	muls.phi55 = 0;
	muls.phi53 = 0;
	muls.phi51 = 0;
	muls.phi66 = 0;
	muls.phi64 = 0;
	muls.phi62 = 0;

	if (readparam("a_33:",buf,1)) {sscanf(buf,"%g",&(muls.a33)); }
	if (readparam("a_31:",buf,1)) {sscanf(buf,"%g",&(muls.a31)); }
	if (readparam("a_44:",buf,1)) {sscanf(buf,"%g",&(muls.a44)); }
	if (readparam("a_42:",buf,1)) {sscanf(buf,"%g",&(muls.a42)); }
	if (readparam("a_55:",buf,1)) {sscanf(buf,"%g",&(muls.a55)); }
	if (readparam("a_53:",buf,1)) {sscanf(buf,"%g",&(muls.a53)); }
	if (readparam("a_51:",buf,1)) {sscanf(buf,"%g",&(muls.a51)); }
	if (readparam("a_66:",buf,1)) {sscanf(buf,"%g",&(muls.a66)); }
	if (readparam("a_64:",buf,1)) {sscanf(buf,"%g",&(muls.a64)); }
	if (readparam("a_62:",buf,1)) {sscanf(buf,"%g",&(muls.a62)); }

	if (readparam("phi_33:",buf,1)) {sscanf(buf,"%g",&(muls.phi33)); }
	if (readparam("phi_31:",buf,1)) {sscanf(buf,"%g",&(muls.phi31)); }
	if (readparam("phi_44:",buf,1)) {sscanf(buf,"%g",&(muls.phi44)); }
	if (readparam("phi_42:",buf,1)) {sscanf(buf,"%g",&(muls.phi42)); }
	if (readparam("phi_55:",buf,1)) {sscanf(buf,"%g",&(muls.phi55)); }
	if (readparam("phi_53:",buf,1)) {sscanf(buf,"%g",&(muls.phi53)); }
	if (readparam("phi_51:",buf,1)) {sscanf(buf,"%g",&(muls.phi51)); }
	if (readparam("phi_66:",buf,1)) {sscanf(buf,"%g",&(muls.phi66)); }
	if (readparam("phi_64:",buf,1)) {sscanf(buf,"%g",&(muls.phi64)); }
	if (readparam("phi_62:",buf,1)) {sscanf(buf,"%g",&(muls.phi62)); }

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


	if (!readparam("alpha:",buf,1)) exit(0); 
	sscanf(buf,"%g",&(muls.alpha)); /* in mrad */

	muls.aAIS = 0;  // initialize AIS aperture to 0 A
	if (readparam("AIS aperture:",buf,1)) 
		sscanf(buf,"%g",&(muls.aAIS)); /* in A */

	///// read beam current and dwell time ///////////////////////////////
	muls.beamCurrent = 1;  // pico Ampere
	muls.dwellTime = 1;    // msec
	if (readparam("beam current:",buf,1)) { 
		sscanf(buf,"%g",&(muls.beamCurrent)); /* in pA */
	}
	if (readparam("dwell time:",buf,1)) { 
		sscanf(buf,"%g",&(muls.dwellTime)); /* in msec */
	}
	muls.electronScale = static_cast<float_tt>(muls.beamCurrent*muls.dwellTime*MILLISEC_PICOAMP);
	//////////////////////////////////////////////////////////////////////

	muls.sourceRadius = 0;
	if (readparam("Source Size (diameter):",buf,1)) 
		muls.sourceRadius = static_cast<float_tt>(atof(buf)/2.0);

	if (readparam("smooth:",buf,1)) sscanf(buf,"%s",answer);
	muls.ismoth = (tolower(answer[0]) == (int)'y');
	muls.gaussScale = 0.05f;
	muls.gaussFlag = 0;
	if (readparam("gaussian:",buf,1)) {
		sscanf(buf,"%s %g",answer,&(muls.gaussScale));
		muls.gaussFlag = (tolower(answer[0]) == (int)'y');
	}

	/**********************************************************************
	* Parameters for image display and directories, etc.
	*/
	muls.imageGamma = 1.0;
	if (readparam("Display Gamma:",buf,1)) {
		muls.imageGamma = static_cast<float_tt>(atof(buf));
	}
	muls.showProbe = 0;
	if (readparam("show Probe:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.showProbe = (tolower(answer[0]) == (int)'y');
	}
	sprintf(muls.folder,"data");
	if (readparam("Folder:",buf,1)) 
		sscanf(buf," %s",muls.folder);

	while ((strptr = strchr(muls.folder,'"')) != NULL) {
		if (strptr == muls.folder) sscanf(strptr+1,"%s",muls.folder);	
		else *strptr = '\0';
	}

	if (muls.folder[strlen(muls.folder)-1] == '/')
		muls.folder[strlen(muls.folder)-1] = 0;
	muls.webUpdate = 0;
	if (readparam("update Web:",buf,1)) {
		sscanf(buf,"%s",answer);
		muls.webUpdate = (tolower(answer[0]) == (int)'y');
	}



	/*  readBeams(parFp); */  
	/************************************************************************/  
	/* read the different detector configurations                           */
	resetParamFile();
	muls.detectorNum = 0;

	if (muls.mode == STEM) 
	{
		int tCount = (int)(ceil((float_tt)((muls.slices * muls.cellDiv) / muls.outputInterval)));

		/* first determine number of detectors */
		while (readparam("detector:",buf,0)) muls.detectorNum++;  
		/* now read in the list of detectors: */
		resetParamFile();

		// loop over thickness planes where we're going to record intermediates
		// TODO: is this too costly in terms of memory?  It simplifies the parallelization to
		//       save each of the thicknesses in memory, then save to disk afterwards.
		for (int islice=0; islice<=tCount; islice++)
		{
			std::vector<DETECTOR> detectors;
			resetParamFile();
			while (readparam("detector:",buf,0)) {
				DETECTOR det;
				det.shiftX = 0.0;
				det.shiftY = 0.0;
				det.error = 0.0;
				det.Navg = 0;
				sscanf(buf,"%g %g %s %g %g",&(det.rInside),
					&(det.rOutside), det.name, &(det.shiftX),&(det.shiftY));  

				det.image = QSfMat(muls.scanXN,muls.scanYN);
				det.image2 = QSfMat(muls.scanXN,muls.scanYN);

				//for (ix=0;ix<muls.scanXN;ix++) for (iy=0;iy < muls.scanYN; iy++)
				//{
					//det.image(iy,ix) = 0.0;
					//det.image2(iy,ix) = 0.0;
				//}

				/* determine v0 specific k^2 values corresponding to the angles */
				det.k2Inside = 
					(float)(sin(det.rInside*0.001)/(wavelength(muls.v0)));
				det.k2Outside = 
					(float)(sin(det.rOutside*0.001)/(wavelength(muls.v0)));
				// printf("Detector %d: %f .. %f, lambda = %f (%f)\n",i,muls.detectors[i].k2Inside,muls.detectors[i].k2Outside,wavelength(muls.v0),muls.v0);
				/* calculate the squares of the ks */
				det.k2Inside *= det.k2Inside;
				det.k2Outside *= det.k2Outside;
				detectors.push_back(det);
			}
			muls.detectors.push_back(detectors);
		}
	}
	/************************************************************************/   

	// in case this file has been written by the tomography function, read the current tilt:
	if (readparam("tomo tilt:",buf,1)) { 
		sscanf(buf,"%lf %s",&(muls.tomoTilt),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.tomoTilt *= static_cast<float_tt>(1000*PI/180.0);
	}
	/************************************************************************
	* Tomography Parameters:
	***********************************************************************/
	if (muls.mode == TOMO) {     
		if (readparam("tomo start:",buf,1)) { 
			sscanf(buf,"%lf %s",&(muls.tomoStart),answer); /* in mrad */
			if (tolower(answer[0]) == 'd')
				muls.tomoStart *= static_cast<float_tt>(1000*PI/180.0);
		}
		if (readparam("tomo step:",buf,1)) {
			sscanf(buf,"%lf %s",&(muls.tomoStep),answer); /* in mrad */
			if (tolower(answer[0]) == 'd')
				muls.tomoStep *= static_cast<float_tt>(1000*PI/180.0);
		}

		if (readparam("tomo count:",buf,1))  
			muls.tomoCount = atoi(buf); 
		if (readparam("zoom factor:",buf,1))  
			sscanf(buf,"%lf",&(muls.zoomFactor));
		if ((muls.tomoStep == 0) && (muls.tomoStep > 1))
			muls.tomoStep = -2.0f*muls.tomoStart/(float_tt)(muls.tomoCount - 1);
	}
	/***********************************************************************/


	/*******************************************************************
	* Read in parameters related to the calculation of the projected
	* Potential
	*******************************************************************/
	muls.atomRadius = 5.0;
	if (readparam("atom radius:",buf,1))  
		sscanf(buf,"%g",&(muls.atomRadius)); /* in A */
	// why ??????  so that number of subdivisions per slice >= number of fitted points!!!
	/*  
	if (muls.atomRadius < muls.sliceThickness)
	muls.atomRadius = muls.sliceThickness;
	*/
	muls.scatFactor = DOYLE_TURNER;
	if (readparam("Structure Factors:",buf,1)) {
		sscanf(buf," %s",answer);
		switch (tolower(answer[0])) {
	case 'w':
		if (tolower(answer[1])=='k') muls.scatFactor = WEICK_KOHL;
		break;
	case 'd':  // DOYLE_TURNER
		muls.scatFactor = DOYLE_TURNER;
		break;
	case 'c':  // CUSTOM - specify k-lookup table and values for all atoms used
		muls.scatFactor = CUSTOM;
		// we already have the kinds of atoms stored in 
		// int *muls.Znums and int muls.atomKinds
		readSFactLUT();
		break;
	default:
		muls.scatFactor = DOYLE_TURNER;
		}
	}




	/***************************************************************
	* We now need to determine the size of the potential array, 
	* and the offset from its edges in A.  We only need to calculate
	* as much potential as we'll be illuminating later with the 
	* electron beam.
	**************************************************************/
	muls.potOffsetX = 0;
	muls.potOffsetY = 0;

	if ((muls.mode == STEM) || (muls.mode == CBED)) {
		/* we are assuming that there is enough atomic position data: */
		muls.potOffsetX = muls.scanXStart - 0.5f*muls.nx*muls.resolutionX;
		muls.potOffsetY = muls.scanYStart - 0.5f*muls.ny*muls.resolutionY;
		/* adjust scanStop so that it coincides with a full pixel: */
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
		muls.potOffsetX = muls.scanXStart - 0.5f*muls.potSizeX;
		muls.potOffsetY = muls.scanYStart - 0.5f*muls.potSizeY;
	}  
	/**************************************************************
	* Check to see if the given scan parameters really fit in cell 
	* dimensions:
	*************************************************************/
	if ((muls.scanXN <=0) ||(muls.scanYN <=0)) {
		printf("The number of scan pixels must be >=1\n");
		exit(0);
	}
	if ((muls.scanXStart<0) || (muls.scanYStart<0) ||
		(muls.scanXStop<0) || (muls.scanYStop<0) ||
		(muls.scanXStart>muls.cellDims[0]) || (muls.scanYStart>muls.cellDims[1]) ||
		(muls.scanXStop>muls.cellDims[0]) || (muls.scanYStop>muls.cellDims[1])) {
			printf("Scanning window is outside model dimensions (%g,%g .. %g,%g) [ax = %g, by = %g]!\n",muls.scanXStart,muls.scanYStart,muls.scanXStop,muls.scanYStop,muls.cellDims[0],muls.cellDims[1]);
			exit(0);
	}
	/*************************************************************
	* read in the beams we want to plot in the pendeloesung plot
	* Only possible if not in STEM or CBED mode 
	*************************************************************/
	muls.lbeams = 0;   /* flag for beam output */
	muls.nbout = 0;    /* number of beams */
	resetParamFile();
	if ((muls.mode != STEM) && (muls.mode != CBED)) {
		if (readparam("Pendelloesung plot:",buf,1)) {
			sscanf(buf,"%s",answer);
			muls.lbeams = (tolower(answer[0]) == (int)'y');
		}
		if (muls.lbeams) {
			while (readparam("beam:",buf,0)) muls.nbout++;  
			printf("will record %d beams\n",muls.nbout);
			muls.hbeam = QSiVec(muls.nbout);//(int*)malloc(muls.nbout*sizeof(int));
			muls.kbeam = QSiVec(muls.nbout);//(int*)malloc(muls.nbout*sizeof(int));
			/* now read in the list of detectors: */
			resetParamFile();
			for (i=0;i<muls.nbout;i++) {
				if (!readparam("beam:",buf,0)) break;
				muls.hbeam[i] = 0;
				muls.kbeam[i] = 0;
				sscanf(buf,"%d %d",muls.hbeam[i],muls.kbeam[i]);
				muls.hbeam[i] *= muls.nCellX;
				muls.kbeam[i] *= muls.nCellY;

				muls.hbeam[i] = (muls.hbeam[i]+muls.nx) % muls.nx;
				muls.kbeam[i] = (muls.kbeam[i]+muls.ny) % muls.ny;
				printf("beam %d [%d %d]\n",i,muls.hbeam[i],muls.kbeam[i]); 			}
		}
	}

	// muls.btiltx = 0;
	// muls.btilty = 0;
	//wave->thickness = 0.0;

	/* TODO: possible breakage here - MCS 2013/04 - made muls.cfgFile be allocated on the struct
	       at runtim - thus this null check doesn't make sense anymore.  Change cfgFile set
	   Old comment:
		if cfgFile != NULL, the program will later write a the atomic config to this file */
	//muls.cfgFile = NULL;
	if (readparam("CFG-file:",buf,1)) 
	{
		sscanf(buf,"%s",muls.cfgFile);
	}

	/* allocate memory for wave function */

	potDimensions[0] = muls.potNx;
	potDimensions[1] = muls.potNy;

	muls.trans = QSVecOfcMat(muls.slices);
	for (int slc=0; slc<muls.slices; slc++)
		muls.trans[slc]=QScMat(muls.potNx,muls.potNy);

#if FLOAT_PRECISION == 1
	// printf("allocated trans %d %d %d\n",muls.slices,muls.potNx,muls.potNy);
	muls.fftPlanPotForw = fftwf_plan_many_dft(2,potDimensions, muls.slices,(fftwf_complex*)muls.trans[0].data(), NULL,
		1, muls.potNx*muls.potNy,(fftwf_complex*)muls.trans[0].data(), NULL,
		1, muls.potNx*muls.potNy, FFTW_FORWARD, fftMeasureFlag);
	muls.fftPlanPotInv = fftwf_plan_many_dft(2,potDimensions, muls.slices,(fftwf_complex*)muls.trans[0].data(), NULL,
		1, muls.potNx*muls.potNy,(fftwf_complex*)muls.trans[0].data(), NULL,
		1, muls.potNx*muls.potNy, FFTW_BACKWARD, fftMeasureFlag);
#else
	muls.fftPlanPotForw = fftw_plan_many_dft(2,potDimensions, (fftw_complex*)muls.slices,muls.trans[0].data(), NULL,
		1, muls.potNx*muls.potNy,(fftw_complex*)muls.trans[0].data(), NULL,
		1, muls.potNx*muls.potNy, FFTW_FORWARD, fftMeasureFlag);
	muls.fftPlanPotInv = fftw_plan_many_dft(2,potDimensions, muls.slices,(fftw_complex*)muls.trans[0].data(), NULL,
		1, muls.potNx*muls.potNy,(fftw_complex*)muls.trans[0].data(), NULL,
		1, muls.potNx*muls.potNy, FFTW_BACKWARD, fftMeasureFlag);
#endif

	////////////////////////////////////
	if (muls.printLevel >= 4) 
		printf("Memory for transmission function (%d x %d x %d) allocated and plans initiated\n",muls.slices,muls.potNx,muls.potNy);


	// printf("%d %d %d %d\n",muls.nx,muls.ny,sizeof(fftw_complex),(int)(&muls.wave(2,2))-(int)(&muls.wave(1,2)));


} /* end of readFile() */



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
void doTOMO() {
	int ix,iy,iz,iTheta,i;
	QSf3Vec u, mpos, boxMin, boxMax;
	QSf3Mat Mm;
	float_tt theta = 0;
	std::vector<atom> atoms;
	char cfgFile[64],stemFile[128],scriptFile[64],diffAnimFile[64];
	FILE *fpScript,*fpDiffAnim;

	Mm = muls.Mm;
	//atoms = (atom *)malloc(muls.natom*sizeof(atom));

	boxMin << muls.cellDims[0]/2.0f, muls.cellDims[1]/2.0f, muls.cellDims[2]/2.0f;;
	boxMax = boxMin;

	// For all tomography tilt angles:
	for (iTheta = 0;iTheta < muls.tomoCount;iTheta++) {
		theta = muls.tomoStart + iTheta*muls.tomoStep; // theta in mrad
		// Try different corners of the box, and see, how far they poke out.
		for (ix=-1;ix<=1;ix++) for (iy=-1;iy<=1;iy++) for (iz=-1;iz<=1;iz++) {
			// Make center of unit cell rotation center
			u << ix*muls.cellDims[0]/2.0f, iy*muls.cellDims[1]/2.0f, iz*muls.cellDims[2]/2.0f;;

			// rotate about y-axis
			rotateVect(u,u,0,theta*1e-3f,0);

			// shift origin back to old (0,0,0):
			u.array() += muls.cellDims/2;
			//u[0]+=muls.cellDims[0]/2; u[1]+=muls.cellDims[1]/2.0f; u[2]+=muls.cellDims[2]/2.0f;

			boxMin[0] = boxMin[0]>u[0] ? u[0] : boxMin[0]; boxMax[0] = boxMax[0]<u[0] ? u[0] : boxMax[0]; 
			boxMin[1] = boxMin[1]>u[1] ? u[1] : boxMin[1]; boxMax[1] = boxMax[1]<u[1] ? u[1] : boxMax[1]; 
			boxMin[2] = boxMin[2]>u[2] ? u[2] : boxMin[2]; boxMax[2] = boxMax[2]<u[2] ? u[2] : boxMax[2]; 

		}
	} /* for iTheta ... */

	// find max. box size:
	boxMax -= boxMin;
	//boxMax[0] -= boxMin[0];
	//boxMax[1] -= boxMin[1];
	//boxMax[2] -= boxMin[2];
	printf("Minimum box size for tomography tilt series: %g x %g x %gA, zoom Factor: %g\n",
		boxMax[0],boxMax[1],boxMax[2],muls.zoomFactor);
	boxMax[0] /= muls.zoomFactor;
	boxMax[1] = boxMax[0]*muls.cellDims[1]/muls.cellDims[0];

	// boxMin will now be boxCenter:
	boxMin = 0.5*boxMax;
	//boxMin[0] = 0.5*boxMax[0];
	//boxMin[1] = 0.5*boxMax[1];
	//boxMin[2] = 0.5*boxMax[2];

	// We have to save the original unit cell dimensions
	mpos<< muls.cellDims[0], muls.cellDims[1], muls.cellDims[2];;
	muls.cellDims[0]=boxMax[0]; muls.cellDims[1]=boxMax[1]; muls.cellDims[2]=boxMax[2];

	printf("Will use box sizes: %g x %g x %gA (kept original aspect ratio). \n"
		"Writing structure files now, please wait ...\n",
		boxMax[0],boxMax[1],boxMax[2]);

	// open the script file and write the stem instructions in there
	sprintf(scriptFile,"%s/run_tomo",muls.folder);
	if ((fpScript=fopen(scriptFile,"w")) == NULL) {
		printf("doTOMO: unable to open scriptFile %s for writing\n",scriptFile);
		exit(0);
	}
	fprintf(fpScript,"#!/bin/bash\n\n");

	// open the diffraction animation file 
	sprintf(diffAnimFile,"%s/diff_anim",muls.folder);
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
			u = (muls.atoms[i].pos - mpos.array()/2.0).matrix(); 
			//u[1] = muls.atoms[i].y - mBy/2.0; 
			//u[2] = muls.atoms[i].z - mCz/2.0; 
			rotateVect(u,u,0,theta*1e-3f,0);
			atoms[i].pos = u+boxMin;
			//atoms[i].x = u[0]+boxMin[0];
			//atoms[i].y = u[1]+boxMin[1]; 
			//atoms[i].z = u[2]+boxMin[2]; 
			atoms[i].Znum = muls.atoms[i].Znum;
			atoms[i].occ = muls.atoms[i].occ;
			atoms[i].dw = muls.atoms[i].dw;
		}

		sprintf(cfgFile,"%s/tomo_%dmrad.cfg",muls.folder,(int)theta);
		sprintf(stemFile,"%s/tomo_%dmrad.dat",muls.folder,(int)theta);
		printf("Writing file %s | ",cfgFile);
		writeCFG(atoms,muls.natom,cfgFile,&muls);
		sprintf(cfgFile,"tomo_%dmrad.cfg",(int)theta);
		writeSTEMinput(stemFile,cfgFile,&muls);

		// add to script files:
		fprintf(fpScript,"stem tomo_%dmrad.dat\n",(int)theta);
		// MCS - commented 07/2010 to limit errors on linux console.
		//     Should be replaced if Christoph's webpage display comes back.
		//fprintf(fpDiffAnim,"showimage tomo_%dmrad/diffAvg_%d.img 6\n",(int)theta,muls.avgRuns);
		//fprintf(fpDiffAnim,"convert -crop 75%%x75%%+250%%+250%% tomo_%dmrad/diffAvg_%d.tif"
		//	" -resize 300x300 -font helvetica -fill white -draw 'text 250,290 \"%4dmrad\"' "
		//	"./diff%03d.jpg\n\n",
		//	(int)theta,muls.avgRuns,(int)theta,iTheta);
	}
	muls.cellDims = mpos;
	//muls.cellDims[0] = mAx; muls.cellDims[1] = mBy; muls.cellDims[2] = mCz;
	sprintf(stemFile,"copy fparams.dat %s/",muls.folder);
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
void doMSCBED() {


}

/************************************************************************
* doCBED performs a CBED calculation
*
***********************************************************************/

void doCBED() {
	int ix,iy,i,pCount,result;
	FILE *avgFp,*fp,*fpPos=0;
	float_tt timer,timerTot;
	float_tt probeCenterX,probeCenterY,probeOffsetX,probeOffsetY;
	char buf[BUF_LEN],avgName[32],systStr[64];
	float_tt t;
	static int oldMulsRepeat1 = 1;
	static int oldMulsRepeat2 = 1;
	static long iseed=0;
	WAVEFUNC *wave = new WAVEFUNC(muls.nx,muls.ny);
	static imageStruct *header = NULL;
	imageStruct *header_read = NULL;

	QSfMat avgArray(muls.nx,muls.ny), diffArray(muls.nx,muls.ny);
	QSfVec chisq(muls.avgRuns);
	QSfMat avgPendelloesung(muls.nbout,
			muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv);

	if (iseed == 0) iseed = -(long) time( NULL );

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

		probeOffsetX = static_cast<float_tt>(muls.sourceRadius*gasdev(&iseed)*SQRT_2);
		probeOffsetY = static_cast<float_tt>(muls.sourceRadius*gasdev(&iseed)*SQRT_2);
		muls.scanXStart = probeCenterX+probeOffsetX;
		muls.scanYStart = probeCenterY+probeOffsetY;
		probe(&muls, wave,muls.scanXStart-muls.potOffsetX,muls.scanYStart-muls.potOffsetY);
		if (muls.saveLevel > 2) {
			if (header == NULL) 
				header = makeNewHeaderCompact(1,muls.nx,muls.ny,wave->thickness,
				muls.resolutionX,muls.resolutionY,
				0,std::vector<float_tt>(),"wave function");
			header->t = 0;
			sprintf(systStr,"%s/wave_probe.img",muls.folder);
			writeComplexImage(wave->wave,header,systStr);
			// writeImage(muls.wave,header,"wave->img");
			// writeImage_old(muls.wave,muls.nx,muls.ny,wave->thickness,"wave->img");
			// system("showimage diff.img 2 &");
		} 	
		// printf("Probe: (%g, %g)\n",muls.scanXStart,muls.scanYStart);
		/*****************************************************************
		* For debugging only!!!
		*
		if (header == NULL) 
		header = makeNewHeaderCompact(1,muls.nx,muls.ny,wave->thickness,
		muls.resolutionX,muls.resolutionY,
		0,NULL,"probe function");
		header->t = wave->thickness;
		writeImage(muls.wave,header,"probe.img");
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
			sprintf(buf,"ee %s/probePlot_0.jpg &",muls.folder);
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
			muls.cin2 = muls.mulsRepeat1;
			//sprintf(muls.cin2,"%d",muls.mulsRepeat1);


			/***********************************************************
			* make sure we have enough memory for the pendelloesung plot
			*/
			if ((muls.lbeams) && 
				((oldMulsRepeat1 !=muls.mulsRepeat1) ||
				(oldMulsRepeat2 !=muls.mulsRepeat2))) {
					oldMulsRepeat1 = muls.mulsRepeat1;
					oldMulsRepeat2 = muls.mulsRepeat2;
					muls.pendelloesung.setZero();
					avgPendelloesung.setZero();
			}
			/*********************************************************/

			// printf("Stacking sequence: %s\n",buf);


			/************************************************
			*make3DSlicesFFT(&muls,muls.slices,atomPosFile,NULL);
			*exit(0);
			************************************************/
			if (muls.equalDivs) {
				if (muls.scatFactor == CUSTOM)
					make3DSlicesFT(&muls);
				else
					make3DSlices(&muls, muls.slices, muls.atomPosFile, NULL);
				initSTEMSlices(&muls,muls.slices);
			}

			muls.saveFlag = 0;
			/****************************************
			* do the (small) loop
			*****************************************/
			for (pCount = 0;pCount<muls.mulsRepeat2*muls.cellDiv;pCount++) {

				// printf("pCount = %d\n",pCount);
				/*******************************************************
				* build the potential slices from atomic configuration
				******************************************************/
				if (!muls.equalDivs) {
					if (muls.scatFactor == CUSTOM)
						make3DSlicesFT(&muls);
					else
						make3DSlices(&muls, muls.slices, muls.atomPosFile, NULL);
					initSTEMSlices(&muls, muls.slices);
				}

				timer = cputim();
				// what probe should runMulsSTEM use here?
				runMulsSTEM(&muls,wave); 

				printf("Thickness: %gA, int.=%g, time: %gsec\n",
					wave->thickness,wave->intIntensity,cputim()-timer);

				/***************** Only if Save level > 2: ****************/
				if ((muls.avgCount == 0) && (muls.saveLevel > 2)) {
					if (header == NULL) 
						header = makeNewHeaderCompact(1,muls.nx,muls.ny,wave->thickness,
						muls.resolutionX,muls.resolutionY,
						0,std::vector<float_tt>(),"wave function");
					header->t = wave->thickness;
					sprintf(systStr,"%s/wave_final.img",muls.folder);
					writeComplexImage(wave->wave,header,systStr);
					// writeImage(muls.wave,header,"wave.img");
					// writeImage_old(muls.wave,muls.nx,muls.ny,wave->thickness,"wave.img");
					// system("showimage diff.img 2 &");
				} 	
#ifdef VIB_IMAGE_TEST_CBED
				if (header == NULL) 
					header = makeNewHeaderCompact(1,muls.nx,muls.ny,wave->thickness,
					muls.resolutionX,muls.resolutionY,
					1,&(muls.tomoTilt),"wave function");
				header->t = wave->thickness;
				sprintf(systStr,"%s/wave_%d.img",muls.folder,muls.avgCount);
				writeImage(muls.wave,header,systStr);
				// writeImage_old(muls.wave,muls.nx,muls.ny,wave->thickness,systStr);
#endif 
				muls.totalSliceCount += muls.slices;

			} // end of for pCount = 0... 
			result = readparam("sequence: ",buf,0);
		}
		/*    printf("Total CPU time = %f sec.\n", cputim()-timerTot ); */

		sprintf(avgName,"%s/diff.img",muls.folder);
		// readRealImage_old(diffArray,muls.nx,muls.ny,&t,avgName);
		header_read = readImage(diffArray,muls.nx,muls.ny,avgName);


		if (muls.avgCount == 0) {
			/* for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++)
			avgArray(iy,ix) = diffArray(iy,ix); */
			avgArray = diffArray;
			// writeRealImage_old(avgArray,muls.nx,muls.ny,wave->thickness,avgName);
			/* move the averaged (raw data) file to the target directory as well */
			sprintf(avgName,"%s/diffAvg_%d.img",muls.folder,muls.avgCount+1);
			sprintf(systStr,"mv %s/diff.img %s",muls.folder,avgName);
			system(systStr);
			if (muls.lbeams) {
				for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
					for (ix=0;ix<muls.nbout;ix++) {
						avgPendelloesung(iy,ix) = muls.pendelloesung(iy,ix);
					}
				}
			}
		} // of if muls.avgCount == 0 ...
		else {
			/*      readRealImage_old(avgArray,muls.nx,muls.ny,&t,"diffAvg.img"); */
			chisq[muls.avgCount-1] = 0.0;
			for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
				t = ((float_tt)muls.avgCount*avgArray(iy,ix)+
					diffArray(iy,ix))/((float_tt)(muls.avgCount+1));
				chisq[muls.avgCount-1] += (avgArray(iy,ix)-t)*(avgArray(iy,ix)-t);
				avgArray(iy,ix) = t;

			}
			chisq[muls.avgCount-1] = chisq[muls.avgCount-1]/(float_tt)(muls.nx*muls.ny);
			sprintf(avgName,"%s/diffAvg_%d.img",muls.folder,muls.avgCount+1);
			// writeRealImage_old(avgArray,muls.nx,muls.ny,wave->thickness,avgName);
			if (header == NULL) 
				header = makeNewHeaderCompact(0,muls.nx,muls.ny,wave->thickness,
				1.0/(muls.nx*muls.resolutionX),1.0/(muls.ny*muls.resolutionY),
				1,std::vector<float_tt>(1,muls.tomoTilt),"Averaged Diffraction pattern, unit: 1/A");
			else {
				header->t = wave->thickness;
				header->dx = 1.0f/(muls.nx*muls.resolutionX);
				header->dy = 1.0f/(muls.ny*muls.resolutionY);
				if (header->paramSize < 1) {
					header->params = std::vector<float_tt>(2);
					header->paramSize = 2;
				}

				header->params[0] = muls.tomoTilt;
				header->params[1] = 1.0f/wavelength(muls.v0);
				setHeaderComment(header,"Averaged Diffraction pattern, unit: 1/A");
			}
			writeRealImage(avgArray,header,avgName);

			/* report the result on the web page */
			// printf("Will write report now\n");
			// add2Webpage(avgName); 

			muls.storeSeries = 1;
			if (muls.saveLevel == 0)	muls.storeSeries = 0;
			else if (muls.avgCount % muls.saveLevel != 0) muls.storeSeries = 0;

			if (muls.storeSeries == 0) {
				// printf("Removing old file \n");
				sprintf(avgName,"%s/diffAvg_%d.img",muls.folder,muls.avgCount);
				sprintf(systStr,"rm %s",avgName);
				system(systStr);
			}

			/* write the data to a file */
			if (muls.saveFlag >-1) {
				sprintf(systStr,"%s/avgresults.dat",muls.folder);
				if ((avgFp = fopen(systStr,"w")) == NULL )
					printf("Sorry, could not open data file for averaging\n");
				else {
					for (ix =0;ix<muls.avgCount;ix++) {
						fprintf(avgFp,"%d %g\n",ix+1,chisq[ix]);
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
						avgPendelloesung(iy,ix) = 
							((float_tt)muls.avgCount*avgPendelloesung(iy,ix)+
							muls.pendelloesung(iy,ix))/(float_tt)(muls.avgCount+1);
					}
				}
			}
		} /* else ... if avgCount was greater than 0 */

		if (muls.lbeams) {
			/**************************************************************
			* The diffraction spot intensities of the selected 
			* diffraction spots are now stored in the 2 dimensional array
			* muls.pendelloesung(slice,beam).
			* We can write the array to a file and display it, just for 
			* demonstration purposes
			*************************************************************/
			sprintf(systStr,"%s/pendelloesung.dat",muls.folder);
			if ((fp=fopen(systStr,"w")) !=NULL) {
				printf("Writing Pendelloesung data\n");
				for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
					/* write the thicknes in the first column of the file */
					fprintf(fp,"%g",iy*muls.cellDims[2]/((float)(muls.slices*muls.cellDiv)));
					/* write the beam intensities in the following columns */
					for (ix=0;ix<muls.nbout;ix++) {
						fprintf(fp,"\t%g",avgPendelloesung(iy,ix));
					}
					/* close the line, and start a new one for the next set of
					* intensities
					*/
					fprintf(fp,"\n");
				}
				fclose(fp);
				/* plot the data in the file */
				/* system("xmgr -nxy pendelloesung.dat &"); */
			}
			else {
				printf("Could not open file for pendelloesung plot\n");
			}  
		} /* end of if lbemas ... */
		displayProgress(1);
	} /* end of for muls.avgCount=0.. */
	delete(wave);
}
/************************************************************************
* End of doCBED()
***********************************************************************/

/************************************************************************
* doTEM performs a TEM calculation (with through focus reconstruction)
*
***********************************************************************/

void doTEM() {
	int ix,iy,i,pCount,result;
	FILE *avgFp,*fp; // *fpPos=0;
	float_tt timer,timerTot;
	float_tt x,y,ktx,kty;
	char buf[BUF_LEN],avgName[256],systStr[512];
	float_tt t;
	QSfMat avgArray(muls.nx,muls.ny), diffArray(muls.nx,muls.ny);
	QSfVec chisq;
	int oldMulsRepeat1 = 1;
	int oldMulsRepeat2 = 1;
	long iseed=0;
	imageStruct *header = NULL;
	imageStruct *header_read = NULL;
	WAVEFUNC *wave = new WAVEFUNC(muls.nx,muls.ny);
	QScMat imageWave(muls.nx,muls.ny);
	QSfMat avgPendelloesung(muls.nbout,
			muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv);

	if (iseed == 0) iseed = -(long) time( NULL );

	muls.chisq = QSfVec(muls.avgRuns);
	// muls.trans = 0;

	if (muls.lbeams) {
		muls.pendelloesung.setZero();
	}

	timerTot = 0; /* cputim();*/
	displayProgress(-1);
	for (muls.avgCount = 0;muls.avgCount < muls.avgRuns;muls.avgCount++) {
		muls.totalSliceCount = 0;

		pCount = 0;
		resetParamFile();

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
			wave->wave.setConstant(std::complex<float_tt>(1,0));
		}
		else {
			// produce a tilted wave function (btiltx,btilty):
			ktx = static_cast<float_tt>(2.0*PI*sin(muls.btiltx)/wavelength(muls.v0));
			kty = static_cast<float_tt>(2.0*PI*sin(muls.btilty)/wavelength(muls.v0));
			for (ix=0;ix<muls.nx;ix++) {
				x = muls.resolutionX*(ix-muls.nx/2);
				for (iy=0;iy<muls.ny;iy++) {
					y = muls.resolutionY*(ix-muls.nx/2);
					wave->wave(iy,ix) = std::complex<float_tt>(cos(ktx*x+kty*y),
															   sin(ktx*x+kty*y));
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
			//sprintf(muls.cin2,"%d",muls.mulsRepeat1);
			muls.cin2 = muls.mulsRepeat1;
			/***********************************************************
			* make sure we have enough memory for the pendelloesung plot
			*/
			if ((muls.lbeams) && 
				((oldMulsRepeat1 !=muls.mulsRepeat1) ||
				(oldMulsRepeat2 !=muls.mulsRepeat2))) {
					oldMulsRepeat1 = muls.mulsRepeat1;
					oldMulsRepeat2 = muls.mulsRepeat2;
					muls.pendelloesung.setZero();
			}
			/*********************************************************/

			// printf("Stacking sequence: %s\n",buf);


			/************************************************
			*make3DSlicesFFT(&muls,muls.slices,atomPosFile,NULL);
			*exit(0);
			************************************************/
			if (muls.equalDivs) {
				if (muls.printLevel > 1) printf("found equal unit cell divisions\n");
				make3DSlices(&muls, muls.slices, muls.atomPosFile, NULL);
				initSTEMSlices(&muls, muls.slices);
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
					make3DSlices(&muls, muls.slices, muls.atomPosFile, NULL);
					initSTEMSlices(&muls,muls.slices);
				}

				timer = cputim();
				runMulsSTEM(&muls,wave); 
				muls.totalSliceCount += muls.slices;

				if (muls.printLevel > 0) {
					printf("t=%gA, int.=%g time: %gsec (avgCount=%d)\n",
						wave->thickness,wave->intIntensity,cputim()-timer,muls.avgCount);
				}

				/***************** FOR DEBUGGING ****************/		
				if ((muls.avgCount == 0) && (muls.saveLevel >=0) && (pCount+1==muls.mulsRepeat2*muls.cellDiv)) {
					// writeImage_old(muls.wave,muls.nx,muls.ny,wave->thickness,"wave.img");
					if (header == NULL) header = makeNewHeader(muls.nx,muls.ny);
					header->t = wave->thickness;
					header->dx = muls.resolutionX;
					header->dy = muls.resolutionY;
					if (muls.tds) setHeaderComment(header,"Test wave function for run 0");
					else setHeaderComment(header,"Exit face wave function for no TDS");
					sprintf(systStr,"%s/wave.img",muls.folder);
					if ((muls.tiltBack) && ((muls.btiltx != 0) || (muls.btilty != 0))) {
						ktx = static_cast<float_tt>(-2.0*PI*sin(muls.btiltx)/wavelength(muls.v0));
						kty = static_cast<float_tt>(-2.0*PI*sin(muls.btilty)/wavelength(muls.v0));
						for (ix=0;ix<muls.nx;ix++) {
							x = muls.resolutionX*(ix-muls.nx/2);
							for (iy=0;iy<muls.ny;iy++) {
								y = muls.resolutionY*(ix-muls.nx/2);
								wave->wave(iy,ix) *= QScf(cos(ktx*x+kty*y), sin(ktx*x+kty*y));
							}
						}
						if (muls.printLevel > 1) printf("** Applied beam tilt compensation **\n");
					}


					writeComplexImage(wave->wave,header,systStr);
					//    system("showimage diff.img 2 &");
				}	
#ifdef VIB_IMAGE_TEST  // doTEM
				if ((muls.tds) && (muls.saveLevel > 2)) {
					sprintf(systStr,"%s/wave_%d.img",muls.folder,muls.avgCount);
					// writeImage_old(muls.wave,muls.nx,muls.ny,wave->thickness,systStr);
					if (header == NULL) header = makeNewHeader(muls.nx,muls.ny);
					header->t = wave->thickness;
					header->dx = muls.resolutionX;
					header->dy = muls.resolutionY;
					header->complexFlag = 1;
					header->paramSize = 9;
					header->params = std::vector<float_tt>(9);
					header->params[0] = muls.v0;  				// high voltage
					header->params[1] = muls.Cs;				// spherical aberration
					header->params[2] = muls.df0;				// defocus
					header->params[3] = muls.astigMag;			// astigmatism
					header->params[4] = muls.astigAngle;	
					header->params[5] = muls.Cc * sqrt(muls.dE_E*muls.dE_E+muls.dV_V*muls.dV_V+muls.dI_I*muls.dI_I);	// focal spread
					// printf("****  Cc = %f, dE_E = %f, Delta = %f ****\n",muls.Cc,muls.dV_V,muls.Cc * muls.dV_V);
					header->params[6] = muls.alpha;				// illumination convergence angle
					// beam tilt:
					header->params[7] = muls.btiltx;			// beam tilt in mrad
					header->params[8] = muls.btilty;			// beam tilt in mrad


					setHeaderComment(header,"complex exit face Wave function");
					writeComplexImage(wave->wave,header,systStr);
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

		// 	sprintf(avgName,"%s/diffAvg.img",muls.folder);
		sprintf(avgName,"%s/diff.img",muls.folder);

		//readRealImage_old(diffArray,muls.nx,muls.ny,&t,avgName);
		header_read = readImage(diffArray,muls.nx,muls.ny,avgName);

		if (muls.avgCount == 0) {
			/***********************************************************
			* Save the diffraction pattern
			**********************************************************/	
			avgArray = diffArray;
			//memcpy((void *)avgArray[0],(void *)diffArray[0],(size_t)(muls.nx*muls.ny*sizeof(float_tt)));
			// writeRealImage_old(avgArray,muls.nx,muls.ny,wave->thickness,avgName);
			/* move the averaged (raw data) file to the target directory as well */
#ifndef WIN32
			sprintf(avgName,"diffAvg_%d.img",muls.avgCount+1);
			sprintf(systStr,"mv %s/diff.img %s/%s",muls.folder,muls.folder,avgName);
			system(systStr);
#else
			sprintf(avgName,"diffAvg_%d.img",muls.avgCount+1);
			sprintf(systStr,"move %s\\diff.img %s/%s",muls.folder,muls.folder,avgName);
			// system(systStr);
#endif
			/***********************************************************
			* Save the Pendelloesung Plot
			**********************************************************/	
			if (muls.lbeams) {
				for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
					for (ix=0;ix<muls.nbout;ix++) {
						avgPendelloesung(iy,ix) = muls.pendelloesung(iy,ix);
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
			// multiply wave (in rec. space) with transfer function and write result to imagewave
			fftwf_execute(wave->fftPlanWaveForw);
			for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
				imageWave(iy,ix) = wave->wave(iy,ix);
			}
			fftwf_execute_dft(wave->fftPlanWaveInv,(fftwf_complex*)imageWave.data(),(fftwf_complex*)imageWave.data());
			// get the amplitude squared:
			for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
				diffArray(iy,ix) = imageWave(iy,ix).real()*imageWave(iy,ix).real()+imageWave(iy,ix).imag()*imageWave(iy,ix).imag();
			}
			if (header == NULL) 
				header = makeNewHeaderCompact(0,muls.nx,muls.ny,wave->thickness,
				muls.resolutionX,muls.resolutionY,
				0,std::vector<float_tt>(),"Image");
			header->t = wave->thickness;
			setHeaderComment(header,"Image intensity");
			sprintf(avgName,"%s/image.img",muls.folder);
			writeRealImage(diffArray,header,avgName);
			// End of Image writing (if avgCount = 0)
			//////////////////////////////////////////////////////////////////////

		} // of if muls.avgCount == 0 ...
		else {
			/* 	 readRealImage_old(avgArray,muls.nx,muls.ny,&t,"diffAvg.img"); */
			chisq[muls.avgCount-1] = 0.0;
			for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
				t = ((float_tt)muls.avgCount*avgArray(iy,ix)+
					diffArray(iy,ix))/((float_tt)(muls.avgCount+1));
				chisq[muls.avgCount-1] += (avgArray(iy,ix)-t)*(avgArray(iy,ix)-t);
				avgArray(iy,ix) = t;
			}
			chisq[muls.avgCount-1] = chisq[muls.avgCount-1]/(float_tt)(muls.nx*muls.ny);
			sprintf(avgName,"%s/diffAvg_%d.img",muls.folder,muls.avgCount+1);
			// writeRealImage_old(avgArray,muls.nx,muls.ny,wave->thickness,avgName);
			if (header == NULL) 
				header = makeNewHeaderCompact(0,muls.nx,muls.ny,wave->thickness,
				muls.resolutionX,muls.resolutionY,
				0,std::vector<float_tt>(),"diffraction pattern");
			header->t = wave->thickness;
			writeRealImage(avgArray,header,avgName);


			/* report the result on the web page */
			// printf("Will write report now\n");
			// add2Webpage(avgName); 

			/* move the averaged (raw data) file to the target directory as well */
			// 	  sprintf(systStr,"mv %s %s/%s",avgName,muls.folder,avgName);
			// system(systStr);

			/* write the data to a file */
			if ((avgFp = fopen("avgresults.dat","w")) == NULL )
				printf("Sorry, could not open data file for averaging\n");
			else {
				for (ix =0;ix<muls.avgCount;ix++) {
					fprintf(avgFp,"%d %g\n",ix+1,chisq[ix]);
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
						avgPendelloesung(iy,ix) = 
							((float_tt)muls.avgCount*avgPendelloesung(iy,ix)+
							muls.pendelloesung(iy,ix))/(float_tt)(muls.avgCount+1);
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
			// multiply wave (in rec. space) with transfer function and write result to imagewave
			fftwf_execute(wave->fftPlanWaveForw);
			for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
				imageWave(iy,ix) = wave->wave(iy,ix);
			}
			fftwf_execute_dft(wave->fftPlanWaveInv,(fftwf_complex*)imageWave.data(),(fftwf_complex*)imageWave.data());

			// save the amplitude squared:
			sprintf(avgName,"%s/image.img",muls.folder); 
			header_read = readImage(diffArray,muls.nx,muls.ny,avgName);
			for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
				t = ((float_tt)muls.avgCount*diffArray(iy,ix)+
					imageWave(iy,ix).real()*imageWave(iy,ix).real()+imageWave(iy,ix).imag()*imageWave(iy,ix).imag())/(float_tt)(muls.avgCount+1);
				diffArray(iy,ix) = t;
			}
			header->t = wave->thickness;
			setHeaderComment(header,"Image intensity");
			writeRealImage(diffArray,header,avgName);
			// End of Image writing (if avgCount > 0)
			//////////////////////////////////////////////////////////////////////

		} /* else ... if avgCount was greater than 0 */


		/////////////////////////////////////////////////////
		// Save the Pendelloesung plot:
		if (muls.lbeams) {
			/**************************************************************
			* The diffraction spot intensities of the selected 
			* diffraction spots are now stored in the 2 dimensional array
			* muls.pendelloesung(slice,beam).
			* We can write the array to a file and display it, just for 
			* demonstration purposes
			*************************************************************/
			sprintf(avgName,"%s/pendelloesung.dat",muls.folder);
			if ((fp=fopen(avgName,"w")) !=NULL) {
				printf("Writing Pendelloesung data\n");
				for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
					/* write the thicknes in the first column of the file */
					fprintf(fp,"%g",iy*muls.cellDims[2]/((float)(muls.slices*muls.cellDiv)));
					/* write the beam intensities in the following columns */
					for (ix=0;ix<muls.nbout;ix++) {
						// store the AMPLITUDE:
						fprintf(fp,"\t%g",sqrt(avgPendelloesung(iy,ix)/(muls.nx*muls.ny)));
					}
					/* close the line, and start a new one for the next set of
					* intensities
					*/
					fprintf(fp,"\n");
				}
				fclose(fp);
				/* plot the data in the file */
				/* system("xmgr -nxy pendelloesung.dat &"); */
			}
			else {
				printf("Could not open file for pendelloesung plot\n");
			}	
		} /* end of if lbeams ... */		 
		displayProgress(1);
	} /* end of for muls.avgCount=0.. */  
	delete(wave);
}
/************************************************************************
* end of doTEM
***********************************************************************/



/************************************************************************
* doSTEM performs a STEM calculation
*
***********************************************************************/

void doSTEM() {
	int ix=0,iy=0,i,pCount,picts,ixa,iya,totalRuns;
	float_tt timer, total_time=0;
	char buf[BUF_LEN];//,avgName[256],systStr[256];
	std::string bufstr;
	//char jpgName[256], tifName[256];
	float_tt t;
	static QSfMat avgArray;
	QSfVec chisq;
	float_tt collectedIntensity;
	static imageStruct *header = NULL;
	static imageStruct *header_read = NULL;
	//float cztot;
	//int islice;

	//waves = (WAVEFUNC *)malloc(muls.scanYN*muls.scanXN*sizeof(WAVEFUNC));
	//WAVEFUNC *wave = *(new WAVEFUNC(muls.nx,muls.ny));
	std::vector<WAVEFUNC *> waves;
	WAVEFUNC *wave;

	for (int th=0; th<omp_get_max_threads(); th++)
	{
		waves.push_back(new WAVEFUNC(muls.nx,muls.ny));
	}

	//chisq = (float_tt *)malloc(muls.avgRuns*sizeof(float_tt));
	// zero-out the chisq array
	//memset(chisq, 0, muls.avgRuns*sizeof(float_tt));
	chisq.setZero();
	muls.chisq = chisq;
	totalRuns = muls.avgRuns;
	timer = cputim();

        //pre-allocate several waves (enough for one row of the scan.  
        //This fixes a memory leak bug.  What was causing it:
        /* 
           - OpenMP required that each thread needs its own FFTW plan
           - To make an FFTW plan for each thread, we had to allocated a wave struct per thread
           - Allocating so many wave structs made the heap too big, causing the leak (I think)
        */
        //for(iy=0;iy<muls.scanYN*muls.scanXN;iy++) waves[iy]=initWave(muls.nx,muls.ny);

	//if (avgArray == NULL)
	//	avgArray = float2D(muls.nx,muls.ny,"avgArray");
	/* create the target folder, if it does not exist yet */
	/*  sprintf(avgName,"diffAvg_%d_%d.img",ix,iy);
	sprintf(tifName,"diffAvg_%d_%d.tif",ix,iy);
	sprintf(jpgName,"diffAvg_%d_%d.jpg",ix,iy);

	sprintf(systStr,"showimage %s 6; convert -crop 17%c -geometry 256x256 %s %s/%s; rm %s",
	avgName,'%',tifName,muls.folder,jpgName,tifName);
	printf("%s\n",systStr);
	system(systStr);
	*/
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
			muls.cin2 = muls.mulsRepeat1;
			//sprintf(muls.cin2,"%d",muls.mulsRepeat1);
			/* if the unit cell is divided into slabs, we need to multiply
			* picts by that number
			*/
			if ((picts > 1)&& (muls.cubex >0) && (muls.cubey >0) && (muls.cubez>0)) {
				printf("Warning: cube size of height %gA has been defined, ignoring sequence\n",muls.cubez);
				picts = 1;
			}
			picts *= muls.cellDiv;

			if (muls.equalDivs) {
				make3DSlices(&muls, muls.slices, muls.atomPosFile, NULL);
				initSTEMSlices(&muls, muls.slices);
				timer = cputim();
			}

			//MCS - initialize the probe wavefunction that will be copied to each thread.
			/* make incident probe wave function with probe exactly in the center */
						/* if the potential array is not big enough, the probe can 
						* then also be adjusted, so that it is off-center
						*/
			//muls.nslic0 = 0;
			//wave->thickness = 0.0;

			/****************************************
			* do the (small) loop over slabs
			*****************************************/
			for (pCount=0;pCount<picts;pCount++) {
				/*******************************************************
				* build the potential slices from atomic configuration
				******************************************************/
				if (!muls.equalDivs) {
					make3DSlices(&muls,muls.slices,muls.atomPosFile,NULL);
					initSTEMSlices(&muls,muls.slices);
					timer = cputim();
				}

				muls.complete_pixels=0;
				/**************************************************
				* scan through the different probe positions
				*************************************************/
				// default(none) forces us to specify all of the variables that are used in the parallel section.  
				//    Otherwise, they are implicitly shared (and this was cause of several bugs.)
#pragma omp parallel firstprivate(header, header_read) \
	private(ix, iy, ixa, iya, wave, t, timer, buf, std::string) \
	shared(pCount, picts, chisq, muls, collectedIntensity, total_time, waves) \
	default(none)
#pragma omp for
				for (i=0; i < (muls.scanXN * muls.scanYN); i++)
				{
					timer=cputim();
					ix = i / muls.scanYN;
					iy = i % muls.scanYN;

					wave = waves[omp_get_thread_num()];
							
					//printf("Scanning: %d %d %d %d\n",ix,iy,pCount,muls.nx);

					//wave = waves[i];
					//xpos = muls.scanXStart+
                                  //ix*(muls.scanXStop-muls.scanXStart)/(float)muls.scanXN;
					//ypos = muls.scanYStart+
                                  //iy*(muls.scanYStop-muls.scanYStart)/(float)muls.scanYN;
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
						sprintf(buf, "%s/mulswav_%d_%d.img", muls.folder, ix, iy);
						wave->fileStart = std::string(buf);
						readStartWave(&muls, wave);  /* this also sets the thickness!!! */
						// TODO: modifying shared value from multiple threads?
						//muls.nslic0 = pCount;
					}
					/* run multislice algorithm
					   and save exit wave function for this position 
					   (done by runMulsSTEM), 
					   but we need to define the file name */
					sprintf(buf,"%s/mulswav_%d_%d.img",muls.folder,ix,iy);
					wave->fileout = std::string(buf);
					// TODO: modifying shared value from multiple threads?
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
					wave->detPosX=ix;
					wave->detPosY=iy;

					runMulsSTEM(&muls,wave); 


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
						sprintf(buf,"%s/diffAvg_%d_%d.img",muls.folder,ix,iy);
						wave->avgName = std::string(buf);
						// printf("Will copy to avgArray %d %d (%d, %d)\n",muls.nx, muls.ny,(int)(muls.diffpat),(int)avgArray);	

						if (muls.saveLevel > 0) 
						{
							if (muls.avgCount == 0)  
							{
								// initialize the avgArray from the diffpat
								for (ixa=0;ixa<muls.nx;ixa++) 
								{
									for (iya=0;iya<muls.ny;iya++)
									{
										wave->avgArray(iya,ixa)=wave->diffpat(iya,ixa);
									}
								}
								/* memcopy((void *)avgArray[0],(void *)muls.diffpat[0],
								(size_t)(muls.nx*muls.ny*sizeof(float_tt)));
								*/
								// printf("Copied to avgArray %d %d\n",muls.nx, muls.ny);	
							}
							else 
							{
								// printf("Will read image %d %d\n",muls.nx, muls.ny);	
								header_read = readImage(wave->avgArray, muls.nx, muls.ny, wave->avgName.c_str());
								for (ixa=0;ixa<muls.nx;ixa++) for (iya=0;iya<muls.ny;iya++) {
									t = ((float_tt)muls.avgCount * wave->avgArray(iya,ixa) +
										wave->diffpat(iya,ixa)) / ((float_tt)(muls.avgCount + 1));
									#pragma omp atomic
									chisq[muls.avgCount-1] += (wave->avgArray(iya,ixa)-t)*
										(wave->avgArray(iya,ixa)-t);
									wave->avgArray(iya,ixa) = t;
								}
							}
							/* Write the array to a file, resize and crop it, 
							* and convert it to jpg format 
							*/
							// writeRealImage_old(avgArray,muls.nx,muls.ny,wave->thickness,avgName);
							if (header == NULL) 
									header = makeNewHeaderCompact(0,muls.nx,muls.ny,wave->thickness,
										muls.resolutionX,muls.resolutionY,
										0,std::vector<float_tt>(),"diffraction pattern");
								// printf("Created header\n");
							header->t = wave->thickness;
							writeRealImage(wave->avgArray, header, wave->avgName.c_str());
							}	
							else {
								if (muls.avgCount > 0)	chisq[muls.avgCount-1] = 0.0;
							}
							/* make file names for tif and jpg files */
							/* crop, resize, convert to jpg, delete tif file */	    
							/* sprintf(systStr,"showimage %s 6; convert -crop 17%c "
							"-geometry 256x256 %s %s/%s; rm %s",
							avgName,'%',tifName,muls.folder,jpgName,tifName);
							*/
	#ifndef WIN32
							// MCS - Commented 07-2010
						
							//sprintf(systStr,"showimage %s 6; convert -crop 70x70+60+60%% "
							//	"-geometry 256x256 %s %s; rm %s",
							//	avgName,tifName,jpgName,tifName);
							/* printf("%s\n",systStr); */
							//system(systStr);
	#endif
					} /* end of if pCount == picts, i.e. conditional code, if this
						  * was the last slice
						  */

					#pragma omp atomic
					muls.complete_pixels+=1;

					if (muls.displayProgInterval > 0) if ((muls.complete_pixels) % muls.displayProgInterval == 0) 
					{
						#pragma omp atomic
						total_time += cputim()-timer;
						printf("Pixels complete: (%d/%d), int.=%.3f, avg time per pixel: %.2fsec\n",
							muls.complete_pixels, muls.scanYN*muls.scanYN, wave->intIntensity,
							(total_time)/muls.complete_pixels);
						timer=cputim();
					}
				} /* end of looping through STEM image pixels */
				/* save STEM images in img files */
				saveSTEMImages(&muls);
				muls.totalSliceCount += muls.slices;
				/*  calculate the total specimen thickness and echo */
				// commented 2013/04 MCS - we're already tracking total 
				//   thickness by reading intermediate mulswav.  
				//   What is this here for?
				/*
				cztot=0.0;
				for( islice=0; islice<muls.slices; islice++) 
				{
					cztot += muls.cz[islice];
				}
				muls.thickness=cztot*(pCount+1);
				*/
			} /* end of loop through thickness (pCount) */
		} /* end of  while (readparam("sequence: ",buf,0)) */
		// printf("Total CPU time = %f sec.\n", cputim()-timerTot ); 

#ifndef WIN32
		//add2STEMWebpage();
		//    makeSTEMPendelloesungPlot();

		/* remove the old diffraction pattern jpg files */
		sprintf(systStr,"rm %s/diffAvg*_%d.jpg",muls.folder,muls.avgCount-1);
		system(systStr);      
#endif
		/*************************************************************/

		chisq[muls.avgCount-1] = chisq[muls.avgCount-1]/(float_tt)(muls.nx*muls.ny);
		muls.intIntensity = collectedIntensity/(muls.scanXN*muls.scanYN);
		displayProgress(1);
	} /* end of for muls.avgCount=0..25 */

	//free(chisq);
	//for (int th=0; th<omp_get_num_threads(); th++)
	//{
	//	delete(waves[th]);
	//}
}

