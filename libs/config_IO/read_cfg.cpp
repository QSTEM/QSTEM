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

#include "readparams.hpp"

/************************************************************************
* readFile() 
*
* reads the parameters from the input file and does some 
* further setup accordingly
*
***********************************************************************/

void CCfgReader::ReadFile(const char *fileName, MULS &muls) {
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
  
  muls.mode = STEM;
  if (readparam("mode:",buf,1)) {
    if (strstr(buf,"STEM")) muls.mode = STEM;
    else if (strstr(buf,"TEM")) muls.mode = TEM;
    else if (strstr(buf,"CBED")) muls.mode = CBED;
    else if (strstr(buf,"TOMO")) muls.mode = TOMO;
    else if (strstr(buf,"REFINE")) muls.mode = REFINE;
  }

  if (readparam("print level:",buf,1)) sscanf(buf,"%d",&(muls.printLevel));
  if (readparam("save level:",buf,1)) sscanf(buf,"%d",&(muls.saveLevel));


  /************************************************************************
   * Basic microscope/simulation parameters, 
   */ 
  if (!readparam("filename:",buf,1)) exit(0); sscanf(buf,"%s",muls.fileBase);
  // search for second '"', in case the filename is in quotation marks:
  if (muls.fileBase[0] == '"') {
    strPtr = strchr(buf,'"');
    strcpy(muls.fileBase,strPtr+1);
    strPtr = strchr(muls.fileBase,'"');
    *strPtr = '\0';
  }

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
			muls.btiltx *= pi/180.0;
	}
	answer[0] = '\0';
	if (readparam("Beam tilt Y:",buf,1)) { 
		sscanf(buf,"%g %s",&(muls.btilty),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.btilty *= pi/180.0;
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
			muls.ctiltx *= pi/180.0;
	}
	answer[0] = '\0';
	if (readparam("Crystal tilt Y:",buf,1)) { 
		sscanf(buf,"%g %s",&(muls.ctilty),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.ctilty *= pi/180.0;
	}  
	answer[0] = '\0';
	if (readparam("Crystal tilt Z:",buf,1)) { 
		sscanf(buf,"%g %s",&(muls.ctiltz),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.ctiltz *= pi/180.0;
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
	sprintf(muls.atomPosFile,muls.fileBase);
	/* remove directory in front of file base: */
	while ((strptr = strchr(muls.fileBase,'\\')) != NULL) strcpy(muls.fileBase,strptr+1);

	/* add a '_' to fileBase, if not existent */
	if (strrchr(muls.fileBase,'_') != muls.fileBase+strlen(muls.fileBase)-1) {
		if ((strPtr = strchr(muls.fileBase,'.')) != NULL) sprintf(strPtr,"_");
		else strcat(muls.fileBase,"_");
	}
	if (strchr(muls.atomPosFile,'.') == NULL) {
		/*   
		strPtr = strrchr(muls.atomPosFile,'_');
		if (strPtr != NULL)
		*(strPtr) = 0;
		*/
		/* take atomPosFile as is, or add an ending to it, if it has none yet
		*/
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
	// We need to initialize a few variables, before reading the atomic 
	// positions for the first time.
	muls.natom = 0;
	muls.atoms = NULL;
	muls.Znums = NULL;
	muls.atomKinds = 0;
	muls.u2 = NULL;
	muls.u2avg = NULL;

	muls.xOffset = 0.0; /* slize z-position offset in cartesian coords */
	if (readparam("xOffset:",buf,1)) sscanf(buf,"%g",&(muls.xOffset));
	muls.yOffset = 0.0; /* slize z-position offset in cartesian coords */
	if (readparam("yOffset:",buf,1)) sscanf(buf,"%g",&(muls.yOffset));
	// printf("Reading Offset: %f, %f\n",muls.xOffset,muls.yOffset);

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
				muls.slices = (int)(muls.c/(muls.cellDiv*muls.sliceThickness)+0.99);
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
						muls.sliceThickness = (muls.centerSlices) ? 2.0*muls.cubez/muls.cellDiv : muls.cubez/muls.cellDiv;
					else
						// muls.sliceThickness = (muls.centerSlices) ? 2.0*muls.c/(muls.cellDiv) : muls.c/(muls.cellDiv);
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
	/* Find out whether we need to recalculate the potential every time, or not
	*/

	muls.equalDivs = ((!muls.tds)  && (muls.nCellZ % muls.cellDiv == 0) && 
		(fabs(muls.slices*muls.sliceThickness-muls.c/muls.cellDiv) < 1e-5));

	// read the output interval:
	muls.outputInterval = muls.slices;
	if (readparam("slices between outputs:",buf,1)) sscanf(buf,"%d",&(muls.outputInterval));
	if (muls.outputInterval < 1) muls.outputInterval= muls.slices;



	initMuls();  
	muls.czOffset = 0.0; /* slize z-position offset in cartesian coords */
	if (readparam("zOffset:",buf,1)) sscanf(buf,"%g",&(muls.czOffset));



	/***********************************************************************
	* Fit the resolution to the wave function array, if not specified different
	*/
	if (muls.resolutionX == 0.0)
		muls.resolutionX = muls.ax / (double)muls.nx;
	if (muls.resolutionY == 0.0)
		muls.resolutionY = muls.by / (double)muls.ny;



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
		if (readparam("scan_x_start:",buf,1)) sscanf(buf,"%g",&(muls.scanXStart));
		if (readparam("scan_y_start:",buf,1)) sscanf(buf,"%g",&(muls.scanYStart));
		muls.scanXStop = muls.scanXStart;
		muls.scanYStop = muls.scanYStart;
		break;

		/////////////////////////////////////////////////////////
		// Read STEM scanning parameters 

	case STEM:
          ReadSTEMParams(fileName, muls);
	}
	muls.displayPotCalcInterval = 1000;
	if (readparam("potential progress interval:",buf,1)) 
		sscanf(buf,"%d",&(muls.displayPotCalcInterval));
	// printf("Potential progress interval: %d\n",muls.displayPotCalcInterval);

	

	/**********************************************************************
	* Parameters for image display and directories, etc.
	*/
	muls.imageGamma = 1.0;
	if (readparam("Display Gamma:",buf,1)) {
		muls.imageGamma = atof(buf);
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
          ReadDetectors(fileName, muls);
	}
	/************************************************************************/   

	// in case this file has been written by the tomography function, read the current tilt:
	if (readparam("tomo tilt:",buf,1)) { 
		sscanf(buf,"%lf %s",&(muls.tomoTilt),answer); /* in mrad */
		if (tolower(answer[0]) == 'd')
			muls.tomoTilt *= 1000*pi/180.0;
	}
	/************************************************************************
	* Tomography Parameters:
	***********************************************************************/
	if (muls.mode == TOMO) {     
          ReadTomoParameters(fileName, muls);
	}
	/***********************************************************************/


	/*******************************************************************
	* Read in parameters related to the calculation of the projected
	* Potential
	*******************************************************************/
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
		muls.potOffsetX = muls.scanXStart - 0.5*muls.nx*muls.resolutionX;
		muls.potOffsetY = muls.scanYStart - 0.5*muls.ny*muls.resolutionY;
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
		muls.potOffsetX = muls.scanXStart - 0.5*muls.potSizeX;
		muls.potOffsetY = muls.scanYStart - 0.5*muls.potSizeY;
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
		(muls.scanXStart>muls.ax) || (muls.scanYStart>muls.by) ||
		(muls.scanXStop>muls.ax) || (muls.scanYStop>muls.by)) {
			printf("Scanning window is outside model dimensions (%g,%g .. %g,%g) [ax = %g, by = %g]!\n",muls.scanXStart,muls.scanYStart,muls.scanXStop,muls.scanYStop,muls.ax,muls.by);
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
			muls.hbeam = (int*)malloc(muls.nbout*sizeof(int));
			muls.kbeam = (int*)malloc(muls.nbout*sizeof(int));
			/* now read in the list of detectors: */
			resetParamFile();
			for (i=0;i<muls.nbout;i++) {
				if (!readparam("beam:",buf,0)) break;
				muls.hbeam[i] = 0;
				muls.kbeam[i] = 0;
				sscanf(buf,"%d %d",muls.hbeam+i,muls.kbeam+i);
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
} /* end of readFile() */

void CCfgReader::ReadSTEMParams(const char *fileName, MULS &muls)
{
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

void CCfgReader::ReadDetectors(const char *fileName, MULS &muls) // 
{
  int tCount = (int)(ceil((double)((muls.slices * muls.cellDiv) / muls.outputInterval)));

  /* first determine number of detectors */
  while (readparam("detector:",buf,0)) muls.detectorNum++;  
  /* now read in the list of detectors: */
  resetParamFile();

  // loop over thickness planes where we're going to record intermediates
  // TODO: is this too costly in terms of memory?  It simplifies the parallelization to
  //       save each of the thicknesses in memory, then save to disk afterwards.
  for (int islice=0; islice<=tCount; islice++)
    {
      std::vector<DetectorPtr> detectors;
      resetParamFile();
      while (readparam("detector:",buf,0)) {
        DetectorPtr det = DetectorPtr(new Detector(muls.scanXN, muls.scanYN, 
                                                   (muls.scanXStop-muls.scanXStart)/(float)muls.scanXN,
                                                   (muls.scanYStop-muls.scanYStart)/(float)muls.scanYN));
        
        sscanf(buf,"%g %g %s %g %g",&(det->rInside),
               &(det->rOutside), det->name, &(det->shiftX),&(det->shiftY));  
        
        /* determine v0 specific k^2 values corresponding to the angles */
        det->k2Inside = (float)(sin(det->rInside*0.001)/(wavelength(muls.v0)));
        det->k2Outside = (float)(sin(det->rOutside*0.001)/(wavelength(muls.v0)));
        /* calculate the squares of the ks */
        det->k2Inside *= det->k2Inside;
        det->k2Outside *= det->k2Outside;
        detectors.push_back(det);
      }
      muls.detectors.push_back(detectors);
    }
}

void CCfgReader::ReadProbeParameters(const char *fileName, MULS &muls)
{
  /**********************************************************************
   * Read STEM/CBED probe parameters 
   */
	muls.dE_E = 0.0;
	muls.dI_I = 0.0;
	muls.dV_V = 0.0;
	muls.Cc = 0.0;
	if (readparam("dE/E:",buf,1))
		muls.dE_E = atof(buf);
	if (readparam("dI/I:",buf,1))
		muls.dI_I = atof(buf);
	if (readparam("dV/V:",buf,1))
		muls.dV_V = atof(buf);
	if (readparam("Cc:",buf,1))
		muls.Cc = 1e7*atof(buf);


	/* memorize dE_E0, and fill the array of well defined energy deviations */
	dE_E0 = sqrt(muls.dE_E*muls.dE_E+
		muls.dI_I*muls.dI_I+
		muls.dV_V*muls.dV_V);
	muls.dE_EArray = (double *)malloc((muls.avgRuns+1)*sizeof(double));
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


	
	

	////////////////////////////////////////////////////////
	// read in more aberrations:
	
        ReadAberrationAmplitudes(muls.Cs, muls.C5, 
                                 muls.df0, muls.astigMag,
                                 muls.a33, muls.a31,
                                 muls.a44, muls.a42,
                                 muls.a55, muls.a53, muls.a51,
                                 muls.a66, muls.a64, muls.a62);

        ReadAberrationAngles(muls.astig, 
                             muls.phi33, muls.phi31,
                             muls.phi44, muls.phi42,
                             muls.phi55, muls.phi53, muls.phi51,
                             muls.phi66, muls.phi64, muls.phi62);

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
	muls.electronScale = muls.beamCurrent*muls.dwellTime*MILLISEC_PICOAMP;
	//////////////////////////////////////////////////////////////////////

	muls.sourceRadius = 0;
	if (readparam("Source Size (diameter):",buf,1)) 
		muls.sourceRadius = atof(buf)/2.0;

	if (readparam("smooth:",buf,1)) sscanf(buf,"%s",answer);
	muls.ismoth = (tolower(answer[0]) == (int)'y');
	muls.gaussScale = 0.05f;
	muls.gaussFlag = 0;
	if (readparam("gaussian:",buf,1)) {
		sscanf(buf,"%s %g",answer,&(muls.gaussScale));
		muls.gaussFlag = (tolower(answer[0]) == (int)'y');
	}
}

void CCfgReader::ReadTomoParameters(float_tt &tomoStart, float_tt &tomoStep, int &tomoCount,
                     float_tt &zoomFactor)
{
  if (readparam("tomo start:",buf,1)) { 
    sscanf(buf,"%lf %s",&(tomoStart),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      tomoStart *= 1000*pi/180.0;
  }
  if (readparam("tomo step:",buf,1)) {
    sscanf(buf,"%lf %s",&(tomoStep),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      tomoStep *= 1000*pi/180.0;
  }
  
  if (readparam("tomo count:",buf,1))  
    tomoCount = atoi(buf); 
  if (readparam("zoom factor:",buf,1))  
    sscanf(buf,"%lf",&(zoomFactor));
  if ((tomoStep == 0) && (tomoStep > 1))
    tomoStep = -2.0*tomoStart/(float_tt)(tomoCount - 1);
}

CCfgReader::ReadAberrationAmplitudes(float_tt &Cs, float_tt &C5,
                           float_tt &df0, float_tt &astig,
                           float_tt &a33, float_tt &a31,
                           float_tt &a44, float_tt &a42,
                           float_tt &a55, float_tt &a53, float_tt &a51,
                           float_tt &a66, float_tt &a64, float_tt &a62)
{
  if (!readparam("Cs:",buf,1))  exit(0); 
  sscanf(buf,"%g",&(Cs)); /* in mm */
  Cs *= 1.0e7; /* convert Cs from mm to Angstroem */

  C5 = 0;
  if (readparam("C5:",buf,1)) { 
    sscanf(buf,"%g",&(C5)); /* in mm */
    C5 *= 1.0e7; /* convert C5 from mm to Angstroem */
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
      muls.df0 = 10.0*muls.df0;       /* convert defocus to A */
      muls.Scherzer = (-(float)sqrt(1.5*muls.Cs*(wavelength(muls.v0)))==muls.df0);
    }
  }
  // Astigmatism:
  astig = 0;
  if (readparam("astigmatism:",buf,1)) sscanf(buf,"%g",&(astig)); 
  // convert to A from nm:
  astig = 10.0*astig;

  a33 = 0;
  a31 = 0;
  a44 = 0;
  a42 = 0;
  a55 = 0;
  a53 = 0;
  a51 = 0;
  a66 = 0;
  a64 = 0;
  a62 = 0;

  if (readparam("a_33:",buf,1)) {sscanf(buf,"%g",&(a33)); }
  if (readparam("a_31:",buf,1)) {sscanf(buf,"%g",&(a31)); }
  if (readparam("a_44:",buf,1)) {sscanf(buf,"%g",&(a44)); }
  if (readparam("a_42:",buf,1)) {sscanf(buf,"%g",&(a42)); }
  if (readparam("a_55:",buf,1)) {sscanf(buf,"%g",&(a55)); }
  if (readparam("a_53:",buf,1)) {sscanf(buf,"%g",&(a53)); }
  if (readparam("a_51:",buf,1)) {sscanf(buf,"%g",&(a51)); }
  if (readparam("a_66:",buf,1)) {sscanf(buf,"%g",&(a66)); }
  if (readparam("a_64:",buf,1)) {sscanf(buf,"%g",&(a64)); }
  if (readparam("a_62:",buf,1)) {sscanf(buf,"%g",&(a62)); }
}

ReadAberrationAngles(float_tt &astig,
                       float_tt &phi33, float_tt &phi31,
                       float_tt &phi44, float_tt &phi42,
                       float_tt &phi55, float_tt &phi53, float_tt &phi51,
                       float_tt &phi66, float_tt &phi64, float_tt &phi62)
{
  astig = 0;
  if (readparam("astigmatism angle:",buf,1)) sscanf(buf,"%g",&(astig)); 
  // convert astigAngle from deg to rad:
  astig *= pi/180.0;

  phi33 = 0;
  phi31 = 0;
  phi44 = 0;
  phi42 = 0;
  phi55 = 0;
  phi53 = 0;
  phi51 = 0;
  phi66 = 0;
  phi64 = 0;
  phi62 = 0;


  if (readparam("phi_33:",buf,1)) {sscanf(buf,"%g",&(phi33)); }
  if (readparam("phi_31:",buf,1)) {sscanf(buf,"%g",&(phi31)); }
  if (readparam("phi_44:",buf,1)) {sscanf(buf,"%g",&(phi44)); }
  if (readparam("phi_42:",buf,1)) {sscanf(buf,"%g",&(phi42)); }
  if (readparam("phi_55:",buf,1)) {sscanf(buf,"%g",&(phi55)); }
  if (readparam("phi_53:",buf,1)) {sscanf(buf,"%g",&(phi53)); }
  if (readparam("phi_51:",buf,1)) {sscanf(buf,"%g",&(phi51)); }
  if (readparam("phi_66:",buf,1)) {sscanf(buf,"%g",&(phi66)); }
  if (readparam("phi_64:",buf,1)) {sscanf(buf,"%g",&(phi64)); }
  if (readparam("phi_62:",buf,1)) {sscanf(buf,"%g",&(phi62)); }

  phi33 /= (float)RAD2DEG;
  phi31 /= (float)RAD2DEG;
  phi44 /= (float)RAD2DEG;
  phi42 /= (float)RAD2DEG;
  phi55 /= (float)RAD2DEG;
  phi53 /= (float)RAD2DEG;
  phi51 /= (float)RAD2DEG;
  phi66 /= (float)RAD2DEG;
  phi64 /= (float)RAD2DEG;
  phi62 /= (float)RAD2DEG;
}
