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

#include "read_qsc.hpp"
#include "readparams.hpp"

#include "string.h"
#include <boost/filesystem.hpp>

CQscReader::CQscReader(std::string &filename) : IConfigReader()
{
  // open the file for reading
  if (parOpen(filename.c_str()) == 0) {
    m_isValid=false;
  }
}
  
CQscReader::~CQscReader()
{
  // make sure the file is closed
  parClose();
}

void CQscReader::ReadMode(int &mode)
{
  mode = STEM;
  if (readparam("mode:",buf,1)) {
    if (strstr(buf,"STEM")) mode = STEM;
    else if (strstr(buf,"TEM")) mode = TEM;
    else if (strstr(buf,"CBED")) mode = CBED;
    else if (strstr(buf,"TOMO")) mode = TOMO;
    else if (strstr(buf,"REFINE")) mode = REFINE;
  }
}

void CQscReader::ReadOutputLevel(int &printLevel, int &saveLevel, 
                                 unsigned &displayPotCalcInterval, unsigned &displayProgInterval)
{
  if (readparam("print level:",buf,1)) sscanf(buf,"%d",&(printLevel));
  if (readparam("save level:",buf,1)) sscanf(buf,"%d",&(saveLevel));
  displayPotCalcInterval = 1000;
  
  if (readparam("potential progress interval:",buf,1)) 
    sscanf(buf,"%d",&(displayPotCalcInterval));

  if (readparam("propagation progress interval:",buf,1)) 
    sscanf(buf,"%d",&(displayProgInterval));
}

// TODO: this should remove quotes from the output directory/filename
void CQscReader::ReadStructureFileName(std::string &directory, std::string &fileOrDirName)
{
  if (!readparam("filename:",buf,1)) exit(0); 

  boost::filesystem::path structure_path = boost::filesystem::path(buf);
  
  // TODO: use boost::filesystem to do any manipulation necessary

  directory = structure_path.parent_path().string();
  fileOrDirName = structure_path.filename().string();
}

void CQscReader::ReadNCells(unsigned &nCellX, unsigned &nCellY, unsigned &nCellZ, int &cellDiv)
{
  char *strPtr;
  if (readparam("NCELLX:",buf,1)) sscanf(buf,"%d",&(nCellX));
  if (readparam("NCELLY:",buf,1)) sscanf(buf,"%d",&(nCellY));

  if (readparam("NCELLZ:",buf,1)) {
    sscanf(buf,"%s",answer);
    if ((strPtr = strchr(answer,'/')) != NULL) {
      strPtr[0] = '\0';
      cellDiv = atoi(strPtr+1);
    }
    nCellZ = atoi(answer);
  }
}
  /*************************************************
   * Read the beam tilt parameters
   */
void CQscReader::ReadBeamTilt(float_tt &btiltx, float_tt &btilty, bool tiltBack)
{
  answer[0] = '\0';
  if (readparam("Beam tilt X:",buf,1)) { 
    sscanf(buf,"%g %s",&(btiltx),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      btiltx *= PI/180.0;
  }
  answer[0] = '\0';
  if (readparam("Beam tilt Y:",buf,1)) { 
    sscanf(buf,"%g %s",&(btilty),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      btilty *= PI/180.0;
  }  
  if (readparam("Tilt back:",buf,1)) { 
    sscanf(buf,"%s",answer);
    tiltBack  = (tolower(answer[0]) == (int)'y');
  }

}
void CQscReader::ReadCrystalCubeAndTilt(float_tt &tiltx, float_tt &tilty, float_tt &tiltz, 
                                     float_tt &cubex, float_tt &cubey, float_tt &cubez, 
                                     bool &adjustCubeSize)
{
  /*************************************************
   * Read the crystal tilt parameters
   */
  answer[0] = '\0';
  if (readparam("Crystal tilt X:",buf,1)) { 
    sscanf(buf,"%g %s",&(tiltx),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      tiltx *= PI/180.0;
  }
  answer[0] = '\0';
  if (readparam("Crystal tilt Y:",buf,1)) { 
    sscanf(buf,"%g %s",&(tilty),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      tilty *= PI/180.0;
  }  
  answer[0] = '\0';
  if (readparam("Crystal tilt Z:",buf,1)) { 
    sscanf(buf,"%g %s",&(tiltz),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      tiltz *= PI/180.0;
  }
  cubex = 0; cubey = 0; cubez = 0;
  if (readparam("Cube:",buf,1)) { 
    sscanf(buf,"%g %g %g",&(cubex),&(cubey),&(cubez)); /* in A */
  }
  
  adjustCubeSize = false;
  if (readparam("Adjust cube size with tilt:",buf,1)) { 
    sscanf(buf,"%s",answer);
    adjustCubeSize  = (tolower(answer[0]) == (int)'y');
  }
}

void CQscReader::ReadTemperatureData(bool &doTDS, float_tt &tdsTemperature, std::string &phononFile, 
                                     bool &useEinstein)
{
  if (readparam("tds:",buf,1)) {
    sscanf(buf,"%s",answer);
    doTDS = (tolower(answer[0]) == (int)'y');
  }
  else doTDS = false;

  if (readparam("temperature:",buf,1)) sscanf(buf,"%g",&(tdsTemperature));
  else tdsTemperature = 300.0;

  useEinstein = true;  //phononFile = NULL;
  if (readparam("phonon-File:",buf,1)) {
    phononFile = buf;
    //sscanf(buf,"%s",phononFile);
    useEinstein = false;
  }
}

void CQscReader::ReadSliceOffset(float_tt &xOffset, float_tt &yOffset)
{
	xOffset = 0.0; /* slize z-position offset in cartesian coords */
	if (readparam("xOffset:",buf,1)) sscanf(buf,"%g",&(xOffset));
	yOffset = 0.0; /* slize z-position offset in cartesian coords */
	if (readparam("yOffset:",buf,1)) sscanf(buf,"%g",&(yOffset));
	// printf("Reading Offset: %f, %f\n",xOffset,yOffset);
}

void CQscReader::ReadProbeArraySize(unsigned &nx, unsigned &ny)
{
	if (!readparam("nx:",buf,1)) exit(0); sscanf(buf,"%d",&(nx));
	if (readparam("ny:",buf,1)) sscanf(buf,"%d",&(ny));
	else ny = nx;
}

void CQscReader::ReadResolution(float_tt &resolutionX, float_tt &resolutionY)
{
	if (readparam("resolutionX:",buf,1)) sscanf(buf,"%g",&(resolutionX));
	if (readparam("resolutionY:",buf,1)) sscanf(buf,"%g",&(resolutionY));
}

void CQscReader::ReadVoltage(float_tt &voltage)
{
	if (!readparam("v0:",buf,1)) exit(0); sscanf(buf,"%g",&(voltage));
}

void CQscReader::ReadSliceParameters(bool &centerSlices, float_tt &sliceThickness, 
                                     unsigned &nslices, unsigned &outputInterval,
                                     float_tt &zOffset)
{
  if (readparam("center slices:",buf,1)) {
    // answer[0] =0;
    sscanf(buf,"%s",answer);
    // printf("center: %s (%s)\n",answer,buf);
    centerSlices = (tolower(answer[0]) == (int)'y');
  }
  // just in case the answer was not exactly 1 or 0:
  // centerSlices = (centerSlices) ? 1 : 0;

  if (readparam("slice-thickness:",buf,1)) sscanf(buf,"%g",&(sliceThickness));
  if (readparam("slices:",buf,1)) sscanf(buf,"%d",&(nslices));
    
  // read the output interval:
  if (readparam("slices between outputs:",buf,1)) sscanf(buf,"%d",&(outputInterval));
  if (readparam("zOffset:",buf,1)) sscanf(buf,"%g",&(zOffset));
}

void CQscReader::ReadPeriodicParameters(bool &periodicXY, bool &periodicZ)
{
  if (readparam("periodicXY:",buf,1)) {
    sscanf(buf,"%s",answer);
    periodicXY = (tolower(answer[0]) == (int)'y');
  }
  if (readparam("periodicZ:",buf,1)) {
    sscanf(buf,"%s",answer);
    periodicZ = (tolower(answer[0]) == (int)'y'); 
  }
}

void CQscReader::ReadBandLimitTrans(bool &bandlimittrans)
{
	if (readparam("bandlimit f_trans:",buf,1)) {
		sscanf(buf,"%s",answer);
		bandlimittrans = (tolower(answer[0]) == (int)'y');
	}
}

void CQscReader::ReadLoadPotential(bool &readPotential)
{    
  if (readparam("read potential:",buf,1)) {
    sscanf(buf," %s",answer);
    readPotential = (tolower(answer[0]) == (int)'y');
  }
}

void CQscReader::ReadPotentialOutputParameters(bool &savePotential, bool &saveTotalPotential, 
                                               bool &plotPotential)
{
	if (readparam("save potential:",buf,1)) {
		sscanf(buf," %s",answer);
		savePotential = (tolower(answer[0]) == (int)'y');
	}  
	if (readparam("save projected potential:",buf,1)) {
		sscanf(buf," %s",answer);
		saveTotalPotential = (tolower(answer[0]) == (int)'y');
	}  
	if (readparam("plot V(r)*r:",buf,1)) {
		sscanf(buf," %s",answer);
		plotPotential = (tolower(answer[0]) == (int)'y');
	}  
}

void CQscReader::ReadPotentialCalculationParameters(bool &fftPotential, bool &potential3D)
{
  if (readparam("one time integration:",buf,1)) {
    sscanf(buf,"%s",answer);
    fftPotential = (tolower(answer[0]) == (int)'y');
  }
  if (readparam("potential3D:",buf,1)) {
    sscanf(buf,"%s",answer);
    potential3D = (tolower(answer[0]) == (int)'y');
  }
}

void CQscReader::ReadAverageParameters(unsigned &avgRuns, bool &storeSeries)
{
	if (readparam("Runs for averaging:",buf,1))
		sscanf(buf,"%d",&(avgRuns));

	if (readparam("Store TDS diffr. patt. series:",buf,1)) {
		sscanf(buf,"%s",answer);
		storeSeries = (tolower(answer[0]) == (int)'y');
	}  
}

void CQscReader::ReadScanParameters(float_tt &scanXStart, float_tt &scanXStop, unsigned &scanXN,
                                    float_tt &scanYStart, float_tt &scanYStop, unsigned &scanYN)
{
  if (!readparam("scan_x_start:",buf,1)) sscanf(buf,"%g",&(scanXStart));
  if (!readparam("scan_x_stop:",buf,1)) sscanf(buf,"%g",&(scanXStop));
  if (!readparam("scan_x_pixels:",buf,1)) sscanf(buf,"%d",&(scanXN));
  if (!readparam("scan_y_start:",buf,1)) sscanf(buf,"%g",&(scanYStart));
  if (!readparam("scan_y_stop:",buf,1)) sscanf(buf,"%g",&(scanYStop));
  if (!readparam("scan_y_pixels:",buf,1)) sscanf(buf,"%d",&(scanYN));
}


void CQscReader::ReadOutputName(std::string &fileOrFolderName)
{
  /**********************************************************************
   * Parameters for image display and directories, etc.
   */
  fileOrFolderName = "data";
  //sprintf(fileOrFolderName,"data");
  if (readparam("Folder:",buf,1)) 
    fileOrFolderName = buf;
  //sscanf(buf," %s",folder);
  
  /*
    // TODO: what is this doing? stripping ""?  do any of this that is required using boost::filesystem.
  while ((strptr = strchr(folder,'"')) != NULL) {
    if (strptr == folder) sscanf(strptr+1,"%s",folder);	
    else *strptr = '\0';
  }
  
  if (folder[strlen(folder)-1] == '/')
  folder[strlen(folder)-1] = 0;
  */
}





void CQscReader::ReadAtomRadius(float_tt &radius)
{
	/*******************************************************************
	* Read in parameters related to the calculation of the projected
	* Potential
	*******************************************************************/
        if (readparam("atom radius:",buf,1))  
		sscanf(buf,"%g",&(radius)); /* in A */
}

void CQscReader::ReadStructureFactorType(int &scatFactor)
{
  if (readparam("Structure Factors:",buf,1)) {
    sscanf(buf," %s",answer);
    switch (tolower(answer[0])) {
    case 'w':
      if (tolower(answer[1])=='k') scatFactor = WEICK_KOHL;
      break;
    case 'd':  // DOYLE_TURNER
      scatFactor = DOYLE_TURNER;
      break;
    case 'c':  // CUSTOM - specify k-lookup table and values for all atoms used
      scatFactor = CUSTOM;
      break;
    default:
      scatFactor = DOYLE_TURNER;
    }
  }
}

void CQscReader::ReadPendelloesungParameters(std::vector<int> &hbeams, std::vector<int> &kbeams, 
                                             bool &lbeams, unsigned &nbout)
{
  /*************************************************************
   * read in the beams we want to plot in the pendeloesung plot
   * Only possible if not in STEM or CBED mode 
   *************************************************************/
  resetParamFile();
  if (readparam("Pendelloesung plot:",buf,1)) {
    sscanf(buf,"%s",answer);
    lbeams = (tolower(answer[0]) == (int)'y');
  }
  if (lbeams) {
    while (readparam("beam:",buf,0)) nbout++;  
    printf("will record %d beams\n",nbout);
    hbeams.resize(nbout);
    kbeams.resize(nbout);
    /* now read in the list of detectors: */
    resetParamFile();
    for (unsigned i=0;i<nbout;i++) {
      if (!readparam("beam:",buf,0)) break;
      hbeams[i] = 0;
      kbeams[i] = 0;
      sscanf(buf,"%d %d",&(hbeams[i]),&(kbeams[i]));
    }
  }
}

void CQscReader::ReadNumberOfDetectors(int &numDetectors)
{
  numDetectors=0;
  resetParamFile();
  while (readparam("detector:",buf, 0)) ++numDetectors;
}

void CQscReader::ReadDetectorParameters(int det_idx, float_tt &rInside, float_tt &rOutside, std::string &name, 
                                        float_tt &shiftX, float_tt &shiftY)
{
  int file_det_idx=0;
  char name_buf[100];
  resetParamFile();
  while (readparam("detector:",buf,0)) {
    if (det_idx==file_det_idx)
      {
        sscanf(buf,"%g %g %s %g %g",&(rInside),
               &(rOutside), name_buf, &(shiftX),&(shiftY));
        name=name_buf;
        return;
      }
    ++file_det_idx;
  }
}

void CQscReader::ReadDoseParameters(float_tt &beamCurrent, float_tt &dwellTimeMs)
{
  ///// read beam current and dwell time ///////////////////////////////
  if (readparam("beam current:",buf,1)) { 
    sscanf(buf,"%g",&(beamCurrent)); /* in pA */
  }
  if (readparam("dwell time:",buf,1)) { 
    sscanf(buf,"%g",&(dwellTimeMs)); /* in msec */
  }
}

void CQscReader::ReadProbeParameters(float_tt &dE_E, float_tt &dI_I, float_tt &dV_V, float_tt &alpha, float_tt &aAIS,
                                     float_tt &sourceRadius, bool &ismoth, float_tt &gaussScale, bool &gaussFlag)
{
  /**********************************************************************
   * Read STEM/CBED probe parameters 
   */
	if (readparam("dE/E:",buf,1))
		dE_E = atof(buf);
	if (readparam("dI/I:",buf,1))
		dI_I = atof(buf);
	if (readparam("dV/V:",buf,1))
		dV_V = atof(buf);

	if (!readparam("alpha:",buf,1)) exit(0); 
	sscanf(buf,"%g",&(alpha)); /* in mrad */

	if (readparam("AIS aperture:",buf,1)) 
		sscanf(buf,"%g",&(aAIS)); /* in A */

	
	//////////////////////////////////////////////////////////////////////

	if (readparam("Source Size (diameter):",buf,1)) 
		sourceRadius = atof(buf)/2.0;

	if (readparam("smooth:",buf,1)) sscanf(buf,"%s",answer);
	ismoth = (tolower(answer[0]) == (int)'y');
	if (readparam("gaussian:",buf,1)) {
		sscanf(buf,"%s %g",answer,&(gaussScale));
		gaussFlag = (tolower(answer[0]) == (int)'y');
	}
}

void CQscReader::ReadTomoParameters(float_tt &tomoTilt, float_tt &tomoStart, float_tt &tomoStep, int &tomoCount,
                     float_tt &zoomFactor)
{
  // in case this file has been written by the tomography function, read the current tilt:
  if (readparam("tomo tilt:",buf,1)) { 
    sscanf(buf,"%g %s",&(tomoTilt),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      tomoTilt *= 1000*PI/180.0;
  }

  if (readparam("tomo start:",buf,1)) { 
    sscanf(buf,"%g %s",&(tomoStart),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      tomoStart *= 1000*PI/180.0;
  }
  if (readparam("tomo step:",buf,1)) {
    sscanf(buf,"%g %s",&(tomoStep),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      tomoStep *= 1000*PI/180.0;
  }
  
  if (readparam("tomo count:",buf,1))  
    tomoCount = atoi(buf); 
  if (readparam("zoom factor:",buf,1))  
    sscanf(buf,"%g",&(zoomFactor));
  if ((tomoStep == 0) && (tomoStep > 1))
    tomoStep = -2.0*tomoStart/(float_tt)(tomoCount - 1);
}

void CQscReader::ReadAberrationAmplitudes(float_tt &Cs, float_tt &C5, float_tt &Cc,
                                          float_tt &df0, int &Scherzer, float_tt &astig,
                                          float_tt &a33, float_tt &a31,
                                          float_tt &a44, float_tt &a42,
                                          float_tt &a55, float_tt &a53, float_tt &a51,
                                          float_tt &a66, float_tt &a64, float_tt &a62)
{
  if (!readparam("Cs:",buf,1))  exit(0); 
  sscanf(buf,"%g",&(Cs)); /* in mm */

  if (readparam("C5:",buf,1)) { 
    sscanf(buf,"%g",&(C5)); /* in mm */
  }
  if (readparam("Cc:",buf,1))
    Cc = atof(buf);

  if (readparam("defocus:",buf,1)) { 
    sscanf(buf,"%s",answer);
    /* if Scherzer defocus */
    if (tolower(answer[0]) == 's') {
      Scherzer = 1;
    }
    else if (tolower(answer[0]) == 'o') {
      Scherzer = 2;
    }
    else {
      sscanf(buf,"%g",&(df0)); /* in nm */
      Scherzer = 0;
    }
  }
  // Astigmatism:
  if (readparam("astigmatism:",buf,1)) sscanf(buf,"%g",&(astig)); 

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

void CQscReader::ReadAberrationAngles(float_tt &astig,
                       float_tt &phi33, float_tt &phi31,
                       float_tt &phi44, float_tt &phi42,
                       float_tt &phi55, float_tt &phi53, float_tt &phi51,
                       float_tt &phi66, float_tt &phi64, float_tt &phi62)
{
  if (readparam("astigmatism angle:",buf,1)) sscanf(buf,"%g",&(astig)); 

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
}
