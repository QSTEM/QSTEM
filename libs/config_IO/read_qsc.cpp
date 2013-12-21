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

CQscReader::CQscReader(boost::filesystem::path &filename) : IConfigReader()
{
  m_fp = fopen(filename.string().c_str(), "r" );
  if( m_fp == NULL ) {
    printf("Cannot open file %s\n",filename.string().c_str());
    m_isValid=false;
  }
}
  
CQscReader::~CQscReader()
{
  // make sure the file is closed
  fclose(m_fp);
}

void CQscReader::ReadMode(int &mode)
{
  mode = STEM;
  if (readparam(m_fp,"mode:",buf,1)) {
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
  if (readparam(m_fp,"print level:",buf,1)) sscanf(buf,"%d",&(printLevel));
  if (readparam(m_fp,"save level:",buf,1)) sscanf(buf,"%d",&(saveLevel));
  displayPotCalcInterval = 1000;
  
  if (readparam(m_fp,"potential progress interval:",buf,1)) 
    sscanf(buf,"%d",&(displayPotCalcInterval));

  if (readparam(m_fp,"propagation progress interval:",buf,1)) 
    sscanf(buf,"%d",&(displayProgInterval));
}

// TODO: this should remove quotes from the output directory/filename
void CQscReader::ReadStructureFileName(boost::filesystem::path &structure_path)
{
  if (!readparam(m_fp,"filename:",buf,1)) exit(0); 
  structure_path = boost::filesystem::path(buf);
}

void CQscReader::ReadNCells(unsigned &nCellX, unsigned &nCellY, unsigned &nCellZ)
{
  if (readparam(m_fp,"NCELLX:",buf,1)) sscanf(buf,"%d",&(nCellX));
  if (readparam(m_fp,"NCELLY:",buf,1)) sscanf(buf,"%d",&(nCellY));

  if (readparam(m_fp,"NCELLZ:",buf,1)) {
    sscanf(buf,"%s",answer);
    nCellZ = atoi(answer);
  }
}

void CQscReader::ReadNSubSlabs(unsigned &cellDiv)
{
  char *strPtr;

  // This is stored on the same line as nCellZ, hence the duplication.
  //    We don't store nCellZ again here.
  if (readparam(m_fp, "NCELLZ:",buf,1)) {
    sscanf(buf,"%s",answer);
    if ((strPtr = strchr(answer,'/')) != NULL) {
      strPtr[0] = '\0';
      cellDiv = atoi(strPtr+1);
    }
  }
}

  /*************************************************
   * Read the beam tilt parameters
   */
void CQscReader::ReadBeamTilt(float_tt &btiltx, float_tt &btilty, bool tiltBack)
{
  answer[0] = '\0';
  if (readparam(m_fp, "Beam tilt X:",buf,1)) { 
    sscanf(buf,"%g %s",&(btiltx),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      btiltx *= PI/180.0;
  }
  answer[0] = '\0';
  if (readparam(m_fp, "Beam tilt Y:",buf,1)) { 
    sscanf(buf,"%g %s",&(btilty),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      btilty *= PI/180.0;
  }  
  if (readparam(m_fp, "Tilt back:",buf,1)) { 
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
  if (readparam(m_fp, "Crystal tilt X:",buf,1)) { 
    sscanf(buf,"%g %s",&(tiltx),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      tiltx *= PI/180.0;
  }
  answer[0] = '\0';
  if (readparam(m_fp, "Crystal tilt Y:",buf,1)) { 
    sscanf(buf,"%g %s",&(tilty),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      tilty *= PI/180.0;
  }  
  answer[0] = '\0';
  if (readparam(m_fp, "Crystal tilt Z:",buf,1)) { 
    sscanf(buf,"%g %s",&(tiltz),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      tiltz *= PI/180.0;
  }
  cubex = 0; cubey = 0; cubez = 0;
  if (readparam(m_fp, "Cube:",buf,1)) { 
    sscanf(buf,"%g %g %g",&(cubex),&(cubey),&(cubez)); /* in A */
  }
  
  adjustCubeSize = false;
  if (readparam(m_fp, "Adjust cube size with tilt:",buf,1)) { 
    sscanf(buf,"%s",answer);
    adjustCubeSize  = (tolower(answer[0]) == (int)'y');
  }
}

void CQscReader::ReadTemperatureData(bool &doTDS, float_tt &tdsTemperature, std::string &phononFile, 
                                     bool &useEinstein)
{
  if (readparam(m_fp, "tds:",buf,1)) {
    sscanf(buf,"%s",answer);
    doTDS = (tolower(answer[0]) == (int)'y');
  }
  else doTDS = false;

  if (readparam(m_fp, "temperature:",buf,1)) sscanf(buf,"%g",&(tdsTemperature));
  else tdsTemperature = 300.0;

  useEinstein = true;  //phononFile = NULL;
  if (readparam(m_fp, "phonon-File:",buf,1)) {
    phononFile = buf;
    //sscanf(buf,"%s",phononFile);
    useEinstein = false;
  }
}

void CQscReader::ReadSliceOffset(float_tt &xOffset, float_tt &yOffset)
{
	xOffset = 0.0; /* slize z-position offset in cartesian coords */
	if (readparam(m_fp, "xOffset:",buf,1)) sscanf(buf,"%g",&(xOffset));
	yOffset = 0.0; /* slize z-position offset in cartesian coords */
	if (readparam(m_fp, "yOffset:",buf,1)) sscanf(buf,"%g",&(yOffset));
	// printf("Reading Offset: %f, %f\n",xOffset,yOffset);
}

void CQscReader::ReadProbeArraySize(unsigned &nx, unsigned &ny)
{
	if (!readparam(m_fp, "nx:",buf,1)) exit(0); sscanf(buf,"%d",&(nx));
	if (readparam(m_fp, "ny:",buf,1)) sscanf(buf,"%d",&(ny));
	else ny = nx;
}

void CQscReader::ReadResolution(float_tt &resolutionX, float_tt &resolutionY)
{
	if (readparam(m_fp, "resolutionX:",buf,1)) sscanf(buf,"%g",&(resolutionX));
	if (readparam(m_fp, "resolutionY:",buf,1)) sscanf(buf,"%g",&(resolutionY));
}

void CQscReader::ReadVoltage(float_tt &voltage)
{
	if (!readparam(m_fp, "v0:",buf,1)) exit(0); sscanf(buf,"%g",&(voltage));
}

void CQscReader::ReadSliceParameters(bool &centerSlices, float_tt &sliceThickness, 
                                     unsigned &nslices, unsigned &outputInterval,
                                     float_tt &zOffset)
{
  if (readparam(m_fp, "center slices:",buf,1)) {
    // answer[0] =0;
    sscanf(buf,"%s",answer);
    // printf("center: %s (%s)\n",answer,buf);
    centerSlices = (tolower(answer[0]) == (int)'y');
  }
  // just in case the answer was not exactly 1 or 0:
  // centerSlices = (centerSlices) ? 1 : 0;

  if (readparam(m_fp, "slice-thickness:",buf,1)) sscanf(buf,"%g",&(sliceThickness));
  if (readparam(m_fp, "slices:",buf,1)) sscanf(buf,"%d",&(nslices));
    
  // read the output interval:
  if (readparam(m_fp, "slices between outputs:",buf,1)) sscanf(buf,"%d",&(outputInterval));
  if (readparam(m_fp, "zOffset:",buf,1)) sscanf(buf,"%g",&(zOffset));
}

void CQscReader::ReadPeriodicParameters(bool &periodicXY, bool &periodicZ)
{
  if (readparam(m_fp, "periodicXY:",buf,1)) {
    sscanf(buf,"%s",answer);
    periodicXY = (tolower(answer[0]) == (int)'y');
  }
  if (readparam(m_fp, "periodicZ:",buf,1)) {
    sscanf(buf,"%s",answer);
    periodicZ = (tolower(answer[0]) == (int)'y'); 
  }
}

void CQscReader::ReadBandLimitTrans(bool &bandlimittrans)
{
	if (readparam(m_fp, "bandlimit f_trans:",buf,1)) {
		sscanf(buf,"%s",answer);
		bandlimittrans = (tolower(answer[0]) == (int)'y');
	}
}

void CQscReader::ReadLoadPotential(bool &readPotential)
{    
  if (readparam(m_fp, "read potential:",buf,1)) {
    sscanf(buf," %s",answer);
    readPotential = (tolower(answer[0]) == (int)'y');
  }
}

void CQscReader::ReadPotentialOutputParameters(bool &savePotential, bool &saveTotalPotential, 
                                               bool &plotPotential)
{
	if (readparam(m_fp, "save potential:",buf,1)) {
		sscanf(buf," %s",answer);
		savePotential = (tolower(answer[0]) == (int)'y');
	}  
	if (readparam(m_fp, "save projected potential:",buf,1)) {
		sscanf(buf," %s",answer);
		saveTotalPotential = (tolower(answer[0]) == (int)'y');
	}  
	if (readparam(m_fp, "plot V(r)*r:",buf,1)) {
		sscanf(buf," %s",answer);
		plotPotential = (tolower(answer[0]) == (int)'y');
	}  
}

void CQscReader::ReadPotentialCalculationParameters(bool &fftPotential, bool &potential3D)
{
  if (readparam(m_fp, "one time integration:",buf,1)) {
    sscanf(buf,"%s",answer);
    fftPotential = (tolower(answer[0]) == (int)'y');
  }
  if (readparam(m_fp, "potential3D:",buf,1)) {
    sscanf(buf,"%s",answer);
    potential3D = (tolower(answer[0]) == (int)'y');
  }
}

void CQscReader::ReadAverageParameters(unsigned &avgRuns, bool &storeSeries)
{
	if (readparam(m_fp, "Runs for averaging:",buf,1))
		sscanf(buf,"%d",&(avgRuns));

	if (readparam(m_fp, "Store TDS diffr. patt. series:",buf,1)) {
		sscanf(buf,"%s",answer);
		storeSeries = (tolower(answer[0]) == (int)'y');
	}  
}

void CQscReader::ReadScanParameters(float_tt &scanXStart, float_tt &scanXStop, unsigned &scanXN,
                                    float_tt &scanYStart, float_tt &scanYStop, unsigned &scanYN)
{
  if (!readparam(m_fp, "scan_x_start:",buf,1)) sscanf(buf,"%g",&(scanXStart));
  if (!readparam(m_fp, "scan_x_stop:",buf,1)) sscanf(buf,"%g",&(scanXStop));
  if (!readparam(m_fp, "scan_x_pixels:",buf,1)) sscanf(buf,"%d",&(scanXN));
  if (!readparam(m_fp, "scan_y_start:",buf,1)) sscanf(buf,"%g",&(scanYStart));
  if (!readparam(m_fp, "scan_y_stop:",buf,1)) sscanf(buf,"%g",&(scanYStop));
  if (!readparam(m_fp, "scan_y_pixels:",buf,1)) sscanf(buf,"%d",&(scanYN));
}


void CQscReader::ReadOutputName(std::string &fileOrFolderName)
{
  /**********************************************************************
   * Parameters for image display and directories, etc.
   */
  fileOrFolderName = "data";
  //sprintf(fileOrFolderName,"data");
  if (readparam(m_fp, "Folder:",buf,1)) 
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
        if (readparam(m_fp, "atom radius:",buf,1))  
		sscanf(buf,"%g",&(radius)); /* in A */
}

void CQscReader::ReadStructureFactorType(int &scatFactor)
{
  if (readparam(m_fp, "Structure Factors:",buf,1)) {
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
  resetParamFile(m_fp);
  if (readparam(m_fp, "Pendelloesung plot:",buf,1)) {
    sscanf(buf,"%s",answer);
    lbeams = (tolower(answer[0]) == (int)'y');
  }
  if (lbeams) {
    while (readparam(m_fp, "beam:",buf,0)) nbout++;  
    printf("will record %d beams\n",nbout);
    hbeams.resize(nbout);
    kbeams.resize(nbout);
    /* now read in the list of detectors: */
    resetParamFile(m_fp);
    for (unsigned i=0;i<nbout;i++) {
      if (!readparam(m_fp, "beam:",buf,0)) break;
      hbeams[i] = 0;
      kbeams[i] = 0;
      sscanf(buf,"%d %d",&(hbeams[i]),&(kbeams[i]));
    }
  }
}

void CQscReader::ReadNumberOfDetectors(int &numDetectors)
{
  numDetectors=0;
  resetParamFile(m_fp);
  while (readparam(m_fp, "detector:",buf, 0)) ++numDetectors;
}

void CQscReader::ReadDetectorParameters(int det_idx, float_tt &rInside, float_tt &rOutside, std::string &name, 
                                        float_tt &shiftX, float_tt &shiftY)
{
  int file_det_idx=0;
  char name_buf[100];
  resetParamFile(m_fp);
  while (readparam(m_fp, "detector:",buf,0)) {
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
  if (readparam(m_fp, "beam current:",buf,1)) { 
    sscanf(buf,"%g",&(beamCurrent)); /* in pA */
  }
  if (readparam(m_fp, "dwell time:",buf,1)) { 
    sscanf(buf,"%g",&(dwellTimeMs)); /* in msec */
  }
}

void CQscReader::ReadProbeParameters(float_tt &dE_E, float_tt &dI_I, float_tt &dV_V, float_tt &alpha, float_tt &aAIS,
                                     float_tt &sourceRadius, bool &ismoth, float_tt &gaussScale, bool &gaussFlag)
{
  /**********************************************************************
   * Read STEM/CBED probe parameters 
   */
	if (readparam(m_fp, "dE/E:",buf,1))
		dE_E = atof(buf);
	if (readparam(m_fp, "dI/I:",buf,1))
		dI_I = atof(buf);
	if (readparam(m_fp, "dV/V:",buf,1))
		dV_V = atof(buf);

	if (!readparam(m_fp, "alpha:",buf,1)) exit(0); 
	sscanf(buf,"%g",&(alpha)); /* in mrad */

	if (readparam(m_fp, "AIS aperture:",buf,1)) 
		sscanf(buf,"%g",&(aAIS)); /* in A */

	
	//////////////////////////////////////////////////////////////////////

	if (readparam(m_fp, "Source Size (diameter):",buf,1)) 
		sourceRadius = atof(buf)/2.0;

	if (readparam(m_fp, "smooth:",buf,1)) sscanf(buf,"%s",answer);
	ismoth = (tolower(answer[0]) == (int)'y');
	if (readparam(m_fp, "gaussian:",buf,1)) {
		sscanf(buf,"%s %g",answer,&(gaussScale));
		gaussFlag = (tolower(answer[0]) == (int)'y');
	}
}

void CQscReader::ReadTomoParameters(float_tt &tomoTilt, float_tt &tomoStart, float_tt &tomoStep, int &tomoCount,
                     float_tt &zoomFactor)
{
  // in case this file has been written by the tomography function, read the current tilt:
  if (readparam(m_fp, "tomo tilt:",buf,1)) { 
    sscanf(buf,"%g %s",&(tomoTilt),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      tomoTilt *= 1000*PI/180.0;
  }

  if (readparam(m_fp, "tomo start:",buf,1)) { 
    sscanf(buf,"%g %s",&(tomoStart),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      tomoStart *= 1000*PI/180.0;
  }
  if (readparam(m_fp, "tomo step:",buf,1)) {
    sscanf(buf,"%g %s",&(tomoStep),answer); /* in mrad */
    if (tolower(answer[0]) == 'd')
      tomoStep *= 1000*PI/180.0;
  }
  
  if (readparam(m_fp, "tomo count:",buf,1))  
    tomoCount = atoi(buf); 
  if (readparam(m_fp, "zoom factor:",buf,1))  
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
  if (!readparam(m_fp, "Cs:",buf,1))  exit(0); 
  sscanf(buf,"%g",&(Cs)); /* in mm */

  if (readparam(m_fp, "C5:",buf,1)) { 
    sscanf(buf,"%g",&(C5)); /* in mm */
  }
  if (readparam(m_fp, "Cc:",buf,1))
    Cc = atof(buf);

  if (readparam(m_fp, "defocus:",buf,1)) { 
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
  if (readparam(m_fp, "astigmatism:",buf,1)) sscanf(buf,"%g",&(astig)); 

  if (readparam(m_fp, "a_33:",buf,1)) {sscanf(buf,"%g",&(a33)); }
  if (readparam(m_fp, "a_31:",buf,1)) {sscanf(buf,"%g",&(a31)); }
  if (readparam(m_fp, "a_44:",buf,1)) {sscanf(buf,"%g",&(a44)); }
  if (readparam(m_fp, "a_42:",buf,1)) {sscanf(buf,"%g",&(a42)); }
  if (readparam(m_fp, "a_55:",buf,1)) {sscanf(buf,"%g",&(a55)); }
  if (readparam(m_fp, "a_53:",buf,1)) {sscanf(buf,"%g",&(a53)); }
  if (readparam(m_fp, "a_51:",buf,1)) {sscanf(buf,"%g",&(a51)); }
  if (readparam(m_fp, "a_66:",buf,1)) {sscanf(buf,"%g",&(a66)); }
  if (readparam(m_fp, "a_64:",buf,1)) {sscanf(buf,"%g",&(a64)); }
  if (readparam(m_fp, "a_62:",buf,1)) {sscanf(buf,"%g",&(a62)); }
}

void CQscReader::ReadAberrationAngles(float_tt &astig,
                       float_tt &phi33, float_tt &phi31,
                       float_tt &phi44, float_tt &phi42,
                       float_tt &phi55, float_tt &phi53, float_tt &phi51,
                       float_tt &phi66, float_tt &phi64, float_tt &phi62)
{
  if (readparam(m_fp, "astigmatism angle:",buf,1)) sscanf(buf,"%g",&(astig)); 

  if (readparam(m_fp, "phi_33:",buf,1)) {sscanf(buf,"%g",&(phi33)); }
  if (readparam(m_fp, "phi_31:",buf,1)) {sscanf(buf,"%g",&(phi31)); }
  if (readparam(m_fp, "phi_44:",buf,1)) {sscanf(buf,"%g",&(phi44)); }
  if (readparam(m_fp, "phi_42:",buf,1)) {sscanf(buf,"%g",&(phi42)); }
  if (readparam(m_fp, "phi_55:",buf,1)) {sscanf(buf,"%g",&(phi55)); }
  if (readparam(m_fp, "phi_53:",buf,1)) {sscanf(buf,"%g",&(phi53)); }
  if (readparam(m_fp, "phi_51:",buf,1)) {sscanf(buf,"%g",&(phi51)); }
  if (readparam(m_fp, "phi_66:",buf,1)) {sscanf(buf,"%g",&(phi66)); }
  if (readparam(m_fp, "phi_64:",buf,1)) {sscanf(buf,"%g",&(phi64)); }
  if (readparam(m_fp, "phi_62:",buf,1)) {sscanf(buf,"%g",&(phi62)); }
}
