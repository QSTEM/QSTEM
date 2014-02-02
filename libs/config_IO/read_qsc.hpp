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

#ifndef READ_QSC_H
#define READ_QSC_H

#include "read_interface.hpp"
#include <boost/filesystem.hpp>

class CQscReader : public IConfigReader
{
public:
  CQscReader(boost::filesystem::path &filename);
  ~CQscReader();

  void ReadMode(std::string &mode);
  void ReadPrintLevel(unsigned &printLevel);
  void ReadSaveLevel(unsigned &saveLevel);
  void ReadPotentialOutputInterval(unsigned &displayPotCalcInterval);
  void ReadSTEMProgressInterval(unsigned &displayProgInterval);
  void ReadOutputName(std::string &fileOrFolderName);
  void ReadNCells(unsigned &nCellX, unsigned &nCellY, unsigned &nCellZ);
  void ReadNSubSlabs(unsigned &cellDiv);
  void ReadBeamTilt(float_tt &btiltx, float_tt &btilty, bool &tiltBack);
  void ReadCrystalCubeAndTilt(float_tt &tiltx, float_tt &tilty, float_tt &tiltz, 
                                      float_tt &cubex, float_tt &cubey, float_tt &cubez,
                                      bool &adjustCubeSize);
  void ReadTemperatureData(bool &doTDS, float_tt &tdsTemperature, 
                                   std::string &phononFile, bool &useEinstein);
  void ReadSliceOffset(float_tt &xOffset, float_tt &yOffset);
  void ReadProbeArraySize(unsigned &nx, unsigned &ny);
  void ReadResolution(float_tt &resolutionX, float_tt &resolutionY);
  void ReadVoltage(float_tt &voltage);
  void ReadSliceParameters(bool &centerSlices, float_tt &sliceThickness, 
                                   unsigned &nslices, unsigned &outputInterval,
                                   float_tt &zOffset);
  void ReadPeriodicParameters(bool &periodicXY, bool &periodicZ);
  void ReadBandLimitTrans(bool &limit);
  void ReadLoadPotential(bool &loadPotential);
  void ReadPotentialOutputParameters(bool &savePotential, bool &saveProjectedPotential, 
                                             bool &plotPotential);
  void ReadPotentialCalculationParameters(bool &fftPotential, bool &potential3D);
  void ReadAverageParameters(unsigned &avgRuns, bool &storeSeries);
  void ReadScanParameters(float_tt &scanXStart, float_tt &scanXStop, unsigned &scanXN,
                                  float_tt &scanYStart, float_tt &scanYStop, unsigned &scanYN);  
  void ReadAtomRadius(float_tt &radius);
  void ReadStructureFactorType(std::string &type);
  void ReadPendelloesungParameters(std::vector<int> &hbeams, std::vector<int> &kbeams,
                                   bool &lbeams, unsigned &nbout);
  void ReadStructureFileName(boost::filesystem::path &structure_file); // 
  
  void ReadNumberOfDetectors(int &numDetectors);
  void ReadDetectorParameters(int det_idx, float_tt &rInside, float_tt &rOutside, std::string &name, 
                              float_tt &shiftX, float_tt &shiftY);
  //void ReadDetectors(std::vector<std::vector<DetectorPtr> > &detectors, std::vector<float_tt> &thicknesses,
  //                           DetectorPtr &detector_to_copy)
  void ReadDoseParameters(float_tt &beamCurrent, float_tt &dwellTimeMs);
  void ReadProbeParameters(float_tt &dE_E, float_tt &dI_I, float_tt &dV_V, float_tt &alpha, float_tt &aAIS);
  void ReadSourceRadius(float_tt &sourceRadius);
  void ReadSmoothingParameters(bool &ismoth, float_tt &gaussScale, bool &gaussFlag);
  void ReadTomoParameters(float_tt &tomoTilt, float_tt &tomoStart, float_tt &tomoStep, int &tomoCount,
                     float_tt &zoomFactor);
  void ReadAberrationAmplitudes(float_tt &Cs, float_tt &C5, float_tt &Cc,
                           float_tt &df0, std::string &Scherzer, float_tt &astig,
                           float_tt &a33, float_tt &a31,
                           float_tt &a44, float_tt &a42,
                           float_tt &a55, float_tt &a53, float_tt &a51,
                           float_tt &a66, float_tt &a64, float_tt &a62);
  void ReadAberrationAngles(float_tt &astig,
                       float_tt &phi33, float_tt &phi31,
                       float_tt &phi44, float_tt &phi42,
                       float_tt &phi55, float_tt &phi53, float_tt &phi51,
                       float_tt &phi66, float_tt &phi64, float_tt &phi62);
protected:
  std::string m_buf;
  char answer[256];
  FILE *m_fp;

  bool IsBufferYes(std::string &buf);
private:
  friend class CConfigReaderFactory;
  // Create an instance of this class, wrapped in a shared ptr
  //     This should not be inherited - any subclass needs its own implementation.
  static ConfigReaderPtr Create(boost::filesystem::path &filename){return ConfigReaderPtr(new CQscReader(filename));}
};

#endif











