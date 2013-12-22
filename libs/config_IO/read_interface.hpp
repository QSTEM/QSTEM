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

#ifndef CONFIG_READER_INTERFACE_H
#define CONFIG_READER_INTERFACE_H

#include "stemtypes_fftw3.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <vector>

class IConfigReader
{
public:
  virtual void ReadMode(int &mode)=0;
  virtual void ReadPrintLevel(unsigned &printLevel)=0;
  virtual void ReadSaveLevel(unsigned &saveLevel)=0;
  virtual void ReadPotentialOutputInterval(unsigned &displayPotCalcInterval)=0;
  virtual void ReadSTEMProgressInterval(unsigned &displayProgInterval)=0;
  virtual void ReadOutputName(std::string &fileOrFolderName)=0;
  virtual void ReadNCells(unsigned &nCellX, unsigned &nCellY, unsigned &nCellZ)=0;
  virtual void ReadNSubSlabs(unsigned &cellDiv)=0;
  virtual void ReadBeamTilt(float_tt &btiltx, float_tt &btilty, bool tiltBack)=0;
  virtual void ReadCrystalCubeAndTilt(float_tt &tiltx, float_tt &tilty, float_tt &tiltz, 
                                      float_tt &cubex, float_tt &cubey, float_tt &cubez,
                                      bool &adjustCubeSize)=0;
  virtual void ReadTemperatureData(bool &doTDS, float_tt &tdsTemperature, 
                                   std::string &phononFile, bool &useEinstein)=0;
  virtual void ReadSliceOffset(float_tt &xOffset, float_tt &yOffset)=0;
  virtual void ReadProbeArraySize(unsigned &nx, unsigned &ny)=0;
  virtual void ReadResolution(float_tt &resolutionX, float_tt &resolutionY)=0;
  virtual void ReadVoltage(float_tt &voltage)=0;
  virtual void ReadSliceParameters(bool &centerSlices, float_tt &sliceThickness, 
                                   unsigned &nslices, unsigned &outputInterval,
                                   float_tt &zOffset)=0;
  virtual void ReadPeriodicParameters(bool &periodicXY, bool &periodicZ)=0;
  virtual void ReadBandLimitTrans(bool &limit)=0;
  virtual void ReadLoadPotential(bool &loadPotential)=0;
  virtual void ReadPotentialOutputParameters(bool &savePotential, bool &saveProjectedPotential, 
                                             bool &plotPotential)=0;
  virtual void ReadPotentialCalculationParameters(bool &fftPotential, bool &potential3D)=0;
  virtual void ReadAtomRadius(float_tt &radius)=0;
  virtual void ReadStructureFactorType(int &type)=0;
  virtual void ReadPendelloesungParameters(std::vector<int> &hbeams, std::vector<int> &kbeams,
                                           bool &lbeams, unsigned &nbout)=0;
  virtual void ReadAverageParameters(unsigned &avgRuns, bool &storeSeries)=0;
  virtual void ReadScanParameters(float_tt &scanXStart, float_tt &scanXStop, unsigned &scanXN,
                                  float_tt &scanYStart, float_tt &scanYStop, unsigned &scanYN)=0;
  // CBED uses this because stop and npixels are irrelevant
  inline void ReadScanParameters(float_tt &scanXStart, float_tt &scanYStart)
  {
    float_tt dummy_f;
    unsigned dummy_u;
    ReadScanParameters(scanXStart, dummy_f, dummy_u, scanYStart, dummy_f, dummy_u);
  }
  virtual void ReadStructureFileName(boost::filesystem::path &structureFile)=0;

  virtual void ReadNumberOfDetectors(int &numDetectors)=0;
  virtual void ReadDetectorParameters(int det_idx, float_tt &rInside, float_tt &rOutside, std::string &name, 
                              float_tt &shiftX, float_tt &shiftY)=0;
  //void ReadDetectors(std::vector<std::vector<DetectorPtr> > &detectors, std::vector<float_tt> &thicknesses,
  //                 DetectorPtr &detector_to_copy)=0;
  virtual void ReadDoseParameters(float_tt &beamCurrent, float_tt &dwellTimeMs)=0;
  virtual void ReadProbeParameters(float_tt &dE_E, float_tt &dI_I, float_tt &dV_V, float_tt &alpha, float_tt &aAIS,
                           float_tt &sourceRadius, bool &ismoth, float_tt &gaussScale, bool &gaussFlag)=0;
  virtual void ReadTomoParameters(float_tt &tomoTilt, float_tt &tomoStart, float_tt &tomoStep, int &tomoCount,
                     float_tt &zoomFactor)=0;
  virtual void ReadAberrationAmplitudes(float_tt &Cs, float_tt &C5, float_tt &Cc,
                           float_tt &df0, int &Scherzer, float_tt &astig,
                           float_tt &a33, float_tt &a31,
                           float_tt &a44, float_tt &a42,
                           float_tt &a55, float_tt &a53, float_tt &a51,
                           float_tt &a66, float_tt &a64, float_tt &a62)=0;
  virtual void ReadAberrationAngles(float_tt &astig,
                       float_tt &phi33, float_tt &phi31,
                       float_tt &phi44, float_tt &phi42,
                       float_tt &phi55, float_tt &phi53, float_tt &phi51,
                       float_tt &phi66, float_tt &phi64, float_tt &phi62)=0;
  inline bool IsValid() {return m_isValid;}
protected:
  bool m_isValid;
};

typedef boost::shared_ptr<IConfigReader> ConfigReaderPtr;

#endif
