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

#include "stemtypes_fftw3.hpp"
#include <boost/shared_ptr>
#include <string>

class IConfigReader
{
public:
  virtual void ReadMode(int &mode)=0;
  virtual void ReadOutputLevel(int &printLevel, int &saveLevel, 
                               unsigned &displayPotCalcInterval, unsigned &displayProgInterval);
  // only STEM mode has the displayProgInterval.  Provide this for other modes...
  inline void ReadOutputLevel(int &printLevel, int &saveLevel, unsigned &displayPotCalcInterval)
  {
    unsigned dummy;
    ReadOutputLevel(printLevel, saveLevel, displayPotCalcInterval, dummy);
  }
  virtual void ReadNCells(unsigned &nCellX, unsigned &nCellY, unsigned &nCellZ, int &cellDiv)=0;
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
  virtual void ReadAverageParmaeters(unsigned &avgRuns, bool &storeSeries)=0;
  virtual void ReadScanParameters(float_tt &scanXStart, float_tt &scanXStop, unsigned &scanXN,
                                  float_tt &scanYStart, float_tt &scanYStop, unsigned &scanYN)=0;
  // CBED uses this because stop and npixels are irrelevant
  inline void ReadScanParameters(float_tt &scanXStart, float_tt &scanYStart)
  {
    float_tt dummy_f;
    unsigned dummy_u;
    ReadScanParameters(scanXStart, dummy_f, dummy_u, scanYStart, dummy_f, dummy_u);
  }
  virtual void ReadTomoParameters(float_tt &tomoStart, float_tt &tomoStep, int &tomoCount,
                     float_tt &zoomFactor)=0;
  virtual void ReadAberrationAmplitudes(float_tt &Cs, float_tt &C5,
                           float_tt &df0, float_tt &astig,
                           float_tt &a33, float_tt &a31,
                           float_tt &a44, float_tt &a42,
                           float_tt &a55, float_tt &a53, float_tt &a51,
                           float_tt &a66, float_tt &a64, float_tt &a62)=0;
  virtual void ReadAberrationAngles(float_tt &astig,
                       float_tt &phi33, float_tt &phi31,
                       float_tt &phi44, float_tt &phi42,
                       float_tt &phi55, float_tt &phi53, float_tt &phi51,
                       float_tt &phi66, float_tt &phi64, float_tt &phi62)=0;
};

typedef boost::shared_ptr<IConfigReader> ConfigReaderPtr






