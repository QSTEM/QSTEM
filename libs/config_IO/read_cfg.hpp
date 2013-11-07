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

#include "read_interface.hpp"

class CCfgReader : public IConfigReader
{
public:
  CCfgReader(std::string filename);
  ~CCfgReader();

  void ReadMode(int &mode);
  void ReadOutputLevel(int &printLevel, int &saveLevel, 
                               unsigned &displayPotCalcInterval, unsigned &displayProgInterval);
  void ReadNCells(unsigned &nCellX, unsigned &nCellY, unsigned &nCellZ, int &cellDiv);
  void ReadBeamTilt(float_tt &btiltx, float_tt &btilty, bool tiltBack);
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
  void ReadAverageParmaeters(unsigned &avgRuns, bool &storeSeries);
  void ReadScanParameters(float_tt &scanXStart, float_tt &scanXStop, unsigned &scanXN,
                                  float_tt &scanYStart, float_tt &scanYStop, unsigned &scanYN);
  void ReadTomoParameters(float_tt &tomoStart, float_tt &tomoStep, int &tomoCount,
                     float_tt &zoomFactor);
  void ReadAberrationAmplitudes(float_tt &Cs, float_tt &C5,
                           float_tt &df0, float_tt &astig,
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
  char buf[1024];
  char answer[256];
  FILE *fpTemp;
};

