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

#ifndef WAVEFUNC_H
#define WAVEFUNC_H

#include "stemtypes_fftw3.hpp"
#include "imagelib_fftw3.hpp"
#include "memory_fftw3.hpp"
#include "config_readers.hpp"

void CreateWaveFunctionDataSets(unsigned x, unsigned y, std::vector<unsigned> positions, std::string output_ext);

std::string waveFilePrefix="mulswav";
std::string dpFilePrefix="diff";
std::string avgFilePrefix="diffAvg";
std::string probeFilePrefix="probe_wave";

// a structure for a probe/parallel beam wavefunction.
// Separate from mulsliceStruct for parallelization.
class WAVEFUNC 
{
  // shared pointer to 
  ImageIOPtr m_imageIO;
public:
  std::string fileStart;
  std::string avgName;
  std::string fileout;
  unsigned detPosX, detPosY;
  unsigned iPosX,iPosY;      /* integer position of probe position array */
  unsigned nx, ny;			/* size of diffpat arrays */
  float_tt **diffpat;
  float_tt **avgArray;
  float_tt thickness;
  float_tt intIntensity;
  std::vector<unsigned> m_position;
  std::map<std::string, double> m_params;

  // These are not used for anything aside from when saving files.
  float_tt resolutionX, resolutionY;

  complex_tt  **wave; /* complex wave function */

#if FLOAT_PRECISION == 1
  fftwf_plan fftPlanWaveForw,fftPlanWaveInv;
#else
  fftw_plan fftPlanWaveForw,fftPlanWaveInv;
#endif

public:
  // initializing constructor:
  WAVEFUNC(unsigned nx, unsigned ny, float_tt resX, float_tt resY, std::string input_ext, std::string output_ext);
  WAVEFUNC(ConfigReaderPtr &configReader);
  // define a copy constructor to create new arrays
  WAVEFUNC( WAVEFUNC& other );

  // For CBED ( &TEM? )
  void SetWavePosition(unsigned navg);
  // For STEM
  void SetWavePosition(unsigned posX, unsigned posY);
  std::vector<unsigned> GetPositionVector();

  void CreateDataSets();

  void WriteProbe();
  void WriteWave(std::string comment="Wavefunction", 
                 std::map<std::string, double>params = std::map<std::string, double>());
  void WriteDiffPat(std::string comment="Diffraction Pattern",
                    std::map<std::string, double>params = std::map<std::string, double>());
  void WriteAvgArray(std::string comment="Average Array",
                     std::map<std::string, double>params = std::map<std::string, double>());

  void ReadWave();
  void ReadWave(unsigned posX, unsigned posY);
  void ReadDiffPat();
  void ReadDiffPat(unsigned posX, unsigned posY);
  void ReadAvgArray();
  void ReadAvgArray(unsigned posX, unsigned posY);
};

typedef boost::shared_ptr<WAVEFUNC> WavePtr;

#endif
