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
#include "imagelib_fftw3.hpp"
#include "memory_fftw3.hpp"
#include "config_readers.hpp"

void CreateWaveFunctionDataSets(unsigned x, unsigned y, std::vector<unsigned> positions, std::string output_ext);

// a structure for a probe/parallel beam wavefunction.
// Separate from mulsliceStruct for parallelization.
class WAVEFUNC 
{
  // shared pointer to 
  ImageIOPtr m_imageIO;
public:
  char fileStart[500];
  char avgName[500];
  char fileout[500];
  unsigned detPosX, detPosY;
  int iPosX,iPosY;      /* integer position of probe position array */
  int nx, ny;			/* size of diffpat arrays */
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
  WAVEFUNC(int nx, int ny, float_tt resX, float_tt resY, std::string input_ext, std::string output_ext);
  WAVEFUNC(ConfigReaderPtr &configReader);
  // define a copy constructor to create new arrays
  WAVEFUNC( WAVEFUNC& other );

  void SetWavePosition(unsigned posX, unsigned posY);
  std::vector<unsigned> GetPositionVector();

  void CreateDataSets();

  void WriteWave(const char *fileName, const char *comment="Wavefunction", 
                 std::map<std::string, double>params = std::map<std::string, double>());
  void WriteDiffPat(const char *fileName, const char *comment="Diffraction Pattern",
                    std::map<std::string, double>params = std::map<std::string, double>());
  void WriteAvgArray(const char *fileName, const char *comment="Average Array",
                     std::map<std::string, double>params = std::map<std::string, double>());

  void ReadWave(const char *fileName);
  void ReadWave(const char *fileName, unsigned posX, unsigned posY);
  void ReadDiffPat(const char *fileName);
  void ReadDiffPat(const char *fileName, unsigned posX, unsigned posY);
  void ReadAvgArray(const char *fileName);
  void ReadAvgArray(const char *fileName, unsigned posX, unsigned posY);
};

typedef boost::shared_ptr<WAVEFUNC> WavePtr;
