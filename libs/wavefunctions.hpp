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
#include "potential.hpp"

void CreateWaveFunctionDataSets(unsigned x, unsigned y, std::vector<unsigned> positions, std::string output_ext);

static std::string waveFilePrefix="mulswav";
static std::string dpFilePrefix="diff";
static std::string avgFilePrefix="diffAvg";
static std::string probeFilePrefix="probe_wave";
static std::string imageFilePrefix="image";
static std::string waveIntensityFilePrefix="waveIntensity";

// a structure for a probe/parallel beam wavefunction.
// Separate from mulsliceStruct for parallelization.
class WAVEFUNC 
{
  // shared pointer to 
  ImageIOPtr m_imageIO;
public:
  std::string m_fileStart;
  std::string m_avgName;
  std::string m_fileout;
  unsigned m_detPosX, m_detPosY; 
  unsigned m_iPosX,m_iPosY;           /* integer position of probe position array */
  unsigned m_nx, m_ny;		      /* size of wavefunc and diffpat arrays */
  float_tt **m_diffpat;
  float_tt **m_avgArray;
  float_tt m_thickness;
  float_tt m_intIntensity;
  float_tt m_electronScale;
  float_tt m_beamCurrent;
  float_tt m_dwellTime;
  float_tt m_v0;
  std::vector<unsigned> m_position;
  std::map<std::string, double> m_params;

  std::vector<float_tt> m_kx2,m_ky2,m_kx,m_ky;
  std::vector<float_tt> m_propxr, m_propxi, m_propyr, m_propyi;
  float_tt m_k2max;

  float_tt m_aAIS, m_rmin, m_rmax, m_aimin, m_aimax;
protected:
  int m_Scherzer;

  float_tt m_a33, m_a31;
  float_tt m_a44, m_a42;
  float_tt m_a55, m_a53, m_a51;
  float_tt m_a66, m_a64, m_a62;
  float_tt m_phi33, m_phi31;
  float_tt m_phi44, m_phi42;
  float_tt m_phi55, m_phi53, m_phi51;
  float_tt m_phi66, m_phi64, m_phi62;

  int m_printLevel;
public:
  float_tt m_C5;
  float_tt m_dE_E;
  float_tt m_dV_V;
  float_tt m_dI_I;
  float_tt m_alpha;
  float_tt m_sourceRadius;
  float_tt m_Cc;
  float_tt m_Cs;
  float_tt m_df0;				/* defocus */
  float_tt m_astigMag;				/* astigmatism*/
  float_tt m_astigAngle;				/* angle of astigmatism */

  bool m_ismoth;                          /* smoothen the probe wave function */
  bool m_gaussFlag;
  float_tt m_gaussScale;

  // These are not used for anything aside from when saving files.
  float_tt m_resolutionX, m_resolutionY;

  complex_tt  **m_wave; /* complex wave function */

#if FLOAT_PRECISION == 1
  fftwf_plan m_fftPlanWaveForw,m_fftPlanWaveInv;
#else
  fftw_plan m_fftPlanWaveForw,m_fftPlanWaveInv;
#endif

public:
  // initializing constructor:
  WAVEFUNC(unsigned nx, unsigned ny, float_tt resX, float_tt resY, std::string input_ext, std::string output_ext);
  WAVEFUNC(ConfigReaderPtr &configReader);
  // define a copy constructor to create new arrays
  WAVEFUNC( WAVEFUNC& other );

  void CreateDataSets();
  void Transmit(PotPtr pot, unsigned sliceIdx);
  void Propagate(float_tt dz);
  void FormProbe();

  void DisplayParams();

  inline void WriteProbe()
  {
    _WriteWave(probeFilePrefix);
  }
  // WriteImage is for TEM mode
  inline void WriteImage()
  {
    _WriteWave(imageFilePrefix, "Image intensity");
  }
  inline void WriteWaveIntensity()
  {
    _WriteDiffPat(waveIntensityFilePrefix, "Wave intensity");
  }

  inline void WriteWave(std::string comment="Wavefunction", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    m_position.clear();
    _WriteWave(waveFilePrefix, comment, params);
  }
  inline void WriteWave(unsigned navg, std::string comment="Wavefunction", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    SetWavePosition(navg);
    _WriteWave(waveFilePrefix, comment, params);
  }
  inline void WriteWave(unsigned posX, unsigned posY, std::string comment="Wavefunction", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    SetWavePosition(posX, posY);
    _WriteWave(waveFilePrefix, comment, params);
  }

  inline void WriteDiffPat(std::string comment="Diffraction Pattern", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    m_position.clear();
    _WriteDiffPat(dpFilePrefix, comment, params);
  }
  inline void WriteDiffPat(unsigned navg, std::string comment="Diffraction Pattern", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    SetWavePosition(navg);
    _WriteDiffPat(dpFilePrefix, comment, params);
  }
  inline void WriteDiffPat(unsigned posX, unsigned posY, std::string comment="Diffraction Pattern", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    SetWavePosition(posX, posY);
    _WriteDiffPat(dpFilePrefix, comment, params);
  }

  inline void WriteAvgArray(std::string comment="Average Array", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    m_position.clear();
    _WriteAvgArray(avgFilePrefix, comment, params);
  }
  inline void WriteAvgArray(unsigned navg, std::string comment="Average Array", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    SetWavePosition(navg);
    _WriteAvgArray(avgFilePrefix, comment, params);
  }
  inline void WriteAvgArray(unsigned posX, unsigned posY, std::string comment="Average Array", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    SetWavePosition(posX, posY);
    _WriteAvgArray(avgFilePrefix, comment, params);
  }

  // ReadImage is for TEM mode
  void ReadImage();
  void ReadWave();
  void ReadWave(unsigned navg);
  void ReadWave(unsigned posX, unsigned posY);
  void ReadDiffPat();
  void ReadDiffPat(unsigned navg);
  void ReadDiffPat(unsigned posX, unsigned posY);
  void ReadAvgArray();
  void ReadAvgArray(unsigned navg);
  void ReadAvgArray(unsigned posX, unsigned posY);

protected:
  void _WriteWave(std::string &prefix, std::string comment="Wavefunction", 
                 std::map<std::string, double>params = std::map<std::string, double>());
  void _WriteDiffPat(std::string &prefix, std::string comment="Diffraction Pattern",
                    std::map<std::string, double>params = std::map<std::string, double>());
  void _WriteAvgArray(std::string &prefix, std::string comment="Average Array",
                     std::map<std::string, double>params = std::map<std::string, double>());

  // For CBED ( &TEM? )
  void SetWavePosition(unsigned navg);
  // For STEM
  void SetWavePosition(unsigned posX, unsigned posY);


  float_tt Wavelength(float_tt keV);
  float_tt m_wavlen;
};

typedef boost::shared_ptr<WAVEFUNC> WavePtr;

#endif
