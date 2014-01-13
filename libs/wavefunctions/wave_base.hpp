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
public:
  // initializing constructor:
  WAVEFUNC(unsigned nx, unsigned ny, float_tt resX, float_tt resY, std::string input_ext, std::string output_ext);
  WAVEFUNC(const ConfigReaderPtr &configReader);
  // define a copy constructor to create new arrays
  WAVEFUNC( const WAVEFUNC& other );

  void CreateDataSets();
  virtual void FormProbe()=0;

  void DisplayParams();

  void ToRealSpace();
  void ToFourierSpace();
  bool IsRealSpace(){return m_realSpace;}

  //void CopyDPToAvgArray(float_tt *avgArray);
  //void AddDPToAvgArray(unsigned avgCount);

  void GetElectronScale(float_tt &electronScale);
  void GetSizePixels(unsigned &x, unsigned &y);
  unsigned GetTotalPixels(){return m_nx*m_ny;}
  void GetResolution(float_tt &x, float_tt &y);
  void GetPositionOffset(unsigned &x, unsigned &y);
  float_tt GetK2(unsigned ix, unsigned iy);
  inline float_tt GetWavelength() {return m_wavlen;}

  float_tt GetPixelIntensity(unsigned i);
  inline float_tt GetPixelIntensity(unsigned x, unsigned y) {return GetPixelIntensity(x+m_nx*y);}
  inline float_tt GetDiffPatPixel(unsigned i) {return m_diffpat[i];}
  inline float_tt GetDiffPatPixel(unsigned x, unsigned y) { return m_diffpat[x+m_nx*y];}
  //inline float_tt GetAvgArrayPixel(unsigned x, unsigned y) {return m_avgArray[x][y];}
  inline void SetDiffPatPixel(unsigned i, float_tt value) {m_diffpat[i]=value;}
  inline void SetDiffPatPixel(unsigned x, unsigned y, float_tt value) {m_diffpat[x+m_nx*y]=value;}
  //inline void SetAvgArrayPixel(unsigned x, unsigned y, float_tt value) {m_avgArray[x][y]=value;}

  void ApplyTransferFunction(complex_tt *wave);

  void WriteBeams(unsigned absoluteSlice);

  inline void WriteProbe()
  {
    _WriteWave(probeFilePrefix);
  }
  // WriteImage is for TEM mode
  inline void WriteImage()
  {
    _WriteWave(imageFilePrefix, "Image intensity");
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

  // People can change the wavefunction - for example, that's what we have to do when we
  //    transmit the wave through the sample's potential.
  complex_tt *GetWavePointer(){return m_wave;}
  // People should not directly change the diffraction pattern, since we'll re-calculate it when 
  //   the wavefunction changes.
  //   They can, however, access it.
  const float_tt *GetDPPointer(){return m_diffpat;}

  float_tt GetIntegratedIntensity();

  // ReadImage is for TEM mode
  void ReadImage();
  void ReadWave();
  void ReadWave(unsigned navg);
  void ReadWave(unsigned posX, unsigned posY);
  void ReadDiffPat();
  void ReadDiffPat(unsigned navg);
  void ReadDiffPat(unsigned posX, unsigned posY);

protected:
  ImageIOPtr m_imageIO;

  bool m_realSpace;  // If true, the m_wave is in real space.  Else, it's in Fourier space.

  std::string m_fileStart;
  std::string m_avgName;
  std::string m_fileout;
  unsigned m_detPosX, m_detPosY; 
  unsigned m_nx, m_ny;		      /* size of wavefunc and diffpat arrays */
  float_tt *m_diffpat;
  //float_tt **m_avgArray;
  float_tt m_thickness;
  //float_tt m_intIntensity;
  //float_tt m_electronScale;
  //float_tt m_beamCurrent;
  //float_tt m_dwellTime;
  float_tt m_v0;
  std::vector<unsigned> m_position;
  std::map<std::string, double> m_params;

  std::vector<float_tt> m_kx2,m_ky2,m_kx,m_ky;
  std::vector<float_tt> m_propxr, m_propxi, m_propyr, m_propyi;
  float_tt m_k2max;

  float_tt m_aAIS, m_rmin, m_rmax, m_aimin, m_aimax;

  // defocus mode: 1 = Scherzer, 2 = ???
  int m_Scherzer;

  // Coefficients to aberration function:
  float_tt m_a33, m_a31;
  float_tt m_a44, m_a42;
  float_tt m_a55, m_a53, m_a51;
  float_tt m_a66, m_a64, m_a62;
  float_tt m_phi33, m_phi31;
  float_tt m_phi44, m_phi42;
  float_tt m_phi55, m_phi53, m_phi51;
  float_tt m_phi66, m_phi64, m_phi62;

  int m_printLevel;

  float_tt m_C5;
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
  float_tt m_dx, m_dy;

  complex_tt  *m_wave; /* complex wave function */

#if FLOAT_PRECISION == 1
  fftwf_plan m_fftPlanWaveForw,m_fftPlanWaveInv;
#else
  fftw_plan m_fftPlanWaveForw,m_fftPlanWaveInv;
#endif

protected:
  void _WriteWave(std::string &prefix, std::string comment="Wavefunction", 
                 std::map<std::string, double>params = std::map<std::string, double>());
  void _WriteDiffPat(std::string &prefix, std::string comment="Diffraction Pattern",
                    std::map<std::string, double>params = std::map<std::string, double>());
  //void _WriteAvgArray(std::string &prefix, std::string comment="Average Array",
  //                   std::map<std::string, double>params = std::map<std::string, double>());

  // For CBED ( &TEM? )
  void SetWavePosition(unsigned navg);
  // For STEM
  void SetWavePosition(unsigned posX, unsigned posY);

  float_tt Wavelength(float_tt keV);
  float_tt m_wavlen;
  void fft_normalize(void **array,int nx, int ny);

  // m_transferFunction  // The transfer function - optionally applied (used by TEM mode)
};

typedef boost::shared_ptr<WAVEFUNC> WavePtr;

#endif
