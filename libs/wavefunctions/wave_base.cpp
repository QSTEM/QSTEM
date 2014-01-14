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

#include "wavefunctions.hpp"

void CreateWaveFunctionDataSets(unsigned x, unsigned y, std::vector<unsigned> positions, std::string output_ext)
{
  CImageIO imageIO(x, y, "", output_ext);
  std::string potDataSetLabel = "Potential";
  std::string mulswavDataSetLabel = "mulswav";
  imageIO.CreateComplexDataSet(potDataSetLabel, positions);
  imageIO.CreateComplexDataSet(mulswavDataSetLabel, positions);
}

WAVEFUNC::WAVEFUNC(unsigned x, unsigned y, float_tt resX, float_tt resY, 
                   std::string input_ext, std::string output_ext) :
  //m_position(std::vector<unsigned>()),
  m_detPosX(0),
  m_detPosY(0),
  m_thickness(0.0),
  m_nx(x),
  m_ny(y),
  m_dx(resX),
  m_dy(resY),
  m_params(std::map<std::string, double>())
{
  Initialize();
}

WAVEFUNC::WAVEFUNC(const ConfigReaderPtr &configReader)
{
  configReader->ReadProbeArraySize(m_nx, m_ny);
  configReader->ReadResolution(m_dx, m_dy);
  configReader->ReadDoseParameters(m_beamCurrent, m_dwellTime);
  configReader->ReadVoltage(m_v0);
  m_electronScale = m_beamCurrent*m_dwellTime*MILLISEC_PICOAMP;

  // TODO: need to figure out how user is going to specify input/output formats
  //WAVEFUNC(m_nx, m_ny, m_dx, m_dy, ".img", ".img");
  Initialize();
}

void WAVEFUNC::Initialize()
{
  m_wavlen = Wavelength(m_v0);
  //m_wavlen = 12.26/ sqrt( m_v0*1.e3 + m_v0*m_v0*0.9788 );

  m_diffpat = float1D(m_nx*m_ny,"diffpat");

  // TODO: need to pass file extension through to this constructor
  m_imageIO=ImageIOPtr(new CImageIO(m_nx, m_ny, input_ext, output_ext));
	
  m_wave = complex1D(m_nx*m_ny, "wave");
#if FLOAT_PRECISION == 1
  m_fftPlanWaveForw = fftwf_plan_dft_2d(m_nx,m_ny,m_wave,m_wave,FFTW_FORWARD, FFTW_ESTIMATE);
  m_fftPlanWaveInv = fftwf_plan_dft_2d(m_nx,m_ny,m_wave,m_wave,FFTW_BACKWARD, FFTW_ESTIMATE);
#else
  m_fftPlanWaveForw = fftw_plan_dft_2d(m_nx,m_ny,m_wave,m_wave,FFTW_FORWARD,
                                     fftMeasureFlag);
  m_fftPlanWaveInv = fftw_plan_dft_2d(m_nx,m_ny,m_wave,m_wave,FFTW_BACKWARD,
                                    fftMeasureFlag);
#endif
}

void WAVEFUNC::DisplayParams()
{
  printf("* Real space res.:      %gA (=%gmrad)\n",
         1.0/m_k2max,GetWavelength()*m_k2max*1000.0);
  printf("* Reciprocal space res: dkx=%g, dky=%g\n",
         1.0/(m_nx*m_dx),1.0/(m_ny*m_dy));

  printf("* Beams:                %d x %d \n",m_nx,m_ny);  

  printf("* Beam tilt:            x=%g deg, y=%g deg (tilt back == %s)\n",m_btiltx*RAD2DEG,m_btilty*RAD2DEG,
         (m_tiltBack == 1 ? "on" : "off"));
  
  printf("* Aperture half angle:  %g mrad\n",m_alpha);
  printf("* AIS aperture:         ");
  if (m_aAIS > 0) printf("%g A\n",m_aAIS);
  else printf("none\n");
  printf("* beam current:         %g pA\n",m_beamCurrent);
  printf("* dwell time:           %g msec (%g electrons)\n",
         m_dwellTime,m_electronScale);

  printf("* Damping dE/E: %g / %g \n",sqrt(m_dE_E*m_dE_E+m_dV_V*m_dV_V+m_dI_I*m_dI_I)*m_v0*1e3,m_v0*1e3);

  printf("* Acc. voltage:         %g (lambda=%gA)\n",m_v0,Wavelength(m_v0));
  printf("* C_3 (C_s):            %g mm\n",m_Cs*1e-7);
  printf("* C_1 (Defocus):        %g nm%s\n",0.1*m_df0,
         (m_Scherzer == 1) ? " (Scherzer)" : (m_Scherzer==2) ? " (opt.)":"");
  printf("* Astigmatism:          %g nm, %g deg\n",0.1*m_astigMag,RAD2DEG*m_astigAngle);

	// more aberrations:
  if (m_a33 > 0)
    printf("* a_3,3:                %g nm, phi=%g deg\n",m_a33*1e-1,m_phi33*RAD2DEG);
  if (m_a31 > 0)
    printf("* a_3,1:                %g nm, phi=%g deg\n",m_a31*1e-1,m_phi31*RAD2DEG);
  
  if (m_a44 > 0)
    printf("* a_4,4:                %g um, phi=%g deg\n",m_a44*1e-4,m_phi44*RAD2DEG);
  if (m_a42 > 0)
    printf("* a_4,2:                %g um, phi=%g deg\n",m_a42*1e-4,m_phi42*RAD2DEG);

  if (m_a55 > 0)
    printf("* a_5,5:                %g um, phi=%g deg\n",m_a55*1e-4,m_phi55*RAD2DEG);
  if (m_a53 > 0)
    printf("* a_5,3:                %g um, phi=%g deg\n",m_a53*1e-4,m_phi53*RAD2DEG);
  if (m_a51 > 0)
    printf("* a_5,1:                %g um, phi=%g deg\n",m_a51*1e-4,m_phi51*RAD2DEG);

  if (m_a66 > 0)
    printf("* a_6,6:                %g um, phi=%g deg\n",m_a66*1e-7,m_phi66*RAD2DEG);
  if (m_a64 > 0)
    printf("* a_6,4:                %g um, phi=%g deg\n",m_a64*1e-7,m_phi64*RAD2DEG);
  if (m_a62 > 0)
    printf("* a_6,2:                %g um, phi=%g deg\n",m_a62*1e-7,m_phi62*RAD2DEG);
  if (m_C5 != 0)
    printf("* C_5:                  %g mm\n",m_C5*1e-7);
  
  printf("* C_c:                  %g mm\n",m_Cc*1e-7);

  if (k_fftMeasureFlag == FFTW_MEASURE)
    printf("* Probe array:          %d x %d pixels (optimized)\n",m_nx,m_ny);
  else
    printf("* Probe array:          %d x %d pixels (estimated)\n",m_nx,m_ny);
  printf("*                       %g x %gA\n",
         m_nx*m_dx,m_ny*m_dy);
}

inline void WAVEFUNC::GetElectronScale(float_tt &electronScale)
{
  electronScale=m_electronScale;
}

inline void WAVEFUNC::GetSizePixels(unsigned &x, unsigned &y)
{
  x=m_nx;
  y=m_ny;
}

inline void WAVEFUNC::GetResolution(float_tt &x, float_tt &y)
{
  x=m_dx;
  y=m_dy;
}

inline void WAVEFUNC::GetPositionOffset(unsigned &x, unsigned &y)
{
  x=m_detPosX;
  y=m_detPosY;
}

inline float_tt WAVEFUNC::GetK2(unsigned ix, unsigned iy)
{
  return m_kx2[ix]+m_ky2[iy];
}

inline float_tt WAVEFUNC::GetPixelIntensity(unsigned x, unsigned y)
{
  return m_wave[x][y][0]*m_wave[x][y][0] + m_wave[x][y][1]*m_wave[x][y][1];
}

float_tt WAVEFUNC::GetIntegratedIntensity()
{
  unsigned px=m_nx*m_ny;
  float_tt intIntensity=0;
  for (unsigned i=0; i<px; i++)
    {
      intIntensity+=m_wave[x][y][0]*m_wave[x][y][0] + m_wave[x][y][1]*m_wave[x][y][1];
    }
  // TODO: divide by px or not?
  return intIntensity/px;
}

inline void WAVEFUNC::SetDiffPatPixel(unsigned x, unsigned y, float_tt value)
{
  m_diffpat[x][y]=value;
}

void WAVEFUNC::ApplyTransferFunction(complex_tt **wave)
{
  // TODO: transfer function should be passed as a 1D vector that is half the size of the wavefunc.
  //       It should be applied by a radial lookup table (with interpolation?)
  //       Alternatively, is it easier to just use a 2D CTF?
  //       Whatever you do, use m_transferFunction as the storage for it.
  if (wave == NULL) wave = complex2D(nx,ny,"imageWave");
      
  // multiply wave (in rec. space) with transfer function and write result to imagewave
#if FLOAT_PRECISION == 1
  fftwf_execute(m_fftPlanWaveForw);
#elif FLOAT_PRECISION == 2
  fftw_execute(m_fftPlanWaveForw);
#endif
  for (unsigned ix=0;ix<m_nx;ix++) for (unsigned iy=0;iy<m_ny;iy++) {
      // here, we apply the CTF:
      // 20140110 - MCS - I think this is where Christoph wanted to apply the CTF - nothing is done ATM.
      wave[ix][iy][0] = m_wave[ix][iy][0];
      wave[ix][iy][1] = m_wave[ix][iy][1];
    }
#if FLOAT_PRECISION == 1
  fftwf_execute_dft(m_fftPlanWaveInv,wave[0],wave[0]);
#elif FLOAT_PRECISION == 2
  fftw_execute_dft(m_fftPlanWaveInv,wave[0],wave[0]);
#endif
}

void WAVEFUNC::_WriteWave(std::string &fileName, std::string comment,
                         std::map<std::string, double>params)
{
  params["dx"]=m_dx;
  params["dy"]=m_dy;
  params["HT"] = m_v0;
  params["Cs"] = m_Cs;
  params["Defocus"] = m_df0;
  params["Astigmatism Magnitude"] = m_astigMag;
  params["Astigmatism Angle"] = m_astigAngle;
  params["Focal Spread"] = m_Cc * sqrt(m_dE_E*m_dE_E+m_dV_V*m_dV_V+m_dI_I*m_dI_I);
  params["Convergence Angle"] = m_alpha;
  params["Beam Tilt X"] = m_btiltx;
  params["Beam Tilt Y"] = m_btilty;
  m_imageIO->WriteComplexImage((void **)m_wave, fileName, params, comment, m_position);
}

void WAVEFUNC::_WriteDiffPat(std::string &fileName, std::string comment,
                            std::map<std::string, double> params)
{
  params["dx"]=1.0/(m_nx*m_dx);
  params["dy"]=1.0/(m_ny*m_dy);
  m_imageIO->WriteRealImage((void **)m_diffpat, fileName, params, comment, m_position);
}



void WAVEFUNC::SetWavePosition(unsigned navg)
{
  m_position.resize(1);
  m_position[0]=navg;
}

void WAVEFUNC::SetWavePosition(unsigned posX, unsigned posY)
{
  m_detPosX=posX;
  m_detPosY=posY;
  m_position.resize(2);
  m_position[0]=posX;
  m_position[1]=posY;
}

void WAVEFUNC::ReadWave()
{
  m_position.clear();
  m_imageIO->ReadImage((void **)m_wave, waveFilePrefix, m_position);
}

void WAVEFUNC::ReadWave(unsigned navg)
{
  SetWavePosition(navg);
  m_imageIO->ReadImage((void **)m_wave, waveFilePrefix, m_position);

}

void WAVEFUNC::ReadWave(unsigned positionx, unsigned positiony)
{
  SetWavePosition(positionx, positiony);
  m_imageIO->ReadImage((void **)m_wave, waveFilePrefix, m_position);
}

void WAVEFUNC::ReadDiffPat()
{
  m_position.clear();
  m_imageIO->ReadImage((void **)m_diffpat, dpFilePrefix, m_position);
}

void WAVEFUNC::ReadDiffPat(unsigned navg)
{
  SetWavePosition(navg);
  m_imageIO->ReadImage((void **)m_diffpat, dpFilePrefix, m_position);
}

void WAVEFUNC::ReadDiffPat(unsigned positionx, unsigned positiony)
{
  SetWavePosition(positionx, positiony);
  m_imageIO->ReadImage((void **)m_diffpat, dpFilePrefix, m_position);
}

/*--------------------- wavelength() -----------------------------------*/
/*
	return the electron wavelength (in Angstroms)
	keep this is one place so I don't have to keep typing in these
	constants (that I can never remember anyhow)

	ref: Physics Vade Mecum, 2nd edit, edit. H. L. Anderson
		(The American Institute of Physics, New York) 1989
		page 4.

	kev = electron energy in keV

*/

float_tt WAVEFUNC::Wavelength(float_tt kev)
{
  double w;
  const double emass=510.99906; /* electron rest mass in keV */
  const double hc=12.3984244; /* Planck's const x speed of light*/
  
  /* electron wavelength in Angstroms */
  return hc/sqrt( kev * ( 2*emass + kev ) );
}  /* end wavelength() */

/*
void WAVEFUNC:WriteBeams(int absolute_slice) {
  static char fileAmpl[32];
  static char filePhase[32];
  static char fileBeam[32];
  static FILE *fp1 = NULL,*fpAmpl = NULL,*fpPhase=NULL;
  int ib;
  static std::vector<int> hbeam,kbeam;
  static float_tt zsum = 0.0f,scale;
  float_tt rPart,iPart,ampl,phase;
  static char systStr[64];
  // static int counter=0;

  if (!muls->lbeams)
    return;  	
}
*/

// FFT to Fourier space, but only if we're current in real space
void WAVEFUNC::ToFourierSpace()
{
  if (IsRealSpace())
    {
#if FLOAT_PRECISION == 1
      fftwf_execute(m_fftPlanWaveForw);
#elif FLOAT_PRECISION == 2
      fftw_execute(m_fftPlanWaveForw);
#endif
    }
}

// FFT back to realspace, but only if we're currently in Fourier space
void WAVEFUNC::ToRealSpace()
{
  if (!IsRealSpace())
    {
#if FLOAT_PRECISION == 1
      fftwf_execute(m_fftPlanWaveInv);
#elif FLOAT_PRECISION == 2
      fftw_execute(m_fftPlanWaveInv);
#endif
    }
}
