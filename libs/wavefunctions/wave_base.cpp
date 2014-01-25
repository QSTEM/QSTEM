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

CBaseWave::CBaseWave(unsigned x, unsigned y, float_tt resX, float_tt resY, 
                   std::string input_ext, std::string output_ext) :
  //m_position(std::vector<unsigned>()),
  m_detPosX(0),
  m_detPosY(0),
  m_nx(x),
  m_ny(y),
  m_dx(resX),
  m_dy(resY),
  m_params(std::map<std::string, double>())
{
  Initialize(".img", ".img");
}

CBaseWave::CBaseWave(const ConfigReaderPtr &configReader)
{
  configReader->ReadProbeArraySize(m_nx, m_ny);
  configReader->ReadResolution(m_dx, m_dy);
  configReader->ReadVoltage(m_v0);
  // TODO: where does this belong?
  //m_electronScale = m_beamCurrent*m_dwellTime*MILLISEC_PICOAMP;

  // TODO: need to figure out how user is going to specify input/output formats
  //CBaseWave(m_nx, m_ny, m_dx, m_dy, ".img", ".img");
  Initialize(".img", ".img");
  printf("Initialized from cfg file");
}

/** Copy constructor - make sure arrays are deep-copied */
CBaseWave::CBaseWave(const WavePtr &other)
{
  // TODO: make sure arrays are deep copied
  other->GetSizePixels(m_nx, m_ny);
  other->GetResolution(m_dx, m_dy);
  m_v0=other->GetVoltage();
  
  Initialize(".img", ".img");
}

CBaseWave::CBaseWave()
{
}


void CBaseWave::Resize(unsigned x, unsigned y)
{
  m_nx=x;
  m_ny=y;
  CreateDataSets();
}

void CBaseWave::CreateDataSets()
{
  m_diffpat.resize(m_nx*m_ny);
  m_wave.resize(m_nx*m_ny);
}

void CBaseWave::Initialize(std::string input_ext, std::string output_ext)
{
  m_wavlen = Wavelength(m_v0);
  //m_wavlen = 12.26/ sqrt( m_v0*1.e3 + m_v0*m_v0*0.9788 );
  // TODO: need to pass file extension through to this constructor
  m_imageIO=ImageIOPtr(new CImageIO(m_nx, m_ny, input_ext, output_ext));
	
  CreateDataSets();

#if FLOAT_PRECISION == 1
  fftwf_complex *ptr = (fftwf_complex *)&m_wave[0];
  m_fftPlanWaveForw = fftwf_plan_dft_2d(m_nx,m_ny,ptr,ptr,FFTW_FORWARD, k_fftMeasureFlag);
  m_fftPlanWaveInv = fftwf_plan_dft_2d(m_nx,m_ny,ptr,ptr,FFTW_BACKWARD, k_fftMeasureFlag);
#else
  fftw_complex *ptr = (fftw_complex *)&m_wave[0];
  m_fftPlanWaveForw = fftw_plan_dft_2d(m_nx,m_ny,&m_wave[0],&m_wave[0],FFTW_FORWARD, k_fftMeasureFlag);
  m_fftPlanWaveInv = fftw_plan_dft_2d(m_nx,m_ny,&m_wave[0],&m_wave[0],FFTW_BACKWARD, k_fftMeasureFlag);
#endif
  InitializeKVectors();
}

 void CBaseWave::InitializeKVectors()
 {
	m_kx.resize(m_nx);
	m_kx2.resize(m_nx);
	m_ky.resize(m_ny);
	m_ky2.resize(m_ny);

  float_tt ax = m_dx*m_nx;
  float_tt by = m_dy*m_ny;
  for(unsigned ixa=0; ixa<m_nx; ixa++) 
    {
  m_kx[ixa] = (ixa>m_nx/2) ? (float_tt)(ixa-m_nx)/ax : 
    (float_tt)ixa/ax;
  m_kx2[ixa] = m_kx[ixa]*m_kx[ixa];
    }
    for(unsigned iya=0; iya<m_ny; iya++) {
      m_ky[iya] = (iya>m_ny/2) ? 
        (float_tt)(iya-m_ny)/by : 
        (float_tt)iya/by;
      m_ky2[iya] = m_ky[iya]*m_ky[iya];
    }
    m_k2max = m_nx/(2.0F*ax);
    if (m_ny/(2.0F*by) < m_k2max ) m_k2max = m_ny/(2.0F*by);
    m_k2max = 2.0/3.0 * m_k2max;
    m_k2max = m_k2max*m_k2max;
}

void CBaseWave::DisplayParams()
{
  printf("* Real space res.:      %gA (=%gmrad)\n",
         1.0/m_k2max,GetWavelength()*m_k2max*1000.0);
  printf("* Reciprocal space res: dkx=%g, dky=%g\n",
         1.0/(m_nx*m_dx),1.0/(m_ny*m_dy));

  printf("* Beams:                %d x %d \n",m_nx,m_ny);  

  printf("* Acc. voltage:         %g (lambda=%gA)\n",m_v0,Wavelength(m_v0));

  if (k_fftMeasureFlag == FFTW_MEASURE)
    printf("* Probe array:          %d x %d pixels (optimized)\n",m_nx,m_ny);
  else
    printf("* Probe array:          %d x %d pixels (estimated)\n",m_nx,m_ny);
  printf("*                       %g x %gA\n",
         m_nx*m_dx,m_ny*m_dy);
}

/*
//TODO: where does this belong?
inline void CBaseWave::GetElectronScale(float_tt &electronScale)
{
  electronScale=m_electronScale;
}
*/

void CBaseWave::GetSizePixels(unsigned &x, unsigned &y) const 
{
  x=m_nx;
  y=m_ny;
}

void CBaseWave::GetResolution(float_tt &x, float_tt &y) const 
{
  x=m_dx;
  y=m_dy;
}

void CBaseWave::GetPositionOffset(unsigned &x, unsigned &y) const 
{
  x=m_detPosX;
  y=m_detPosY;
}

float_tt CBaseWave::GetK2(unsigned ix, unsigned iy) const 
{
  return m_kx2[ix]+m_ky2[iy];
}

float_tt CBaseWave::GetIntegratedIntensity() const 
{
  unsigned px=m_nx*m_ny;
  float_tt intIntensity=0;
  for (unsigned i=0; i<px; i++)
    {
      intIntensity+=m_wave[i][0]*m_wave[i][0] + m_wave[i][1]*m_wave[i][1];
    }
  // TODO: divide by px or not?
  return intIntensity/px;
}

void CBaseWave::ApplyTransferFunction(boost::shared_array<complex_tt> &wave)
{
  // TODO: transfer function should be passed as a 1D vector that is half the size of the wavefunc.
  //       It should be applied by a radial lookup table (with interpolation?)
  //       Alternatively, is it easier to just use a 2D CTF?
  //       Whatever you do, use m_transferFunction as the storage for it.
  if (wave == boost::shared_array<complex_tt>()) wave = complex1D(m_nx*m_ny,"imageWave");
  unsigned px=GetTotalPixels();
      
  // multiply wave (in rec. space) with transfer function and write result to imagewave
  ToFourierSpace();
  for (unsigned i=0;i<px;i++)
    {
      // here, we apply the CTF:
      // 20140110 - MCS - I think this is where Christoph wanted to apply the CTF - nothing is done ATM.

      // TODO: use these for calculating a radius (to get the CTF value from)
      //ix=i%m_nx;
      //iy=i/m_ny;

      wave[i][0] = m_wave[i][0];
      wave[i][1] = m_wave[i][1];
    }
  ToRealSpace();
}

void CBaseWave::_WriteWave(std::string &fileName, std::string comment,
                         std::map<std::string, double>params)
{
  params["dx"]=m_dx;
  params["dy"]=m_dy;
  //params["HT"] = m_v0;
  //params["Cs"] = m_Cs;
  //params["Defocus"] = m_df0;
  //params["Astigmatism Magnitude"] = m_astigMag;
  //params["Astigmatism Angle"] = m_astigAngle;
  //params["Focal Spread"] = m_Cc * sqrt(m_dE_E*m_dE_E+m_dV_V*m_dV_V+m_dI_I*m_dI_I);
  //params["Convergence Angle"] = m_alpha;
  //params["Beam Tilt X"] = m_btiltx;
  //params["Beam Tilt Y"] = m_btilty;
  m_imageIO->WriteComplexImage((void **)&m_wave[0], fileName, params, comment, m_position);
}

void CBaseWave::_WriteDiffPat(std::string &fileName, std::string comment,
                            std::map<std::string, double> params)
{
  params["dx"]=1.0/(m_nx*m_dx);
  params["dy"]=1.0/(m_ny*m_dy);
  m_imageIO->WriteRealImage((void **)&m_diffpat[0], fileName, params, comment, m_position);
}



void CBaseWave::SetWavePosition(unsigned navg)
{
  m_position.resize(1);
  m_position[0]=navg;
}

void CBaseWave::SetWavePosition(unsigned posX, unsigned posY)
{
  m_detPosX=posX;
  m_detPosY=posY;
  m_position.resize(2);
  m_position[0]=posX;
  m_position[1]=posY;
}

void CBaseWave::SetWavePosition(unsigned posX, unsigned posY, unsigned posZ)
{
  m_detPosX=posX;
  m_detPosY=posY;
  m_position.resize(3);
  m_position[0]=posX;
  m_position[1]=posY;
  m_position[2]=posZ;
}

void CBaseWave::ReadWave()
{
  m_position.clear();
  m_imageIO->ReadImage((void **)&m_wave[0], waveFilePrefix, m_position);
}

void CBaseWave::ReadWave(unsigned navg)
{
  SetWavePosition(navg);
  m_imageIO->ReadImage((void **)&m_wave[0], waveFilePrefix, m_position);

}

void CBaseWave::ReadWave(unsigned positionx, unsigned positiony)
{
  SetWavePosition(positionx, positiony);
  m_imageIO->ReadImage((void **)&m_wave[0], waveFilePrefix, m_position);
}

void CBaseWave::ReadDiffPat()
{
  m_position.clear();
  m_imageIO->ReadImage((void **)&m_diffpat[0], dpFilePrefix, m_position);
}

void CBaseWave::ReadDiffPat(unsigned navg)
{
  SetWavePosition(navg);
  m_imageIO->ReadImage((void **)&m_diffpat[0], dpFilePrefix, m_position);
}

void CBaseWave::ReadDiffPat(unsigned positionx, unsigned positiony)
{
  SetWavePosition(positionx, positiony);
  m_imageIO->ReadImage((void **)&m_diffpat[0], dpFilePrefix, m_position);
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

float_tt CBaseWave::Wavelength(float_tt kev)
{
  double w;
  const double emass=510.99906; /* electron rest mass in keV */
  const double hc=12.3984244; /* Planck's const x speed of light*/
  
  /* electron wavelength in Angstroms */
  return hc/sqrt( kev * ( 2*emass + kev ) );
}  /* end wavelength() */

/*
void CBaseWave:WriteBeams(int absolute_slice) {
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
void CBaseWave::ToFourierSpace()
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
void CBaseWave::ToRealSpace()
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
