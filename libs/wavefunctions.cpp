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
  detPosX(0),
  detPosY(0),
  iPosX(0),
  iPosY(0),
  thickness(0.0),
  nx(x),
  ny(y),
  resolutionX(resX),
  resolutionY(resY),
  m_params(std::map<std::string, double>())
{
  diffpat = float2D(nx,ny,"diffpat");
  avgArray = float2D(nx,ny,"avgArray");

  // TODO: need to pass file extension through to this constructor
  m_imageIO=ImageIOPtr(new CImageIO(nx, ny, input_ext, output_ext));
	
  wave = complex2D(nx, ny, "wave");
#if FLOAT_PRECISION == 1
  fftPlanWaveForw = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD, FFTW_ESTIMATE);
  fftPlanWaveInv = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD, FFTW_ESTIMATE);
#else
  fftPlanWaveForw = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD,
                                     fftMeasureFlag);
  fftPlanWaveInv = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD,
                                    fftMeasureFlag);
#endif
}

WAVEFUNC::WAVEFUNC(ConfigReaderPtr &configReader)
  : detPosX(0)
  , detPosY(0)
  , iPosX(0)
  , iPosY(0)
  , thickness(0)
{
  configReader->ReadProbeArraySize(nx, ny);
  configReader->ReadResolution(resolutionX, resolutionY);
  // TODO: need to figure out how user is going to specify input/output formats
  WAVEFUNC(nx, ny, resolutionX, resolutionY, ".img", ".img");
  
}

void WAVEFUNC::WriteWave(std::string comment,
                         std::map<std::string, double>params)
{
  params["dx"]=resolutionX;
  params["dy"]=resolutionY;
  params["Thickness"]=thickness;
  m_imageIO->WriteComplexImage((void **)wave, waveFilePrefix, params, comment, m_position);
}

void WAVEFUNC::WriteDiffPat(std::string comment,
                            std::map<std::string, double> params)
{
  params["dx"]=1.0/(nx*resolutionX);
  params["dy"]=1.0/(ny*resolutionY);
  params["Thickness"]=thickness;
  m_imageIO->WriteRealImage((void **)diffpat, dpFilePrefix, params, comment, m_position);
}

  void WAVEFUNC::WriteAvgArray(std::string comment, 
                               std::map<std::string, double> params)
{
  params["dx"]=1.0/(nx*resolutionX);
  params["dy"]=1.0/(ny*resolutionY);
  params["Thickness"]=thickness;
  m_imageIO->WriteRealImage((void **)avgArray, avgFilePrefix, params, comment, m_position);
}

void WAVEFUNC::SetWavePosition(unsigned navg)
{
  m_position.resize(1);
  m_position[0]=navg;
}

void WAVEFUNC::SetWavePosition(unsigned posX, unsigned posY)
{
  detPosX=posX;
  detPosY=posY;
  m_position.resize(2);
  m_position[0]=posX;
  m_position[1]=posY;
}

void WAVEFUNC::ReadWave()
{
  m_imageIO->ReadImage((void **)wave, waveFilePrefix, m_position);
}

void WAVEFUNC::ReadWave(unsigned positionx, unsigned positiony)
{
  SetWavePosition(positionx, positiony);
  return ReadWave();
}

void WAVEFUNC::ReadDiffPat()
{
  m_imageIO->ReadImage((void **)diffpat, dpFilePrefix, m_position);
}

void WAVEFUNC::ReadDiffPat(unsigned positionx, unsigned positiony)
{
  SetWavePosition(positionx, positiony);
  return ReadDiffPat();
}

void WAVEFUNC::ReadAvgArray()
{
  m_imageIO->ReadImage((void **)avgArray, avgFilePrefix, m_position);
}

void WAVEFUNC::ReadAvgArray(unsigned positionx, unsigned positiony)
{
  SetWavePosition(positionx, positiony);
  return ReadAvgArray();
}
