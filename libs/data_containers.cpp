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

#include "stdio.h"
#include <string.h>
#include "data_containers.hpp"
#include "memory_fftw3.hpp"

WAVEFUNC::WAVEFUNC(int x, int y, float_tt resX, float_tt resY) :
  m_position(std::vector<unsigned>()),
  iPosX(0),
  iPosY(0),
  thickness(0.0),
  nx(x),
  ny(y),
  resolutionX(resX),
  resolutionY(resY),
  m_params(std::map<std::string, double>())
{
	char waveFile[256];
	const char *waveFileBase = "mulswav";

	diffpat = float2D(nx,ny,"diffpat");
	avgArray = float2D(nx,ny,"avgArray");

        // TODO: need to pass file extension through to this constructor
	m_imageIO=ImageIOPtr(new CImageIO(nx, ny));
	
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

	sprintf(waveFile,"%s",waveFileBase);
	strcpy(fileout,waveFile);
	sprintf(fileStart,"mulswav");
}

std::vector<unsigned> WAVEFUNC::GetPositionVector()
{
  return m_position;
}

void WAVEFUNC::WriteWave(const char *fileName, const char *comment,
                         std::map<std::string, double>params)
{
  params["dx"]=resolutionX;
  params["dy"]=resolutionY;
  params["Thickness"]=thickness;
  m_imageIO->WriteComplexImage((void **)wave, fileName, params, std::string(comment), GetPositionVector());
}

void WAVEFUNC::WriteDiffPat(const char *fileName, const char *comment,
                            std::map<std::string, double> params)
{
  params["dx"]=1.0/(nx*resolutionX);
  params["dy"]=1.0/(ny*resolutionY);
  params["Thickness"]=thickness;
  m_imageIO->WriteRealImage((void **)diffpat, fileName, params, std::string(comment), GetPositionVector());
}

  void WAVEFUNC::WriteAvgArray(const char *fileName, const char *comment, 
                               std::map<std::string, double> params)
{
  params["dx"]=1.0/(nx*resolutionX);
  params["dy"]=1.0/(ny*resolutionY);
  params["Thickness"]=thickness;
  m_imageIO->WriteRealImage((void **)avgArray, fileName, params, std::string(comment), GetPositionVector());
}

void WAVEFUNC::SetWavePosition(unsigned posX, unsigned posY)
{
  if (m_position==std::vector<unsigned>())
    {
      m_position.push_back(posX);
      m_position.push_back(posY);
    }
  else
    {
      m_position[0]=posX;
      m_position[1]=posY;
    }
}

void WAVEFUNC::ReadWave(const char *fileName)
{
  m_imageIO->ReadImage((void **)wave, fileName, GetPositionVector());
}

void WAVEFUNC::ReadWave(const char *fileName, unsigned positionx, unsigned positiony)
{
  SetWavePosition(positionx, positiony);
  ReadWave(fileName);
}

void WAVEFUNC::ReadDiffPat(const char *fileName)
{
  m_imageIO->ReadImage((void **)diffpat, fileName, GetPositionVector());
}

void WAVEFUNC::ReadDiffPat(const char *fileName, unsigned positionx, unsigned positiony)
{
  SetWavePosition(positionx, positiony);
  ReadDiffPat(fileName);
}

void WAVEFUNC::ReadAvgArray(const char *fileName)
{
  m_imageIO->ReadImage((void **)avgArray, fileName, GetPositionVector());
}

void WAVEFUNC::ReadAvgArray(const char *fileName, unsigned positionx, unsigned positiony)
{
  SetWavePosition(positionx, positiony);
  ReadAvgArray(fileName);
}

Detector::Detector(int nx, int ny, float_tt resX, float_tt resY) :
  error(0),
  shiftX(0),
  shiftY(0),
  Navg(0),
  thickness(0)
{
  image = float2D(nx,ny,"ADFimag");
  image2 = float2D(nx,ny,"ADFimag");
  
  // TODO: need way of passing file formats into constructor here
  m_imageIO=ImageIOPtr(new CImageIO(nx, ny));
}

    void Detector::WriteImage(const char *fileName, const char *comment, std::map<std::string, double> &params)
{
  
  m_imageIO->WriteRealImage((void **)image, fileName, params, std::string(comment));
}

void Detector::SetThickness(float_tt t)
{
  thickness=t;
}


MULS::MULS():
  cubex(0), cubey(0), cubez(0),
  printLevel(2),
  saveLevel(0),
  lpartl(0),
  atomRadius(5.0),
  saveFlag(0),
  sigmaf(0),
  dfdelt(0),
  acmax(0),
  acmin(0),
  aobj(0),
  Cs(0),
  aAIS(0),
  tomoTilt(0),
  tomoStart(0),
  tomoStep(0),
  tomoCount(0),
  trans(NULL),
  cz(NULL),
  normHolog(0),
  gaussianProp(0),
  kx(NULL), kx2(NULL), ky(NULL), ky2(NULL),
  pendelloesung(NULL),
  onlyFresnel(NULL),
  showPhaseplate(NULL),
  czOffset(NULL)
{

}
