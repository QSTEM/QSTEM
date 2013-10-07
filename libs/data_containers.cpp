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
detPosX(0),
detPosY(0),
iPosX(0),
iPosY(0),
thickness(0.0),
nx(x),
ny(y),
resolutionX(resX),
resolutionY(resY)
{
	char waveFile[256];
	const char *waveFileBase = "mulswav";

	diffpat = float2D(nx,ny,"diffpat");
	avgArray = float2D(nx,ny,"avgArray");

	m_imageIO=ImageIOPtr(new CImageIO(nx, ny, thickness, resolutionX, resolutionY));
	
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

	sprintf(waveFile,"%s.img",waveFileBase);
	strcpy(fileout,waveFile);
	sprintf(fileStart,"mulswav.img");
}

std::vector<unsigned long> WAVEFUNC::GetPositionVector()
{
  std::vector<unsigned long> position(2, unsigned long());
  position[0]=detPosX;
  position[1]=detPosY;
  return position;
}

void WAVEFUNC::WriteWave(const char *fileName, const char *comment,
                         std::map<std::string, double>params)
{
	m_imageIO->SetComment(comment);
	m_imageIO->SetResolution(resolutionX, resolutionY);
	m_imageIO->SetParams(params);
	m_imageIO->SetThickness(thickness);
	m_imageIO->WriteComplexImage((void **)wave, fileName, GetPositionVector());
}

void WAVEFUNC::WriteDiffPat(const char *fileName, const char *comment,
                            std::map<std::string, double>params)
{
	m_imageIO->SetComment(comment);
	m_imageIO->SetResolution(1.0/(nx*resolutionX), 1.0/(ny*resolutionY));
	m_imageIO->SetParams(params);
	m_imageIO->SetThickness(thickness);
	m_imageIO->WriteRealImage((void**)diffpat, fileName, GetPositionVector());
}

void WAVEFUNC::WriteAvgArray(const char *fileName, const char *comment,
                             std::map<std::string, double>params)
{
	m_imageIO->SetComment(comment);
	m_imageIO->SetResolution(1.0/(nx*resolutionX), 1.0/(ny*resolutionY));
	m_imageIO->SetParams(params);
	m_imageIO->SetThickness(thickness);
	m_imageIO->WriteRealImage((void **)avgArray, fileName, GetPositionVector());
}

void WAVEFUNC::SetWavePosition(unsigned long posX, unsigned long posY)
{
  detPosX=posX;
  detPosY=posY;
}

void WAVEFUNC::ReadWave(const char *fileName)
{
  m_imageIO->ReadImage((void **)wave, nx, ny, fileName);
}

void WAVEFUNC::ReadWave(const char *fileName, unsigned long positionx, unsigned long positiony)
{
  ReadWave(fileName);
  SetWavePosition(positionx, positiony);
}

void WAVEFUNC::ReadDiffPat(const char *fileName)
{
  m_imageIO->ReadImage((void **)diffpat, nx, ny, fileName);
}

void WAVEFUNC::ReadDiffPat(const char *fileName, unsigned long positionx, unsigned long positiony)
{
  ReadDiffPat(fileName);
  SetWavePosition(positionx, positiony);
}

void WAVEFUNC::ReadAvgArray(const char *fileName)
{
  m_imageIO->ReadImage((void **)avgArray, nx, ny, fileName);
}

void WAVEFUNC::ReadAvgArray(const char *fileName, unsigned long positionx, unsigned long positiony)
{
  ReadAvgArray(fileName);
  SetWavePosition(positionx, positiony);
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

  m_imageIO=ImageIOPtr(new CImageIO(nx, ny, thickness, resX, resY, std::map<std::string, double>(), "STEM image"));
}

void Detector::WriteImage(const char *fileName)
{
	m_imageIO->SetThickness(thickness);
	m_imageIO->WriteRealImage((void **)image, fileName);
}

void Detector::SetThickness(float_tt t)
{
	thickness=t;
}

void Detector::SetParameter(std::string key, double value)
{
	m_imageIO->SetParameter(key, value);
}

 void Detector::SetParams(std::map<std::string, double> &params)
{
	m_imageIO->SetParams(params);
}

void Detector::SetComment(const char *comment)
{
	m_imageIO->SetComment(comment);
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
