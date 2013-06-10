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
#include "data_containers.h"

WAVEFUNC::WAVEFUNC(int x, int y, double resX, double resY) :
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
	

#if FLOAT_PRECISION == 1
	wave = complex2Df(nx, ny, "wave");
	fftPlanWaveForw = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD, FFTW_ESTIMATE);
	fftPlanWaveInv = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD, FFTW_ESTIMATE);
#else
	wave = complex2D(nx, ny, "wave");
	fftPlanWaveForw = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD,
		fftMeasureFlag);
	fftPlanWaveInv = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD,
		fftMeasureFlag);
#endif

	sprintf(waveFile,"%s.img",waveFileBase);
	strcpy(fileout,waveFile);
	sprintf(fileStart,"mulswav.img");
}

/*
WAVEFUNC::WAVEFUNC( WAVEFUNC& other )
{
	iPosX = other.iPosX;
	iPosY = other.iPosY;
	detPosX = other.detPosX;
	detPosY = other.detPosY;

	strcpy(fileStart, other.fileStart);
	strcpy(fileout, other.fileout);
	strcpy(avgName, other.avgName);

	nx = other.nx;
	ny = other.ny;

	resolutionX = other.resolutionX
	resolutionY = other.resolutionY

	diffpat = float2D(nx,ny,"diffpat");
	avgArray = float2D(nx,ny,"avgArray");

#if FLOAT_PRECISION == 1
	wave = complex2Df(nx, ny, "wave");
	fftPlanWaveForw = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD, FFTW_ESTIMATE);
	fftPlanWaveInv = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD, FFTW_ESTIMATE);
#else
	wave = complex2D(nx, ny, "wave");
	fftPlanWaveForw = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD,
		fftMeasureFlag);
	fftPlanWaveInv = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD,
		fftMeasureFlag);
#endif
}
*/

// reset the wave's thickness
void WAVEFUNC::ZeroWave(void)
{

}

void WAVEFUNC::WriteWave(const char *fileName, const char *comment,
	std::vector<double>params)
{
	m_imageIO->SetComment(comment);
	m_imageIO->SetResolution(resolutionX, resolutionY);
	m_imageIO->SetParams(params);
	m_imageIO->SetThickness(thickness);
	m_imageIO->WriteComplexImage((const void**)wave, fileName);
}

void WAVEFUNC::WriteDiffPat(const char *fileName, const char *comment,
	std::vector<double>params)
{
	m_imageIO->SetComment(comment);
	m_imageIO->SetResolution(1.0/(nx*resolutionX), 1.0/(ny*resolutionY));
	m_imageIO->SetParams(params);
	m_imageIO->SetThickness(thickness);
	m_imageIO->WriteRealImage((const void**)diffpat, fileName);
}

void WAVEFUNC::WriteAvgArray(const char *fileName, const char *comment,
	std::vector<double>params)
{
	m_imageIO->SetComment(comment);
	m_imageIO->SetResolution(1.0/(nx*resolutionX), 1.0/(ny*resolutionY));
	m_imageIO->SetParams(params);
	m_imageIO->SetThickness(thickness);
	m_imageIO->WriteRealImage((const void**)avgArray, fileName);
}