#include "stdio.h"
#include <string.h>
#include "data_containers.h"

WAVEFUNC::WAVEFUNC(int x, int y) :
detPosX(0),
detPosY(0),
iPosX(0),
iPosY(0),
thickness(0.0),
nx(0),
ny(0)
{

  std::string waveFile;
  std::string waveFileBase = "mulswav";
  nx = x;
  ny = y;
  diffpat = QSfMat::Zero(nx, ny);
  avgArray = QSfMat::Zero(nx, ny);
  wave = QScMat::Zero(nx, ny);
  
  /*
#if FLOAT_PRECISION == 1
	fftPlanWaveForw = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD, FFTW_ESTIMATE);
	fftPlanWaveInv = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD, FFTW_ESTIMATE);
#else
	fftPlanWaveForw = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD,
		fftMeasureFlag);
	fftPlanWaveInv = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD,
		fftMeasureFlag);
#endif
  */
        waveFile = waveFileBase+".img";
        fileout = waveFile;
        fileStart = "mulswav.img";
}

WAVEFUNC::WAVEFUNC( WAVEFUNC& other )
{
	iPosX = other.iPosX;
	iPosY = other.iPosY;
	detPosX = other.detPosX;
	detPosY = other.detPosY;

        fileStart = other.fileStart;
        fileout = other.fileout;
        avgName = other.avgName;

	nx = other.nx;
	ny = other.ny;

	diffpat = QSfMat::Zero(nx, ny);
	avgArray = QSfMat::Zero(nx, ny);

        wave = QScMat::Zero(nx,ny);

#if FLOAT_PRECISION == 1
        //	fftPlanWaveForw = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD, FFTW_ESTIMATE);
        //	fftPlanWaveInv = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD, FFTW_ESTIMATE);
#else
        //	fftPlanWaveForw = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD,
        //		fftMeasureFlag);
        //	fftPlanWaveInv = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD,
        //		fftMeasureFlag);
#endif
}

// reset the wave's thickness
void WAVEFUNC::ZeroWave(void)
{

}
