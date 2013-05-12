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

MULS::MULS (int slices) : 
atomRadius(5.0), /* radius in A for making the potential boxes */
nlayer(slices),
saveFlag(0),
sigmaf(0), dfdelt(0), acmax(0), acmin(0), aobj(0), Cs(0), aAIS(0),
// Tomography stuff - tomoCount = 0 indicates no Tomo simulation.
tomoTilt(0), tomoStart(0), tomoStep(0), tomoCount(0),
onlyFresnel(0),
czOffset(0), /* defines the offset for the first slice in 
						fractional coordinates        */
normHolog(0),
gaussianProp(0)
{
	int sCount;
	/* general setup: */
	// cin2 should be a vector of strings - with each element being "a" followed by the slice index.  
	//     Is a a character (integer), and this is some kind of offset based on that?
	for (sCount =0;sCount<slices;sCount++)
		cin2[sCount] = 'a'+sCount;
	for (sCount = slices;sCount < NCINMAX;sCount++)
		cin2[sCount] = 0;
	// muls.areaAIS = 1.0;

	/* make multislice read the inout files and assign transr and transi: */
	//muls.trans = NULL;
	//muls.cz = NULL;  // (real *)malloc(muls.slices*sizeof(real));

	sparam = QSfVec::Zero(NPARAM);

	//muls.kx = NULL;
	//muls.kx2= NULL;
	//muls.ky = NULL;
	//muls.ky2= NULL;

	/****************************************************/
	/* copied from slicecell.c                          */
	//muls.pendelloesung = NULL;
}