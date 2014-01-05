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

#include "stem.hpp"

void CExperimentSTEM::run()
{
  int ix=0,iy=0,i,pCount,picts,ixa,iya,totalRuns;
	double timer, total_time=0;
	char buf[BUF_LEN];
	float_tt t;
	static float_tt **avgArray=NULL;
	double collectedIntensity;

	std::vector<WavePtr> waves;
	WavePtr wave;

	//pre-allocate several waves (enough for one row of the scan.  
	for (int th=0; th<omp_get_max_threads(); th++)
	{
		waves.push_back(WavePtr(new WAVEFUNC(muls.nx, muls.ny, muls.resolutionX, muls.resolutionY, muls.input_ext, muls.output_ext)));
	}

	muls.chisq = std::vector<double>(muls.avgRuns);
	totalRuns = muls.avgRuns;
	timer = cputim();

	/* average over several runs of for TDS */
	displayProgress(-1);

	for (muls.avgCount = 0;muls.avgCount < totalRuns; muls.avgCount++) {
		total_time = 0;
		collectedIntensity = 0;
		muls.totalSliceCount = 0;
		muls.dE_E = muls.dE_EArray[muls.avgCount];


		/****************************************
		* do the (big) loop
		*****************************************/
		pCount = 0;
		/* make sure we start at the beginning of the file 
		so we won't miss any line that contains a sequence,
		because we will not do any EOF wrapping
		*/
		resetParamFile();
		while (readparam("sequence: ",buf,0)) {
			if (((buf[0] < 'a') || (buf[0] > 'z')) && 
				((buf[0] < '1') || (buf[0] > '9')) &&
				((buf[0] < 'A') || (buf[0] > 'Z'))) {
					printf("Stacking sequence: %s\n",buf);
					printf("Can only work with old stacking sequence\n");
					break;
			}

			// printf("Stacking sequence: %s\n",buf);

			picts = 0;
			/* for the dislocation models picts will be 1, because the atomcoordinates
			* are expressed explicitly for every atom in the whole specimen
			* For perfect Si samples mulsRepeat1 will remain 1, but picts will give
			* the number of unit cells in Z-direction, where the unit cell is defined by 
			* the unit cell in the cssr file multiplied by NCELLZ.  
			* cellDiv will usually be 1 in that case.
			*/
			sscanf(buf,"%d %d",&muls.mulsRepeat1,&picts);
			for (i=0;i<(int)strlen(buf);i++) buf[i] = 0;
			if (picts < 1) picts = 1;
			muls.mulsRepeat2 = picts;
			sprintf(muls.cin2,"%d",muls.mulsRepeat1);
			/* if the unit cell is divided into slabs, we need to multiply
			* picts by that number
			*/
			if ((picts > 1)&& (muls.cubex >0) && (muls.cubey >0) && (muls.cubez>0)) {
				printf("Warning: cube size of height %gA has been defined, ignoring sequence\n",muls.cubez);
				picts = 1;
			}
			picts *= muls.cellDiv;

			if (muls.equalDivs) {
                          pot->Refresh();
                          timer = cputim();
			}

			/****************************************
			* do the (small) loop over slabs
			*****************************************/
			for (pCount=0;pCount<picts;pCount++) {
				/*******************************************************
				* build the potential slices from atomic configuration
				******************************************************/
				if (!muls.equalDivs) {
                                  pot->Refresh();
                                  timer = cputim();
				}

				muls.complete_pixels=0;
				/**************************************************
				* scan through the different probe positions
				*************************************************/
				// default(none) forces us to specify all of the variables that are used in the parallel section.  
				//    Otherwise, they are implicitly shared (and this was cause of several bugs.)
#pragma omp parallel \
	private(ix, iy, ixa, iya, wave, t, timer) \
  shared(pot, pCount, picts, muls, collectedIntensity, total_time, waves) \
	default(none)
#pragma omp for
				for (i=0; i < (muls.scanXN * muls.scanYN); i++)
				{
					timer=cputim();
					ix = i / muls.scanYN;
					iy = i % muls.scanYN;

					wave = waves[omp_get_thread_num()];

					//printf("Scanning: %d %d %d %d\n",ix,iy,pCount,muls.nx);

					/* if this is run=0, create the inc. probe wave function */
					if (pCount == 0) 
					{
						probe(&muls, wave, muls.nx/2*muls.resolutionX, muls.ny/2*muls.resolutionY);

						// TODO: modifying shared value from multiple threads?
						//muls.nslic0 = 0;
						//wave->thickness = 0.0;
					}
                                          
					else 
					{
        					/* load incident wave function and then propagate it */
                                          
                                          wave->ReadWave(ix, iy); /* this also sets the thickness!!! */
						// TODO: modifying shared value from multiple threads?
						//muls.nslic0 = pCount;
					}
					/* run multislice algorithm
					   and save exit wave function for this position 
					   (done by runMulsSTEM), 
					   but we need to define the file name */
					muls.saveFlag = 1;

					wave->iPosX =(int)(ix*(muls.scanXStop-muls.scanXStart)/
									  ((float)muls.scanXN*muls.resolutionX));
					wave->iPosY = (int)(iy*(muls.scanYStop-muls.scanYStart)/
									   ((float)muls.scanYN*muls.resolutionY));
					if (wave->iPosX > muls.potNx-muls.nx)
					{
						wave->iPosX = muls.potNx-muls.nx;  
					}
					if (wave->iPosY > muls.potNy-muls.ny)
					{
						wave->iPosY = muls.potNy-muls.ny;
					}

					// MCS - update the probe wavefunction with its position

					runMulsSTEM(&muls,wave, pot); 


					/***************************************************************
					* In order to save some disk space we will add the diffraction 
					* patterns to their averages now.  The diffraction pattern 
					* should be stored in wave->diffpat (which each thread has independently), 
					* if collectIntensity() has been executed correctly.
					***************************************************************/

					#pragma omp atomic
					collectedIntensity += wave->intIntensity;

					if (pCount == picts-1)  /* if this is the last slice ... */
					{
                                          if (muls.saveLevel > 0) 
						{
							if (muls.avgCount == 0)  
							{
								// initialize the avgArray from the diffpat
								for (ixa=0;ixa<muls.nx;ixa++) 
								{
									for (iya=0;iya<muls.ny;iya++)
									{
										wave->avgArray[ixa][iya]=wave->diffpat[ixa][iya];
									}
								}
							}
							else 
							{
								// printf("Will read image %d %d\n",muls.nx, muls.ny);	
                                                          wave->ReadAvgArray(ix, iy);
								for (ixa=0;ixa<muls.nx;ixa++) for (iya=0;iya<muls.ny;iya++) {
									t = ((float_tt)muls.avgCount * wave->avgArray[ixa][iya] +
										wave->diffpat[ixa][iya]) / ((float_tt)(muls.avgCount + 1));
									if (muls.avgCount>1)
									{
										#pragma omp atomic
										muls.chisq[muls.avgCount-1] += (wave->avgArray[ixa][iya]-t)*
											(wave->avgArray[ixa][iya]-t);
									}
									wave->avgArray[ixa][iya] = t;
								}
							}
							// Write the array to a file, resize and crop it, 
							wave->WriteAvgArray(ix, iy);
							}	
							else {
								if (muls.avgCount > 0)	muls.chisq[muls.avgCount-1] = 0.0;
							}
					} /* end of if pCount == picts, i.e. conditional code, if this
						  * was the last slice
						  */

					#pragma omp atomic
					++muls.complete_pixels;

					if (muls.displayProgInterval > 0) if ((muls.complete_pixels) % muls.displayProgInterval == 0) 
					{
						#pragma omp atomic
						total_time += cputim()-timer;
						printf("Pixels complete: (%d/%d), int.=%.3f, avg time per pixel: %.2fsec\n",
							muls.complete_pixels, muls.scanXN*muls.scanYN, wave->intIntensity,
							(total_time)/muls.complete_pixels);
						timer=cputim();
					}
				} /* end of looping through STEM image pixels */
				/* save STEM images in img files */
				saveSTEMImages(&muls);
				muls.totalSliceCount += muls.slices;
			} /* end of loop through thickness (pCount) */
		} /* end of  while (readparam("sequence: ",buf,0)) */
		// printf("Total CPU time = %f sec.\n", cputim()-timerTot ); 

		/*************************************************************/
		if (muls.avgCount>1)
			muls.chisq[muls.avgCount-1] = muls.chisq[muls.avgCount-1]/(double)(muls.nx*muls.ny);
		muls.intIntensity = collectedIntensity/(muls.scanXN*muls.scanYN);
		displayProgress(1);
	} /* end of loop over muls.avgCount */
}

void CExperimentSTEM::displayParams()
{
    printf("*\n"
           "* STEM parameters:\n");
    printf("* Maximum scattering angle:  %.0f mrad\n",
           0.5*2.0/3.0*wavelength(muls.v0)/muls.resolutionX*1000);    
    muls.detectors->PrintDetectors();
    
    printf("* Scan window:          (%g,%g) to (%g,%g)A, %d x %d = %d pixels\n",
           muls.scanXStart,muls.scanYStart,muls.scanXStop,muls.scanYStop,
           muls.scanXN,muls.scanYN,muls.scanXN*muls.scanYN);
}
