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

void CExperimentTEM::run()
{
  const double pi=3.1415926535897;
  int ix,iy,i,pCount,result;
  FILE *avgFp,*fp; // *fpPos=0;
  double timer,timerTot;
  double x,y,ktx,kty;
  char buf[BUF_LEN];//,avgName[256],systStr[512];
  std::string comment;
  float_tt t;
  float_tt **avgPendelloesung = NULL;
  int oldMulsRepeat1 = 1;
  int oldMulsRepeat2 = 1;
  long iseed=0;
  std::map<std::string, double> params;
	fftwf_complex **imageWave = NULL;

	if (iseed == 0) iseed = -(long) time( NULL );

	muls.chisq=std::vector<double>(muls.avgRuns);

	if (muls.lbeams) {
          muls.pendelloesung = NULL;
          if (avgPendelloesung == NULL) {
            avgPendelloesung = float2D(muls.nbout,
                                       muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv,
                                       "pendelloesung");
          }	  
	}

	timerTot = 0; /* cputim();*/
	displayProgress(-1);
	for (muls.avgCount = 0;muls.avgCount < muls.avgRuns;muls.avgCount++) {
          muls.totalSliceCount = 0;
          
          pCount = 0;

          /* make incident probe wave function with probe exactly in the center */
          /* if the potential array is not big enough, the probe can 
           * then also be adjusted, so that it is off-center
           */

          // muls.scanXStart = muls.nx/2*muls.resolutionX+muls.sourceRadius*gasdev(&iseed)*sqrt(2);
          // muls.scanYStart = muls.ny/2*muls.resolutionY+muls.sourceRadius*gasdev(&iseed)*sqrt(2);
          // probe(&muls,muls.scanXStart,muls.scanYStart);

          //muls.nslic0 = 0;
          // produce an incident plane wave:
          if ((muls.btiltx == 0) && (muls.btilty == 0)) {
            for (ix=0;ix<initialWave->m_nx;ix++) for (iy=0;iy<initialWave->m_ny;iy++) {
                initialWave->m_wave[ix][iy][0] = 1;	initialWave->m_wave[ix][iy][1] = 0;
              }
          }
          else {
            // produce a tilted wave function (btiltx,btilty):
            ktx = 2.0*pi*sin(muls.btiltx)/wavelength(initialWave->m_v0);
            kty = 2.0*pi*sin(muls.btilty)/wavelength(initialWave->m_v0);
            for (ix=0;ix<initialWave->m_nx;ix++) {
              x = initialWave->m_dx*(ix-initialWave->m_nx/2);
              for (iy=0;iy<initialWave->m_ny;iy++) {
                y = initialWave->m_dx*(ix-initialWave->m_nx/2);
                initialWave->m_wave[ix][iy][0] = (float)cos(ktx*x+kty*y);	
                initialWave->m_wave[ix][iy][1] = (float)sin(ktx*x+kty*y);
              }
            }
          }

          result = readparam(fp, "sequence: ",buf,0);
          while (result) {
            if (((buf[0] < 'a') || (buf[0] > 'z')) && 
                ((buf[0] < '1') || (buf[0] > '9')) &&
                ((buf[0] < 'A') || (buf[0] > 'Z'))) {
              // printf("Stacking sequence: %s\n",buf);
              printf("Can only work with old stacking sequence\n");
              break;
            }

            muls.mulsRepeat1 = 1;
            muls.mulsRepeat2 = 1;
            sscanf(buf,"%d %d",&muls.mulsRepeat1,&muls.mulsRepeat2);
            for (i=0;i<(int)strlen(buf);i++)
              buf[i] = 0;
            if (muls.mulsRepeat2 < 1) muls.mulsRepeat2 = 1;
            sprintf(muls.cin2,"%d",muls.mulsRepeat1);
            /***********************************************************
             * make sure we have enough memory for the pendelloesung plot
             */
            if ((muls.lbeams) && 
                ((oldMulsRepeat1 !=muls.mulsRepeat1) ||
                 (oldMulsRepeat2 !=muls.mulsRepeat2))) {
              oldMulsRepeat1 = muls.mulsRepeat1;
              oldMulsRepeat2 = muls.mulsRepeat2;
              if (muls.pendelloesung != NULL)  free(muls.pendelloesung[0]);
              free(avgPendelloesung[0]);
              muls.pendelloesung = NULL;
              avgPendelloesung = float2D(muls.nbout,
                                         muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv,
                                         "pendelloesung");
            }
            /*********************************************************/
            
            // printf("Stacking sequence: %s\n",buf);

            if (muls.equalDivs) {
              if (muls.printLevel > 1) printf("found equal unit cell divisions\n");
              pot->Refresh();
            }

            muls.saveFlag = 0;
            /****************************************
             * do the (small) loop through the slabs
             *****************************************/
            for (pCount = 0;pCount<muls.mulsRepeat2*muls.cellDiv;pCount++) {

              /*******************************************************
               * build the potential slices from atomic configuration
               ******************************************************/
              // if ((muls.tds) || (muls.nCellZ % muls.cellDiv != 0)) {
              if (!muls.equalDivs) {
                pot->Refresh();
              }

              timer = cputim();
              runMulsSTEM(&muls,wave, pot); 
              muls.totalSliceCount += muls.slices;

              if (muls.printLevel > 0) {
                printf("t=%gA, int.=%g time: %gsec (avgCount=%d)\n",
                       wave->thickness,wave->intIntensity,cputim()-timer,muls.avgCount);
              }

              /***************** FOR DEBUGGING ****************/		
              if ((muls.avgCount == 0) && (muls.saveLevel >=0) && (pCount+1==muls.mulsRepeat2*muls.cellDiv)) {
                if (muls.tds) comment = "Test wave function for run 0";
                else comment = "Exit face wave function for no TDS";
                if ((muls.tiltBack) && ((muls.btiltx != 0) || (muls.btilty != 0))) {
                  ktx = -2.0*pi*sin(muls.btiltx)/wavelength(muls.v0);
                  kty = -2.0*pi*sin(muls.btilty)/wavelength(muls.v0);
                  for (ix=0;ix<muls.nx;ix++) {
                    x = muls.resolutionX*(ix-muls.nx/2);
                    for (iy=0;iy<muls.ny;iy++) {
                      y = muls.resolutionY*(ix-muls.nx/2);
                      wave->wave[ix][iy][0] *= cos(ktx*x+kty*y);	
                      wave->wave[ix][iy][1] *= sin(ktx*x+kty*y);
                    }
                  }
                  if (muls.printLevel > 1) printf("** Applied beam tilt compensation **\n");
                }
                
                wave->WriteWave();
              }	
#ifdef VIB_IMAGE_TEST  // doTEM
              if ((muls.tds) && (muls.saveLevel > 2)) {
                params["HT"] = muls.v0;
                params["Cs"] = muls.Cs;
                params["Defocus"] = muls.df0;
                params["Astigmatism Magnitude"] = muls.astigMag;
                params["Astigmatism Angle"] = muls.astigAngle;
                params["Focal Spread"] = muls.Cc * sqrt(muls.dE_E*muls.dE_E+muls.dV_V*muls.dV_V+muls.dI_I*muls.dI_I);
                params["Convergence Angle"] = muls.alpha;
                params["Beam Tilt X"] = muls.btiltx;
                params["Beam Tilt Y"] = muls.btilty;
                comment = "complex exit face Wave function";

                wave->WriteWave(muls.avgCount, comment, params);
              }
#endif 

            } 
            result = readparam("sequence: ",buf,0);
          } 
          /////////////////////////////////////////////////////////////////////////////
          // finished propagating through whole sample, we're at the exit surface now.
          // This means the wave function is used for nothing else than producing image(s)
          // and diffraction patterns.
          //////////////////////////////////////////////////////////////////////////////

          wave->ReadDiffPat();

          if (muls.avgCount == 0) {
            /***********************************************************
             * Save the diffraction pattern
             **********************************************************/	
            memcpy((void *)wave->avgArray[0],(void *)wave->diffpat[0],(size_t)(muls.nx*muls.ny*sizeof(float_tt)));
            /* move the averaged (raw data) file to the target directory as well */
            wave->WriteAvgArray(muls.avgCount+1);

            /***********************************************************
             * Save the Pendelloesung Plot
             **********************************************************/	
            if (muls.lbeams) {
              for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
                for (ix=0;ix<muls.nbout;ix++) {
                  avgPendelloesung[ix][iy] = muls.pendelloesung[ix][iy];
                }
              }
            }
            /***********************************************************
             * Save the defocused image, we can do with the wave what 
             * we want, since it is not used after this anymore. 
             * We will therefore multiply with the transfer function for
             * all the different defoci, inverse FFT and save each image.
             * diffArray will be overwritten with the image.
             **********************************************************/ 
            if (imageWave == NULL) imageWave = complex2D(muls.nx,muls.ny,"imageWave");
            // multiply wave (in rec. space) with transfer function and write result to imagewave
            fftwf_execute(wave->fftPlanWaveForw);
            for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
                // here, we apply the CTF:
                imageWave[ix][iy][0] = wave->wave[ix][iy][0];
                imageWave[ix][iy][1] = wave->wave[ix][iy][1];
              }
            fftwf_execute_dft(wave->fftPlanWaveInv,imageWave[0],imageWave[0]);
            // get the amplitude squared:
            for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
                wave->diffpat[ix][iy] = imageWave[ix][iy][0]*imageWave[ix][iy][0]+imageWave[ix][iy][1]*imageWave[ix][iy][1];
              }
            wave->WriteWaveIntensity();
            // End of Image writing (if avgCount = 0)
            //////////////////////////////////////////////////////////////////////
            
          } // of if muls.avgCount == 0 ...
          else {
            /* 	 readRealImage_old(avgArray,muls.nx,muls.ny,&t,"diffAvg.img"); */
            muls.chisq[muls.avgCount-1] = 0.0;
            for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
                t = ((float_tt)muls.avgCount*wave->avgArray[ix][iy]+
                     wave->diffpat[ix][iy])/((float_tt)(muls.avgCount+1));
                muls.chisq[muls.avgCount-1] += (wave->avgArray[ix][iy]-t)*(wave->avgArray[ix][iy]-t);
                wave->avgArray[ix][iy] = t;
              }
            muls.chisq[muls.avgCount-1] = muls.chisq[muls.avgCount-1]/(double)(muls.nx*muls.ny);
            wave->WriteAvgArray(muls.avgCount+1);

            /* write the data to a file */
            if ((avgFp = fopen("avgresults.dat","w")) == NULL )
              printf("Sorry, could not open data file for averaging\n");
            else {
              for (ix =0;ix<muls.avgCount;ix++) {
                fprintf(avgFp,"%d %g\n",ix+1,muls.chisq[ix]);
              }
              fclose(avgFp);
            }
            /*************************************************************/
            
            /***********************************************************
             * Average over the pendelloesung plot as well
             */
            if (muls.lbeams) {
              for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
                for (ix=0;ix<muls.nbout;ix++) {
                  avgPendelloesung[ix][iy] = 
                    ((float_tt)muls.avgCount*avgPendelloesung[ix][iy]+
                     muls.pendelloesung[ix][iy])/(float_tt)(muls.avgCount+1);
                }
              }
            }
            /***********************************************************
             * Save the defocused image, we can do with the wave what 
             * we want, since it is not used after this anymore. 
             * We will therefore multiply with the transfer function for
             * all the different defoci, inverse FFT and save each image.
             * diffArray will be overwritten with the image.
             **********************************************************/ 
            if (imageWave == NULL) imageWave = complex2D(muls.nx,muls.ny,"imageWave");
            // multiply wave (in rec. space) with transfer function and write result to imagewave
#if FLOAT_PRECISION == 1
            fftwf_execute(wave->fftPlanWaveForw);
#elif FLOAT_PRECISION == 2
            fftw_execute(wave->fftPlanWaveForw);
#endif
            for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
                imageWave[ix][iy][0] = wave->wave[ix][iy][0];
                imageWave[ix][iy][1] = wave->wave[ix][iy][1];
              }
#if FLOAT_PRECISION == 1
            fftwf_execute_dft(wave->fftPlanWaveInv,imageWave[0],imageWave[0]);
#elif FLOAT_PRECISION == 2
            fftw_execute_dft(wave->fftPlanWaveInv,imageWave[0],imageWave[0]);
#endif

            // save the amplitude squared:
            wave->ReadImage();
            for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
                t = ((float_tt)muls.avgCount*wave->diffpat[ix][iy]+
                     imageWave[ix][iy][0]*imageWave[ix][iy][0]+imageWave[ix][iy][1]*imageWave[ix][iy][1])/(float_tt)(muls.avgCount+1);
                wave->diffpat[ix][iy] = t;
              }
            wave->WriteImage();
            // End of Image writing (if avgCount > 0)
            //////////////////////////////////////////////////////////////////////

          } /* else ... if avgCount was greater than 0 */


          /////////////////////////////////////////////////////
          // Save the Pendelloesung plot:
          if (muls.lbeams) {
            /**************************************************************
             * The diffraction spot intensities of the selected 
             * diffraction spots are now stored in the 2 dimensional array
             * muls.pendelloesung[beam][slice].
             * We can write the array to a file and display it, just for 
             * demonstration purposes
             *************************************************************/
            char avgName[255];
            sprintf(avgName,"%s/pendelloesung.dat",muls.folder.c_str());
            if ((fp=fopen(avgName,"w")) !=NULL) {
              printf("Writing Pendelloesung data\n");
              for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
                /* write the thicknes in the first column of the file */
                fprintf(fp,"%g",iy*muls.c/((float)(muls.slices*muls.cellDiv)));
                /* write the beam intensities in the following columns */
                for (ix=0;ix<muls.nbout;ix++) {
                  // store the AMPLITUDE:
                  fprintf(fp,"\t%g",sqrt(avgPendelloesung[ix][iy]/(muls.nx*muls.ny)));
                }
                /* close the line, and start a new one for the next set of
                 * intensities
                 */
                fprintf(fp,"\n");
              }
              fclose(fp);
            }
            else {
              printf("Could not open file for pendelloesung plot\n");
            }	
          } /* end of if lbeams ... */		 
          displayProgress(1);
	} /* end of for muls.avgCount=0.. */  
}
