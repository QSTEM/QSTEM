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

#include "tem.hpp"

CExperimentTEM::CExperimentTEM(const ConfigReaderPtr &configReader) : CExperimentBase(configReader)
	, m_mode("TEM")
{
}

void CExperimentTEM::Run()
{
  const double pi=3.1415926535897;
  int ix,iy,i,pCount,result;
  FILE *avgFp,*fp; // *fpPos=0;
  double timer,timerTot;
  double x,y,ktx,kty;
  char buf[256];//,avgName[256],systStr[512];
  std::string comment;
  float_tt t;
  float_tt **avgPendelloesung = NULL;
  int oldMulsRepeat1 = 1;
  int oldMulsRepeat2 = 1;
  long iseed=0;
  std::map<std::string, double> params;
  fftwf_complex **imageWave = NULL;

  unsigned nx, ny;
  m_wave->GetSizePixels(nx, ny);

  float_tt dx, dy;
  m_wave->GetResolution(dx, dy);

  if (iseed == 0) iseed = -(long) time( NULL );

  m_chisq.resize(m_avgRuns);

  if (m_lbeams) {
    m_pendelloesung = NULL;
    if (avgPendelloesung == NULL) {
      avgPendelloesung = float2D(m_nbout,
                                 m_slices*oldMulsRepeat1*oldMulsRepeat2*m_cellDiv,
                                 "pendelloesung");
    }	  
  }

  timerTot = 0; /* cputim();*/
  DisplayProgress(-1);
  for (m_avgCount = 0;m_avgCount < m_avgRuns;m_avgCount++) {
    m_totalSliceCount = 0;
    
    pCount = 0;

    /* make incident probe wave function with probe exactly in the center */
    /* if the potential array is not big enough, the probe can 
     * then also be adjusted, so that it is off-center
     */

    // m_scanXStart = m_nx/2*m_resolutionX+m_sourceRadius*gasdev(&iseed)*sqrt(2);
    // m_scanYStart = m_ny/2*m_resolutionY+m_sourceRadius*gasdev(&iseed)*sqrt(2);
    // probe(&muls,m_scanXStart,m_scanYStart);

    //m_nslic0 = 0;
    // produce an incident plane wave:
    if ((m_btiltx == 0) && (m_btilty == 0)) {
      for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
          m_wave->m_wave[ix][iy][0] = 1; 
          m_wave->m_wave[ix][iy][1] = 0;
        }
    }
    else {
      // produce a tilted wave function (btiltx,btilty):
      ktx = 2.0*pi*sin(m_btiltx)/m_wave->GetWavelength();
      kty = 2.0*pi*sin(m_btilty)/m_wave->GetWavelength();
      for (ix=0;ix<nx;ix++) {
        x = dx*(ix-nx/2);
        for (iy=0;iy<ny;iy++) {
          y = dy*(iy-ny/2);
          m_wave->m_wave[ix][iy][0] = (float)cos(ktx*x+kty*y);	
          m_wave->m_wave[ix][iy][1] = (float)sin(ktx*x+kty*y);
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

      m_mulsRepeat1 = 1;
      m_mulsRepeat2 = 1;
      sscanf(buf,"%d %d",&m_mulsRepeat1,&m_mulsRepeat2);
      for (i=0;i<(int)strlen(buf);i++)
        buf[i] = 0;
      if (m_mulsRepeat2 < 1) m_mulsRepeat2 = 1;
      sprintf(m_cin2,"%d",m_mulsRepeat1);
      /***********************************************************
       * make sure we have enough memory for the pendelloesung plot
       */
      if ((m_lbeams) && 
          ((oldMulsRepeat1 !=m_mulsRepeat1) ||
           (oldMulsRepeat2 !=m_mulsRepeat2))) {
        oldMulsRepeat1 = m_mulsRepeat1;
        oldMulsRepeat2 = m_mulsRepeat2;
        if (m_pendelloesung != NULL)  free(m_pendelloesung[0]);
        free(avgPendelloesung[0]);
        m_pendelloesung = NULL;
        avgPendelloesung = float2D(m_nbout,
                                   m_slices*oldMulsRepeat1*oldMulsRepeat2*m_cellDiv,
                                   "pendelloesung");
      }
      /*********************************************************/
      
      // printf("Stacking sequence: %s\n",buf);

      if (m_equalDivs) {
        if (m_printLevel > 1) printf("found equal unit cell divisions\n");
        m_pot->Refresh();
      }

      m_saveFlag = 0;
      /****************************************
       * do the (small) loop through the slabs
       *****************************************/
      for (pCount = 0;pCount<m_mulsRepeat2*m_cellDiv;pCount++) {

        /*******************************************************
         * build the potential slices from atomic configuration
         ******************************************************/
        // if ((m_tds) || (m_nCellZ % m_cellDiv != 0)) {
        if (!m_equalDivs) {
          m_pot->Refresh();
        }

        timer = cputim();
        RunMuls(); 
        m_totalSliceCount += m_slices;

        if (m_printLevel > 0) {
          printf("t=%gA, int.=%g time: %gsec (avgCount=%d)\n",
                 m_wave->thickness,m_wave->intIntensity,cputim()-timer,m_avgCount);
        }

        /***************** FOR DEBUGGING ****************/		
        if ((m_avgCount == 0) && (m_saveLevel >=0) && (pCount+1==m_mulsRepeat2*m_cellDiv)) {
          if (m_tds) comment = "Test wave function for run 0";
          else comment = "Exit face wave function for no TDS";
          if ((m_tiltBack) && ((m_btiltx != 0) || (m_btilty != 0))) {
            ktx = -2.0*pi*sin(m_btiltx)/m_wave->GetWavelength();
            kty = -2.0*pi*sin(m_btilty)/m_wave->GetWavelength();
            for (ix=0;ix<nx;ix++) {
              x = dx*(ix-nx/2);
              for (iy=0;iy<ny;iy++) {
                y = dy*(iy-ny/2);
                m_wave->wave[ix][iy][0] *= cos(ktx*x+kty*y);	
                m_wave->wave[ix][iy][1] *= sin(ktx*x+kty*y);
              }
            }
            if (m_printLevel > 1) printf("** Applied beam tilt compensation **\n");
          }
          
          wave->WriteWave();
        }	
#ifdef VIB_IMAGE_TEST  // doTEM
        if ((m_tds) && (m_saveLevel > 2)) {
          params["HT"] = m_v0;
          params["Cs"] = m_Cs;
          params["Defocus"] = m_df0;
          params["Astigmatism Magnitude"] = m_astigMag;
          params["Astigmatism Angle"] = m_astigAngle;
          params["Focal Spread"] = m_Cc * sqrt(m_dE_E*m_dE_E+m_dV_V*m_dV_V+m_dI_I*m_dI_I);
          params["Convergence Angle"] = m_alpha;
          params["Beam Tilt X"] = m_btiltx;
          params["Beam Tilt Y"] = m_btilty;
          comment = "complex exit face Wave function";

          wave->WriteWave(m_avgCount, comment, params);
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
    
    if (m_avgCount == 0) {
      /***********************************************************
       * Save the diffraction pattern
       **********************************************************/	
      memcpy((void *)wave->avgArray[0],(void *)wave->diffpat[0],(size_t)(nx*ny*sizeof(float_tt)));
      /* move the averaged (raw data) file to the target directory as well */
      wave->WriteAvgArray(m_avgCount+1);

      /***********************************************************
       * Save the Pendelloesung Plot
       **********************************************************/	
      if (m_lbeams) {
        for (iy=0;iy<m_slices*m_mulsRepeat1*m_mulsRepeat2*m_cellDiv;iy++) {
          for (ix=0;ix<m_nbout;ix++) {
            avgPendelloesung[ix][iy] = m_pendelloesung[ix][iy];
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
      if (imageWave == NULL) imageWave = complex2D(nx,ny,"imageWave");
      // multiply wave (in rec. space) with transfer function and write result to imagewave
      fftwf_execute(wave->fftPlanWaveForw);
      for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
          // here, we apply the CTF:
          imageWave[ix][iy][0] = wave->wave[ix][iy][0];
          imageWave[ix][iy][1] = wave->wave[ix][iy][1];
        }
      fftwf_execute_dft(wave->fftPlanWaveInv,imageWave[0],imageWave[0]);
      // get the amplitude squared:
      for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
          wave->diffpat[ix][iy] = imageWave[ix][iy][0]*imageWave[ix][iy][0]+imageWave[ix][iy][1]*imageWave[ix][iy][1];
        }
      wave->WriteWaveIntensity();
      // End of Image writing (if avgCount = 0)
      //////////////////////////////////////////////////////////////////////
      
    } // of if m_avgCount == 0 ...
    else {
      /* 	 readRealImage_old(avgArray,m_nx,m_ny,&t,"diffAvg.img"); */
      m_chisq[m_avgCount-1] = 0.0;
      for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
          t = ((float_tt)m_avgCount*wave->avgArray[ix][iy]+
               wave->diffpat[ix][iy])/((float_tt)(m_avgCount+1));
          m_chisq[m_avgCount-1] += (wave->avgArray[ix][iy]-t)*(wave->avgArray[ix][iy]-t);
          wave->avgArray[ix][iy] = t;
        }
      m_chisq[m_avgCount-1] = m_chisq[m_avgCount-1]/(double)(nx*ny);
      wave->WriteAvgArray(m_avgCount+1);

      /* write the data to a file */
      if ((avgFp = fopen("avgresults.dat","w")) == NULL )
        printf("Sorry, could not open data file for averaging\n");
      else {
        for (ix =0;ix<m_avgCount;ix++) {
          fprintf(avgFp,"%d %g\n",ix+1,m_chisq[ix]);
        }
        fclose(avgFp);
      }
      /*************************************************************/
      
      /***********************************************************
       * Average over the pendelloesung plot as well
       */
      if (m_lbeams) {
        for (iy=0;iy<m_slices*m_mulsRepeat1*m_mulsRepeat2*m_cellDiv;iy++) {
          for (ix=0;ix<m_nbout;ix++) {
            avgPendelloesung[ix][iy] = 
              ((float_tt)m_avgCount*avgPendelloesung[ix][iy]+
               m_pendelloesung[ix][iy])/(float_tt)(m_avgCount+1);
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
      if (imageWave == NULL) imageWave = complex2D(nx,ny,"imageWave");
      // multiply wave (in rec. space) with transfer function and write result to imagewave
#if FLOAT_PRECISION == 1
      fftwf_execute(wave->fftPlanWaveForw);
#elif FLOAT_PRECISION == 2
      fftw_execute(wave->fftPlanWaveForw);
#endif
      for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
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
      for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
          t = ((float_tt)m_avgCount*wave->diffpat[ix][iy]+
               imageWave[ix][iy][0]*imageWave[ix][iy][0]+imageWave[ix][iy][1]*imageWave[ix][iy][1])/(float_tt)(m_avgCount+1);
          wave->diffpat[ix][iy] = t;
        }
      wave->WriteImage();
      // End of Image writing (if avgCount > 0)
      //////////////////////////////////////////////////////////////////////

    } /* else ... if avgCount was greater than 0 */


    /////////////////////////////////////////////////////
    // Save the Pendelloesung plot:
    if (m_lbeams) {
      /**************************************************************
       * The diffraction spot intensities of the selected 
       * diffraction spots are now stored in the 2 dimensional array
       * m_pendelloesung[beam][slice].
       * We can write the array to a file and display it, just for 
       * demonstration purposes
       *************************************************************/
      char avgName[255];
      sprintf(avgName,"%s/pendelloesung.dat",m_folder.c_str());
      if ((fp=fopen(avgName,"w")) !=NULL) {
        printf("Writing Pendelloesung data\n");
        for (iy=0;iy<m_slices*m_mulsRepeat1*m_mulsRepeat2*m_cellDiv;iy++) {
          /* write the thicknes in the first column of the file */
          fprintf(fp,"%g",iy*m_c/((float)(m_slices*m_cellDiv)));
          /* write the beam intensities in the following columns */
          for (ix=0;ix<m_nbout;ix++) {
            // store the AMPLITUDE:
            fprintf(fp,"\t%g",sqrt(avgPendelloesung[ix][iy]/(nx*ny)));
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
  } /* end of for m_avgCount=0.. */  
}

void CExperimentTEM::CollectIntensity(unsigned absoluteSlice)
{
	writeBeams(muls,m_wave,islice, absolute_slice);
}