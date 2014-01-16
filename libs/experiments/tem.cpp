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
{
  m_mode="TEM";
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
  fftwf_complex *imageWave = NULL;

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
                                 m_potential->GetNSlices()*oldMulsRepeat1*oldMulsRepeat2*m_cellDiv,
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
    m_wave->FormProbe();

    /***********************************************************
     * make sure we have enough memory for the pendelloesung plot
     */
    if (m_lbeams)
      {
        if (m_pendelloesung != NULL)  free(m_pendelloesung[0]);
        free(avgPendelloesung[0]);
        m_pendelloesung = NULL;
        avgPendelloesung = float2D(m_nbout, m_potential->GetNSlices()*m_cellDiv,
                                   "pendelloesung");
      }
      /*********************************************************/
      
      // printf("Stacking sequence: %s\n",buf);

      if (m_equalDivs) {
        if (m_printLevel > 1) printf("found equal unit cell divisions\n");
        m_potential->Refresh();
      }

      /****************************************
       * do the (small) loop through the slabs
       *****************************************/
      for (pCount = 0;pCount<m_cellDiv;pCount++) {

        /*******************************************************
         * build the potential slices from atomic configuration
         ******************************************************/
        // if ((m_tds) || (m_nCellZ % m_cellDiv != 0)) {
        if (!m_equalDivs) {
          m_potential->Refresh();
        }

        RunMuls(m_wave); 
        m_totalSliceCount += m_potential->GetNSlices();

        if (m_printLevel > 0) {
          printf("t=%gA, int.=%g (avgCount=%d)\n",
                 m_thickness,m_intIntensity,m_avgCount);
        }

        /***************** FOR DEBUGGING ****************/		
        if ((m_avgCount == 0) && (m_saveLevel >=0) && (pCount+1==m_cellDiv)) {
          if (m_tds) comment = "Test wave function for run 0";
          else comment = "Exit face wave function for no TDS";
          if (m_tiltBack) m_wave->TiltBack();
          if (m_printLevel > 1) printf("** Applied beam tilt compensation **\n");
          m_wave->WriteWave(comment);
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

          m_wave->WriteWave(m_avgCount, comment, params);
        }
#endif 
  } 
    /////////////////////////////////////////////////////////////////////////////
    // finished propagating through whole sample, we're at the exit surface now.
    // This means the wave function is used for nothing else than producing image(s)
    // and diffraction patterns.
    //////////////////////////////////////////////////////////////////////////////

    m_wave->ReadDiffPat();
    
    if (m_avgCount == 0) {
      /***********************************************************
       * Save the diffraction pattern
       **********************************************************/
      AddDPToAvgArray(m_wave);
      /* move the averaged (raw data) file to the target directory as well */
      WriteAvgArray(m_avgCount+1);

      /***********************************************************
       * Save the Pendelloesung Plot
       **********************************************************/	
      if (m_lbeams) {
        for (iy=0;iy<m_potential->GetNSlices()*m_cellDiv;iy++) {
          for (ix=0;ix<m_nbout;ix++) {
            avgPendelloesung[ix][iy] = m_pendelloesung[ix][iy];
          }
        }
      }
      // End of Image writing (if avgCount = 0)
      //////////////////////////////////////////////////////////////////////
      
    } // of if m_avgCount == 0 ...
    else {
      /* 	 readRealImage_old(avgArray,m_nx,m_ny,&t,"diffAvg.img"); */
      m_chisq[m_avgCount-1] = 0.0;
      AddDPToAvgArray(m_wave);
      WriteAvgArray(m_avgCount+1);

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
      for (iy=0;iy<m_potential->GetNSlices()*m_cellDiv;iy++) {
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
      m_wave->ApplyTransferFunction(imageWave);

      // save the amplitude squared:
      m_wave->ReadImage();
      unsigned px=nx*ny;
      for (unsigned i=0;i<px;i++){
        t = ((float_tt)m_avgCount*m_wave->GetDiffPatPixel(i)+
               imageWave[i][0]*imageWave[i][0]+imageWave[i][1]*imageWave[i][1])/(float_tt)(m_avgCount+1);
          m_wave->SetDiffPatPixel(i,t);
        }
      m_wave->WriteImage();
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
      sprintf(avgName,"%s/pendelloesung.dat",m_outputLocation.c_str());
      if ((fp=fopen(avgName,"w")) !=NULL) {
        printf("Writing Pendelloesung data\n");
        for (iy=0;iy<m_potential->GetNSlices()*m_cellDiv;iy++) {
          /* write the thicknes in the first column of the file */
          fprintf(fp,"%g",iy*m_potential->GetSliceThickness()/((float)(m_potential->GetNSlices()*m_cellDiv)));
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
    DisplayProgress(1);
  } /* end of for m_avgCount=0.. */  
}

void CExperimentTEM::CollectIntensity(unsigned absoluteSlice)
{
  WriteBeams(absoluteSlice);
}

void CExperimentTEM::WriteBeams(unsigned int absoluteSlice)
{
  unsigned nx, ny;
  m_wave->GetSizePixels(nx, ny);
  float_tt scale = 1.0/(nx*ny); 

  if (m_pendelloesung == NULL) 
    {
      m_pendelloesung = float2D(m_nbout, m_potential->GetNSlices()*m_cellDiv,
                              "pendelloesung");
      printf("Allocated memory for pendelloesung plot (%d x %d)\n",
             m_nbout,m_potential->GetNSlices());
    }
    for(unsigned ib=0; ib<m_nbout; ib++) {
      m_pendelloesung[ib][absoluteSlice] = scale*m_wave->GetPixelIntensity(m_hbeams[ib],m_kbeams[ib]);
      // printf("slice: %d beam: %d [%d,%d], intensity: %g\n",muls->nslic0,ib,muls->hbeam[ib],muls->kbeam[ib],muls->pendelloesung[ib][muls->nslic0]);			
    } // end of ib=0 ... 

    // TODO: This isn't actually saving anything...
}


void CExperimentTEM::PostSliceProcess(unsigned absoluteSlice)
{
  InterimWave(absoluteSlice);
}
