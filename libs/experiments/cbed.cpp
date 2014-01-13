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

#include "cbed.hpp"
#include "random.hpp"

CExperimentCBED::CExperimentCBED(const ConfigReaderPtr &configReader) : CExperimentBase(configReader)
{
  m_mode="CBED";
}

void CExperimentCBED::Run()
{
  int ix,iy,i,pCount,result;
  FILE *avgFp,*fp,*fpPos=0;
  double timer,timerTot;
  double probeCenterX,probeCenterY,probeOffsetX,probeOffsetY;
  float_tt t=0;
  float_tt **avgPendelloesung = NULL;

  unsigned nx, ny;
  m_wave->GetSizePixels(nx, ny);

  std::map<std::string, double> params;

  std::vector<unsigned> position(1);         // Used to indicate the number of averages

  m_chisq.resize(m_avgRuns);

  if (m_lbeams) {
    m_pendelloesung = NULL;
    if (avgPendelloesung == NULL) {
      avgPendelloesung = float2D(m_nbout,
                                 m_potential->GetNSlices()*m_cellDiv,
                                 "pendelloesung");
    }    
  }
  probeCenterX = m_scanXStart;
  probeCenterY = m_scanYStart;

  timerTot = 0; /* cputim();*/
  DisplayProgress(-1);

  for (m_avgCount = 0;m_avgCount < m_avgRuns;m_avgCount++) {
    m_totalSliceCount = 0;
    pCount = 0;
    
    /* probe(&muls,xpos,ypos); */
    /* make incident probe wave function with probe exactly in the center */
    /* if the potential array is not big enough, the probe can 
     * then also be adjusted, so that it is off-center
     */

    probeOffsetX = m_sourceRadius*gasdev()*SQRT_2;
    probeOffsetY = m_sourceRadius*gasdev()*SQRT_2;
    m_scanXStart = probeCenterX+probeOffsetX;
    m_scanYStart = probeCenterY+probeOffsetY;
    m_wave->FormProbe();
    //probe(&muls, wave,m_scanXStart-m_potOffsetX,m_scanYStart-m_potOffsetY);
    if (m_saveLevel > 2) {
      m_wave->WriteProbe();
    } 	
    // printf("Probe: (%g, %g)\n",m_scanXStart,m_scanYStart);
    /*****************************************************************
     * For debugging only!!!
     *
		muls->WriteWave("probe.img")
    *****************************************************************/


    if (m_sourceRadius > 0) {
      if (m_avgCount == 0) fpPos = fopen("probepos.dat","w");
      else fpPos = fopen("probepos.dat","a");
      if (fpPos == NULL) {
        printf("Was unable to open file probepos.dat for writing\n");
      }
      else {
        fprintf(fpPos,"%g %g\n",m_scanXStart,m_scanYStart);
        fclose(fpPos);
      }
    }

    /***********************************************************
     * make sure we have enough memory for the pendelloesung plot
       */
      if (m_lbeams)
        {
        if (m_pendelloesung != NULL)
          free(m_pendelloesung[0]);
        free(avgPendelloesung[0]);
        m_pendelloesung = NULL;
        avgPendelloesung = float2D(m_nbout,
                                   m_potential->GetNSlices()*m_cellDiv,
                                   "pendelloesung");
      }
      /*********************************************************/
      
      // printf("Stacking sequence: %s\n",buf);

      //m_saveFlag = 0;
      /****************************************
       * do the (small) loop
       *****************************************/
      for (pCount = 0;pCount<m_cellDiv;pCount++) {
        
        m_potential->Refresh();
        
        // what probe should runMulsSTEM use here?
        RunMuls(); 
        
        //printf("Thickness: %gA, int.=%g\n",
        //       m_wave->thickness,m_wave->intIntensity);

        /***************** Only if Save level > 2: ****************/
        if ((m_avgCount == 0) && (m_saveLevel > 2)) {
          m_wave->WriteWave();
        } 	
#ifdef VIB_IMAGE_TEST_CBED
        m_wave->WriteWave()
#endif 
          m_totalSliceCount += m_potential->GetNSlices();
        
      } // end of for pCount = 0... 

    /*    printf("Total CPU time = %f sec.\n", cputim()-timerTot ); */
    
    // TODO: Why are we reading in a DP at this point?  Do we have one yet?  
    //     What happens if it isn't there?
    m_wave->ReadDiffPat();
    AddDPToAvgArray(m_wave);

    
    if (m_avgCount == 0) {
      /* move the averaged (raw data) file to the target directory as well */
      // TODO: make sure that DP average gets created properly
      //sprintf(avgName,"%s/diffAvg_%d.img",m_folder.c_str(),m_avgCount+1);
      //sprintf(systStr,"mv %s/diff.img %s",m_folder.c_str(),avgName);
      //system(systStr);
      if (m_lbeams) {
        for (iy=0;iy<m_potential->GetNSlices()*m_cellDiv;iy++) {
          for (ix=0;ix<m_nbout;ix++) {
            avgPendelloesung[ix][iy] = m_pendelloesung[ix][iy];
          }
        }
      }
    } // of if m_avgCount == 0 ...
    else {
                        
      m_storeSeries = 1;
      if (m_saveLevel == 0)	m_storeSeries = 0;
      else if (m_avgCount % m_saveLevel != 0) m_storeSeries = 0;

      if (m_storeSeries) 
        {
          params["1/Wavelength"] = 1.0/m_wave->GetWavelength();
          WriteAvgArray(m_avgCount+1, "Averaged Diffraction pattern, unit: 1/A", params);
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
    } /* else ... if avgCount was greater than 0 */
    
    if (m_lbeams) {
      /**************************************************************
       * The diffraction spot intensities of the selected 
       * diffraction spots are now stored in the 2 dimensional array
       * m_pendelloesung[beam][slice].
       * We can write the array to a file and display it, just for 
       * demonstration purposes
       *************************************************************/
      char systStr[255];
      sprintf(systStr,"%s/pendelloesung.dat",m_outputLocation.c_str());
      if ((fp=fopen(systStr,"w")) !=NULL) {
        printf("Writing Pendelloesung data\n");
        for (iy=0;iy<m_potential->GetNSlices()*m_cellDiv;iy++) {
          /* write the thicknes in the first column of the file */
          fprintf(fp,"%g",iy*m_potential->GetSliceThickness());//((float)(m_potential->GetNSlices()*m_cellDiv)));
          /* write the beam intensities in the following columns */
          for (ix=0;ix<m_nbout;ix++) {
            fprintf(fp,"\t%g",avgPendelloesung[ix][iy]);
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
  //delete(wave);
}

void CExperimentCBED::CollectIntensity(unsigned absoluteSlice)
{
	WriteBeams(absoluteSlice);
}

void CExperimentCBED::WriteBeams(unsigned int absoluteSlice)
{
  // TODO: this needs to be reconsidered in terms of not using static variables

  /*
  if ((fp1 == NULL) || (fpAmpl == NULL) || (fpPhase == NULL)) {
    scale = 1.0F / ( ((float_tt)m_nx) * ((float_tt)m_ny) );
    hbeam = (*muls).hbeams;
    kbeam = (*muls).kbeams;
    if ((hbeam.empty()) || (kbeam.empty())) {
      printf("ERROR: hbeam or kbeam == NULL!\n");
      exit(0);
    }
      
    sprintf(fileAmpl,"%s/beams_amp.dat",(*muls).folder.c_str());
    sprintf(filePhase,"%s/beams_phase.dat",(*muls).folder.c_str());
    sprintf(fileBeam,"%s/beams_all.dat",(*muls).folder.c_str());
    fp1 = fopen(fileBeam, "w" );
    fpAmpl = fopen( fileAmpl, "w" );
    fpPhase = fopen( filePhase, "w" );
    if(fp1==NULL) {
      printf("can't open file %s\n", fileBeam);
      exit(0);
    }
    if(fpAmpl==NULL) {
      printf("can't open amplitude file %s\n",fileAmpl);
      exit(0);
    }
    if(fpPhase==NULL) {
      printf("can't open phase file %s\n", filePhase);
      exit(0);
    }
    fprintf(fp1, " (h,k) = ");
    for(ib=0; ib<(*muls).nbout; ib++) {
      fprintf(fp1," (%d,%d)", muls->hbeams[ib],  muls->kbeams[ib]);
    }
    fprintf( fp1, "\n" );
    fprintf( fp1, "nslice, (real,imag) (real,imag) ...\n\n");
    for( ib=0; ib<muls->nbout; ib++)
      {
        // printf("beam: %d [%d,%d]",ib,hbeam[ib],kbeam[ib]);			
        if(hbeam[ib] < 0 ) hbeam[ib] = muls->nx + hbeam[ib];
        if(kbeam[ib] < 0 ) kbeam[ib] = muls->ny + kbeam[ib];
        if(hbeam[ib] < 0 ) hbeam[ib] = 0;
        if(kbeam[ib] < 0 ) kbeam[ib] = 0;
        if(hbeam[ib] > muls->nx-1 ) hbeam[ib] = muls->nx-1;
        if(kbeam[ib] > muls->ny-1 ) kbeam[ib] = muls->ny-1;
        // printf(" => [%d,%d] %d %d\n",hbeam[ib],kbeam[ib],muls->nx,muls->ny);			
      }

    // setup of beam files, include the t=0 information 
    fprintf( fpAmpl, "%g",0.0);
    fprintf( fpPhase, "%g",0.0);
    for( ib=0; ib<muls->nbout; ib++) {
      ampl = 0.0;
      if ((hbeam[ib] == 0) && (kbeam[ib]==0))
        ampl = 1.0;
      fprintf(fpAmpl,"\t%g",ampl);
      fprintf(fpPhase,"\t%g",0.0);
    }
    fprintf( fpAmpl, "\n");
    fprintf( fpPhase, "\n");
  } // end of if fp1 == NULL ... i.e. setup 


  zsum += (*muls).cz[ilayer];

  fprintf( fp1, "%g", zsum);
  fprintf( fpAmpl, "%g",zsum);
  fprintf( fpPhase, "%g",zsum);
  for( ib=0; ib<(*muls).nbout; ib++) {
    fprintf(fp1, "\t%g\t%g",
            rPart = scale*(*wave).wave[hbeam[ib]][kbeam[ib]][0],
            iPart = scale*(*wave).wave[hbeam[ib]][kbeam[ib]][1]);
    ampl = (float_tt)sqrt(rPart*rPart+iPart*iPart);
    phase = (float_tt)atan2(iPart,rPart);	
    fprintf(fpAmpl,"\t%g",ampl);
    fprintf(fpPhase,"\t%g",phase);
  }
  fprintf( fp1, "\n");
  fprintf( fpAmpl, "\n");
  fprintf( fpPhase, "\n");
  */
}

void CExperimentCBED::PostSliceProcess(unsigned absoluteSlice)
{
  if (m_saveLevel>1)
    {
      InterimWave(absoluteSlice); 
      // TODO: does CBED actually have detectors?
      //m_detectors->CollectIntensity(m_wave, absoluteSlice);
    }
}
