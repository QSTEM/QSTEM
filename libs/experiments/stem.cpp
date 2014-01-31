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
#include "wavefunctions/wave_convergent.hpp"

#include <omp.h>

CExperimentSTEM::CExperimentSTEM(const ConfigReaderPtr &configReader) : CExperimentBase(configReader)
{
  m_mode="STEM";
}

void CExperimentSTEM::Run()
{
  int ix=0,iy=0,i,pCount,picts,ixa,iya,totalRuns;
  double timer, total_time=0;
  float_tt t;
  double collectedIntensity;

  std::vector<WavePtr> waves;
  WavePtr wave;

  std::vector<std::vector<float_tt> > avgArrays; //temporary buffers for each thread to have its own avg array

  unsigned nx, ny;
  m_wave->GetSizePixels(nx, ny);

  unsigned potNX, potNY;
  m_potential->GetSizePixels(potNX, potNY);

  float_tt dx, dy;
  m_wave->GetResolution(dx, dy);

  //pre-allocate several waves (enough for one row of the scan.  
  for (int th=0; th<omp_get_max_threads(); th++)
    {
      waves.push_back(WavePtr(new CConvergentWave(m_wave)));
      avgArrays.push_back(std::vector<float_tt>(nx*ny));
    }

  m_chisq.resize(m_avgRuns);
  totalRuns = m_avgRuns;

  /* average over several runs of for TDS */
  DisplayProgress(-1);

  for (m_avgCount = 0;m_avgCount < totalRuns; m_avgCount++) {
    total_time = 0;
    collectedIntensity = 0;
    m_totalSliceCount = 0;
    //m_dE_E = m_dE_EArray[m_avgCount];


    if (m_equalDivs) {
      m_potential->Refresh();
    }

    /****************************************
     * do the (small) loop over slabs
     *****************************************/
    for (pCount=0;pCount<m_cellDiv;pCount++) {
      /*******************************************************
       * build the potential slices from atomic configuration
       ******************************************************/
      if (!m_equalDivs) {
        m_potential->Refresh();
      }

      m_completePixels=0;
      /**************************************************
       * scan through the different probe positions
       *************************************************/
      // default(none) forces us to specify all of the variables that are used in the parallel section.  
      //    Otherwise, they are implicitly shared (and this was cause of several bugs.)
#pragma omp parallel \
  private(ix, iy, ixa, iya, wave, t, timer)                   \
  shared(pCount, picts, collectedIntensity, total_time, waves, avgArrays, dx, dy, nx, ny, potNX, potNY) \
  default(none)
#pragma omp for
      for (i=0; i < (m_scanXN * m_scanYN); i++)
        {
          ix = i / m_scanYN;
          iy = i % m_scanYN;
            
          wave = waves[omp_get_thread_num()];
          m_avgArray = &avgArrays[omp_get_thread_num()][0];
            
            //printf("Scanning: %d %d %d %d\n",ix,iy,pCount,m_nx);
            
            /* if this is run=0, create the inc. probe wave function */
            if (pCount == 0) 
              {
                wave->FormProbe();
              }
                                          
            else 
              {
                /* load incident wave function and then propagate it */
                
                wave->ReadWave(ix, iy); /* this also sets the thickness!!! */
                // TODO: modifying shared value from multiple threads?
                //m_nslic0 = pCount;
              }
            /* run multislice algorithm
               and save exit wave function for this position 
               (done by runMulsSTEM), 
               but we need to define the file name */
            
            m_iPosX =(int)(ix*(m_scanXStop-m_scanXStart)/
                               ((float)m_scanXN*dx));
            m_iPosY = (int)(iy*(m_scanYStop-m_scanYStart)/
                                ((float)m_scanYN*dy));
            if (m_iPosX > potNX-nx)
              {
                m_iPosX = potNX-nx;  
              }
            if (m_iPosY > potNY-ny)
              {
                m_iPosY = potNY-ny;
              }

            // MCS - update the probe wavefunction with its position

            RunMuls(wave); 


            /***************************************************************
             * In order to save some disk space we will add the diffraction 
             * patterns to their averages now.  The diffraction pattern 
             * should be stored in wave->diffpat (which each thread has independently), 
             * if collectIntensity() has been executed correctly.
             ***************************************************************/
            
#pragma omp atomic
            collectedIntensity += wave->GetIntegratedIntensity();
            
            if (pCount == picts-1)  /* if this is the last slice ... */
              {
                if (m_saveLevel > 0) 
                  {
                    ReadAvgArray(ix, iy);
                    AddDPToAvgArray(m_wave);
                    // Write the array to a file, resize and crop it, 
                    WriteAvgArray(ix, iy);
                  }	
                else {
                  if (m_avgCount > 0)	m_chisq[m_avgCount-1] = 0.0;
                }
              } /* end of if pCount == picts, i.e. conditional code, if this
                 * was the last slice
                 */

#pragma omp atomic
            ++m_completePixels;
            
            if (m_displayProgInterval > 0) if ((m_completePixels) % m_displayProgInterval == 0) 
               {
                 //#pragma omp atomic
                 //total_time += cputim()-timer;
                 printf("Pixels complete: (%d/%d), int.=%.3f, avg time per pixel: %.2fsec\n",
                        m_completePixels, m_scanXN*m_scanYN, m_intIntensity,
                        (total_time)/m_completePixels);
                 //timer=cputim();
               }
          } /* end of looping through STEM image pixels */
        /* save STEM images in img files */
        SaveImages();
        m_totalSliceCount += m_potential->GetNSlices();
      } /* end of loop through thickness (pCount) */
    // printf("Total CPU time = %f sec.\n", cputim()-timerTot ); 

    /*************************************************************/
    m_intIntensity = collectedIntensity/(m_scanXN*m_scanYN);
    DisplayProgress(1);
  } /* end of loop over m_avgCount */
}

void CExperimentSTEM::DisplayParams()
{
  float_tt dx;
  m_wave->GetResolution(dx, dx);
  printf("*\n"
         "* STEM parameters:\n");
  printf("* Maximum scattering angle:  %.0f mrad\n",
         0.5*2.0/3.0*m_wave->GetWavelength()/dx*1000);    
  m_detectors->PrintDetectors();
    
  printf("* Scan window:          (%g,%g) to (%g,%g)A, %d x %d = %d pixels\n",
         m_scanXStart,m_scanYStart,m_scanXStop,m_scanYStop,
         m_scanXN,m_scanYN,m_scanXN*m_scanYN);
}


/*****  saveSTEMImages *******/
// Saves all detector images (STEM images) that are defined in m_
//   When saving intermediate STEM images is enabled, this also saves
//   the intermediate STEM images for each detector.
void CExperimentSTEM::SaveImages()
{
  std::map<std::string, double> params;
  params["Runs Averaged"]=(double)m_avgCount+1;
  m_detectors->SaveDetectors(params);  
}

void CExperimentSTEM::CollectIntensity(unsigned absoluteSlice)
{
	m_detectors->CollectIntensity(m_wave, absoluteSlice);//m_totalSliceCount+islice*(1+mRepeat));
}
