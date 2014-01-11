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

#include "base.hpp"

// This file defines a base class that should cover most multislice simulations.
// Things like displaying progress after a multislice run, handling multiple runs, etc are covered here.

CExperimentBase::CExperimentBase(const ConfigReaderPtr &configReader) : IExperiment()
	, m_mode("Undefined")
{
  // Read potential parameters and initialize a pot object
  m_wave = WavePtr(new WAVEFUNC(configReader));
  m_potential = GetPotential(configReader);
  readFile(configReader);
  DisplayParams();
}

void CExperimentBase::DisplayParams() {
  FILE *fpDir;
  char systStr[64];
  double k2max,temp;
  int i,j;
  static char Date[16],Time[16];
  time_t caltime;
  struct tm *mytime;
  const double pi=3.1415926535897;

  /*
  if (wave->printLevel < 1) {
    if ((fpDir = fopen(muls.folder.c_str(),"r"))) {
      fclose(fpDir);
      // printf(" (already exists)\n");
    }
    else {
      sprintf(systStr,"mkdir %s",muls.folder.c_str());
      system(systStr);
      // printf(" (created)\n");
    }	  
    return;
  }
  */
  caltime = time( NULL );
  mytime = localtime( &caltime );
  strftime( Date, 12, "%Y:%m:%d", mytime );
  strftime( Time, 9, "%H:%M:%S", mytime );
  
  printf("\n*****************************************************\n");
  printf("* Running program STEM3 (version %.2f) in %s mode\n",VERSION, m_mode);
  printf("* Date: %s, Time: %s\n",Date,Time);
  
  // create the data folder ... 
  printf("* Output file/folder:          ./%s/ ",m_outputLocation.c_str()); 
	
  printf("* Super cell divisions: %d (in z direction) %s\n",m_cellDiv, m_equalDivs ? "equal" : "non-equal");
  printf("* Output every:         %d slices\n",m_outputInterval);
  

  /* 
     if (muls.ismoth) printf("Type 1 (=smooth aperture), ");
     if (muls.gaussFlag) printf("will apply gaussian smoothing"); 
     printf("\n");
  */

  /***************************************************/
  /*  printf("Optimizing fftw plans according to probe array (%d x %dpixels = %g x %gA) ...\n",
      muls.nx,muls.ny,muls.nx*muls.resolutionX,muls.ny*muls.resolutionY);
  */
  
  printf("* TDS:                  %d runs)\n",m_avgRuns);

  printf("*\n*****************************************************\n");
}

void CExperimentBase::DisplayProgress(int flag)
{
  // static double timer;
  static double timeAvg = 0;
  static double intensityAvg = 0;
  static time_t time0,time1;
  double curTime;
  int jz;

  if (flag < 0) {
    time(&time0);
    // timer = cputim();
    return;
  }
  time(&time1);  
  curTime = difftime(time1,time0);
  /*   curTime = cputim()-timer;
       if (curTime < 0) {
       printf("timer: %g, curr. time: %g, diff: %g\n",timer,cputim(),curTime);
       }
  */
  if (m_printLevel > 0) {
    if (m_crystal->GetTDS()) {
      timeAvg = ((m_avgCount)*timeAvg+curTime)/(m_avgCount+1);
      intensityAvg = ((m_avgCount)*intensityAvg+m_intIntensity)/(m_avgCount+1);
      printf("\n********************** run %3d ************************\n",m_avgCount+1);
      // if (muls.avgCount < 1) {

	  std::vector<unsigned> atomTypes(m_crystal->GetAtomTypes());
	  std::vector<unsigned>::iterator atom=atomTypes.begin(), end=atomTypes.end();

      printf("* <u>: %3d |",(*atom++));
	  while(atom!=end) printf(" %8d |",(*atom++));  

      printf(" intensity | time(sec) |    chi^2  |\n");
      // }
      /*
        printf("* %9g | %9g | %9g \n",muls.u2,muls.intIntensity,curTime);  
        }
        else {
      */
      printf("*");

      atom = atomTypes.begin();
      while (atom!=end) printf(" %8f |",(float)(m_crystal->GetU2((*atom++))));  
      printf(" %9f | %9f | %9f |\n",m_intIntensity,curTime,m_avgCount > 0 ? m_chisq[m_avgCount-1] : 0);
      printf("*");

      atom = atomTypes.begin();
      while (atom!=end) printf(" %8f |",(float)(m_crystal->GetU2avg((*atom++))));  
      printf(" %9f | %9f \n",intensityAvg,timeAvg);
    }
    else {
      printf("\n**************** finished after %.1f sec ******************\n",curTime);
    }
  }  // end of printLevel check.

  time(&time0);
  //  timer = cputim();
}


////////////////////////////////////////////////////////////////
// save the current wave function at this intermediate thickness:
void CExperimentBase::InterimWave(int slice) {
  int t;
  char fileName[256]; 
  std::map<std::string, double> params;

  if ((slice < m_potential->GetNSlices()*m_cellDiv-1) && ((slice+1) % m_outputInterval != 0)) return;
  
  t = (int)((slice)/m_outputInterval);
	
  // produce the following filename:
  // wave_avgCount_thicknessIndex.img or
  // wave_thicknessIndex.img if tds is turned off
  if (m_tds) m_wave->WriteWave(m_avgCount, t, "Wave Function", params);
  else m_wave->WriteWave(t, "Wave Function", params);
}

/******************************************************************
* runMulsSTEM() - do the multislice propagation in STEM/CBED mode
* 
*    Each probe position is running this function.  Each CPU is thus
*      running a separate instance of the function.  It is nested in
*      the main OpenMP parallel region - specifying critical, single, and
*      barrier OpenMP pragmas should be OK.
*
* waver, wavei are expected to contain incident wave function 
* they will be updated at return
*****************************************************************/
int CExperimentBase::RunMuls() 
{
  int printFlag = 0; 
  int showEverySlice=1;
  int islice,i,ix,iy,mRepeat;
  float_tt cztot=0.0;
  float_tt wavlen,sum=0.0; //,zsum=0.0
  // static int *layer=NULL;
  float_tt x,y;
  int absolute_slice;

  char outStr[64];
  double fftScale;

  unsigned nx, ny;

  m_wave->GetSizePixels(nx, ny);

  printFlag = (m_printLevel > 3);
  fftScale = 1.0/(nx*ny);

  wavlen = m_wave->GetWavelength();

  /*  calculate the total specimen thickness and echo */
  cztot=0.0;
  for( islice=0; islice<m_potential->GetNSlices(); islice++) {
    cztot += (*muls).cz[islice];
  }
  if (printFlag)
    printf("Specimen thickness: %g Angstroms\n", cztot);

  for (mRepeat = 0; mRepeat < muls->mulsRepeat1; mRepeat++) 
    {
      for( islice=0; islice < muls->slices; islice++ ) 
        {
          absolute_slice = (muls->totalSliceCount+islice);
          
          /***********************************************************************
           * Transmit is a simple multiplication of wave with trans in real space
           **********************************************************************/
          m_wave->Transmit(m_potential, islice);   
          /***************************************************** 
           * remember: prop must be here to anti-alias
           * propagate is a simple multiplication of wave with prop
           * but it also takes care of the bandwidth limiting
           *******************************************************/
#if FLOAT_PRECISION == 1
          fftwf_execute(m_wave->m_fftPlanWaveForw);
#else
          fftw_execute(m_wave->m_fftPlanWaveForw);
#endif
          m_wave->Propagate();
          //propagate_slow(wave, muls->nx, muls->ny, muls);

		  CollectIntensity(absolute_slice);
          
          if (muls->mode != STEM) {
            /* write pendelloesung plots, if this is not STEM */
            writeBeams(muls,m_wave,islice, absolute_slice);
          }

          // go back to real space:
#if FLOAT_PRECISION == 1
          fftwf_execute(m_wave->fftPlanWaveInv);
#else
          fftw_execute(wave->fftPlanWaveInv);
#endif
          // old code: fftwnd_one((*muls).fftPlanInv,(complex_tt *)wave[0][0], NULL);
          fft_normalize((void **)m_wave->wave,nx,ny);
          
          // write the intermediate TEM wave function:
          
          /********************************************************************
           * show progress:
           ********************************************************************/
          m_wave->thickness = (absolute_slice+1)*muls->sliceThickness;
          if ((printFlag)) {
            sum = 0.0;
            for( ix=0; ix<nx; ix++)  for( iy=0; iy<ny; iy++) {
                sum +=  m_wave->GetPixelIntensity(ix,iy);
              }
            sum *= fftScale;
            
            sprintf(outStr,"position (%3d, %3d), slice %4d (%.2f), int. = %f", 
                    m_wave->detPosX, m_wave->detPosY,
                    muls->totalSliceCount+islice,m_wave->thickness,sum );
            if (showEverySlice)
              printf("%s\n",outStr);
            else {
              printf("%s",outStr);
              for (i=0;i<(int)strlen(outStr);i++) printf("\b");
            }
          }
			
          if ((muls->mode == TEM) || ((muls->mode == CBED)&&(muls->saveLevel > 1))) 
            {
              // TODO (MCS 2013/04): this restructure probably broke this file saving - 
              //   need to rewrite a function to save things for TEM/CBED?
              // This used to call interimWave(muls,wave,muls->totalSliceCount+islice*(1+mRepeat));
              InterimWave(absolute_slice*(1+mRepeat)); 
              muls->detectors->CollectIntensity(m_wave, absolute_slice*(1+mRepeat));
              
              //collectIntensity(muls,wave,absolute_slice*(1+mRepeat));
            }
        } /* end for(islice...) */
      // collect intensity at the final slice
      //collectIntensity(muls, wave, muls->totalSliceCount+muls->slices*(1+mRepeat));
    } /* end of mRepeat = 0 ... */
  if (printFlag) printf("\n***************************************\n");

  // TODO: modifying shared value from multiple threads?
  muls->rmin  = m_wave->wave[0][0][0];
  //#pragma omp single
  muls->rmax  = (*muls).rmin;
  //#pragma omp single
  muls->aimin = m_wave->wave[0][0][1];
  //#pragma omp single
  muls->aimax = (*muls).aimin;

  sum = 0.0;
  for( ix=0; ix<muls->nx; ix++)  
    {
    for( iy=0; iy<muls->ny; iy++) 
      {
        x =  m_wave->wave[ix][iy][0];
        y =  m_wave->wave[ix][iy][1];
        if( x < (*muls).rmin ) (*muls).rmin = x;
        if( x > (*muls).rmax ) (*muls).rmax = x;
        if( y < (*muls).aimin ) (*muls).aimin = y;
        if( y > (*muls).aimax ) (*muls).aimax = y;
        sum += x*x+y*y;
      }
    }
  // TODO: modifying shared value from multiple threads?
  //  Is this sum supposed to be across multiple pixels?
  //#pragma omp critical
  m_wave->intIntensity = sum*fftScale;

  if (printFlag) {
    printf( "pix range %g to %g real,\n"
            "          %g to %g imag\n",  
            (*muls).rmin,(*muls).rmax,(*muls).aimin,(*muls).aimax);
    
  }
  if (muls->saveFlag) {
    if ((muls->saveLevel > 1) || (muls->cellDiv > 1)) {
      m_wave->WriteWave();
      if (printFlag)
        printf("Created complex image file %s\n",(*wave).fileout.c_str());
    }
  }
  return 0;
}  // end of runMulsSTEM
