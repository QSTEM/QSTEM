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

// This file defines a base class that should cover most multislice simulations.
// Things like displaying progress after a multislice run, handling multiple runs, etc are covered here.

CExperimentBase::CExperimentBase(const ConfigReaderPtr &reader) : IExperiment()
{
  // Read potential parameters and initialize a pot object
  m_wave = WavePtr(new WAVEFUNC(configReader));
  m_pot = GetPotential(configReader);
  readFile(configReader);
  displayParams(initialWave, potential);
}

void CExperimentBase::displayParams() {
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
  printf("* Running program STEM3 (version %.2f) in %s mode\n",VERSION,
         (muls.mode == STEM) ? "STEM" : (muls.mode==TEM) ? "TEM" : 
         (muls.mode == CBED) ? "CBED" : (muls.mode==TOMO)? "TOMO" : 
         "???"); 
  printf("* Date: %s, Time: %s\n",Date,Time);
	
  printf("* Input file:           %s\n",muls.atomPosFile);
  
  // create the data folder ... 
  printf("* Data folder:          ./%s/ ",muls.folder.c_str()); 
	
  printf("* Super cell divisions: %d (in z direction) %s\n",muls.cellDiv,muls.equalDivs ? "equal" : "non-equal");
  printf("* Slices per division:  %d (%gA thick slices [%scentered])\n",
         muls.slices,muls.sliceThickness,(muls.centerSlices) ? "" : "not ");
  printf("* Output every:         %d slices\n",muls.outputInterval);
  

  /* 
     if (muls.ismoth) printf("Type 1 (=smooth aperture), ");
     if (muls.gaussFlag) printf("will apply gaussian smoothing"); 
     printf("\n");
  */

  /***************************************************/
  /*  printf("Optimizing fftw plans according to probe array (%d x %dpixels = %g x %gA) ...\n",
      muls.nx,muls.ny,muls.nx*muls.resolutionX,muls.ny*muls.resolutionY);
  */
  
  printf("* TDS:                  %d runs)\n",muls.avgRuns);

  printf("*\n*****************************************************\n");
}

void CExperimentBase::DisplayProgress(int flag, MULS &muls, WavePtr &wave, StructurePtr &crystal)
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
  if (muls.m_printLevel > 0) {
    if (crystal->GetTDS()) {
      timeAvg = ((muls.avgCount)*timeAvg+curTime)/(muls.avgCount+1);
      intensityAvg = ((muls.avgCount)*intensityAvg+muls.intIntensity)/(muls.avgCount+1);
      printf("\n********************** run %3d ************************\n",muls.avgCount+1);
      // if (muls.avgCount < 1) {

	  std::vector<unsigned> atomTypes(crystal->GetAtomTypes());
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
      while (atom!=end) printf(" %8f |",(float)(crystal->GetU2((*atom++))));  
      printf(" %9f | %9f | %9f |\n",muls.intIntensity,curTime,muls.avgCount > 0 ? muls.chisq[muls.avgCount-1] : 0);
      printf("*");

      atom = atomTypes.begin();
      while (atom!=end) printf(" %8f |",(float)(crystal->GetU2avg((*atom++))));  
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

  if ((slice < muls->slices*muls->cellDiv-1) && ((slice+1) % muls->outputInterval != 0)) return;
  
  t = (int)((slice)/muls->outputInterval);
	
  // produce the following filename:
  // wave_avgCount_thicknessIndex.img or
  // wave_thicknessIndex.img if tds is turned off
  if (muls->tds) wave->WriteWave(muls->avgCount, t, "Wave Function", params);
  else wave->WriteWave(t, "Wave Function", params);
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
int CExperimentBase::runMuls() {
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

  wave->GetSizePixels(nx, ny);

  printFlag = (muls->printLevel > 3);
  fftScale = 1.0/(nx*ny);

  wavlen = wave->GetWavelength();

  /*  calculate the total specimen thickness and echo */
  cztot=0.0;
  for( islice=0; islice<(*muls).slices; islice++) {
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
          wave->Transmit(pot, islice);   
          /***************************************************** 
           * remember: prop must be here to anti-alias
           * propagate is a simple multiplication of wave with prop
           * but it also takes care of the bandwidth limiting
           *******************************************************/
#if FLOAT_PRECISION == 1
          fftwf_execute(wave->m_fftPlanWaveForw);
#else
          fftw_execute(wave->m_fftPlanWaveForw);
#endif
          wave->Propagate();
          //propagate_slow(wave, muls->nx, muls->ny, muls);

          muls->detectors->CollectIntensity(wave, muls->totalSliceCount+islice*(1+mRepeat));
          //collectIntensity(muls, wave, muls->totalSliceCount+islice*(1+mRepeat));
          
          if (muls->mode != STEM) {
            /* write pendelloesung plots, if this is not STEM */
            writeBeams(muls,wave,islice, absolute_slice);
          }

          // go back to real space:
#if FLOAT_PRECISION == 1
          fftwf_execute(wave->fftPlanWaveInv);
#else
          fftw_execute(wave->fftPlanWaveInv);
#endif
          // old code: fftwnd_one((*muls).fftPlanInv,(complex_tt *)wave[0][0], NULL);
          fft_normalize((void **)wave->wave,nx,ny);
          
          // write the intermediate TEM wave function:
          
          /********************************************************************
           * show progress:
           ********************************************************************/
          wave->thickness = (absolute_slice+1)*muls->sliceThickness;
          if ((printFlag)) {
            sum = 0.0;
            for( ix=0; ix<nx; ix++)  for( iy=0; iy<ny; iy++) {
                sum +=  wave->GetIntensity(ix,iy);
              }
            sum *= fftScale;
            
            sprintf(outStr,"position (%3d, %3d), slice %4d (%.2f), int. = %f", 
                    wave->detPosX, wave->detPosY,
                    muls->totalSliceCount+islice,wave->thickness,sum );
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
              interimWave(muls,wave,absolute_slice*(1+mRepeat)); 
              muls->detectors->CollectIntensity(wave, absolute_slice*(1+mRepeat));
              
              //collectIntensity(muls,wave,absolute_slice*(1+mRepeat));
            }
        } /* end for(islice...) */
      // collect intensity at the final slice
      //collectIntensity(muls, wave, muls->totalSliceCount+muls->slices*(1+mRepeat));
    } /* end of mRepeat = 0 ... */
  if (printFlag) printf("\n***************************************\n");

  // TODO: modifying shared value from multiple threads?
  muls->rmin  = wave->wave[0][0][0];
  //#pragma omp single
  muls->rmax  = (*muls).rmin;
  //#pragma omp single
  muls->aimin = wave->wave[0][0][1];
  //#pragma omp single
  muls->aimax = (*muls).aimin;

  sum = 0.0;
  for( ix=0; ix<muls->nx; ix++)  
    {
    for( iy=0; iy<muls->ny; iy++) 
      {
        x =  wave->wave[ix][iy][0];
        y =  wave->wave[ix][iy][1];
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
  wave->intIntensity = sum*fftScale;

  if (printFlag) {
    printf( "pix range %g to %g real,\n"
            "          %g to %g imag\n",  
            (*muls).rmin,(*muls).rmax,(*muls).aimin,(*muls).aimax);
    
  }
  if (muls->saveFlag) {
    if ((muls->saveLevel > 1) || (muls->cellDiv > 1)) {
      wave->WriteWave();
      if (printFlag)
        printf("Created complex image file %s\n",(*wave).fileout.c_str());
    }
  }
  return 0;
}  // end of runMulsSTEM
