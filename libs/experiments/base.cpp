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



CExperimentBase::CExperimentBase(ConfigReaderPtr &reader)
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
