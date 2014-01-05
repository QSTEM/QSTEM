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
