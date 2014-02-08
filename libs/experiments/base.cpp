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
  // Subclasses should choose their wavefunction type appropriately - the config file does not
  //    provide enough information to determine this ATM.

  // Read potential parameters and initialize a pot object
  m_potential = CPotFactory::Get()->GetPotential(configReader);
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
  printf("* Running program STEM3 (version %.2f) in %s mode\n",VERSION, m_mode.c_str());
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

      std::map<unsigned, float_tt> displacements(m_crystal->GetU2());
      std::map<unsigned, float_tt>::iterator disp=displacements.begin(), end=displacements.end();

      printf("* <u>: %3d |",(*disp++).first);
	  while(disp!=end) printf(" %8d |",(*disp++).first);  

      printf(" intensity | time(sec) |    chi^2  |\n");
      // }
      /*
        printf("* %9g | %9g | %9g \n",muls.u2,muls.intIntensity,curTime);  
        }
        else {
      */
      printf("*");

      //ComputeAverageU2();

      disp = displacements.begin();
      while (disp!=end) printf(" %8f |",(*disp++).second);  
      printf(" %9f | %9f | %9f |\n",m_intIntensity,curTime,m_avgCount > 0 ? m_chisq[m_avgCount-1] : 0);
      printf("*");

      /*
        // TODO: averaging should be handled on this class, not on lower level crystal class.
      atom = atomTypes.begin();
      while (atom!=end) printf(" %8f |",(float)(m_crystal->GetU2avg((*atom++))));  
      */
      printf(" %9f | %9f \n",intensityAvg,timeAvg);
    }
    else {
      printf("\n**************** finished after %.1f sec ******************\n",curTime);
    }
  }  // end of printLevel check.

  time(&time0);
  //  timer = cputim();
}

/*
// TODO: this is very broken.  displaced here from Crystal because I want to handle averages outside of the 
//     lower level classes.
void CExperimentBase::ComputeAverageU2()
{
  
  (*z)->second /= u2Count[(*z)->first];
  if (runCount > 0) 
    m_u2avg[(*z)] = sqrt(((runCount-1)*(m_u2avg[(*z)]*m_u2avg[(*z)])+u2[(*z)])/runCount);
  else
    m_u2avg[(*z)] = sqrt(m_u2[(*z)]);
}
*/

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
int CExperimentBase::RunMuls(WavePtr wave) 
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

  wave->GetSizePixels(nx, ny);

  printFlag = (m_printLevel > 3);
  fftScale = 1.0/(nx*ny);

  wavlen = wave->GetWavelength();

  /*  calculate the total specimen thickness and echo */
  cztot=0.0;
  for( islice=0; islice<m_potential->GetNSlices(); islice++) {
    cztot += m_potential->GetSliceThickness(islice);
  }
  if (printFlag)
    printf("Specimen thickness: %g Angstroms\n", cztot);

  for( islice=0; islice < m_potential->GetNSlices(); islice++ ) 
    {
      absolute_slice = (m_totalSliceCount+islice);
          
      /***********************************************************************
       * Transmit is a simple multiplication of wave with trans in real space
       **********************************************************************/
      Transmit(wave, islice);   
      /***************************************************** 
       * remember: prop must be here to anti-alias
       * propagate is a simple multiplication of wave with prop
       * but it also takes care of the bandwidth limiting
       *******************************************************/
      wave->ToFourierSpace();
      Propagate(wave, islice);
      //propagate_slow(wave, m_nx, m_ny, muls);

      CollectIntensity(absolute_slice);

      // go back to real space:
      wave->ToRealSpace();
      // TODO: Is this normalization necessary?
      fft_normalize(wave);
      
      // write the intermediate TEM wave function:
      
      /********************************************************************
       * show progress:
       ********************************************************************/
      //m_wave->thickness = (absolute_slice+1)*m_sliceThickness;
      /*
      // TODO: if we want this, move it to STEM class

      if ((printFlag)) {
        sum=wave->GetIntegratedIntensity()*fftScale;
        
        sprintf(outStr,"position (%3d, %3d), slice %4d (%.2f), int. = %f", 
                m_wave->detPosX, m_wave->detPosY,
                m_totalSliceCount+islice,m_wave->thickness,sum );
        if (showEverySlice)
          printf("%s\n",outStr);
        else {
          printf("%s",outStr);
          for (i=0;i<(int)strlen(outStr);i++) printf("\b");
        }
      }
      */
	
      // Call any additional saving/post-processing that should occur on a per-slice basis
      PostSliceProcess(absolute_slice);
    } /* end for(islice...) */
      // collect intensity at the final slice
      //collectIntensity(muls, wave, m_totalSliceCount+m_slices*(1+mRepeat));
  if (printFlag) printf("\n***************************************\n");

  /*
  // TODO: modifying shared value from multiple threads?
  m_rmin  = m_wave->wave[0][0][0];
  //#pragma omp single
  m_rmax  = (*muls).rmin;
  //#pragma omp single
  m_aimin = m_wave->wave[0][0][1];
  //#pragma omp single
  m_aimax = (*muls).aimin;

  sum = 0.0;
  for( ix=0; ix<nx; ix++)  
    {
    for( iy=0; iy<ny; iy++) 
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
  */
  if ((m_saveLevel > 1) || (m_cellDiv > 1)) {
    wave->WriteWave();
  }
  return 0;
}  // end of runMulsSTEM

/******************************************************************
* propagate_slow() 
* Propagates a wave
*****************************************************************/
void CExperimentBase::Propagate(WavePtr wave, float_tt dz)
{
  int ixa, iya;
  float_tt wr, wi, tr, ti;
  float_tt scale,t; 
  float_tt dzs=0;

  float_tt dx, dy;
  unsigned nx, ny, px;
  
  wave->GetResolution(dx, dy);
  wave->GetSizePixels(nx, ny);

  complex_tt *w=wave->GetWavePointer();

  px=nx*ny;

  // TODO: this will not be thread safe, since m_propxr will be shared amongst threads
  if (dz != dzs) {
    dzs = dz;
    scale = dz*PI;

    for( ixa=0; ixa<nx; ixa++) {
      t = scale * (wave->GetKX2(ixa)*wave->GetWavelength());
      m_propxr[ixa] = (float_tt)  cos(t);
      m_propxi[ixa] = (float_tt) -sin(t);
    }
    for( iya=0; iya<ny; iya++) {
      
      t = scale * (wave->GetKY2(iya)*wave->GetWavelength());
      m_propyr[iya] = (float_tt)  cos(t);
      m_propyi[iya] = (float_tt) -sin(t);
    }
    
  } 
  /* end of: if dz != dzs */
  /*************************************************************/
  
  /*************************************************************
   * Propagation
   ************************************************************/
  for (unsigned i=0; i<px; i++)
    {
      ixa=i%nx;
      iya=i/nx;
      if( wave->GetKX2(ixa) < wave->GetK2Max() ) {
        if( (wave->GetKX2(ixa) + wave->GetKY2(iya)) < wave->GetK2Max() ) {
          wr = w[i][0];
          wi = w[i][1];
          tr = wr*m_propyr[iya] - wi*m_propyi[iya];
          ti = wr*m_propyi[iya] + wi*m_propyr[iya];
          w[i][0] = tr*m_propxr[ixa] - ti*m_propxi[ixa];
          w[i][1] = tr*m_propxi[ixa] + ti*m_propxr[ixa];
        } else
          w[i][0] = w[i][1] = 0.0F;
      } /* end for(iy..) */
      
      else w[i][0] = w[i][1] = 0.0F;
    } /* end for(ix..) */
} /* end propagate */

/*------------------------ transmit() ------------------------*/
/*
transmit the wavefunction thru one layer 
(simply multiply wave by transmission function)

waver,i[ix][iy]  = real and imaginary parts of wavefunction
transr,i[ix][iy] = real and imag parts of transmission functions

nx, ny = size of array

on entrance waver,i and transr,i are in real space

only waver,i will be changed by this routine
*/
void CExperimentBase::Transmit(WavePtr wave, unsigned sliceIdx) {
  double wr, wi, tr, ti;
  
  complex_tt *w;
  unsigned nx, ny;
  
  w = wave->GetWavePointer();
  wave->GetSizePixels(nx, ny);

  /*  trans += posx; */
  for(unsigned ix=0; ix<nx; ix++) for(unsigned iy=0; iy<ny; iy++) {
      unsigned offset=ix+nx*iy;
      complex_tt t = m_potential->GetSlicePixel(sliceIdx, ix+m_iPosX, iy+m_iPosY);

      wr = w[offset][0];
      wi = w[offset][1];
      tr = t[0];
      ti = t[1];
      w[offset][0] = wr*tr - wi*ti;
      w[offset][1] = wr*ti + wi*tr;
    } /* end for(iy.. ix .) */
} /* end transmit() */

void CExperimentBase::AddDPToAvgArray(const WavePtr &wave)
{
  unsigned px=wave->GetTotalPixels();
  // get the pointer to the first data element, and do 1D addressing (it's faster)
  float_tt chisq;

  const float_tt *dp = wave->GetDPPointer();

  for (unsigned i=0; i<px; i++)
    {
      float_tt t=m_avgArray[i]*m_avgCount+dp[i]/(m_avgCount+1);
      chisq+=(m_avgArray[i]-t)*(m_avgArray[i]-t);
      m_avgArray[i]=t;
    }
  #pragma omp atomic
  m_chisq[m_avgCount]+=chisq/px;
}

void CExperimentBase::_WriteAvgArray(std::string &fileName, std::string &comment, 
                                      std::map<std::string, double> &params,
                                      std::vector<unsigned> &position)
{
  //params["dx"]=1.0/(m_nx*m_dx);
  //params["dy"]=1.0/(m_ny*m_dy);
  params["Thickness"]=m_thickness;
  m_imageIO->WriteRealImage(m_avgArray, fileName, params, comment, position);
}

void CExperimentBase::ReadAvgArray()
{
  std::vector<unsigned> position;
  m_imageIO->ReadImage(m_avgArray, avgFilePrefix, position);
}

void CExperimentBase::ReadAvgArray(unsigned navg)
{
  std::vector<unsigned> position(1);
  position[0]=navg;
  m_imageIO->ReadImage(m_avgArray, avgFilePrefix, position);
}

void CExperimentBase::ReadAvgArray(unsigned positionx, unsigned positiony)
{
  std::vector<unsigned>position(2);
  position[0]=positionx;
  position[1]=positiony;
  m_imageIO->ReadImage(m_avgArray, avgFilePrefix, position);
}

void CExperimentBase::fft_normalize(WavePtr wave) 
{
  complex_tt *w;
  unsigned px = wave->GetTotalPixels();

  float_tt fftScale = 1.0/px;
  for (unsigned i=0; i<px; i++)
    {
      w[i][0] *= fftScale;
      w[i][1] *= fftScale;
    }
}










