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

CExperimentCBED::CExperimentCBED(const ConfigReaderPtr &configReader) : CExperimentBase(configReader)
{
}

void CExperimentCBED::Run()
{
  int ix,iy,i,pCount,result;
  FILE *avgFp,*fp,*fpPos=0;
  double timer,timerTot;
  double probeCenterX,probeCenterY,probeOffsetX,probeOffsetY;
  char buf[BUF_LEN];
  float_tt t=0;
  float_tt **avgPendelloesung = NULL;
  int oldMulsRepeat1 = 1;
  int oldMulsRepeat2 = 1;
  long iseed=0;
  std::map<std::string, double> params;

  std::vector<unsigned> position(1);         // Used to indicate the number of averages

  muls.chisq = std::vector<double>(muls.avgRuns);

  if (iseed == 0) iseed = -(long) time( NULL );

  if (muls.lbeams) {
    muls.pendelloesung = NULL;
    if (avgPendelloesung == NULL) {
      avgPendelloesung = float2D(muls.nbout,
                                 muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv,
                                 "pendelloesung");
    }    
  }
  probeCenterX = muls.scanXStart;
  probeCenterY = muls.scanYStart;

  timerTot = 0; /* cputim();*/
  displayProgress(-1);

  for (muls.avgCount = 0;muls.avgCount < muls.avgRuns;muls.avgCount++) {
    muls.totalSliceCount = 0;
    pCount = 0;
    /* make sure we start at the beginning of the file 
       so we won't miss any line that contains a sequence,
       because we will not do any EOF wrapping
    */
    resetParamFile();
    
    /* probe(&muls,xpos,ypos); */
    /* make incident probe wave function with probe exactly in the center */
    /* if the potential array is not big enough, the probe can 
     * then also be adjusted, so that it is off-center
     */

    probeOffsetX = muls.sourceRadius*gasdev(&iseed)*SQRT_2;
    probeOffsetY = muls.sourceRadius*gasdev(&iseed)*SQRT_2;
    muls.scanXStart = probeCenterX+probeOffsetX;
    muls.scanYStart = probeCenterY+probeOffsetY;
    wave->FormProbe();
    //probe(&muls, wave,muls.scanXStart-muls.potOffsetX,muls.scanYStart-muls.potOffsetY);
    if (muls.saveLevel > 2) {
      wave->WriteProbe();
    } 	
    // printf("Probe: (%g, %g)\n",muls.scanXStart,muls.scanYStart);
    /*****************************************************************
     * For debugging only!!!
     *
		muls->WriteWave("probe.img")
    *****************************************************************/


    if (muls.sourceRadius > 0) {
      if (muls.avgCount == 0) fpPos = fopen("probepos.dat","w");
      else fpPos = fopen("probepos.dat","a");
      if (fpPos == NULL) {
        printf("Was unable to open file probepos.dat for writing\n");
      }
      else {
        fprintf(fpPos,"%g %g\n",muls.scanXStart,muls.scanYStart);
        fclose(fpPos);
      }
    }

    if ((muls.showProbe) && (muls.avgCount == 0)) {
#ifndef WIN32
      //probePlot(&muls);
      sprintf(buf,"ee %s/probePlot_0.jpg &",muls.folder.c_str());
      system(buf);
#endif
    }
    //muls.nslic0 = 0;

    result = readparam("sequence: ",buf,0);
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
      for (i=0;i<(int)strlen(buf);i++) buf[i] = 0;
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
        if (muls.pendelloesung != NULL)
          free(muls.pendelloesung[0]);
        free(avgPendelloesung[0]);
        muls.pendelloesung = NULL;
        avgPendelloesung = float2D(muls.nbout,
                                   muls.slices*oldMulsRepeat1*oldMulsRepeat2*muls.cellDiv,
                                   "pendelloesung");
      }
      /*********************************************************/
      
      // printf("Stacking sequence: %s\n",buf);

      muls.saveFlag = 0;
      /****************************************
       * do the (small) loop
       *****************************************/
      for (pCount = 0;pCount<muls.mulsRepeat2*muls.cellDiv;pCount++) {
        
        pot->Refresh();
        
        timer = cputim();
        // what probe should runMulsSTEM use here?
        runMulsSTEM(&muls, wave, pot); 
        
        printf("Thickness: %gA, int.=%g, time: %gsec\n",
               wave->thickness,wave->intIntensity,cputim()-timer);

        /***************** Only if Save level > 2: ****************/
        if ((muls.avgCount == 0) && (muls.saveLevel > 2)) {
          wave->WriteWave();
        } 	
#ifdef VIB_IMAGE_TEST_CBED
        wave->WriteWave()
#endif 
          muls.totalSliceCount += muls.slices;
        
      } // end of for pCount = 0... 
      result = readparam("sequence: ",buf,0);
    }
    /*    printf("Total CPU time = %f sec.\n", cputim()-timerTot ); */
    
    // TODO: Why are we reading in a DP at this point?  Do we have one yet?  
    //     What happens if it isn't there?
    wave->ReadDiffPat();
    
    if (muls.avgCount == 0) {
      memcpy((void *)wave->avgArray[0],(void *)wave->diffpat[0],
             (size_t)(muls.nx*muls.ny*sizeof(float_tt)));
      /* move the averaged (raw data) file to the target directory as well */
      // TODO: make sure that DP average gets created properly
      //sprintf(avgName,"%s/diffAvg_%d.img",muls.folder.c_str(),muls.avgCount+1);
      //sprintf(systStr,"mv %s/diff.img %s",muls.folder.c_str(),avgName);
      //system(systStr);
      if (muls.lbeams) {
        for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
          for (ix=0;ix<muls.nbout;ix++) {
            avgPendelloesung[ix][iy] = muls.pendelloesung[ix][iy];
          }
        }
      }
    } // of if muls.avgCount == 0 ...
    else {
      muls.chisq[muls.avgCount-1] = 0.0;
      for (ix=0;ix<muls.nx;ix++) for (iy=0;iy<muls.ny;iy++) {
          t = ((float_tt)muls.avgCount*wave->avgArray[ix][iy]+
               wave->diffpat[ix][iy])/((float_tt)(muls.avgCount+1));
          muls.chisq[muls.avgCount-1] += (wave->avgArray[ix][iy]-t)*(wave->avgArray[ix][iy]-t);
          wave->avgArray[ix][iy] = t;

        }
      muls.chisq[muls.avgCount-1] = muls.chisq[muls.avgCount-1]/(double)(muls.nx*muls.ny);
      params["Tilt"] = muls.tomoTilt;
      params["1/Wavelength"] = 1.0/wavelength(muls.v0);
      wave->WriteDiffPat("Averaged Diffraction pattern, unit: 1/A", params);
                        
      muls.storeSeries = 1;
      if (muls.saveLevel == 0)	muls.storeSeries = 0;
      else if (muls.avgCount % muls.saveLevel != 0) muls.storeSeries = 0;

      if (muls.storeSeries) 
        wave->WriteAvgArray(muls.avgCount+1, "Averaged Diffraction pattern, unit: 1/A", params);


      /* write the data to a file */
      // TODO: vestigial code?  does anyone use this?
      if (muls.saveFlag >-1) {
        char systStr[255];
        sprintf(systStr,"%s/avgresults.dat",muls.folder.c_str());
        if ((avgFp = fopen(systStr,"w")) == NULL )
          printf("Sorry, could not open data file for averaging\n");
        else {
          for (ix =0;ix<muls.avgCount;ix++) {
            fprintf(avgFp,"%d %g\n",ix+1,muls.chisq[ix]);
          }
          fclose(avgFp);
        }
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
    } /* else ... if avgCount was greater than 0 */
    
    if (muls.lbeams) {
      /**************************************************************
       * The diffraction spot intensities of the selected 
       * diffraction spots are now stored in the 2 dimensional array
       * muls.pendelloesung[beam][slice].
       * We can write the array to a file and display it, just for 
       * demonstration purposes
       *************************************************************/
      char systStr[255];
      sprintf(systStr,"%s/pendelloesung.dat",muls.folder.c_str());
      if ((fp=fopen(systStr,"w")) !=NULL) {
        printf("Writing Pendelloesung data\n");
        for (iy=0;iy<muls.slices*muls.mulsRepeat1*muls.mulsRepeat2*muls.cellDiv;iy++) {
          /* write the thicknes in the first column of the file */
          fprintf(fp,"%g",iy*muls.cz/((float)(muls.slices*muls.cellDiv)));
          /* write the beam intensities in the following columns */
          for (ix=0;ix<muls.nbout;ix++) {
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
    displayProgress(1);
  } /* end of for muls.avgCount=0.. */
  //delete(wave);
}
