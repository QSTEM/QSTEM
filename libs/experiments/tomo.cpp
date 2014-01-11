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

#include "tomo.hpp"

/************************************************************************
* doTOMO performs a Diffraction Tomography simulation
*
* This routine creates a script file which can then be run by a separate 
* command.
* For now this routine will only allow tomography about the y-axis.
* To do anything else one can start with a previously rotated super-cell.
*
* Important parameters: tomoStart, tomoStep, tomoCount, zoomFactor
***********************************************************************/
void CExperimentTomo::run()
{
  double boxXmin=0,boxXmax=0,boxYmin=0,boxYmax=0,boxZmin=0,boxZmax=0;
  double mAx,mBy,mCz;
  int ix,iy,iz,iTheta,i;
  double u[3],**Mm = NULL;
  double theta = 0;
  atom *atoms = NULL;
  char cfgFile[64],stemFile[128],scriptFile[64],diffAnimFile[64];
  FILE *fpScript,*fpDiffAnim;

  float_tt ax = pot->GetCellAX();
  float_tt by = pot->GetCellBY();
  float_tt cz = pot->GetCellCZ();

  Mm = muls.Mm;
  atoms = (atom *)malloc(muls.natom*sizeof(atom));

  boxXmin = boxXmax = ax/2.0;
  boxYmin = boxYmax = by/2.0;
  boxZmin = boxZmax = cz/2.0;

  // For all tomography tilt angles:
  for (iTheta = 0;iTheta < muls.tomoCount;iTheta++) {
    theta = muls.tomoStart + iTheta*muls.tomoStep; // theta in mrad
    // Try different corners of the box, and see, how far they poke out.
    for (ix=-1;ix<=1;ix++) for (iy=-1;iy<=1;iy++) for (iz=-1;iz<=1;iz++) {
          // Make center of unit cell rotation center
          u[0]=ix*ax/2; u[1]=iy*by/2.0; u[2]=iz*cz/2.0;
          
          // rotate about y-axis
          rotateVect(u,u,0,theta*1e-3,0);
          
          // shift origin back to old (0,0,0):
          u[0]+=ax/2; u[1]+=by/2.0; u[2]+=cz/2.0;
          
          boxXmin = boxXmin>u[0] ? u[0] : boxXmin; boxXmax = boxXmax<u[0] ? u[0] : boxXmax; 
          boxYmin = boxYmin>u[1] ? u[1] : boxYmin; boxYmax = boxYmax<u[1] ? u[1] : boxYmax; 
          boxZmin = boxZmin>u[2] ? u[2] : boxZmin; boxZmax = boxZmax<u[2] ? u[2] : boxZmax; 
          
        }
  } // for iTheta ... 

  // find max. box size:
  boxXmax -= boxXmin;
  boxYmax -= boxYmin;
  boxZmax -= boxZmin;
  printf("Minimum box size for tomography tilt series: %g x %g x %gA, zoom Factor: %g\n",
         boxXmax,boxYmax,boxZmax,muls.zoomFactor);
  boxXmax /= muls.zoomFactor;
  boxYmax = boxXmax*by/ax;

  // boxMin will now be boxCenter:
  boxXmin = 0.5*boxXmax;
  boxYmin = 0.5*boxYmax;
  boxZmin = 0.5*boxZmax;

  // We have to save the original unit cell dimensions
  mAx = ax; mBy = by; mCz = cz;
  ax=boxXmax; by=boxYmax; cz=boxZmax;

  printf("Will use box sizes: %g x %g x %gA (kept original aspect ratio). \n"
         "Writing structure files now, please wait ...\n",
         boxXmax,boxYmax,boxZmax);

  // open the script file and write the stem instructions in there
  sprintf(scriptFile,"%s/run_tomo",muls.folder.c_str());
  if ((fpScript=fopen(scriptFile,"w")) == NULL) {
    printf("doTOMO: unable to open scriptFile %s for writing\n",scriptFile);
    exit(0);
  }
  fprintf(fpScript,"#!/bin/bash\n\n");

  // open the diffraction animation file 
  sprintf(diffAnimFile,"%s/diff_anim",muls.folder.c_str());
  if ((fpDiffAnim=fopen(diffAnimFile,"w")) == NULL) {
    printf("doTOMO: unable to open diffraction animation %s for writing\n",diffAnimFile);
    exit(0);
  }
  fprintf(fpDiffAnim,"#!/bin/bash\n\n");


  // We will now rotate all the unit cell and write the answer to
  // separate structure output files.
  for (iTheta = 0;iTheta < muls.tomoCount;iTheta++) {
    theta = muls.tomoStart + iTheta*muls.tomoStep; // theta in mrad
    muls.tomoTilt = theta;

    // rotate the structure and write result to local atom array
    for(i=0;i<(muls.natom);i++) {	
      u[0] = pot->GetAtom(i).x - mAx/2.0; 
      u[1] = pot->GetAtom(i).y - mBy/2.0; 
      u[2] = pot->GetAtom(i).z - mCz/2.0; 
      rotateVect(u,u,0,theta*1e-3,0);
      atoms[i].x = u[0]+boxXmin;
      atoms[i].y = u[1]+boxYmin; 
      atoms[i].z = u[2]+boxZmin; 
      atoms[i].Znum = pot->GetAtom(i).Znum;
      atoms[i].occ = pot->GetAtom(i).occ;
      atoms[i].dw = pot->GetAtom(i).dw;
    }

    sprintf(cfgFile,"%s/tomo_%dmrad.cfg",muls.folder.c_str(),(int)theta);
    sprintf(stemFile,"%s/tomo_%dmrad.dat",muls.folder.c_str(),(int)theta);
    printf("Writing file %s | ",cfgFile);
    writeCFG(atoms,cfgFile,&muls);
    sprintf(cfgFile,"tomo_%dmrad.cfg",(int)theta);
    writeSTEMinput(stemFile,cfgFile,&muls);

    // add to script files:
    fprintf(fpScript,"stem tomo_%dmrad.dat\n",(int)theta);
  }
  // TODO: this should probably set cell params on potential's crystal object
  ax = mAx; by = mBy; cz = mCz;
  sprintf(stemFile,"copy fparams.dat %s/",muls.folder.c_str());
  system(stemFile);

  // close script files again
  fprintf(fpDiffAnim,"convert -delay 20 diff*.jpg diff.gif\n");  
  fclose(fpScript);
  fclose(fpDiffAnim);
  sprintf(stemFile,"chmod +x %s",scriptFile);
  system(stemFile);
  sprintf(stemFile,"chmod +x %s",diffAnimFile);
  system(stemFile);
  
  exit(0);
}


void CExperimentTomo::DisplayParams()
{
    printf("*\n"
           "* TOMO parameters:\n");
    printf("* Starting angle:       %g mrad (%g deg)\n",
           m_tomoStart,m_tomoStart*0.18/pi);
    printf("* Angular increase:     %g mrad (%g deg)\n",
           m_tomoStep,m_tomoStep*0.180/pi);
    printf("* Number of dp's:       %d\n",m_tomoCount);
    printf("* Zoom factor:          %g\n",m_zoomFactor);
}
