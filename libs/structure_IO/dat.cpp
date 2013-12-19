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

#include "dat.hpp"

/***********************************************************************
* The following function returns the number of atoms in the specified
* CFG file and updates the cell parameters in the muls struct
***********************************************************************/
int CDatReader::ReadCellParams(float_tt **Mm) {
  int ncoord;
  char buf[256];
  float_tt a,b,c,alpha=90.0,beta=90.0,gamma=90.0;

  // printf("Paramter File pointer (1): %d\n",(int)getFp());
  parFpPush(); /* push the current parameter file file pointer 
                * on the stack to make this function totaly 
                * transparent 
                */
  if (!parOpen(fileName)) {
    printf("Could not open CFG input file %s\n",fileName);
    parFpPull();  /* restore old parameter file pointer */
    return 0;
  }
  // printf("Paramter File pointer (2): %d\n",(int)getFp());

  resetParamFile();  
  setComment('#');  
  if (readparam("Number of atoms =",buf,1)) sscanf(buf,"%d",&ncoord);

  if (readparam("a =",buf,1)) sscanf(buf,"%lf",&a);
  if (readparam("b =",buf,1)) sscanf(buf,"%lf",&b);
  if (readparam("c =",buf,1)) sscanf(buf,"%lf",&c);

  if (readparam("alpha =",buf,1)) sscanf(buf,"%lf",&alpha);
  if (readparam("beta =",buf,1)) sscanf(buf,"%lf",&beta);
  if (readparam("gamma =",buf,1)) sscanf(buf,"%lf",&gamma);

  setComment('%');
  parClose();   
  parFpPull();  /* restore old parameter file pointer */


  muls->ax = a;
  muls->by = b;
  muls->c  = c;
  muls->cGamma = alpha;
  muls->cBeta  = beta;
  muls->cAlpha = gamma;
  // construct the unit cell metric from the lattice parameters and angles:
  makeCellVectMuls(muls, Mm[0], Mm[1], Mm[2]);   
  if (ncoord < 1) {
    printf("Number of atoms in CFG file not specified!\n");
    ncoord = 0;
  }
  return ncoord;
}

/*******************************************************************************
* This function reads the atomic position and element data for a single atom
* from a .dat file.  The atomic positions are given in reduced coordinates.
* Calling this function several times will read one atom after another from the 
* file.
* This function will return -1, if the end of file is reached prematurely, 0 otherwise.
*******************************************************************************/

int CDatReader::ReadNextAtom(atom *newAtom, int flag) {
  int printFlag = 0;
  static FILE *fp=NULL;
  static int noVelocityFlag = 1,entryCount = 3,element = 1;
  static char buf[NCMAX];
  static float_tt *atomData = NULL;
  char *str,elementStr[3];
  int j;

  if (flag < 0) {
    if (fp != NULL) {
      parClose();   
      parFpPull();  /* restore old parameter file pointer */      
      fp = NULL;
      setComment('%');
    }
    return -1;
  }

  if (fp == NULL) {
    parFpPush();  /* save old parameter file pointer */      
    if (!parOpen(fileName)) {
      printf("Could not open DAT input file %s\n",fileName);
      parFpPull();  /* restore old parameter file pointer */
      return -1;
    }
    resetParamFile();  
    // advance file pointer to last (known) header line
    readparam("gamma =",buf,1);
    fp = getFp();  /* get the file pointer from the parameter file routines */
    atomData = (float_tt *)malloc(entryCount*sizeof(float_tt));
  }

  element = 0;
  do {
    if (fgets(buf,NCMAX,fp) == NULL) {
      if (printFlag) printf("Found end of file!\n");
      return -1;
    }
    /* check, if this is a new element name */	
    // str = strnext(buf," \t");
    if (buf != NULL) {
      memcpy(elementStr,buf,2);
      element = getZNumber(elementStr);
    }
  } while (element == 0);
  if (printFlag) printf("Found Z=%d\n",element);
  if (element > 0) {
    str = buf+2;
    // skip leading spaces:
    while (strchr(" \t",*str) != NULL) str++; 
    for (j=0;j<entryCount;j++) {
      if (str==NULL) {
        printf("CDatReader::ReadNextAtom: Error: incomplete data line: >%s<\n",buf);
        return -1;
      }
      atomData[j] = atof(str); str=strnext(str," \t");
    }		
  }

  newAtom->Znum = element;
  newAtom->x    = atomData[0];
  newAtom->y    = atomData[1];
  newAtom->z    = atomData[2];
  newAtom->dw   = 0.45*28.0/(float_tt)(2.0*element);	
  newAtom->occ  = 1.0;
  newAtom->q    = 0.0;
  printf("Atom: %d (%g %g %g), occ=%g, q=%g\n",newAtom->Znum,newAtom->x,newAtom->y,newAtom->z,newAtom->occ,newAtom->q);	  

  return 0;
}
