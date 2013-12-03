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

#include "cssr.hpp"


/***********************************************************************
* The following function returns the number of atoms in the specified
* CFG file and updates the cell parameters in the muls struct
***********************************************************************/
int CStructureCSSR::ReadCellParams(MULS *muls, double **Mm, char *fileName) 
{
  FILE *fp;
  char s1[32],s2[32],s3[32],buf[NCMAX];
  int spaceGrp = 1,ncoord;

  fp = fopen(fileName, "r" );
  if( fp == NULL ) {
    printf("Cannot open file %s\n",fileName);
    return 0;
  }
  /* read the cell parameters from first and secons line */
  ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
  sscanf( buf, " %s %s %s",s1,s2,s3);
  muls->ax = atof(s1); muls->by = atof(s2); muls->c = atof(s3);
  
  ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
  sscanf( buf, " %s %s %s",s1,s2,s3);
  muls->cAlpha = atof(s1); muls->cBeta  = atof(s2); muls->cGamma = atof(s3);
  makeCellVectMuls(muls, Mm[0], Mm[1], Mm[2]);   

  /* check the space group: */
  spaceGrp = atoi(strstr(buf,"SPGR =")+strlen("SPGR ="));
  if (spaceGrp != 1) {
    printf("cannot interpret space group %d\n",spaceGrp);
    exit(0);
  }
  ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
  ncoord = atoi(buf);
  fclose(fp);
  return ncoord;
}


/*******************************************************************************
* This function reads the atomic position and element data for a single atom
* from a .cssr file.  The atomic positions are given in reduced coordinates.
* Calling this function several times will read one atom after another from the 
* file.
* A ReadLine error will occur, if the end of file is reached prematurely.
*******************************************************************************/
int CStructureCSSR::ReadNextAtom(atom *newAtom,int flag, char *fileName) {
  static FILE *fp=NULL;
  static char buf[NCMAX];
  int count;
  static char element[32],s1[32],s2[32],s3[32];
  double dw;

  if (flag < 0) {
    if (fp != NULL) fclose(fp);
    fp = NULL;
    return 0;
  }

  if (fp == NULL) {
    fp = fopen(fileName, "r" );
    if( fp == NULL ) {
      printf("Cannot open file %s\n",fileName);
      exit( 0 );
    }
    // skip the header:
    ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
    ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
    ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
    ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
  }

  ReadLine( fp, buf, NCMAX, "in ReadXYZcoord()" );
  /* for Si */
  /*    dw = 0.444; */
  sscanf(buf,"%d %s %s %s %s %d %d %d %d %d %d %d %d %lf",
         &count,element,s1,s2,s3,&count,&count,
         &count,&count,&count,&count,&count,&count,&dw);
  
  newAtom->x = atof(s1);
  newAtom->y = atof(s2);
  newAtom->z = atof(s3);
  newAtom->occ = 1.0;
  newAtom->Znum = getZNumber(element); 
  newAtom->dw = dw;

  /*
    switch (newAtom->Znum) {
    case 14:  newAtom->dw = 0.45f;
    break;
    case 29:  newAtom->dw = 0.21f;
    break;
    default: newAtom->dw = 0.21f;
    }
  */
  return 0;
}
