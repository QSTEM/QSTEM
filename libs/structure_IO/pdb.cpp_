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

#include "pdb.hpp"

/********************************************************
* writePDB(atoms,natoms,fileName)
* This function will write the atomic coordinates in the 
* atoms array to the file <fileName> in pdb format, which
* is readable by AtomicEye
* The return value is 1 for success, or 0 for failure.
********************************************************/

int CStructurePDB::Write(atom *atoms,int natoms,char *fileName,MULS *muls) 
{
  FILE *fp;
  int j,i;
  static char *elTable = {
    "H HeLiBeB C N O F NeNaMgAlSiP S Cl"
    "ArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBr"
    "KrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTe"
    "I XeCsBaLaCePrNdPmSmEuGdTbDyHoErTm"
    "YbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn"
    "FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLr"};
  char elem[16];
  double ax,by,cz;
  
  if (natoms < 1) {
    printf("Atom array empty - no file written\n");
    return 1;
  }
  
  fp = fopen(fileName, "w" );
  if( fp == NULL ) {
    printf("Cannot open file %s\n",fileName);
    return 0;
  }
  
  ax = muls->ax;
  by=muls->by;
  cz=muls->c;


  fprintf(fp,"HEADER    libAtoms:Config_save_as_pdb; %d atoms\n",natoms);
  fprintf(fp,"CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1\n",
          muls->ax,muls->by,muls->c);
  elem[2] = '\0';
  for (j=0;j<natoms;j++) {
    elem[0] = elTable[2*atoms[j].Znum-2];
    elem[1] = elTable[2*atoms[j].Znum-1];
    if (elem[1] == ' ') elem[1] = '\0';
    fprintf(fp,"ATOM   %4d %s",j+1,elem);
    for (i=strlen(elem);i<13;i++)
      fprintf(fp," ");
    fprintf(fp,"1   ");
    fprintf(fp," %8.3f%8.3f%8.3f\n",atoms[j].x,atoms[j].y,atoms[j].z);
    /*
      fprintf(fp," %8.3f%8.3f%8.3f\n",atoms[j].x/muls->ax,
      atoms[j].y/muls->by,atoms[j].z/muls->c);
    */
  } 

  fclose(fp);
  return 1;
}
