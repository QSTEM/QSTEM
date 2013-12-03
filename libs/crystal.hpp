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

#ifndef CRYSTAL_H
#define CRYSTAL_H

#include "stemtypes_fftw3.hpp"

class CCrystal
{
public:
  CCrystal(ConfigReaderPtr &configReader);
  ~CCrystal();
  
  void ReadUnitCell(unsigned &natom, char *fileName, int handleVacancies);
  void TiltBoxed(int ncoord,&atoms,int handleVacancies);
  void PhononDisplacement(double *u,MULS *muls,int id,int icx,int icy,
                          int icz,int atomCount,double dw,int maxAtom,int ZnumIndex);
  void ReplicateUnitCell(int handleVacancies);
  
protected:
  std::vector<atom> m_atoms;
}

#endif










