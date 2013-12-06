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
  
  void ReadUnitCell(char *fileName, int handleVacancies);
  void TiltBoxed(int ncoord,int handleVacancies);
  void PhononDisplacement(double *u,int id,int icx,int icy,
                          int icz,double dw,int maxAtom,int ZnumIndex);
  void ReplicateUnitCell(int handleVacancies);

  float_tt GetCZ(){return m_cz;}
  
protected:
  unsigned m_natoms;
  std::vector<atom> m_atoms;
  float_tt **m_Mm;
  float_tt m_ax, m_by, m_cz;
  float_tt m_cAlpha, m_cBeta, m_cGamma;
};

typedef boost::shared_ptr<CCrystal> StructurePtr;

#endif










