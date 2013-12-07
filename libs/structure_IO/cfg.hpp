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

#ifndef CFG_STRUCTURE_H
#define CFG_STRUCTURE_H

#include "structureInterface.hpp"
#include <boost/filesystem.hpp>

class CStructureCfg : public IStructureIO
{
public:
  CStructureCfg(boost::filesystem::path &structure_file);
  ~CStructureCfg();
  int Write(unsigned run_number);
  int WriteFractCubic(double *pos,int *Znum,double *dw,int natoms,char *fileName,
                  double a,double b,double c);
  int ReadCellParams(float_tt **Mm);
  int ReadNextAtom(atom *newAtom, int flag);
};


#endif
