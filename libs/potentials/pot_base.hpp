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

#include "../stemtypes_fftw3.hpp"
#include "../imagelib_fftw3.hpp"
#include "../config_readers.hpp"
#include "../memory_fftw3.hpp"
#include <map>

#ifndef POTENTIAL_BASE_H
#define POTENTIAL_BASE_H

#define OVERSAMPLING 3
#define OVERSAMPLINGZ 3*OVERSAMPLING

class CPotential
{
public:
  CPotential(unsigned nx, unsigned ny, unsigned nz, float_tt dx, float_tt dy, float_tt dz, float_tt atomRadius, float_tt v0);
  CPotential(ConfigReaderPtr &configReader);
  ~CPotential();

  virtual void atomBoxLookUp(complex_tt &val, int Znum, float_tt x, float_tt y, float_tt z, float_tt B);
  virtual void make3DSlices(int nlayer,char *fileName,atom *center);
  virtual void initSTEMSlices();
  // encapsulates make slices and initSTEMslices - used to refresh the potential with a new structure (after a random
  //    shake)
  virtual void Refresh();
protected:
  void Initialize();

  ImageIOPtr m_imageIO;
  complex_tt ***m_trans;
  unsigned m_nx, m_ny;
  // resolutions
  float_tt m_dx, m_dy, m_dz;
  // oversampled resolutions
  float_tt m_ddx, m_ddy, m_ddz;
  //
  int m_boxNx, m_boxNy, m_boxNz;
  // Atom radius
  float_tt m_radius, m_radius2;
  // voltage
  float_tt m_v0;
  std::map<unsigned, atomBoxPtr> m_atomBoxes;

  int m_printLevel;
  bool m_savePotential, m_saveProjectedPotential, m_plotPotential;
};

typedef boost::shared_ptr<CPotential> PotPtr;

#endif
