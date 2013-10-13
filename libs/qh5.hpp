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

#ifndef QH5_H
#define QH5_H

#include "stemtypes_fftw3.hpp"

#include "H5Cpp.h"

#include <string>

// This class provides utilities for opening hdf5 files and for creating default data sets.

class QH5
{
public:
  void OpenFile(const std::string &filename);

  void GetConfigGroup();

  void CreatePotentialDataSet(unsigned size_x, unsigned size_y, unsigned size_z);
  void CreateWaveDataSet(unsigned size_x, unsigned size_y, std::vector<unsigned> positions);
  void CreateDPDataSet(unsigned size_x, unsigned size_y, std::vector<unsigned> positions);
  void CreateDetectorDataSet(std::string name, unsigned size_x, unsigned size_y, unsigned nslices);
  
private:
  void CreateRealDataSet(const std::string &filename, std::vector<unsigned> size);
  void CreateComplexDataSet(const std::string &filename, std::vector<unsigned> size);

  H5::H5File *m_file;
  hid_t m_config;
  std::string m_run_id;
  bool save_individual_wave;
  bool save_individual_dp;
  bool save_individual_pot_slice;
  bool save_individual_pot_volume;
  bool save_average_pot_volume;
};

#if FLOAT_PRECISION == 1
typedef H5::PredType::NATIVE_FLOAT QH5_NATIVE_FLOAT
#else
typedef H5::PredType::NATIVE_DOUBLE QH5_NATIVE_FLOAT
#endif

// Create the complex data type for complex image files
hid_t QH5_NATIVE_COMPLEX = H5Tcreate (H5T_COMPOUND, sizeof (complex_tt));
H5Tinsert (complex_id, “real”, HOFFSET(complex_tt,re),
           QH5_NATIVE_FLOAT);
H5Tinsert (complex_id, “imag”, HOFFSET(complex_tt,im),
           QH5_NATIVE_FLOAT); 


typedef boost::shared_ptr<QH5> QH5ptr;

















