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

  inline void CreatePotentialDataSet(unsigned size_x, unsigned size_y, unsigned size_z)
  { CreateComplexDataSet(GetDataSetPath("Potential"), size_x, size_y, size_z); }
  inline void CreateWaveDataSet(unsigned size_x, unsigned size_y, std::vector<unsigned> &positions)
  { CreateComplexDataSet(GetDataSetPath("WaveFunctions"), size_x, size_y, positions); }
  inline void CreateDPDataSet(unsigned size_x, unsigned size_y, std::vector<unsigned> &positions)
  { CreateRealDataSet(GetDataSetPath("DiffractionPatterns"), size_x, size_y, positions); }
  inline void CreateDetectorDataSet(std::string name, unsigned size_x, unsigned size_y, unsigned nslices)
  { CreateRealDataSet(GetDetectorPath(name), size_x, size_y, nslices); }


  inline void WritePotentialSlice(complex_tt *pot, unsigned size_x, unsigned size_y, unsigned slice)
  { WriteComplexDataSlab(GetDataSetPath("Potential"), pot, slice); }
  inline void WriteWave(complex_tt *wave, unsigned size_x, unsigned size_y, std::vector<unsigned> &position)
  { WriteComplexDataSlab(GetDataSetPath("WaveFunctions"), wave, position); }
  inline void WriteDiffractionPattern(float_tt *dp, unsigned size_x, unsigned size_y, std::vector<unsigned> &position)
  { WriteRealDataSlab(GetDataSetPath("DiffractionPatterns"), dp, position); }
  inline void WriteDetector(float_tt *det, std::string name, unsigned size_x, unsigned size_y, unsigned slice)
  { WriteRealDataSlab(GetDetectorPath(name), det, slice); }

  inline void ReadPotentialSlice(complex_tt *pot, unsigned size_x, unsigned size_y, unsigned slice);
  inline void ReadWave(complex_tt *wave, unsigned size_x, unsigned size_y, std::vector<unsigned> position);
  inline void ReadDiffractionPattern(float_tt *dp, unsigned size_x, unsigned size_y, std::vector<unsigned> &position);
  inline void ReadDetector(float_tt *det, std::string name, unsigned size_x, unsigned size_y, unsigned slice);
  
private:
  std::string GetDataSetPath(std::string &dataSetName);
  std::string GetDetectorPath(std::string &detectorName);

  void CreateComplexDataSet(std::string &path, unsigned size_x, unsigned size_y, std::vector<unsigned> &positions);
  inline void CreateComplexDataSet(std::string &path, unsigned size_x, unsigned size_y, unsigned position)
  {
    std::vector<unsigned> pos(1);
    pos[0]=position;
    CreateComplexDataSet(path, size_x, size_y, pos);
  }
  void CreateRealDataSet(std::string &path, unsigned size_x, unsigned size_y, std::vector<unsigned> &positions);
  inline void CreateRealDataSet(std::string &path, unsigned size_x, unsigned size_y, unsigned position)
  {
    std::vector<unsigned> pos(1);
    pos[0]=position;
    CreateRealDataSet(path, size_x, size_y, pos);
  }
  void WriteRealDataSlab(float_tt *pix, const std::string &path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters);
  void WriteComplexDataSlab(complex_tt *pix, const std::string &path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters);

  void ReadRealDataSlab(float_tt *pix, const std::string &path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters);
  void ReadComplexDataSlab(complex_tt *pix, const std::string &path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters);

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

















