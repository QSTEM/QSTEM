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

#include "hdf5.h"
//#include "H5Cpp.h"

#include <string>
#include <vector>
#include <map>

// This class provides utilities for opening hdf5 files and for creating default data sets.

class QH5
{
public:
  QH5();
  QH5(std::string fileName);
  ~QH5();

  void OpenFile(std::string fileName);

  void SetRunID(std::string run_id);

  void GetConfigGroup();

  inline void CreatePotentialDataSet(unsigned size_x, unsigned size_y, unsigned size_z)
  { return CreateComplexDataSet(GetDataSetPath("Potential"), size_x, size_y, size_z); }
  inline void CreateWaveDataSet(unsigned size_x, unsigned size_y, std::vector<unsigned> &positions)
  { return CreateComplexDataSet(GetDataSetPath("WaveFunctions"), size_x, size_y, positions); }
  inline void CreateDPDataSet(unsigned size_x, unsigned size_y, std::vector<unsigned> &positions)
  { return CreateRealDataSet(GetDataSetPath("DiffractionPatterns"), size_x, size_y, positions); }
  inline void CreateDetectorDataSet(std::string name, unsigned size_x, unsigned size_y, unsigned nslices)
  { return CreateRealDataSet(GetDetectorPath(name), size_x, size_y, nslices); }

  /************* Write image methods ***************/

  inline void WritePotentialSlice(complex_tt *pot, unsigned size_x, unsigned size_y, unsigned slice, 
                                  std::map<std::string, double> &parameters)
  { 
    return WriteComplexDataSlab(pot, GetDataSetPath("Potential"), size_x, size_y, slice, parameters); 
  }
  inline void WriteWave(complex_tt *wave, unsigned size_x, unsigned size_y, std::vector<unsigned> &position,
                        std::map<std::string, double> &parameters)
  { 
    return WriteComplexDataSlab(wave, GetDataSetPath("WaveFunctions"), size_x, size_y, position, parameters); 
  }
  inline void WriteDiffractionPattern(float_tt *dp, unsigned size_x, unsigned size_y, std::vector<unsigned> &position,
                                      std::map<std::string, double> &parameters)
  { 
    return WriteRealDataSlab(dp, GetDataSetPath("DiffractionPatterns"), size_x, size_y, position, parameters); 
  }
  inline void WriteDetector(float_tt *det, std::string name, unsigned size_x, unsigned size_y, unsigned slice,
                            std::map<std::string, double> &parameters)
  { 
    return WriteRealDataSlab(det, GetDetectorPath(name), size_x, size_y, slice, parameters); 
  }

  /*************  Read image methods ************/

  inline void ReadPotentialSlice(complex_tt *pot, unsigned size_x, unsigned size_y, unsigned slice, 
                                 std::map<std::string, double> &parameters)
  { 
    return ReadComplexDataSlab(pot, GetDataSetPath("Potential"), size_x, size_y, slice, parameters);
  }
  inline void ReadWave(complex_tt *wave, unsigned size_x, unsigned size_y, std::vector<unsigned> position,
                       std::map<std::string, double> &parameters)
  {
    return ReadComplexDataSlab(wave, GetDataSetPath("WaveFunctions"), size_x, size_y, position, parameters);
  }
  inline void ReadDiffractionPattern(float_tt *dp, unsigned size_x, unsigned size_y, std::vector<unsigned> &position,
                                     std::map<std::string, double> &parameters)
  { 
    return ReadRealDataSlab(dp, GetDataSetPath("DiffractionPatterns"), size_x, size_y, position, parameters);
  }
  inline void ReadDetector(float_tt *det, std::string name, unsigned size_x, unsigned size_y, unsigned slice,
                           std::map<std::string, double> &parameters)
  { 
    return ReadRealDataSlab(det, GetDetectorPath(name), size_x, size_y, slice, parameters);
  }
  
private:
  std::string GetRunPath();

  std::string GetDataSetPath(std::string dataSetName);
  inline std::string GetDataSetPath(const char *dataSetName)
  {return GetDataSetPath(std::string(dataSetName));}
  std::string GetDetectorPath(std::string detectorName);
  inline std::string GetDetectorPath(const char *detectorName)
  {return GetDetectorPath(std::string(detectorName));}

  void CreateComplexDataSet(std::string path, unsigned size_x, unsigned size_y, std::vector<unsigned> &positions);
  inline void CreateComplexDataSet(std::string path, unsigned size_x, unsigned size_y, unsigned position)
  {
    std::vector<unsigned> pos(1);
    pos[0]=position;
    CreateComplexDataSet(path, size_x, size_y, pos);
  }

  void CreateRealDataSet(std::string path, unsigned size_x, unsigned size_y, std::vector<unsigned> &positions);
  inline void CreateRealDataSet(std::string path, unsigned size_x, unsigned size_y, unsigned position)
  {
    std::vector<unsigned> pos(1);
    pos[0]=position;
    CreateRealDataSet(path, size_x, size_y, pos);
  }

  void WriteRealDataSlab(float_tt *pix, std::string path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters);
  inline void WriteRealDataSlab(float_tt *pix, std::string path, unsigned size_x, unsigned size_y,
                                unsigned slice, std::map<std::string, double> &parameters)
  {
    std::vector<unsigned> slice_vec(1);
    slice_vec[0]=slice;
    WriteRealDataSlab(pix, path, size_x, size_y, slice_vec, parameters);
  }

  void WriteComplexDataSlab(complex_tt *pix, std::string path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters);
  inline void WriteComplexDataSlab(complex_tt *pix, std::string path, unsigned size_x, unsigned size_y,
                                unsigned slice, std::map<std::string, double> &parameters)
  {
    std::vector<unsigned> slice_vec(1);
    slice_vec[0]=slice;
    WriteComplexDataSlab(pix, path, size_x, size_y, slice_vec, parameters);
  }

  void ReadRealDataSlab(float_tt *pix, std::string path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters);
  inline void ReadRealDataSlab(float_tt *pix, std::string path, unsigned size_x, unsigned size_y, 
                        unsigned slice, std::map<std::string, double> &parameters)
  {
    std::vector<unsigned> slice_vec(1);
    slice_vec[0]=slice;
    ReadRealDataSlab(pix, path, size_x, size_y, slice_vec, parameters);
  }

  void ReadComplexDataSlab(complex_tt *pix, std::string path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters);
  inline void ReadComplexDataSlab(complex_tt *pix, std::string path, unsigned size_x, unsigned size_y, 
                        unsigned slice, std::map<std::string, double> &parameters)
  {
    std::vector<unsigned> slice_vec(1);
    slice_vec[0]=slice;
    ReadComplexDataSlab(pix, path, size_x, size_y, slice_vec, parameters);
  }

  static int ReadParameterOp(hid_t location_id, const char *attr_name, const H5A_info_t *ainfo, void *op_data);
  void ReadParameters(hid_t dataset, std::map<std::string, double> &parameters);
  void WriteParameters(hid_t dataset, std::map<std::string, double> &parameters);

  hid_t m_file;
  hid_t m_config;
  // This is the complex datatype for storage purposes
  hid_t QH5_NATIVE_FLOAT, QH5_NATIVE_COMPLEX;
  std::string m_run_id;
  bool save_individual_wave;
  bool save_individual_dp;
  bool save_individual_pot_slice;
  
};

typedef boost::shared_ptr<QH5> QH5ptr;

#endif



