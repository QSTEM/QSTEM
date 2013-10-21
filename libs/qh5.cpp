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

#include "qh5.hpp"

#include <sstream>
#include <stdexcept>

QH5::QH5() :
  m_file(-1)
{
#if FLOAT_PRECISION == 1
  QH5_NATIVE_FLOAT = H5Tcopy(H5T_NATIVE_FLOAT);
#else
  QH5_NATIVE_FLOAT = H5Tcopy(H5T_NATIVE_DOUBLE);
#endif
  H5Tlock(QH5_NATIVE_FLOAT);

  // Create the complex data type for complex image files
  QH5_NATIVE_COMPLEX = H5Tcreate (H5T_COMPOUND, sizeof(complex_tt));
  H5Tinsert(QH5_NATIVE_COMPLEX, "real", 0, 
             QH5_NATIVE_FLOAT);
  H5Tinsert(QH5_NATIVE_COMPLEX, "imag", sizeof(float_tt),
             QH5_NATIVE_FLOAT);
  H5Tlock(QH5_NATIVE_COMPLEX);
}

QH5::QH5(std::string fileName)
{
  QH5();
  OpenFile(fileName);
}

QH5::~QH5()
{
  H5Fclose(m_file);
}

void QH5::OpenFile(std::string fileName)
{
  int file = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file<0)
    {
      file = H5Fcreate(fileName.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    }
  if (file < 0)
    {
      throw std::runtime_error("Error opening HDF5 file");
    }

  int group = H5Gcreate(file, "Runs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group);
  // TODO: add group for config
  //group = H5Gcreate(file, "");
  //H5Gclose(group);
  
  m_file = file;
}

/**
Sets the group under which datasets will go.  Creates the group as necessary.
 */
void QH5::SetRunID(std::string run_id)
{
  m_run_id = run_id;
  int run_group;
  // Make sure the group exists
  run_group = H5Gopen(m_file, GetRunPath().c_str(), H5P_DEFAULT);
  if (run_group<0)
    {
      run_group = H5Gcreate(m_file, GetRunPath().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
  H5Gclose(run_group);
}

std::string QH5::GetRunPath()
{
  std::stringstream path;
  path << "/Runs/" << m_run_id;
  return path.str();
}

std::string QH5::GetDataSetPath(std::string dataSetName)
{
  std::stringstream path;
  path << GetRunPath() << "/" << dataSetName;
  return path.str();
}

std::string QH5::GetDetectorPath(std::string detectorName)
{
  std::stringstream path;
  path << GetDataSetPath("Detectors") << "/" << detectorName;
  return path.str();
}

void QH5::CreateRealDataSet(std::string path, unsigned int size_x, unsigned int size_y, 
                            std::vector<unsigned> &position)
{
  std::vector<unsigned>::iterator pos;
  std::vector<hsize_t> dims(2+position.size(),0);
  dims[0]=size_x;
  dims[1]=size_y;
  
  for (unsigned pos=0; pos<position.size(); pos++)
    {
      dims[2+pos] = position[pos];
    }
  
  hsize_t space = H5Screate_simple(2+position.size(), &dims[0], NULL);
  // Try opening it first
  hsize_t data = H5Dopen(m_file, path.c_str(), H5P_DEFAULT);

  // If not there, create it
  if (data < 0)
    {
      data = H5Dcreate(m_file, path.c_str(), QH5_NATIVE_FLOAT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

  H5Sclose(space);

  if (data <0)
    {
      throw std::runtime_error("Error creating dataset");
    }

  // The dataset is created - close it.  When anyone else wants it, they should open it, then close 
  //     it when they're done.
  H5Dclose(data);
}

void QH5::CreateComplexDataSet(std::string path, unsigned int size_x, unsigned int size_y, 
                               std::vector<unsigned> &position)
{
  std::vector<unsigned>::iterator pos;
  std::vector<hsize_t> dims(2+position.size(),0); 
  dims[0]=size_x;
  dims[1]=size_y;
  
  for (unsigned pos=0; pos<position.size(); pos++)
    {
      dims[2+pos] = position[pos];
    }
  
  hsize_t space = H5Screate_simple(2+position.size(), &dims[0], NULL);
  // Try opening it first
  hsize_t data = H5Dopen(m_file, path.c_str(), H5P_DEFAULT);

  // If not there, create it
  if (data < 0)
    {
      data = H5Dcreate(m_file, path.c_str(), QH5_NATIVE_FLOAT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

  H5Sclose(space);

  if (data <0)
    {
      throw std::runtime_error("Error creating dataset");
    }

  // The dataset is created - close it.  When anyone else wants it, they should open it, then close 
  //     it when they're done.
  H5Dclose(data);
}

void QH5::WriteRealDataSlab(float_tt *pix, std::string path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters)
{
  hsize_t memory_dims[2]={size_x, size_y};
  // Position offset is where we sit in the n-dimensional dataset.  In other words,
  //    which position, or which slice, are we going to write?
  std::vector<hsize_t> position_offset(2+position.size(), 0);
  
  for (size_t idx=0; idx<position.size(); idx++)
    {
      position_offset[2+idx]=position[idx];
    }

  hid_t data = H5Dopen(m_file, path.c_str(), H5P_DEFAULT);
  hid_t dataspace = H5Dget_space(data);
  hid_t memoryspace = H5Screate_simple(2, memory_dims, NULL);
  int status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &position_offset[0], NULL, memory_dims, NULL);
  
  H5Dwrite(data, QH5_NATIVE_FLOAT, memoryspace, dataspace, H5P_DEFAULT, (void *)pix);
  H5Sclose(dataspace);
  H5Sclose(memoryspace);
  H5Dclose(data);
}

void QH5::WriteComplexDataSlab(complex_tt *pix, std::string path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters)
{
  hsize_t memory_dims[2]={size_x, size_y};
  // Position offset is where we sit in the n-dimensional dataset.  In other words,
  //    which position, or which slice, are we going to write?
  std::vector<hsize_t> position_offset(2+position.size(), 0);
  
  for (size_t idx=0; idx<position.size(); idx++)
    {
      position_offset[2+idx]=position[idx];
    }

  hid_t data = H5Dopen(m_file, path.c_str(), H5P_DEFAULT);
  hid_t dataspace = H5Dget_space(data);
  hid_t memoryspace = H5Screate_simple(2, memory_dims, NULL);
  int status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &position_offset[0], NULL, memory_dims, NULL);
  
  H5Dwrite(data, QH5_NATIVE_COMPLEX, memoryspace, dataspace, H5P_DEFAULT, (void *)pix);
  H5Sclose(dataspace);
  H5Sclose(memoryspace);
  H5Dclose(data);
}

void QH5::ReadRealDataSlab(float_tt *pix, std::string path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters)
{
  hsize_t memory_dims[2]={size_x, size_y};
  // Position offset is where we sit in the n-dimensional dataset.  In other words,
  //    which position, or which slice, are we going to write?
  std::vector<hsize_t> position_offset(2+position.size(), 0);
  
  for (size_t idx=0; idx<position.size(); idx++)
    {
      position_offset[2+idx]=position[idx];
    }

  hid_t data = H5Dopen(m_file, path.c_str(), H5P_DEFAULT);
  hid_t dataspace = H5Dget_space(data);
  hid_t memoryspace = H5Screate_simple(2, memory_dims, NULL);
  int status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &position_offset[0], NULL, memory_dims, NULL);
  
  H5Dread(data, QH5_NATIVE_FLOAT, memoryspace, dataspace, H5P_DEFAULT, (void *)pix);
  H5Sclose(dataspace);
  H5Sclose(memoryspace);
  H5Dclose(data);
}

void QH5::ReadComplexDataSlab(complex_tt *pix, std::string path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters)
{
  hsize_t memory_dims[2]={size_x, size_y};
  // Position offset is where we sit in the n-dimensional dataset.  In other words,
  //    which position, or which slice, are we going to write?
  std::vector<hsize_t> position_offset(2+position.size(), 0);
  
  for (size_t idx=0; idx<position.size(); idx++)
    {
      position_offset[2+idx]=position[idx];
    }

  hid_t data = H5Dopen(m_file, path.c_str(), H5P_DEFAULT);
  hid_t dataspace = H5Dget_space(data);
  hid_t memoryspace = H5Screate_simple(2, memory_dims, NULL);
  int status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &position_offset[0], NULL, memory_dims, NULL);
  
  H5Dread(data, QH5_NATIVE_COMPLEX, memoryspace, dataspace, H5P_DEFAULT, (void *)pix);
  H5Sclose(dataspace);
  H5Sclose(memoryspace);
  H5Dclose(data);
}

