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

// create or open a file
H5FilePtr QH5::OpenFile(const std::string& fname)
{
    H5::Exception::dontPrint();

    try {
        m_file = new H5::H5File(fname.c_str(), H5F_ACC_RDWR);
    } catch(const H5::FileIException&) {
        m_file = new H5::H5File(fname.c_str(), H5F_ACC_TRUNC);
    }
}

/**
Sets the group under which datasets will go.  Creates the group as necessary.
 */
void QH5::SetRunID(std::string run_id)
{
  m_run_id = run_id;
}

std::string QH5::GetRunPath()
{
  std::stringstream path;
  path << "/Runs/" << m_run_id;
  return path.str();
}

std::string QH5::GetDataSetPath(std::string &dataSetName)
{
  std::stringstream path;
  path << GetRunPath() << "/" << dataSetName;
  return path.str();
}

std::string QH5::GetDetectorPath(std::string &detectorName)
{
  std::stringstream path;
  path << GetDataSetPath("Detectors") << "/" << detectorName;
  return path.str();
}

void QH5::CreateRealDataSet(std::string &path, unsigned int size_x, unsigned int size_y, 
                            std::vector<unsigned> position)
{
  hsize_t dims = new hsize_t(2+position.size());
  DataSpace ds (2+position.size(), dims);

  DataSet data = m_file.createDataSet(path, QH5_NATIVE_FLOAT, ds);  
  
  delete dims;
}

void QH5::CreateComplexDataSet(std::string &path, unsigned int size_x, unsigned int size_y, 
                               std::vector<unsigned> position)
{
  hsize_t dims = new hsize_t(2+position.size());
  DataSpace ds (2+position.size(), dims);

  DataSet data = m_file.createDataSet(path, QH5_NATIVE_COMPLEX, ds);  
  
  delete dims;
}

void WriteRealDataSlab(float_tt *pix, const std::string &path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters)
{
  hsize_t memory_dims[2]={size_x, size_y};
  DataSet ds = m_file.openDataSet(path);
  DataSpace memspace(2, memory_dims);
  DataSpace filespace = ds.getSpace();
  // Position offset is where we sit in the n-dimensional dataset.  In other words,
  //    which position, or which slice, are we going to write?
  hsize_t position_offset = new hsize_t(2+position.size(),0);
  
  for (size_t idx=0; idx<position.size(); idx++)
    {
      position_offset[2+idx]=position[idx];
    }

  filespace.selectHyperslab( H5S_SELECT_SET, memory_dims, position_offset );
  
  ds.write((void *)pix, QH5_NATIVE_FLOAT, memspace, filespace);
  
  delete position_offset;
}

void WriteComplexDataSlab(complex_tt *pix, const std::string &path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters)
{
  hsize_t memory_dims[2]={size_x, size_y};
  DataSet ds = m_file.openDataSet(path);
  DataSpace memspace(2, memory_dims);
  DataSpace filespace = ds.getSpace();
  // Position offset is where we sit in the n-dimensional dataset.  In other words,
  //    which position, or which slice, are we going to write?
  hsize_t position_offset = new hsize_t(2+position.size(),0);
  
  for (size_t idx=0; idx<position.size(); idx++)
    {
      position_offset[2+idx]=position[idx];
    }

  filespace.selectHyperslab( H5S_SELECT_SET, memory_dims, position_offset );
  
  ds.write((void *)pix, QH5_NATIVE_COMPLEX, memspace, filespace);
  
  delete position_offset;
}

void ReadRealDataSlab(float_tt *pix, const std::string &path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters)
{
  hsize_t memory_dims[2]={size_x, size_y};
  DataSet ds = m_file.openDataSet(path);
  DataSpace memspace(2, memory_dims);
  DataSpace filespace = ds.getSpace();
  // Position offset is where we sit in the n-dimensional dataset.  In other words,
  //    which position, or which slice, are we going to write?
  hsize_t position_offset = new hsize_t(2+position.size(),0);
  
  for (size_t idx=0; idx<position.size(); idx++)
    {
      position_offset[2+idx]=position[idx];
    }

  filespace.selectHyperslab( H5S_SELECT_SET, memory_dims, position_offset );
  
  ds.read((void *)pix, QH5_NATIVE_FLOAT, memspace, filespace);
  
  delete position_offset;
}

void ReadComplexDataSlab(complex_tt *pix, const std::string &path, unsigned size_x, unsigned size_y, 
                        std::vector<unsigned> &position, std::map<std::string, double> &parameters)
{
  hsize_t memory_dims[2]={size_x, size_y};
  DataSet ds = m_file.openDataSet(path);
  DataSpace memspace(2, memory_dims);
  DataSpace filespace = ds.getSpace();
  // Position offset is where we sit in the n-dimensional dataset.  In other words,
  //    which position, or which slice, are we going to write?
  hsize_t position_offset = new hsize_t(2+position.size(),0);
  
  for (size_t idx=0; idx<position.size(); idx++)
    {
      position_offset[2+idx]=position[idx];
    }

  filespace.selectHyperslab( H5S_SELECT_SET, memory_dims, position_offset );
  
  ds.read((void *)pix, QH5_NATIVE_COMPLEX, memspace, filespace);
  
  delete position_offset;
}
