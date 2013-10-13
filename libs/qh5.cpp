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

std::string QH5::GetDataSetPath(std::string dataSetName)
{
  std::stringstream path;
  path << GetRunPath() << "/" << dataSetName;
  return path.str();
}

/**
Create a dataset for the 3D potential - either the 3D volume, or slices.
 */
void QH5::CreatePotentialVolumeDataSet(unsigned int size_x, unsigned int size_y, unsigned int size_z)
{
  hsize_t dims[3] = {size_x, size_y, size_z};
  DataSpace ds (3, dims);

  DataSet data = m_file.createDataSet(GetDataSetPath("Potential volume"), QH5_NATIVE_COMPLEX, ds);
}

void QH5::CreatePotentialDataSet(unsigned int size_x, unsigned int size_y, unsigned int size_z)
{
  hsize_t dims[3] = {size_x, size_y, size_z};
  DataSpace ds (3, dims);

  DataSet data = m_file.createDataSet(GetDataSetPath("Potential slices"), QH5_NATIVE_COMPLEX, ds);
}


void QH5::CreateWaveDataSet(unsigned int size_x, unsigned int size_y, std::vector<unsigned> positions)
{
  hsize_t dims = new hsize_t(2+positions.size());
  DataSpace ds (2+positions.size(), dims);

  DataSet data = m_file.createDataSet(GetDataSetPath("Potential volume"), QH5_NATIVE_COMPLEX, ds);  
  
  delete dims;
}

void QH5::CreateDPDataSet(unsigned int size_x, unsigned int size_y, std::vector<unsigned> positions)
{
  hsize_t dims = new hsize_t(2+positions.size());
  DataSpace ds (2+positions.size(), dims);

  DataSet data = m_file.createDataSet(GetDataSetPath("Potential volume"), QH5_NATIVE_FLOAT, ds);  
  
  delete dims;
}

void QH5::CreateDetectorDataSet(std::string name, unsigned int size_x, unsigned int size_y, unsigned int nslices)
{
  hsize_t dims[3]={size_x, size_y, size_z};
  DataSpace ds (3, dims);

  // TODO: need to make sure that detectors group has been created before creating this dataset!

  std::stringstream path;
  path << "Detectors/" << name;

  DataSet data = m_file.createDataSet(GetDataSetPath(path.str()), QH5_NATIVE_FLOAT, ds);  
}











