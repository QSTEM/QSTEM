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

#include "qh5_input.hpp"

CQH5Input::CQH5Input()
{}

CQH5Input::~CQH5Input()
{
status = H5Fclose (m_file);
}

void CQH5Input::ReadImage(void **pix, std::string label, 
                 std::map<std::string, double> &params,
                 std::string &comment, 
                 std::vector<unsigned> position=std::vector<unsigned>())
{
  DataSet ds = m_file.openDataSet(label);
  DataSpace space = ds.getSpace();
  int ndim = dataspace.getSimpleExtentNdims();
  hsize_t offset[2];
  hsize_t img_size[2];

  space.getSimpleExtentDims(img_size, NULL);

  DataSpace memspace (2, img_size);


}

void CQH5Input::SetFile(std::string filename)
{
  m_file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY);
}
















