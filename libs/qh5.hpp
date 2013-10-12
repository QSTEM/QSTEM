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

#include "H5Cpp.h"

// This class provides utilities for opening hdf5 files and for creating default data sets.

class QH5
{
public:
  void OpenFile(const std::string &filename);
  void CreateRealDataSet(const std::string &filename, std::vector<unsigned> size);
  void CreateComplexDataSet(const std::string &filename, std::vector<unsigned> size);
  
private:
  H5::H5File *m_file;
};

typedef boost::shared_ptr<QH5> QH5ptr;














