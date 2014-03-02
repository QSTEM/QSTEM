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

#ifndef STRUCTURE_READERS_H
#define STRUCTURE_READERS_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "structure_IO/cfg.hpp"

namespace QSTEM
{

inline StructureReaderPtr GetStructureReader(boost::filesystem::path &structure_file)
{
  std::string extension = structure_file.extension().string();
  boost::algorithm::to_lower(extension);
  if (extension == ".cfg")
    return StructureReaderPtr(new CCfgReader(structure_file));
}

}
#endif
