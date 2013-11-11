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

#include "config_readers.hpp"
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

ConfigReaderPtr GetConfigReader(std::string &filename)
{
  boost::filesystem::path filepath( filename );
  std::string extension = filepath.extension().string();
  boost::algorithm::to_lower(extension);
  // check file extension, instantiate appropriate reader
  if (extension == ".qsc")
    return ConfigReaderPtr(new CQscReader(filename));
  else if (extension == ".qh5")
    {
      // TODO: flesh out qh5 as config option
    }
}















