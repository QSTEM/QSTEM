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

// TODO: implement some kind of dynamic registration here, like maybe:
// http://stackoverflow.com/questions/9975672/c-automatic-factory-registration-of-derived-types

#ifndef CONFIG_READER_H
#define CONFIG_READER_H

#include "config_IO/read_interface.hpp"
#include "config_IO/read_qsc.hpp"

ConfigReaderPtr QSTEM_HELPER_DLL_EXPORT GetConfigReader(std::string filename);

#endif
