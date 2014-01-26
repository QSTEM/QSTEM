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

#ifndef CONFIG_READER_FACTORY_H
#define CONFIG_READER_FACTORY_H

#include <map>

// Include this so that people know what a ConfigReaderPtr is and how to create them
#include "read_interface.hpp"

// Factory for creating instances of IPotential
class QSTEM_HELPER_DLL_EXPORT CConfigReaderFactory
{
public:
    ~CConfigReaderFactory() { m_FactoryMap.clear(); }

    static CConfigReaderFactory *Get()
    {
        static CConfigReaderFactory instance;
        return &instance;
    }

  void Register(const std::string &extension, CreateReaderFn pfnCreate);
  // Looks up which reader to get based on string mapping of registered readers
  ConfigReaderPtr GetReader(const std::string &animalName);

private:
  CConfigReaderFactory();
  CConfigReaderFactory(const CConfigReaderFactory &) { }
  CConfigReaderFactory &operator=(const CConfigReaderFactory &) { return *this; }

  typedef std::map<std::string, CreateReaderFn> FactoryMap;
  FactoryMap m_FactoryMap;
};

#endif







