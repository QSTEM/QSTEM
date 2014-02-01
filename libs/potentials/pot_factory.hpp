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

#ifndef POT_FACTORY_H
#define POT_FACTORY_H

#include <map>

// Include this so that people know what a PotPtr is and how to create them
#include "pot_interface.hpp"

// Include this because we use ConfigReaderPtr below.
#include "config_IO/config_reader_factory.hpp"


// Factory for creating instances of IPotential
class QSTEM_HELPER_DLL_EXPORT CPotFactory
{
private:
  CPotFactory();
  CPotFactory(const CPotFactory &) { }
  CPotFactory &operator=(const CPotFactory &) { return *this; }

  typedef std::map<std::string, CreatePotentialFn> FactoryMap;
  FactoryMap m_FactoryMap;
public:
    ~CPotFactory() { m_FactoryMap.clear(); }

    static CPotFactory *Get()
    {
        static CPotFactory instance;
        return &instance;
    }

    void Register(const std::string &potentialName, CreatePotentialFn pfnCreate);
  // Looks up which potential to get based on string mapping of registered potentials
  PotPtr GetPotential(const std::string &animalName);
  // ultimately uses the string based method, but builds the string based on the standard 2d/3d fft/real-space options
  PotPtr GetPotential(bool _3D, bool fft);
  // ultimately uses the string-based method, but parses the config reader for you to make that string.  
  //    Has side effect of initializing pot automatically for you using the passed in configReader.
  PotPtr GetPotential(const ConfigReaderPtr &configReader);
};

#endif
