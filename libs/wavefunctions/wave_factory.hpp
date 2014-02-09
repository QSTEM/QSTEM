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

#ifndef WAVE_FACTORY_H
#define WAVE_FACTORY_H

#include <map>

#include "wave_interface.hpp"
#include "config_IO/config_reader_factory.hpp"
#include <boost/filesystem.hpp>

// Factory for creating instances of IWave
class QSTEM_HELPER_DLL_EXPORT CWaveFactory
{
public:
  ~CWaveFactory() { m_FactoryMap.clear(); }

  static CWaveFactory *Get();

  void Register(const std::string &type, CreateWaveFn pfnCreate);
  // Looks up which reader to get based on string mapping of registered readers
  WavePtr GetWave(const std::string &type, const ConfigReaderPtr &reader);

private:
  CWaveFactory();
  CWaveFactory(const CWaveFactory &) { }
  CWaveFactory &operator=(const CWaveFactory &) { return *this; }

  typedef std::map<std::string, CreateWaveFn> FactoryMap;
  FactoryMap m_FactoryMap;
};

#endif