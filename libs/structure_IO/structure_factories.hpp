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

#ifndef STRUCTURE_FACTORIES_H
#define STRUCTURE_FACTORIES_H

#include <map>

#include "structureInterface.hpp"
#include <boost/filesystem.hpp>

// Factory for creating instances of IStructureReader
class QSTEM_HELPER_DLL_EXPORT CStructureReaderFactory
{
public:
  ~CStructureReaderFactory() { m_FactoryMap.clear(); }

  static CStructureReaderFactory *Get();

  void Register(const std::string &extension, CreateStructureReaderFn pfnCreate);
  // Looks up which reader to get based on string mapping of registered readers
  StructureReaderPtr GetReader(const boost::filesystem::path &filename);

private:
  CStructureReaderFactory();
  CStructureReaderFactory(const CStructureReaderFactory &) { }
  CStructureReaderFactory &operator=(const CStructureReaderFactory &) { return *this; }

  typedef std::map<std::string, CreateStructureReaderFn> FactoryMap;
  FactoryMap m_FactoryMap;
};


// Factory for creating instances of IStructureWriter
class QSTEM_HELPER_DLL_EXPORT CStructureWriterFactory
{
public:
  ~CStructureWriterFactory() { m_FactoryMap.clear(); }

  static CStructureWriterFactory *Get();  

  void Register(const std::string &extension, CreateStructureWriterFn pfnCreate);
  // Looks up which reader to get based on string mapping of registered readers
  StructureWriterPtr GetWriter(const boost::filesystem::path &filename, float_tt ax, float_tt by, float_tt cz);

private:
  CStructureWriterFactory();
  CStructureWriterFactory(const CStructureWriterFactory &) { }
  CStructureWriterFactory &operator=(const CStructureWriterFactory &) { return *this; }

  typedef std::map<std::string, CreateStructureWriterFn> FactoryMap;
  FactoryMap m_FactoryMap;
};


#endif
