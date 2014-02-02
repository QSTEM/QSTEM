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

#include <boost/algorithm/string.hpp>

#include "structure_factories.hpp"
#include "cfg.hpp"

/************** Structure readers *************/

CStructureReaderFactory::CStructureReaderFactory()
{
  Register(".cfg",    &CCfgReader::Create);
}

CStructureReaderFactory *CStructureReaderFactory::Get()
{
static CStructureReaderFactory instance;
return &instance;
}

void CStructureReaderFactory::Register(const std::string &extension, CreateStructureReaderFn pfnCreate)
{
	m_FactoryMap[extension] = pfnCreate;
}

StructureReaderPtr CStructureReaderFactory::GetReader(const boost::filesystem::path &filepath)
{
  std::string extension = filepath.extension().string();
  boost::algorithm::to_lower(extension);
  FactoryMap::iterator it = m_FactoryMap.find(extension);
  if( it != m_FactoryMap.end() )
    return it->second(filepath);
  return StructureReaderPtr();
}


/********** Structure writers **************/

CStructureWriterFactory::CStructureWriterFactory()
{
  Register(".cfg",    &CCfgWriter::Create);
}

CStructureWriterFactory *CStructureWriterFactory::Get()
{
static CStructureWriterFactory instance;
return &instance;
}

void CStructureWriterFactory::Register(const std::string &extension, CreateStructureWriterFn pfnCreate)
{
	m_FactoryMap[extension] = pfnCreate;
}

StructureWriterPtr CStructureWriterFactory::GetWriter(const boost::filesystem::path &filename, 
                                                        float_tt ax, float_tt by, float_tt cz)
{
std::string extension = filename.extension().string();
  boost::algorithm::to_lower(extension);
  FactoryMap::iterator it = m_FactoryMap.find(extension);
  if( it != m_FactoryMap.end() )
    return it->second(filename, ax, by, cz);
  return StructureWriterPtr();
}
