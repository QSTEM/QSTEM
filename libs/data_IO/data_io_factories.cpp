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

#include "data_io_factories.hpp"


/************** data readers *************/

#include "img_input.hpp"

CDataReaderFactory::CDataReaderFactory()
{
  Register(".img",    &CImgReader::Create);
}

void CDataReaderFactory::Register(const std::string &extension, CreateDataReaderFn pfnCreate)
{
	m_FactoryMap[extension] = pfnCreate;
}

DataReaderPtr CDataReaderFactory::GetReader(const std::string &extension)
{
  FactoryMap::iterator it = m_FactoryMap.find(extension);
  if( it != m_FactoryMap.end() )
    return it->second();
  return DataReaderPtr();
}


/********** data writers **************/

#include "img_output.hpp"

CDataWriterFactory::CDataWriterFactory()
{
  Register(".img",    &CImgWriter::Create);
}

void CDataWriterFactory::Register(const std::string &extension, CreateDataWriterFn pfnCreate)
{
	m_FactoryMap[extension] = pfnCreate;
}

DataWriterPtr CDataWriterFactory::GetWriter(const std::string &extension)
{
  FactoryMap::iterator it = m_FactoryMap.find(extension);
  if( it != m_FactoryMap.end() )
    return it->second();
  return DataWriterPtr();
}

