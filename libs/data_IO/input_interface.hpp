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

#ifndef INPUT_INTERFACE_H
#define INPUT_INTERFACE_H

#include <map>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

class IDataReader;
typedef boost::shared_ptr<IDataReader> DataReaderPtr;
typedef DataReaderPtr (*CreateDataReaderFn)();

class IDataReader
{
public:
  virtual void ReadImageData(const std::string &filename, void *pix)=0;
  virtual void ReadComment(const std::string &filename, std::string &comment)=0;
  virtual void ReadParameters(const std::string &filename, std::map<std::string, double> &parameters)=0;
  virtual void ReadSize(const std::string &filename, unsigned &nx, unsigned &ny)=0;
  virtual void ReadElementByteSize(const std::string &filename, unsigned &elementByteSize)=0;
  virtual void ReadComplex(const std::string &filename, bool &complex)=0;
  // One-shot catch-all
  virtual void ReadImage(const std::string &filename, void *pix, std::map<std::string, double> &params,
                         std::string &comment)=0;
};


#endif













