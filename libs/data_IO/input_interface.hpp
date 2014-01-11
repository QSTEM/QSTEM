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


class IDataInput
{
public:
  virtual void ReadImage(void **pix, const std::string label, std::map<std::string, double> &params,
                         std::string &comment, const std::vector<unsigned> position=std::vector<unsigned>())=0;
};

typedef boost::shared_ptr<IDataInput> DataReaderPtr;

#endif
