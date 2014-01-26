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

#ifndef QH5_INPUT_H
#define QH5_INPUT_H

#include "input_interface.hpp"
#include "qh5.hpp"

class CQH5Input : public IDataReader
{
public:
  CQH5Input();
  ~CQH5Input();
  virtual void ReadImage(void **pix, std::string label, 
                         std::map<std::string, double> &params,
                         std::string &comment, 
                         std::vector<unsigned> &position);
  void SetFile(std::string filename);
};

#endif
