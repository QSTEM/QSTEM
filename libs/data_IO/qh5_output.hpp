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

#ifndef QH5_OUTPUT_H
#define QH5_OUTPUT_H

#include "output_interface.hpp"

class CQH5Output : public IDataWriter
{
  public:
  virtual void WriteRealVolume(float_tt *data, std::vector<unsigned> shape, std::string label, 
                               std::vector<unsigned> position=std::vector<unsigned>(), 
                               std::string comment=std::string(),
                               std::map<std::string, double> parameters=std::map<std::string, double>())=0;
  virtual void WriteComplexVolume(complex_tt *data, std::vector<unsigned> shape, std::string label, 
                                  std::vector<unsigned> position=std::vector<unsigned>(), 
                                  std::string comment=std::string(),
                                  std::map<std::string, double> parameters=std::map<std::string, double>())=0;
  virtual void WriteRealImage(float_tt **data, std::vector<unsigned> shape, std::string label, 
                              std::vector<unsigned> position=std::vector<unsigned>(), 
                              std::string comment=std::string(),
                              std::map<std::string, double> parameters=std::map<std::string, double>())=0;
  virtual void WriteComplexImage(complex_tt **data, std::vector<unsigned> shape, std::string label, 
                                 std::vector<unsigned> position=std::vector<unsigned>(), 
                                 std::string comment=std::string(),
                                 std::map<std::string, double> parameters=std::map<std::string, double>())=0;
};

