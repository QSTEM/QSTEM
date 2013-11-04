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


#ifndef OUTPUT_INTERFACE_H
#define OUTPUT_INTERFACE_H

#include <string>
#include <vector>
#include <map>
#include "stemtypes_fftw3.hpp"
#include "boost/shared_ptr.hpp"

class IDataWriter
{
public:
  virtual void Initialize(std::string dirOrFileName, std::string run_id="avg")=0;
  // io plugins don't have to do anything to create a data set, but they might want to (qh5 uses this
  //   to create data sets which are then filled later, rather than saving individual files.)
  virtual void CreateRealDataSet(std::string name, unsigned nx, unsigned ny, std::vector<unsigned> &positions){};
  virtual void CreateComplexDataSet(std::string name, unsigned nx, unsigned ny, std::vector<unsigned> &positions){};
  virtual void WriteRealImage(float_tt **data, std::vector<unsigned> shape, std::string label, 
                              std::vector<unsigned> position=std::vector<unsigned>(), 
                              std::string comment=std::string(),
                              std::map<std::string, double> parameters=std::map<std::string, double>())=0;
  virtual void WriteComplexImage(complex_tt **data, std::vector<unsigned> shape, std::string label, 
                                 std::vector<unsigned> position=std::vector<unsigned>(), 
                                 std::string comment=std::string(),
                                 std::map<std::string, double> parameters=std::map<std::string, double>())=0;
};

typedef boost::shared_ptr<IDataWriter> DataWriterPtr;

#endif










