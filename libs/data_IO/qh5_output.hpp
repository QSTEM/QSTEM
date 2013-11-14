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

#include "../qh5.hpp"
#include "output_interface.hpp"

class CQH5Output : public IDataWriter
{
public:
  CQH5Output(std::string fileName, std::string run_id);
  ~CQH5Output();
  virtual void Initialize(std::string fileName, std::string run_id="avg");
  virtual void CreateRealDataSet(std::string name, unsigned nx, unsigned ny, std::vector<unsigned> &positions);
  virtual void CreateComplexDataSet(std::string name, unsigned nx, unsigned ny, std::vector<unsigned> &positions);
  virtual void WriteRealImage(float_tt **data, std::vector<unsigned> &shape, std::string label, 
                              std::vector<unsigned> &position, std::string &comment,
                              std::map<std::string, double> &parameters);
  virtual void WriteComplexImage(complex_tt **data, std::vector<unsigned> &shape, std::string label, 
                                 std::vector<unsigned> &position, std::string &comment,
                                 std::map<std::string, double> &parameters);
private:
  QH5ptr m_qh5;
};

#endif
