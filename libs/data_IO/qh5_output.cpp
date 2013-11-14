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

#include "qh5_output.hpp"

CQH5Output::CQH5Output(std::string fileName, std::string run_id)
{
  Initialize(fileName, run_id);
}

CQH5Output::~CQH5Output()
{}

void CQH5Output::Initialize(std::string fileName, std::string run_id)
{
  m_qh5 = QH5ptr(new QH5(fileName));
  m_qh5->SetRunID(run_id);
}

void CQH5Output::CreateRealDataSet(std::string name, unsigned size_x, unsigned size_y, 
                                   std::vector<unsigned int> &positions)
{
  m_qh5->CreateRealDataSet(name, size_x, size_y, positions);
}

void CQH5Output::CreateComplexDataSet(std::string name, unsigned size_x, unsigned size_y,
                                      std::vector<unsigned int> &positions)
{
  m_qh5->CreateComplexDataSet(name, size_x, size_y, positions);
}

void CQH5Output::WriteRealImage(float_tt **data, std::vector<unsigned int> &shape, std::string &label, 
                           std::vector<unsigned int> &position, std::string &comment, 
                           std::map<std::string, double> &parameters)
{  
  m_qh5->WriteRealDataSlab(*data, label, shape[0], shape[1], position, parameters);
}

void CQH5Output::WriteComplexImage(complex_tt **data, std::vector<unsigned> &shape, std::string &label,
                              std::vector<unsigned> &position, std::string &comment, 
                              std::map<std::string, double> &parameters)
{
  m_qh5->WriteComplexDataSlab(*data, label, shape[0], shape[1], position, parameters);
}







