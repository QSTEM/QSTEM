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


#include "img_output.hpp"
#include "img_VERSION.hpp"

CImgOutput::CImgOutput() : CBinaryOutput(),
                           m_headerSize(56),
                           m_version(IMG_VERSION)
{
}

void CImgOutput::Initialize(std::string dirOrFileName, std::string run_id)
{
  // TODO: create the output folder if it doesn't exist
}

void CImgOutput::WriteComplexImage(complex_tt **data, const std::vector<unsigned> &shape, const std::string &label, 
                                 const std::vector<unsigned> &position, const std::string &comment,
                                 std::map<std::string, double> &parameters)
{
  unsigned dataSize = 2*sizeof(float_tt);
  WriteData((void **)data, true, dataSize, shape, label, position, comment, parameters);
}

void CImgOutput::WriteRealImage(float_tt **data, const std::vector<unsigned> &shape, const std::string &label, 
                                  const std::vector<unsigned> &position, const std::string &comment,
                                  std::map<std::string, double> &parameters)
{
  unsigned dataSize = sizeof(float_tt);
    WriteData((void **)data, false, dataSize, shape, label, position, comment, parameters);
}

void CImgOutput::WriteData(void **pix, bool is_complex, unsigned dataSize, std::vector<unsigned> shape, 
                             std::string label, std::vector<unsigned> position, std::string comment,
                             std::map<std::string, double> parameters)
{
  //FILE *fp;
  std::stringstream filename;
  char buf[200];
  filename<<label;
  for (unsigned idx=0; idx<position.size(); idx++)
    {
      filename<<"_"<<position[idx];
    }
  filename<<".img";

  std::fstream file(filename.str().c_str(), std::ios::out|std::ios::binary);

  // Sychronize lengths of comments and parameters
  size_t paramSize = parameters.size();
  std::map<std::string, double>::iterator param;
  
  // TODO: need to parse parameters to pull out thickness, complexFlag, and dataSize.
  //    Hardcode version?  Things aren't going to change anymore...

  size_t commentSize = comment.size();
  double thickness = parameters["Thickness"];
  double resX = parameters["dx"];
  double resY = parameters["dy"];
  int complexFlag = (int) is_complex;

  if(!file.is_open()) 
    {
      sprintf(buf,"WriteData: Could open file %s for writing\n",filename.str().c_str());
      throw std::runtime_error(buf);
    }

  file.write(reinterpret_cast <const char*> (&m_headerSize), 4);
  file.write(reinterpret_cast <const char*> (&paramSize), 4);
  file.write(reinterpret_cast <const char*> (&commentSize), 4);
  file.write(reinterpret_cast <const char*> (&shape[0]), 4);
  file.write(reinterpret_cast <const char*> (&shape[1]), 4);
  file.write(reinterpret_cast <const char*> (&complexFlag), 4);
  file.write(reinterpret_cast <const char*> (&dataSize), 4);
  file.write(reinterpret_cast <const char*> (&m_version), 4);
  file.write(reinterpret_cast <const char*> (&thickness), 8);
  file.write(reinterpret_cast <const char*> (&resX), 8);
  file.write(reinterpret_cast <const char*> (&resY), 8);

  for (param=parameters.begin(); param!=parameters.end(); param++)
    {
      file.write(reinterpret_cast<const char*>(&(param->second)), sizeof(double));
    }
  
  file.write(comment.c_str(), commentSize);
  file.write(reinterpret_cast<char*>(pix[0]), shape[0]*shape[1]*dataSize);
  file.close();
}





