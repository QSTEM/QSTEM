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

CImgWriter::CImgWriter() : CBinaryOutput(),
                           m_headerSize(56),
                           m_version(IMG_VERSION)
{
}

void CImgWriter::Initialize(std::string dirOrFileName, std::string run_id)
{
  // TODO: create the output folder if it doesn't exist
}

void CImgWriter::WriteComplexImage(const ComplexVector &data, const std::vector<unsigned> &shape, const std::string &label, 
                                  const std::vector<unsigned> &position, const std::string &comment, 
								  const std::map<std::string, double> &parameters)
{
  unsigned dataSize = 2*sizeof(float_tt);
  bool is_complex=true;
  WriteData((void *)&data[0], is_complex, dataSize, shape, label, position, comment, parameters);
}

void CImgWriter::WriteRealImage(const RealVector &data, const std::vector<unsigned> &shape, const std::string &label, 
                                  const std::vector<unsigned> &position, const std::string &comment, 
								  const std::map<std::string, double> &parameters)
{
  unsigned dataSize = sizeof(float_tt);
  bool is_complex=false;
  WriteData((void *)&data[0], is_complex, dataSize, shape, label, position, comment, parameters);
}

void CImgWriter::WriteData(void *pix, bool is_complex, unsigned dataSize, const std::vector<unsigned> &shape, 
                             const std::string &filebase, const std::vector<unsigned> &position, 
							 const std::string &comment, const std::map<std::string, double> &parameters)
{
  //FILE *fp;
  std::stringstream filename;
  char buf[200];
  filename<<filebase;
  
  for (unsigned idx=0; idx<position.size(); idx++)
    {
      filename<<"_"<<position[idx];
    }
  filename<<".img";

  std::fstream file(filename.str().c_str(), std::ios::out|std::ios::binary);

  // Sychronize lengths of comments and parameters
  size_t paramSize = parameters.size();
  std::map<std::string, double>::const_iterator param;

  // TODO: need to parse parameters to pull out thickness, complexFlag, and dataSize.
  //    Hardcode version?  Things aren't going to change anymore...

  size_t commentSize = comment.size();
  double thickness = parameters.at("Thickness");
  double resX = parameters.at("dx");
  double resY = parameters.at("dy");
  int complexFlag = (int) is_complex;
  std::string usedParams = "Thickness dx dy";

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

  // TODO: this will probably break the old code, since order of parameters matters, and map
  //     has its own order (need to fix Matlab code, or indicate a different file version?)
  for (param=parameters.begin(); param!=parameters.end(); param++)
    {
		// if this parameter isn't already recorded in the header, record it
		if (usedParams.find(param->first)==std::string::npos)
			file.write(reinterpret_cast<const char*>(&(param->second)), sizeof(double));
    }
  
  file.write(comment.c_str(), commentSize);
  file.write(reinterpret_cast<char*>(pix), shape[0]*shape[1]*dataSize);
  file.close();
}





