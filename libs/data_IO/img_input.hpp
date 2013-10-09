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

#include "input_interface.hpp"

class CImgInput : public IDataInput 
{
public:
  CImgInput();
  virtual void ReadImage(void **pix, std::string label, std::map<std::string, double> &params,
                         std::string &comment, std::vector<unsigned> position=std::vector<unsigned>());
private:
  int m_headerSize;  // first byte of image will be size of image header (in bytes)
                   // This is the size without the data, parameters, and comment!!!
  int m_paramSize;   // number of additional parameters
  int m_commentSize; // length of comment string
  int m_nx,m_ny;
  int m_complexFlag;
  int m_dataSize;    // size of one data element in bytes (e.g. complex double: 16)
  int m_version;     // The version flag will later help to find out how to 
                   // distinguish between images produced by different versions of stem
  double m_t;        // thickness
  double m_dx,m_dy;    // size of one pixel
  std::vector<double> m_params;  // array for additional parameters
  std::string m_comment;   // comment of prev. specified length
  char m_buf[200];  // General purpose temporary text buffer
  virtual void ReadHeader(const char *fileName);
};
