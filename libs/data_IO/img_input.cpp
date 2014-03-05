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

// TODO: C++-ify this?
#include "stdio.h"

#include <stdexcept>
#include <sstream>
#include <iostream>
#include <fstream>

#include "img_input.hpp"
#include "img_VERSION.hpp"
#include "filename_utilities.hpp"

namespace QSTEM
{

CImgReader::CImgReader() :
m_headerSize(56),
m_version(IMG_VERSION)
{
}

void CImgReader::ReadHeader(const std::string &filebase, const std::vector<unsigned> &position)
{
  FILE *fp;
  std::vector<double> params;
  std::string filename=BuildFilenameString(filebase, position);
  if ((fp = fopen(filename.c_str(),"rb"))==NULL)
  {
    sprintf(m_buf, "Could not open file %s for reading header.\n",filename.c_str());
    throw std::runtime_error(m_buf);
  }
  fread((void*)&m_headerSize, 4, 1,fp);
  fread((void*)&m_paramSize, 4, 1, fp);
  fread((void*)&m_commentSize, 4, 1, fp);
  fread((void*)&m_nx, 4, 1, fp);
  fread((void*)&m_ny, 4, 1, fp);
  fread((void*)&m_complexFlag, 4, 1, fp);
  fread((void*)&m_dataSize, 4, 1, fp);
  fread((void*)&m_version, 4, 1, fp);
  fread((void*)&m_t, 8, 1, fp);  // thickness
  fread((void*)&m_dx, 8, 1, fp); // resolution in X
  fread((void*)&m_dy, 8, 1, fp); // resolution in Y

  if (m_paramSize>0)
    {
      m_params=std::vector<double>(m_paramSize);
      fread((void *)&m_params[0],sizeof(double),m_paramSize,fp);
    }
  if (m_commentSize>0)
    {
      fread((void*)m_buf, 1, m_commentSize, fp);
      m_comment = std::string(m_buf);
	  m_comment.resize(m_commentSize);
    }

  if (fp != NULL) fclose(fp);
}

void CImgReader::ReadParameters(const std::string &filebase, const std::vector<unsigned> &position,
                                  std::map<std::string, double> &params)
{
  ReadHeader(filebase, position);
  _ReadParameters(params);
}

void CImgReader::_ReadParameters(std::map<std::string, double> &params)
{
  std::stringstream paramname;
  // Get a clear map for the parameters
  params=std::map<std::string, double>();
  params["Thickness"]=m_t;
  params["dx"]=m_dx;
  params["dy"]=m_dy;

  // Loop over any additional parameters
  for (size_t param=0; param<m_paramSize; param++)
    {
      // clear any existing parameter name
      paramname.str(std::string());
      paramname<<"parameter "<<param;
      params[paramname.str().c_str()]=m_params[param];
    }
}

void CImgReader::ReadComment(const std::string &filebase, const std::vector<unsigned> &position,
                               std::string &comment)
{
  ReadHeader(filebase,position);
  _ReadComment(comment);
}

void CImgReader::_ReadComment(std::string &comment)
{
  comment=m_comment;
}

/** stores the image size, as read from the header. Reads the header first. */
void CImgReader::ReadSize(const std::string &filebase,  const std::vector<unsigned> &position,
                          unsigned int &nx, unsigned int &ny)
{
  ReadHeader(filebase,position);
  _ReadSize(nx, ny);
}

/** stores the image size, as read from the header. */
void CImgReader::_ReadSize(unsigned &nx, unsigned &ny)
{
  nx=m_nx;
  ny=m_ny;
}

void CImgReader::ReadComplex(const std::string &filebase, const std::vector<unsigned> &position, 
                             bool &complex)
{
  ReadHeader(filebase,position);
  _ReadComplex(complex);
}

void CImgReader::_ReadComplex(bool &complex)
{
  complex=(bool)m_complexFlag;
}

void CImgReader::ReadElementByteSize(const std::string &filebase, const std::vector<unsigned> &position, 
                                     unsigned &elementSizeBytes)
{
  ReadHeader(filebase,position);
  _ReadElementByteSize(elementSizeBytes);
}

void CImgReader::_ReadElementByteSize(unsigned &elementByteSize)
{
	elementByteSize=m_dataSize;
}

/**
Reads data from file into buffer given by pointer *pix.  Requires that you've already read the header using the ReadHeader function,
   so that the image size and complexFlag is set properly.
*/
void CImgReader::_ReadImageData(const std::string &filebase, const std::vector<unsigned> &position, void *pix)
{
  std::ifstream inputfile;

  std::string filename=BuildFilenameString(filebase, position);
  inputfile.open (filename, std::ios::in | std::ios::binary);

	if (inputfile.good()) 
    {
      // Seek to the location of the actual data
      //fseek(fp,56+(m_commentSize)+(2*m_paramSize),SEEK_SET);
	  inputfile.seekg(56+(m_commentSize)+(2*m_paramSize));
      
      // this is type-agnostic - the type interpretation is done by the
      //   function sending in the pointer.  It casts it as char for the reading,
      //   but it then "knows" that it is double, complex, whatever, based on the
      //   type of the data that it passed into this function.
      //   
      //   Complex data is determined/communicated by the m_complexFlag, which is read in the header.  
	  //      The data size is doubled for complex data.
	  inputfile.read(reinterpret_cast <char*> (pix), m_nx*m_ny*m_dataSize);
    }
	else {
      printf("Could not open file %s for reading\n",filename.c_str());
      /* wait a short while */
      //while (nRead < 1e5) nRead++;
    }
  
  inputfile.close();
}

std::string CImgReader::BuildFilenameString(const std::string &label, const std::vector<unsigned> &position)
{
  std::stringstream filename;
  filename<<AddPositionToFilename(label,position);
  filename<<".img";
  return filename.str();
}

void CImgReader::ReadImageData(const std::string &filebase, const std::vector<unsigned> &position, RealVector &data)
{
	ReadHeader(filebase, position);
	_ReadImageData(filebase, position, (void *)&data[0]);
}

void CImgReader::ReadImageData(const std::string &filebase, const std::vector<unsigned> &position, ComplexVector &data)
{
	ReadHeader(filebase, position);
	_ReadImageData(filebase, position, (void *)&data[0]);
}

void CImgReader::ReadImage(const std::string &filebase, const std::vector<unsigned> &position, RealVector &data,
				std::map<std::string, double> &params, std::string &comment)
{
	ReadHeader(filebase, position);
	_ReadImageData(filebase, position, (void *)&data[0]);
	_ReadComment(comment);
	_ReadParameters(params);
}

void CImgReader::ReadImage(const std::string &filebase, const std::vector<unsigned> &position, ComplexVector &data,
				std::map<std::string, double> &params, std::string &comment)
{
	ReadHeader(filebase, position);
	_ReadImageData(filebase, position, (void *)&data[0]);
	_ReadComment(comment);
	_ReadParameters(params);
}

} // end namespace QSTEM