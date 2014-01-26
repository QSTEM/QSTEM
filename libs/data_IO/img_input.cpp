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

#include "img_input.hpp"
#include "img_VERSION.hpp"

// TODO: C++-ify this?
#include "stdio.h"

#include <stdexcept>
#include <sstream>

CImgReader::CImgReader() :
m_headerSize(56),
m_version(IMG_VERSION)
{
}

void CImgReader::ReadHeader(const char *fileName)
{
  FILE *fp;
  std::vector<double> params;
  if ((fp = fopen(fileName,"rb"))==NULL)
  {
      sprintf(m_buf, "Could not open file %s for reading header.\n",fileName);
	  throw std::runtime_error(m_buf);
  }
  // This sets many of the m_ variables!
  fread((void*)this, 1, 56, fp);
  if (m_paramSize>0)
    {
      m_params=std::vector<double>(m_paramSize);
      fread((void *)&m_params[0],sizeof(double),m_paramSize,fp);
    }
  if (m_commentSize>0)
    {
      fread((void*)m_buf, 1, m_commentSize, fp);
      m_comment = std::string(m_buf);
    }

  if (fp != NULL) fclose(fp);
}

void CImgReader::ReadParameters(const std::string &filename, std::map<std::string, double> &params)
{
  ReadHeader(BuildFilenameString(filename).c_str());
  _ReadParameters(params);
}

void CImgReader::ReadComment(const std::string &filename, std::string &comment)
{
  ReadHeader(BuildFilenameString(filename).c_str());
  _ReadComment(comment);
}

void CImgReader::_ReadComment(std::string &comment)
{
  comment=m_comment;
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

/** stores the image size, as read from the header. Reads the header first. */
void CImgReader::ReadSize(const std::string &filename, unsigned int &nx, unsigned int &ny)
{
  ReadHeader(BuildFilenameString(filename).c_str());
  _ReadSize(nx, ny);
}

/** stores the image size, as read from the header. */
void CImgReader::_ReadSize(unsigned &nx, unsigned &ny)
{
  nx=m_nx;
  ny=m_ny;
}

void CImgReader::ReadImageData(const std::string &filebase, void *pix)
{
  std::string filename=BuildFilenameString(filebase);
  ReadHeader(filename.c_str());
  _ReadImageData(filename, pix);
}

void CImgReader::_ReadImageData(const std::string &filename, void *pix)
{
  FILE *fp;
  size_t nRead=0;
  int trial=0,maxTrial=3,freadError=0;

  do {
    if ((fp = fopen(filename.c_str(),"rb"))==NULL) {
      printf("Could not open file %s for reading\n",filename.c_str());
      /* wait a short while */
      while (nRead < 1e5) nRead++;
    }
    else 
    {
      // Seek to the location of the actual data
      fseek(fp,56+(m_commentSize)+(2*m_paramSize),SEEK_SET);
      
      // this is type-agnostic - the type interpretation is done by the
      //   function sending in the pointer.  It casts it as void for the reading,
      //   but it then "knows" that it is double, complex, whatever, based on the
      //   type of the data that it passed into this function.
      //   
      //   Complex data is determined/communicated by the m_complexFlag, which is read in the header.
      nRead = fread(pix, sizeof(m_dataSize),(size_t)(m_nx*m_ny),fp);
      if (nRead != m_nx*m_ny) 
      {
        freadError = 1;
        sprintf(m_buf, "Error while reading data from file %s:"
                " %d (of %d specified) elements read\n"
                "EOF: %d, Ferror: %d, dataSize: %d\n",
                filename.c_str(),nRead,m_nx*m_ny,
                feof(fp),ferror(fp),m_dataSize);
        fclose(fp);
        fp = NULL;
        throw std::runtime_error(std::string(m_buf));
      }
    }
    /* we will try three times to read this file. */
  }  while ((freadError > 0) && (++trial < maxTrial));
  
  if (fp != NULL) fclose(fp);
}

std::string CImgReader::BuildFilenameString(const std::string &label)
{
  std::stringstream filename, paramname;
  filename<<label;
  /*
  for (unsigned idx=0; idx<position.size(); idx++)
    {
      filename<<"_"<<idx;
    }
  */
  filename<<".img";
  return filename.str();
}

void CImgReader::ReadImage(const std::string &filebase, void *pix, std::map<std::string, double> &params,
                         std::string &comment)
{
  // sets the important info from the header - most importantly, where to start reading the image.
  //   This sets the member variables that are used in the read below - that's why you don't have to
  //   specify them.
  std::string filename=BuildFilenameString(filebase);
  ReadHeader(filename.c_str());

  // These essentially dump information read from the header into the desired variables
  _ReadComment(comment);
  _ReadParameters(params);
  _ReadImageData(filename, pix);

  
}
