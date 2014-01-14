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

CImgInput::CImgInput() :
m_headerSize(56),
m_version(IMG_VERSION)
{
}

void CImgInput::ReadHeader(const char *fileName)
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

void CImgInput::ReadImage(void **pix, std::string label, std::map<std::string, double> &params,
                          std::string &comment, std::vector<unsigned> position)
{
  FILE *fp;
  size_t nRead=0;
  int trial=0,maxTrial=3,freadError=0;
 
  std::stringstream filename, paramname;
  filename<<label;
  for (unsigned idx=0; idx<position.size(); idx++)
    {
      filename<<"_"<<idx;
    }
  filename<<".img";

  // sets the important info from the header - most importantly, where to start reading the image.
  //   This sets the member variables that are used in the read below - that's why you don't have to
  //   specify them.
  ReadHeader(filename.str().c_str());


  comment=m_comment;

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
  
  do {
    if ((fp = fopen(filename.str().c_str(),"rb"))==NULL) {
      printf("Could not open file %s for reading\n",filename.str().c_str());
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
      nRead = fread(pix[0], sizeof(m_dataSize),(size_t)(m_nx*m_ny),fp);
      if (nRead != m_nx*m_ny) 
      {
        freadError = 1;
        sprintf(m_buf, "Error while reading data from file %s:"
                " %d (of %d specified) elements read\n"
                "EOF: %d, Ferror: %d, dataSize: %d\n",
                filename.str().c_str(),nRead,m_nx*m_ny,
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
