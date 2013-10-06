#include <stdio.h>	/*  ANSI-C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <fstream>

#include <stdexcept>
#include <map>

//#include "boost/shared_ptr.hpp"

#include "stemtypes_fftw3.hpp"
#include "imagelib_fftw3.hpp"
#include "memory_fftw3.hpp"	/* memory allocation routines */

#include "file_IO/data_writers.hpp"

#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#ifdef WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif

#define VERSION 1  // please update this number, if anything in the image header 
                   // format changes.

CImageIO::CImageIO(int nx, int ny, std::string extension) :
  m_headerSize(56),
  m_params(std::map<std::string, double>()),
  m_paramSize(0),
  m_nx(nx),
  m_ny(ny),
  m_version(VERSION),
  m_t(0.0),
  m_dx(1.0),
  m_dy(1.0),
  m_comment("")
{
  m_imageWriter = GetDataWriter(extension);
};

CImageIO::CImageIO(int nx, int ny, double t, double dx, double dy,
                   std::map<std::string, double> params, std::string comment, std::string extension) :
m_headerSize(56),
m_params(params),
m_nx(nx),
m_ny(ny),
m_version(VERSION),
m_t(t),
m_dx(dx),
m_dy(dy),
m_comment("")
{
  m_imageWriter=GetDataWriter(extension);
};

void CImageIO::WriteComplexImage(void **pix, const char *fileName) {
  m_imageWriter->WriteComplexImage((complex_tt **)pix, GetShapeVector(), std::string(fileName),
                                   std::vector<ulong>(), m_comment, m_params, GetResolutionVector());
}

void CImageIO::WriteRealImage(void **pix, const char *fileName) {
  m_imageWriter->WriteRealImage((float_tt **)pix, GetShapeVector(), std::string(fileName),
                                std::vector<ulong>(), m_comment, m_params, GetResolutionVector());
}

void CImageIO::ReadHeader(const char *fileName)
{
  FILE *fp;
  std::vector<double> params;
  if ((fp = fopen(fileName,"rb"))==NULL)
  {
      sprintf(m_buf, "Could not open file %s for reading header.\n",fileName);
	  throw std::runtime_error(m_buf);
  }
  fread((void*)this, 1, 56, fp);
  if (m_paramSize>0)
    {
      params=std::vector<double>(m_paramSize);
      fread((void *)&params[0],sizeof(double),m_paramSize,fp);
    }
  if (m_commentSize>0)
    {
      fread((void*)m_buf, 1, m_commentSize, fp);
      m_comment = std::string(m_buf);
    }

  if (fp != NULL) fclose(fp);
}

void CImageIO::ReadImage(void **pix, int nx, int ny, const char *fileName) 
{
  FILE *fp;
  size_t nRead=0;
  int trial=0,maxTrial=3,freadError=0;
 
  // sets the important info from the header - most importantly, where to start reading the image.
  ReadHeader(fileName);

  do {
    if ((fp = fopen(fileName,"rb"))==NULL) {
      printf("Could not open file %s for reading\n",fileName);
      /* wait a short while */
      while (nRead < 1e5) nRead++;
    }
    else 
    {
      if ((m_nx != nx)||(m_ny != nx)) {
        sprintf(m_buf, "readImage: image size mismatch nx = %d (%d), ny = %d (%d)\n", m_nx,nx,m_ny,ny);
        throw std::runtime_error(std::string(m_buf));
      }
      
      // Seek to the location of the actual data
      fseek(fp,56+(m_commentSize)+(2*m_paramSize),SEEK_SET);
      
      // this is type-agnostic - the type interpretation is done by the
      //   function sending in the pointer.  It casts it as void for the reading,
      //   but it then "knows" that it is double, complex, whatever, based on the
      //   type of the data that it passed into this function.
      //   
      //   Complex data is determined/communicated by the m_complexFlag, which is read in the header.
      nRead = fread(pix[0], sizeof(m_dataSize),(size_t)(nx*ny),fp);
      if (nRead != nx*ny) 
      {
        freadError = 1;
        sprintf(m_buf, "Error while reading data from file %s:"
                " %d (of %d specified) elements read\n"
                "EOF: %d, Ferror: %d, dataSize: %d\n",
                fileName,nRead,nx*ny,feof(fp),ferror(fp),m_dataSize);
        fclose(fp);
        fp = NULL;
        throw std::runtime_error(std::string(m_buf));
      }
    }
    /* we will try three times to read this file. */
  }  while ((freadError > 0) && (++trial < maxTrial));
  
  if (fp != NULL) fclose(fp);
}

/*****************************************************************
 * Image header routines
 ****************************************************************/

void CImageIO::SetComment(std::string comment) 
{
  m_comment = comment;
}

void CImageIO::SetThickness(double thickness)
{
  m_t = thickness;
}

void CImageIO::SetResolution(double resX, double resY)
{
	m_dx = resX;
	m_dy = resY;
}

void CImageIO::SetParams(std::map<std::string, double> &params)
{
  m_params=params;
}

void CImageIO::SetParameter(std::string key, double value)
{
  m_params[key] = value;
}

/*
void CImageIO::SetParameter(int index, double value)
{
	if (index < m_params.size())
		m_params[index] = value;
	else
		throw std::runtime_error("Tried to set out of bounds parameter.");
}
*/

std::vector<ulong> CImageIO::GetShapeVector()
{
  std::vector<ulong> size(2, ulong());
  size[0]=m_nx;
  size[1]=m_ny;
  return size;
}

std::vector<float_tt> CImageIO::GetResolutionVector()
{
  std::vector<float_tt> size(2, float_tt());
  size[0]=m_dx;
  size[1]=m_dy;
  return size;
};


















