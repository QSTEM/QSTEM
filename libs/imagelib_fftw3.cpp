#include <stdio.h>	/*  ANSI-C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <fstream>

#include <stdexcept>

//#include "boost/shared_ptr.hpp"

#include "stemtypes_fftw3.h"
#include "imagelib_fftw3.h"
#include "memory_fftw3.h"	/* memory allocation routines */

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

CImageIO::CImageIO(int nx, int ny) :
  m_headerSize(56),
  m_params(std::vector<double>()),
  m_paramSize(0),
  m_nx(nx),
  m_ny(ny),
  m_version(VERSION),
  m_t(0.0),
  m_dx(1.0),
  m_dy(1.0),
  m_comment("")
{
};

CImageIO::CImageIO(int nx, int ny, double t, double dx, double dy,
			  std::vector<double> params, std::string comment) :
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
};

void CImageIO::WriteComplexImage(const void *pix, const char *fileName) {
  m_dataSize = 2*sizeof(float_tt);
  m_complexFlag = 1;
  
  WriteData((const void *)pix, fileName);
}

void CImageIO::WriteRealImage(const void *pix, const char *fileName) {
  m_dataSize = sizeof(float_tt);
  m_complexFlag = 0;
  
  WriteData((const void *)pix, fileName);
}

void CImageIO::WriteData(const void *pix, const char *fileName)
{
	//FILE *fp;
	std::fstream file(fileName, std::ios::out|std::ios::binary);

	// Sychronize lengths of comments and parameters
	m_paramSize = m_params.size();
	m_commentSize = m_comment.size();

	if(!file.is_open()) 
	{
		sprintf(m_buf,"WriteData: Could open file %s for writing\n",fileName);
		throw std::runtime_error(m_buf);
	}

	// TODO: should we write each element individually for clarity?
	// write the 56-byte header
	file.write(reinterpret_cast<const char*>(this), 56);
	//fwrite((void *)this,m_headerSize,1,fp);
	/*
	fwrite((void *)&m_headerSize, 4, 1, fp);
	fwrite((void *)&m_paramSize, 4, 1, fp);
	fwrite((void *)&m_commentSize, 4, 1, fp);
	fwrite((void *)&m_nx, 4, 1, fp);
	fwrite((void *)&m_ny, 4, 1, fp);
	fwrite((void *)&m_complexFlag, 4, 1, fp);
	fwrite((void *)&m_dataSize, 4, 1, fp);
	fwrite((void *)&m_version, 4, 1, fp);
	fwrite((void *)&m_t, 8, 1, fp);
	fwrite((void *)&m_dx, 8, 1, fp);
	 fwrite((void *)&m_dy, 8, 1, fp);
	*/
	if (m_paramSize>0)
	{
		//fwrite((void *)(&m_params[0]), sizeof(double), m_paramSize, fp);
		file.write(reinterpret_cast<const char*>(&m_params[0]), m_paramSize*sizeof(double));
	}
	//fwrite((void *)(m_comment.c_str()), 1, m_commentSize, fp);
	file.write(m_comment.c_str(), m_commentSize);
	file.write(reinterpret_cast<const char*>(pix), m_nx*m_ny*m_dataSize);
	//size_t pixWritten = fwrite(pix, m_dataSize,(size_t)(m_nx*m_ny),fp);
	//if (fwrite(pix,m_dataSize,(size_t)(m_nx*m_ny),fp) != m_nx*m_ny) {
	//if (pixWritten != m_nx*m_ny) {
	//	sprintf(m_buf,"writeRealImage: Error while writing data to file %s\n",fileName);
		//fclose(fp);
		//file.close();
		//throw std::runtime_error(m_buf);
	//}
	//fclose(fp);
	file.close();
}

void CImageIO::ReadHeader(const char *fileName)
{
  FILE *fp;
  if ((fp = fopen(fileName,"rb"))==NULL)
  {
      sprintf(m_buf, "Could not open file %s for reading header.\n",fileName);
	  throw std::runtime_error(m_buf);
  }
  fread((void*)this, 1, 56, fp);
  if (m_paramSize>0)
    {
      m_params=std::vector<double>(m_paramSize);
      //for (int i=0; i<m_paramSize; i++)
      //{
          fread((void *)&m_params[0],sizeof(double),m_paramSize,fp);
          //}
    }
  if (m_commentSize>0)
    {
      fread((void*)m_buf, 1, m_commentSize, fp);
      m_comment = std::string(m_buf);
    }

  if (fp != NULL) fclose(fp);
}

/***********************************************************
 * This function will read an image.  It reuses the same header 
 * struct over and over.  Therefore, values must be copied from 
 * the header members before calling this function again.
 *
 * The image pointer may also be NULL, in which case memory will be
 * allocated for it, and its size will be returned in the header struct
 * members nx, and ny.
 ***********************************************************/
void CImageIO::ReadImage(void *pix, int nx, int ny, const char *fileName) 
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
      nRead = fread(pix, sizeof(m_dataSize),(size_t)(nx*ny),fp);
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

void CImageIO::SetParams(std::vector<double> params)
{
  m_params=params;
}

void CImageIO::SetParameter(int index, double value)
{
	if (index < m_params.size())
		m_params[index] = value;
	else
		throw std::runtime_error("Tried to set out of bounds parameter.");
}

/* The pointers params and comment must either be initialized to NULL,
 * or point to a valid memory region already, because the will be 
 * REALLOCed
 */
/*
void CImageIO::ReadImageHeader(FILE *fp) {
  int hSize=sizeof(imageStruct);
  char buf[200];

  // reset the file pointer, just to make sure we are at the beginning
  clearerr(fp);
  fseek(fp,0,SEEK_SET);
  
  // first we will read in the fixed header of headerSize bytes
  fread((void *)&hSize,sizeof(int),1,fp);
  fseek(fp,0,SEEK_SET);
  fread((void *)header,1,hSize,fp);
  // read the additional parameters:
  if (m_paramSize >0) {
	  m_params = std::vector<float_tt>(m_paramSize);//(double *)realloc(header->params,header->paramSize*sizeof(double));
	  fread((void *)(&m_params[0]),sizeof(float_tt),m_paramSize,fp);
  }
  // read the comment
  if (m_commentSize>0) {
    fread((void *)(buf),1,m_commentSize,fp);
	m_comment = std::string(buf);
  }
  // printf("DataSize: %d, complex: %d\n",header->dataSize,header->complexFlag);
}
*/

