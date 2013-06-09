#include <stdio.h>	/*  ANSI-C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

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
  m_params(std::vector<float_tt>()),
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

CImageIO::CImageIO(int nx, int ny, float_tt t, float_tt dx, float_tt dy,
			  int paramSize, std::vector<float_tt> params, std::string comment) :
m_headerSize(56),
m_params(params),
m_paramSize(paramSize),
m_nx(nx),
m_ny(ny),
m_version(VERSION),
m_t(t),
m_dx(dx),
m_dy(dy),
m_comment("")
{
};

void CImageIO::WriteComplexImage(const void **pix, const char *fileName) {
  m_dataSize = 2*sizeof(float_tt);
  m_complexFlag = 1;
  
  WriteData(pix, fileName);
}

void CImageIO::WriteRealImage(const void **pix, const char *fileName) {
  m_dataSize = sizeof(float_tt);
  m_complexFlag = 0;
  
  WriteData(pix, fileName);
}

void CImageIO::WriteData(void **pix, char *fileName)
{
  FILE *fp;

  // Sychronize lengths of comments and parameters
  m_paramSize = m_params.size();
  m_commentSize = m_comment.size();

  if ((fp = fopen(fileName,"wb"))==NULL) {
    printf("writeRealImage: Could open file %s for writing\n",fileName);
    exit(0);
  }

  // TODO: should we write each element individually for clarity?
  fwrite((void *)this,m_headerSize,1,fp);
  fwrite((void *)(&m_params[0]), sizeof(double), m_paramSize, fp);
  fwrite((void *)(m_comment.c_str()), 1, m_commentSize, fp);
  if (fwrite(pix.data(),m_dataSize,(size_t)(m_nx*m_ny),fp) != m_nx*m_ny) {
    printf("writeRealImage: Error while writing data to file %s\n",fileName);
    fclose(fp);
    exit(0);
  }
  fclose(fp);
}

void CImageIO::ReadHeader(const char *fileName)
{
  FILE *fp;
  if ((fp = fopen(fileName,"rb"))==NULL)
  {
      printf("Could not open file %s for reading header.\n",fileName);
  }
  fread((void*)this, 1, 56, fp);
  if (m_paramSize>0)
    {
      m_params=std::vector<double>(m_paramSize)
      for (int i=0; i<m_paramSize; i++)
        {
          fread((void *)&m_params[0],sizeof(double),m_paramSize,fp);
        }
    }
  if (this->commentSize>0)
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
        throw std::exception(m_buf);
      }
      
      // Seek to the location of the actual data
      fseek(fp,56+(m_commentSize)+(2*m_paramSize),SEEK_SET);
      
      // this is type-agnostic - the type interpretation is done by the
      //   function sending in the pointer.  It casts it as void for the reading,
      //   but it then "knows" that it is double, complex, whatever, based on the
      //   type of the data that it passed into this function.
      //   
      //   Complex data is determined/communicated by the m_complexFlag, which is read in the header.
      nRead = fread((void *)pix, sizeof(m_dataSize),(size_t)(nx*ny),fp);
      if (nRead != nx*ny) 
      {
        freadError = 1;
        sprintf(m_buf, "Error while reading data from file %s:"
                " %d (of %d specified) elements read\n"
                "EOF: %d, Ferror: %d, dataSize: %d\n",
                fileName,nRead,nx*ny,feof(fp),ferror(fp),m_dataSize);
        fclose(fp);
        fp = NULL;
        throw std::exception(m_buf);
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

void SetParameters(std::vector<double> params)
{
  m_params=params;
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

