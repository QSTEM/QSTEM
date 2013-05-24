

#include <stdio.h>	/*  ANSI-C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "boost/shared_ptr.hpp"

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
  m_headerSize(sizeof(imageStruct)-sizeof(float_tt *)-sizeof(char *));
  m_params(std::vector<float_tt>()),
  m_paramSize(0),
  m_nx(nx),
  m_ny(ny),
  m_version(VERSION),
  m_t(0.0),
  m_dx(1.0),
  m_dy(1.0),
  m_complexFlag(0),
  m_comment("")
{
}

CImageIO::CImageIO(int nx, int ny, float_tt t, float_tt dx, float_tt dy,
			  int paramSize, std::vector<float_tt> params, std::string comment) :
m_headerSize(sizeof(CImageIO)-sizeof(float_tt *)-sizeof(char *));
m_params(params);
m_paramSize(paramSize);
m_m_nx(nx);
m_ny(ny);
m_version(VERSION);
m_t(t);
m_dx(dx);
m_dy(dy);
m_comment("");
m_complexFlag(_cFlag);
{
}

void CImageIO::WriteComplexImage(const QScMat pix, const char *fileName) {
  FILE *fp;

  if ((fp = fopen(fileName,"wb"))==NULL) {
    printf("Could not write to file %s\n",fileName);
    exit(0);
  }
  header->headerSize = sizeof(imageStruct)-sizeof(double *)-sizeof(char *);
  m_complexFlag = 1;
  m_datasize = 2*sizeof(float_tt);

  // fwrite((void *)header,header->headerSize,1,fp);
  fwrite((void *)this,m_headerSize,1,fp);
  fwrite((void *)(&m_params[0]), sizeof(float_tt), m_paramSize, fp);
  fwrite((void *)(m_comment.c_str()), 1, m_commentSize, fp);
  
  if (fwrite(pix.data(),m_dataSize,(size_t)(pix.cols()*pix.rows()),fp) != pix.cols()*pix.rows()) {
    printf("Error while writing %d x %d data to file %s\n",pix.cols(),pix.row(),fileName);
    fclose(fp);
    exit(0);
  }
  fclose(fp);
}

void CImageIO::WriteRealImage(const QSfMat pix, const char *fileName) {
  FILE *fp;
  // double rmin,rmax;
  // int flag=2,ix,iy;

  if ((fp = fopen(fileName,"wb"))==NULL) {
    printf("writeRealImage: Could open file %s for writing\n",fileName);
    exit(0);
  }
  //header->headerSize = sizeof(imageStruct)-sizeof(float_tt *)-sizeof(char *);
  m_dataSize = sizeof(float_tt);
  m_complexFlag = 0;

  // fwrite((void *)header,header->headerSize,1,fp);
  fwrite((void *)this,m_headerSize,1,fp);
  fwrite((void *)(&m_params[0]), sizeof(float_tt), m_paramSize, fp);
  fwrite((void *)(m_comment.c_str()), 1, m_commentSize, fp);
  if (fwrite(pix.data(),m_dataSize,(size_t)(pix.cols()*pix.rows()),fp) != pix.cols()*pix.rows()) {
    printf("writeRealImage: Error while writing data to file %s\n",fileName);
    fclose(fp);
    exit(0);
  }
  fclose(fp);
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
void CImageIO::ReadRealImage(QSfMat &pix, const char *fileName) {
  FILE *fp;
  size_t nRead=0;
  int trial=0,maxTrial=3,freadError=0;
  char error_string[300];
 
  do {
    if ((fp = fopen(fileName,"rb"))==NULL) {
      printf("Could not open file %s for reading\n",fileName);
      /* wait a short while */
      while (nRead < 1e5) nRead++;
    }
    else {
	  // sets m_nx, m_ny, m_datasize, m_complex
      getImageHeader(header,fp);

	  // TODO: check if image to be read is complex - if so, they should be using the ReadComplexImage function instead.

      if ((m_nx != pix.cols())||(m_ny != pix.rows())) {
		sprintf(error_string, "readImage: image size mismatch nx = %d (%d), ny = %d (%d)\n", header->nx,nx,header->ny,ny);
		throw std::exception(error_string);
      }
      nRead = fread((void *)pix.data(),sizeof(float_tt),(size_t)(pix.cols()*pix.rows()),fp);
      if (nRead != pix.cols()*pix.rows()) {
		freadError = 1;
		sprintf(error_string, "Error while reading data from file %s:"
	       " %d (of %d) elements read\n"
		   "EOF: %d, Ferror: %d, dataSize: %d\n",
		   fileName,nRead,(pix.cols())*(pix.rows()),feof(fp),ferror(fp),m_dataSize);
		fclose(fp);
		fp = NULL;
		throw std::exception(error_string);
      }
    }
    /* we will try three times to read this file. */
  }  while ((freadError > 0) && (++trial < maxTrial));
  
  if (fp != NULL) fclose(fp);

  return header;
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
void CImageIO::ReadComplexImage(QScMat &pix, int nx, int ny, const char *fileName) {
  FILE *fp;
  size_t nRead=0;
  int trial=0,maxTrial=3,freadError=0;
  char error_string[300];
 
  do {
    if ((fp = fopen(fileName,"rb"))==NULL) {
      printf("Could not open file %s for reading\n",fileName);
      /* wait a short while */
      while (nRead < 1e5) nRead++;
    }
    else {
      getImageHeader(header,fp);

      if ((nx != m_nx)||(ny != m_ny)) {
		sprintf(error_string, "readImage: image size mismatch nx = %d (%d), ny = %d (%d)\n", m_nx,nx,m_ny,ny);
		throw std::exception(error_string);
      }
      nRead = fread((void *)pix.data(),2*sizeof(float_tt),(size_t)(nx*ny),fp);
      if (nRead != nx*ny) {
		freadError = 1;
		sprintf(error_string, "Error while reading data from file %s:"
	       " %d (of %d) elements read\n"
		   "EOF: %d, Ferror: %d, dataSize: %d\n",
		   fileName,nRead,(nx)*(ny),feof(fp),ferror(fp),m_dataSize);
		fclose(fp);
		fp = NULL;
		throw std::exception(error_string);
      }
    }
    /* we will try three times to read this file. */
  }  while ((freadError > 0) && (++trial < maxTrial));
  
  if (fp != NULL) fclose(fp);

  return header;
}

/*****************************************************************
 * Image header routines
 ****************************************************************/

void ImageIO::SetHeaderComment(std::string comment) {

  if (!comment.empty()) {
	m_comment = comment;
    m_commentSize = (int)comment.length();
  }
  else {
    m_commentSize = 0;
  }
}

/* The pointers params and comment must either be initialized to NULL,
 * or point to a valid memory region already, because the will be 
 * REALLOCed
 */
void ImageIO::ReadImageHeader(FILE *fp) {
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


