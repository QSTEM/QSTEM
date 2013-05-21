

#include <stdio.h>	/*  ANSI-C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

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


void writeImage(void **pix, imageStruct *header, const char *fileName) {
  FILE *fp;
  int nx,ny;
  // double rmin,rmax;
  // int flag=2,ix,iy;

  if ((fp = fopen(fileName,"wb"))==NULL) {
    printf("Could not write to file %s\n",fileName);
    exit(0);
  }
  header->version = VERSION;
  header->headerSize = sizeof(imageStruct)-sizeof(double *)-sizeof(char *);
  header->complexFlag = 1;
  if (header->complexFlag == 1) {
#if FLOAT_PRECISION == 1
    header->dataSize = sizeof(fftwf_complex);
#else
    header->dataSize = sizeof(fftw_complex);
#endif
  }
  else header->dataSize = sizeof(float_tt);
  nx = header->nx;  ny = header->ny;

  /*
  rmin = pix[0][0].re;
  rmax = rmin;
  for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
    if (rmin > pix[ix][iy].re) rmin = pix[ix][iy].re;
    if (rmax < pix[ix][iy].re) rmax = pix[ix][iy].re;    
  }
  printf("value: %g .. %g\n",rmin,rmax); 
  */

  fwrite((void *)header,header->headerSize,1,fp);
  fwrite((void *)(&header->params[0]),sizeof(float_tt),header->paramSize,fp);
  fwrite((void *)(header->comment.c_str()),1,header->commentSize,fp);
  
  if (fwrite(pix[0],header->dataSize,(size_t)(nx*ny),fp) != nx*ny) {
    printf("Error while writing %d x %d data to file %s\n",nx,ny,fileName);
    fclose(fp);
    exit(0);
  }
  fclose(fp);
}


void writeRealImage(QSfMat pix, imageStruct *header, const char *fileName, int dataSize) {
  FILE *fp;
  int nx,ny;
  // double rmin,rmax;
  // int flag=2,ix,iy;

  if ((fp = fopen(fileName,"wb"))==NULL) {
    printf("writeRealImage: Could open file %s for writing\n",fileName);
    exit(0);
  }
  header->version = VERSION;
  header->headerSize = sizeof(imageStruct)-sizeof(double *)-sizeof(char *);
  header->dataSize = dataSize;
  header->complexFlag = 0;
  nx = header->nx;  ny = header->ny;

  /*
  rmin = ((float **)pix)[0][0];
  rmax = rmin;
  for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
    if (rmin > ((float **)pix)[ix][iy]) rmin = ((float **)pix)[ix][iy];
    if (rmax < ((float **)pix)[ix][iy]) rmax = ((float **)pix)[ix][iy];    
  }
  printf("value: %g .. %g\n",rmin,rmax); 
  */  
  fwrite((void *)header,header->headerSize,1,fp);
  fwrite((void *)(&header->params[0]),sizeof(float_tt),header->paramSize,fp);
  fwrite((void *)(header->comment.c_str()),1,header->commentSize,fp);
  if (fwrite((void *)pix.data(),header->dataSize,(size_t)(nx*ny),fp) != nx*ny) {
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
imageStruct *readImage(QSfMat pix, int nx, int ny, const char *fileName) {
  FILE *fp;
  size_t nRead=0;
  int trial=0,maxTrial=3,freadError=0;
  static imageStruct *header = NULL;
  char error_string[300];

  if (header == NULL) header = makeNewHeader(1,1);
 
  do {
    if ((fp = fopen(fileName,"rb"))==NULL) {
      printf("Could not open file %s for reading\n",fileName);
      /* wait a short while */
      while (nRead < 1e5) nRead++;
    }
    else {
      getImageHeader(header,fp);


      if ((nx != header->nx)||(ny != header->ny)) {
		sprintf(error_string, "readImage: image size mismatch nx = %d (%d), ny = %d (%d)\n", header->nx,nx,header->ny,ny);
		throw std::exception(error_string);
      }
      nRead = fread((void *)pix.data(),sizeof(float_tt),(size_t)(nx*ny),fp);
      if (nRead != nx*ny) {
		freadError = 1;
		sprintf(error_string, "Error while reading data from file %s:"
	       " %d (of %d) elements read\n"
		   "EOF: %d, Ferror: %d, dataSize: %d\n",
		   fileName,nRead,(nx)*(ny),feof(fp),ferror(fp),header->dataSize);
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

imageStruct *makeNewHeader(int nx,int ny) {  
  imageStruct *header = NULL;
  
  
  header = (imageStruct *)malloc(sizeof(imageStruct));
  header->headerSize = sizeof(imageStruct)-sizeof(double *)-sizeof(char *);
  header->params = std::vector<float_tt>();
  header->paramSize = 0;
  header->nx = nx;
  header->ny = ny;
  header->version = VERSION;
  header->t = 0.0;
  header->dx = 1.0;
  header->dy = 1.0;
  header->complexFlag = 0;
  header->comment = "";
  setHeaderComment(header,"");
  return header;
}

imageStruct *makeNewHeaderCompact(int cFlag,int nx,int ny,float_tt t,float_tt dx,float_tt dy,
								  int paramSize, std::vector<float_tt> params, std::string comment) {  
  imageStruct *header;
  
  header = (imageStruct *)malloc(sizeof(imageStruct));
  header->headerSize = sizeof(imageStruct)-sizeof(double *)-sizeof(char *);
  header->params = params;
  header->paramSize = paramSize;
  header->nx = nx;
  header->ny = ny;
  header->version = VERSION;
  header->t = t;
  header->dx = dx;
  header->dy = dy;
  header->complexFlag = cFlag;

  header->comment = "";
  header->commentSize = 0;
  setHeaderComment(header,comment);	  
  return header;
}

void setHeaderComment(imageStruct *header, std::string comment) {

  if (!comment.empty()) {
    //header->comment = (char *)malloc(strlen(comment)+1);
	header->comment = comment;
    header->commentSize = (int)comment.length();
  }
  else {
    //header->comment = (char *)malloc(1);
    //*(header->comment) = '\0';  
    header->commentSize = 0;
  }
}

/* The pointers params and comment must either be initialized to NULL,
 * or point to a valid memory region already, because the will be 
 * REALLOCed
 */
void getImageHeader(imageStruct *header, FILE *fp) {
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
  if (header->paramSize >0) {
	  header->params = std::vector<float_tt>(header->paramSize);//(double *)realloc(header->params,header->paramSize*sizeof(double));
	  fread((void *)(&header->params[0]),sizeof(float_tt),header->paramSize,fp);
  }
  // read the comment
  if (header->commentSize>0) {
    fread((void *)(buf),1,header->commentSize,fp);
	header->comment = std::string(buf);
  }
  // printf("DataSize: %d, complex: %d\n",header->dataSize,header->complexFlag);
}


