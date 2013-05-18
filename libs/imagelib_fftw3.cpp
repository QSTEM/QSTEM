

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


void writeImage(void **pix, imageStruct *header, char *fileName) {
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
  fwrite((void *)(header->params),sizeof(double),header->paramSize,fp);
  fwrite((void *)(header->comment),1,header->commentSize,fp);
  
  if (fwrite(pix[0],header->dataSize,(size_t)(nx*ny),fp) != nx*ny) {
    printf("Error while writing %d x %d data to file %s\n",nx,ny,fileName);
    fclose(fp);
    exit(0);
  }
  fclose(fp);
}


void writeRealImage(void **pix, imageStruct *header, char *fileName, int dataSize) {
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
  fwrite((void *)(header->params),sizeof(double),header->paramSize,fp);
  fwrite((void *)(header->comment),1,header->commentSize,fp);
  if (fwrite((void *)pix[0],header->dataSize,(size_t)(nx*ny),fp) != nx*ny) {
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
imageStruct *readImage(QSfMat pix,int nx,int ny,char *fileName) {
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
  header->params = NULL;
  header->paramSize = 0;
  header->nx = nx;
  header->ny = ny;
  header->version = VERSION;
  header->t = 0.0;
  header->dx = 1.0;
  header->dy = 1.0;
  header->complexFlag = 0;
  header->comment = NULL;
  setHeaderComment(header,NULL);
  return header;
}

imageStruct *makeNewHeaderCompact(int cFlag,int nx,int ny,double t,double dx,double dy,
				  int paramSize, double *params,char *comment) {  
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

  header->comment = NULL;
  header->commentSize = 0;
  setHeaderComment(header,comment);	  
  return header;
}

void setHeaderComment(imageStruct *header, char *comment) {

  if (header->comment != NULL) {
    free(header->comment);
    // header->comment = NULL;
  }
  if (comment != NULL) {
    header->comment = (char *)malloc(strlen(comment)+1);
    strcpy(header->comment,comment);
    header->commentSize = (int)strlen(comment);
  }
  else {
    header->comment = (char *)malloc(1);
    *(header->comment) = '\0';  
    header->commentSize = 0;
  }
}

/* The pointers params and comment must either be initialized to NULL,
 * or point to a valid memory region already, because the will be 
 * REALLOCed
 */
void getImageHeader(imageStruct *header, FILE *fp) {
  int hSize=sizeof(imageStruct);

  // reset the file pointer, just to make sure we are at the beginning
  clearerr(fp);
  fseek(fp,0,SEEK_SET);
  
  // first we will read in the fixed header of headerSize bytes
  fread((void *)&hSize,sizeof(int),1,fp);
  fseek(fp,0,SEEK_SET);
  fread((void *)header,1,hSize,fp);
  // read the additional parameters:
  if (header->paramSize >0) {
    header->params = (double *)realloc(header->params,header->paramSize*sizeof(double));
    fread((void *)(header->params),sizeof(double),header->paramSize,fp);
  }
  // read the comment
  if (header->commentSize>0) {
    header->comment = (char *)realloc(header->comment,(header->commentSize)+1);
    header->comment[header->commentSize] = '\0';
    fread((void *)(header->comment),1,header->commentSize,fp);  	
  }
  // printf("DataSize: %d, complex: %d\n",header->dataSize,header->complexFlag);
}



/*****************************************************************
 * Old Image routines:
 ****************************************************************/


void writeImage_old(fftw_complex **pix, int nx, int ny,float t,char *fileName) {
  FILE *fp;
  int size[2],flag=2,ix,iy;
  double rmin,rmax;

  size[0] = nx;
  size[1] = ny;
  if ((fp = fopen(fileName,"w"))==NULL) {
    printf("Could not write file %s\n",fileName);
    exit(0);
  }

  rmin = pix[0][0][0];
  rmax = rmin;
  for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
    if (rmin > pix[ix][iy][0]) rmin = pix[ix][iy][0];
    if (rmax < pix[ix][iy][0]) rmax = pix[ix][iy][0];    
  }
  //  printf("value: %g .. %g\n",rmin,rmax);

  fwrite((void *)size,sizeof(int),2,fp);
  fwrite((void *)(&t),sizeof(float),1,fp);  
  fwrite((void *)(&flag),sizeof(int),1,fp);
  fseek(fp,2*sizeof(fftw_complex),SEEK_SET);
  if (fwrite((void *)pix[0],sizeof(fftw_complex), 
	     (size_t)(nx*ny),fp) != nx*ny)
    {
      printf("Error while writing data to file %s\n",fileName);
      fclose(fp);
      exit(0);
    }
  fclose(fp);
}


void writeRealImage_old(fftw_real **pix, int nx, int ny,float t,char *fileName) {
  FILE *fp;
  int size[2],flag=0;

  size[0] = nx;
  size[1] = ny;
  if ((fp = fopen(fileName,"w"))==NULL) {
    printf("Could not write file %s\n",fileName);
    return;
  }

  fwrite((void *)size,sizeof(int),2,fp);
  fwrite((void *)(&t),sizeof(float),1,fp);  
  fwrite((void *)(&flag),sizeof(int),1,fp);  
  fseek(fp,4*sizeof(fftw_real),SEEK_SET);
  if (fwrite((void *)pix[0],sizeof(fftw_real), 
	     (size_t)(nx*ny),fp) != nx*ny)
    {
      printf("Error while writing data to file %s\n",fileName);
      fclose(fp);
      exit(0);
    }
  fclose(fp);
}

void readImage_old(fftw_complex **pix, int nx, int ny,float *t, char *fileName) {
  FILE *fp;
  size_t nRead=0;
  int size[2];
  int trial=0,maxTrial=3,freadError=0;

  do {
    if ((fp = fopen(fileName,"r"))==NULL) {
      printf("Could not open file %s for reading\n",fileName);
      // wait a short while 
      while (nRead < 1e5) nRead++;
    }
    else {
      clearerr(fp);
      
      fread((void *)size,sizeof(int),2,fp);
      fread((void *)t,sizeof(float),1,fp);
      fseek(fp,2*sizeof(fftw_complex),SEEK_SET);
      nRead = fread((void *)pix[0],sizeof(fftw_complex), (size_t)(nx*ny),fp);
      if (nRead != nx*ny) {
	freadError = 1;
	printf("Error while reading data from file %s:"
	       " %d (of %d) elements read\n",
	       fileName,nRead,(nx)*(ny));
	printf("EOF: %d, Ferror: %d\n",feof(fp),ferror(fp));
	fclose(fp);
      }
    }
    // we will try three times to read this file. 
  }  while ((freadError > 0) && (++trial < maxTrial));
  
  if ((nx != size[0]) || (ny != size[1])) {
    printf("Stored image %s has wrong size!\n",fileName);
    fclose(fp);
    exit(0);
  }
  fclose(fp);
}

void readRealImage_old(fftw_real **pix, int nx, int ny,float *t, char *fileName) {
  FILE *fp;
  size_t nRead,ix,iy;
  int size[2];

  if ((fp = fopen(fileName,"r"))==NULL) {
    printf("Could not open file %s for reading\n",fileName);
    for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) pix[ix][iy] = 0.0;
    return;
  }
  clearerr(fp);

  fread((void *)size,sizeof(int),2,fp);
  fread((void *)t,sizeof(float),1,fp);
  fseek(fp,2*sizeof(fftw_complex),SEEK_SET);
  nRead = fread((void *)pix[0],sizeof(fftw_real), (size_t)(nx*ny),fp);
  if (nRead != nx*ny) {
    printf("Error while reading data from file %s:"
	   " %d (of %d) elements read ",
	   fileName,nRead,(nx)*(ny));
    printf("(EOF: %d, Ferror: %d)\n",feof(fp),ferror(fp));
    fclose(fp);
    exit(0);
  }
  
  if ((nx != size[0]) || (ny != size[1])) {
    printf("Stored image %s has wrong size!\n",fileName);
    fclose(fp);
    exit(0);
  }
  fclose(fp);
}


