#ifndef IMAGELIB_H
#define IMAGELIB_H

#include "stemtypes_fftw3.h"
#include <vector>
#include <string>

/**************************************************************
 * Here is how to use the new image writing routines
 *
 * static imageStruct *header = NULL;
 *
 * if (header == NULL) header = makeNewHeaderCompact(cFlag,Nx,Ny,t,dx,dy,0,NULL,comment);
 * writeImage(cimage,header,filename);
 * or : writeRealImage(rimage,header,filename,sizeof(float));
 *
 **************************************************************
 * Reading an image works like this:
 * 
 * imageStruct *header;
 * header = readImage((void ***)(&pix),nx,ny,fileName);
 *
 * [This function will read an image.  It reuses the same header 
 * struct over and over.  Therefore, values must be copied from 
 * the header members before calling this function again.
 *
 * The image pointer may also be NULL, in which case memory will be
 * allocated for it, and its size will be returned in the header struct
 * members nx, and ny.]
 **************************************************************/

class CImageIO {
  int m_headerSize;  // first byte of image will be size of image header (in bytes)
                   // This is the size without the data and comment pointers!!!
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
public:
  CImageIO(int nx, int ny);
  CImageIO(int nx, int ny, double t, double dx, double dy,
           int paramSize, std::vector<double> params, std::string comment);

  void WriteRealImage(const void **pix, const char *fileName);
  void WriteComplexImage(const void **pix, const char *fileName);
  // reads in the header; returns the byte offset at which we should start reading image data.
  void ReadHeader(const char *fileName);
  void ReadImage(void **pix, int nx, int ny, const char *fileName);
  
  void WriteImage( std::string fileName);
        
  void SetComment(std::string comment);
  void SetThickness(double thickness);
  void SetParameters(std::vector<double> params);
private:
  WriteData(void **pix, char *fileName);
};

#endif
