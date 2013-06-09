#ifndef IMAGELIB_H
#define IMAGELIB_H

#include "stemtypes_fftw3.h"

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
  float_tt m_t;        // thickness
  float_tt m_dx,m_dy;    // size of one pixel
  std::vector<float_tt> m_params;  // array for additional parameters
  std::string m_comment;   // comment of prev. specified length
public:
	CImageIO(int nx, int ny);
	CImageIO(int nx, int ny, float_tt t, float_tt dx, float_tt dy,
			  int paramSize, std::vector<float_tt> params, std::string comment);

	// Following functions write .img files, which will be deprecated in favor of HDF5 files.
	void WriteRealImage(const QSfMat pix, const char *fileName);
	void WriteComplexImage(const QScMat pix, const char *fileName);
	void ReadRealImage(QSfMat &pix, const char *fileName);
	void ReadComplexImage(QScMat &pix, const char *fileName);

	void WriteImage( std::string fileName);

	void SetComment(std::string comment);

	// Following functions are placeholders for reading/writing done with HDF5.
	template <typename Derived>
	void WriteImage(const EigenBase<Derived> m, std::string name) {
		
	}

	template <typename Derived>
	void ReadImage(EigenBase<Derived>& m, std::string name) {
		
	}
	
};

#endif
