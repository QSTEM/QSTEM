#ifndef IMAGELIB_H
#define IMAGELIB_H

#include "stemtypes_fftw3.hpp"
#include <vector>
#include <map>
#include <string>
#include "boost/shared_ptr.hpp"
#include "data_IO/data_io_factories.hpp"

/**************************************************************
 * Here is how to use the new image writing routines
 *
 * ImageIOPtr imageio = ImageIOPtr(new CImageIO(nx,ny))
 *
 * imageIO->WriteRealImage((void**)real_image,filename);
 * imageIO->WriteComplexImage((void**)complex_image,filename);
 *
 **************************************************************
 * Reading an image works like this:
 * 
 * ImageIOPtr imageio = ImageIOPtr(new CImageIO(nx,ny))
 * imageio->ReadImage((void **)pix,nx,ny,fileName);
 *
 * Note that header parameters are persistent on any object.  You
 *   should use the various Set* functions to set parameters as
 *   necessary.  You should not need to read values from this class - 
 *   only set them.  They will be recorded to any file saved from this
 *   this object.
 **************************************************************/

class CImageIO {
  int m_nx,m_ny;
  char m_buf[200];  // General purpose temporary text buffer
  DataWriterPtr m_imageWriter;
  DataReaderPtr m_imageReader;
public:
  CImageIO(int nx, int ny, std::string input_extension=".img", 
           std::string output_extension=".img");

  void WriteRealImage(void **pix, const char *fileName, 
                      std::map<std::string, double> &params,
                      std::string comment, 
                      std::vector<unsigned> position=std::vector<unsigned>());
  // If you want to add a comment, but no parameters
  inline void WriteRealImage(void **pix, const char *fileName, std::string comment,
                             std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::map<std::string, double> params;
    WriteRealImage(pix, fileName, params, comment, position);
  }
  // If you want to add parameters, but no comment
  inline void WriteRealImage(void **pix, const char *fileName, std::map<std::string, double> &params,
                               std::vector<unsigned> position=std::vector<unsigned>())
  {
    WriteRealImage(pix, fileName, params, "", position);
  }  
  // If you don't care about parameters or comment, use this simplified overload:
  inline void WriteRealImage(void **pix, const char *fileName, 
                        std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::map<std::string, double> params;
    WriteRealImage(pix, fileName, params, "", position);
  }


  void WriteComplexImage(void **pix, const char *fileName,
                         std::map<std::string, double> &params,
                         std::string comment,
                         std::vector<unsigned> position=std::vector<unsigned>());
  
  // If you want to add a comment, but no parameters
  inline void WriteComplexImage(void **pix, const char *fileName, std::string comment,
                             std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::map<std::string, double> params;
    WriteComplexImage(pix, fileName, params, comment, position);
  }
  // If you want to add parameters, but no comment
  inline void WriteComplexImage(void **pix, const char *fileName, std::map<std::string, double> &params,
                               std::vector<unsigned> position=std::vector<unsigned>())
  {
    WriteComplexImage(pix, fileName, params, "", position);
  }  
  // If you don't care about parameters or comment, use this simplified overload:
  inline void WriteComplexImage(void **pix, const char *fileName, 
                        std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::map<std::string, double> params;
    WriteComplexImage(pix, fileName, params, "", position);
  }


  void ReadImage(void **pix, const char *fileName, std::map<std::string, double> &params,
                 std::string &comment,
                 std::vector<unsigned> position=std::vector<unsigned>());
  // If you don't care about parameters, this reads the data without you passing them in.
  inline void ReadImage(void **pix, const char *fileName, std::string &comment,
                        std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::map<std::string, double> params;
    ReadImage(pix, fileName, params, comment, position);
  }
  // If you want the parameters, but don't care about the comment
  inline void ReadImage(void **pix, const char *fileName, std::map<std::string, double> &params,
                        std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::string comment;
    ReadImage(pix, fileName, params, comment, position);
  }
  // If you don't care about comment or parameters
  inline void ReadImage(void **pix, const char *fileName,
                        std::vector<unsigned> position=std::vector<unsigned>())
  {
    std::map<std::string, double> params;
    std::string comment;
    ReadImage(pix, fileName, params, comment, position);
  }
  
private:
  std::vector<unsigned> GetShapeVector();
};

typedef boost::shared_ptr<CImageIO> ImageIOPtr;

#endif







