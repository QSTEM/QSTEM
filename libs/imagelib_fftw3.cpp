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

#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#ifdef WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif

CImageIO::CImageIO(int nx, int ny, std::string input_extension, std::string output_extension) :
m_nx(nx),
m_ny(ny)
{
  m_imageReader=GetDataReader(input_extension);
  m_imageWriter=GetDataWriter(output_extension);
};

void CImageIO::CreateRealDataSet(std::string &name, std::vector<unsigned int> &positions)
{
  m_imageWriter->CreateRealDataSet(name, m_nx, m_ny, positions);
}

void CImageIO::CreateComplexDataSet(std::string &name, std::vector<unsigned int> &positions)
{
  m_imageWriter->CreateComplexDataSet(name, m_nx, m_ny, positions);
}

void CImageIO::WriteRealImage(void **pix, std::string &fileName, std::map<std::string, double> &params,
                              std::string &comment, std::vector<unsigned> &position) {
  std::vector<unsigned> shape=GetShapeVector();
  m_imageWriter->WriteRealImage((float_tt **)pix, shape, fileName, position, comment, params);
}

void CImageIO::WriteComplexImage(void **pix, std::string &fileName, std::map<std::string, double> &params,
                                 std::string &comment, std::vector<unsigned> &position) {
  std::vector<unsigned> shape=GetShapeVector();
  m_imageWriter->WriteComplexImage((complex_tt **)pix, shape, fileName,
                                   position, comment, params);
}

void CImageIO::ReadImage(void **pix, std::string &fileName, std::map<std::string, double> &params,
                         std::string &comment, std::vector<unsigned> position)
{
  
  m_imageReader->ReadImage(pix, fileName, params, comment, position);
}

/*****************************************************************
 * Image header routines
 ****************************************************************/

std::vector<unsigned> CImageIO::GetShapeVector()
{
  std::vector<unsigned> size(2);
  size[0]=m_nx;
  size[1]=m_ny;
  return size;
}
