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
  m_imageReader=CDataReaderFactory::Get()->GetReader(input_extension);
  m_imageWriter=CDataWriterFactory::Get()->GetWriter(output_extension);
};

void CImageIO::CreateRealDataSet(const std::string &name, const std::vector<unsigned int> &positions)
{
  m_imageWriter->CreateRealDataSet(name, m_nx, m_ny, positions);
}

void CImageIO::CreateComplexDataSet(const std::string &name, const std::vector<unsigned int> &positions)
{
  m_imageWriter->CreateComplexDataSet(name, m_nx, m_ny, positions);
}

void CImageIO::WriteRealImage(const RealVector &pix, const std::string &fileName, std::map<std::string, double> &params,
                              const std::string &comment, const std::vector<unsigned> &position) {
  std::vector<unsigned> shape=GetShapeVector();
  m_imageWriter->WriteRealImage(pix, shape, fileName, position, comment, params);
}

void CImageIO::WriteComplexImage(const ComplexVector &pix, const std::string &fileName, std::map<std::string, double> &params,
                                 const std::string &comment, const std::vector<unsigned> &position) {
  std::vector<unsigned> shape=GetShapeVector();
  m_imageWriter->WriteComplexImage(pix, shape, fileName,
                                   position, comment, params);
}

void CImageIO::ReadImage(RealVector &pix, std::string &fileName, std::map<std::string, double> &params,
                         std::string &comment, std::vector<unsigned> position)
{
	unsigned byteSize;
	unsigned nx, ny;
	m_imageReader->ReadElementByteSize(fileName, position, byteSize);
	m_imageReader->ReadSize(fileName, position, nx, ny);
	// TODO: we should be able to handle reading in files of a different floating point precision at some point
	assert(byteSize == sizeof(float_tt));
	assert(pix.size() == nx*ny);
  m_imageReader->ReadImage(fileName, (void *)&pix[0], params, comment);
}

void CImageIO::ReadImage(ComplexVector &pix, std::string &fileName, std::map<std::string, double> &params,
                         std::string &comment, std::vector<unsigned> position)
{
	bool complex;
	unsigned byteSize;
	unsigned nx, ny;
  m_imageReader->ReadComplex(fileName, position, complex);
  m_imageReader->ReadElementByteSize(fileName, position, byteSize);
	m_imageReader->ReadSize(fileName, position, nx, ny);
  assert(complex);
  assert(byteSize == sizeof(float_tt));
	assert(pix.size() == nx*ny);
  m_imageReader->ReadImage(fileName, (void *)&pix[0], params, comment);
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
