#include "../stemtypes_fftw3.hpp"
#include "../imagelib_fftw3.hpp"
#include <map>

#ifndef POTENTIAL_BASE_H
#define POTENTIAL_BASE_H

class CPotential
{
public:
	CPotential(unsigned nx, unsigned ny, unsigned nz, float_tt dx, float_tt dy, float_tt dz, float_tt atomRadius, float_tt v0);
	CPotential(std::string parameter_file);
	~CPotential();

	virtual void atomBoxLookUp(complex_tt &val, int Znum, float_tt x, float_tt y, float_tt z, float_tt B);
	virtual void make3DSlices(int nlayer,char *fileName,atom *center);
protected:
  void Initialize();


  ImageIOPtr m_imageIO;
  complex_tt ***m_trans;
  unsigned m_nx, m_ny;
  // resolutions
  float_tt m_dx, m_dy, m_dz;
  // oversampled resolutions
  float_tt m_ddx, m_ddy, m_ddz;
  //
  int m_boxNx, m_boxNy, m_boxNz;
  // Atom radius
  float_tt m_radius, m_radius2;
  // voltage
  float_tt m_v0;
  std::map<unsigned, atomBoxPtr> m_atomBoxes;

  int m_printLevel;
};

typedef boost::shared_ptr<CPotential> PotPtr;

#endif