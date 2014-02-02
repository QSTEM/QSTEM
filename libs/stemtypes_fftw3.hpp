/*
QSTEM - image simulation for TEM/STEM/CBED
    Copyright (C) 2000-2010  Christoph Koch
	Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef STEMTYPES_H
#define STEMTYPES_H

#include <boost/shared_ptr.hpp>

#include "fftw_allocator.hpp"

#include "elTable.hpp"

////////////////////////////////////////////////////////////////////////
// define whether to use single or double precision
///////////////////////////////////////////////////////////////////////
#define FLOAT_PRECISION 1


#define BW (2.0F/3.0F)	/* bandwidth limit */

// Mode definitions
#define STEM    1
#define CBED    2
#define TEM     3
#define REFINE  4
#define MSCBED  5
#define TOMO    6

// Helpful mathematical constants
#define RAD2DEG 57.2958
#define SQRT_2 1.4142135
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

// Scattering factor types
#define DOYLE_TURNER 0
#define WEICK_KOHL 1
#define CUSTOM 2

////////////////////////////////////////////////////////////////////////
// Define physical constants
////////////////////////////////////////////////////////////////////////
#define ELECTRON_CHARGE (1.6021773e-19)
#define PICO_AMPERE (1e-12/ELECTRON_CHARGE)
#define MILLISEC_PICOAMP (1e-3*PICO_AMPERE)

// #include "floatdef.hpp"
#include "fftw3.h"
#include <cmath>

////////////////////////////////////////////////////////////////
#if FLOAT_PRECISION == 1
typedef FFTWComplex<float> complex_tt;
typedef float float_tt;
#else  // FLOAT_PRECISION
typedef FFTWComplex<double> complex_tt;
typedef double float_tt;
#endif  // FLOAT_PRECISION
typedef std::vector<float_tt, FFTWAllocator<float_tt> > RealVector;
typedef std::vector<complex_tt, FFTWAllocator<complex_tt > > ComplexVector;
////////////////////////////////////////////////////////////////

const float_tt PI = 2*acos(0.0);

// FFTW constants
const int k_fftMeasureFlag = FFTW_ESTIMATE;

typedef struct atomStruct {
  float_tt z,y,x;
  // float dx,dy,dz;  // thermal displacements
  float_tt dw;      // Debye-Waller factor
  float_tt occ;     // occupancy
  float_tt q;       // charge 
  unsigned Znum;
  float_tt mass;
  atomStruct();
  // constructor (used for testing)
  atomStruct(float_tt _mass, std::string _symbol, 
             float_tt _x, float_tt _y, float_tt _z, 
             float_tt _dw, float_tt _occ, float_tt _charge);
} atom;

inline atomStruct::atomStruct()
  : x(0)
  , y(0)
  , z(0)
  , dw(0)
  , occ(1)
  , q(0)
  , Znum(0)
  , mass(0)
{
}

inline atomStruct::atomStruct(float_tt _mass, std::string _symbol, 
                              float_tt _x, float_tt _y, float_tt _z, 
                              float_tt _dw, float_tt _occ, float_tt _charge)
  : x(_x)
  , y(_y)
  , z(_z)
  , dw(_dw)
  , occ(_occ)
  , q(_charge)
  , mass(_mass)
{
  Znum = getZNumber(_symbol.c_str());
}

typedef boost::shared_ptr<atom> atomPtr;

/* Planes will be defined by the standard equation for a plane, i.e.
 * a point (point) and 2 vectors (vect1, vect2)
 */
typedef struct planeStruct {
  float_tt normX,normY,normZ;
  float_tt vect1X,vect1Y,vect1Z;
  float_tt vect2X,vect2Y,vect2Z;
  float_tt pointX,pointY,pointZ;
} plane;

typedef boost::shared_ptr<plane> planePtr;

typedef struct grainBoxStruct {
  bool amorphFlag;
  float_tt density,rmin, rFactor;  /* density, atomic distance, reduced atomic distance
				  * for amorphous material.  Red. r is for making a hex.
				  * closed packed structure, which will fill all space, 
				  * but will only be sparsely filled, and later relaxed.
				  */ 
  char *name;
  atom *unitCell; /* definition of unit cell */
  unsigned natoms;     /* number of atoms in unit cell */
  float_tt ax,by,cz; /* unit cell parameters */
  float_tt alpha, beta, gamma; /* unit cell parameters */
  float_tt tiltx,tilty,tiltz;
  float_tt shiftx,shifty,shiftz;
  plane *planes;   /* pointer to array of bounding planes */
  float_tt sphereRadius, sphereX,sphereY,sphereZ; /* defines a sphere instead of a grain with straight edges */
  int nplanes; /* number of planes in array planes */
} grainBox;

typedef boost::shared_ptr<grainBox> grainBoxPtr;

typedef struct superCellBoxStruct {
  float_tt cmx,cmy,cmz;  /* fractional center of mass coordinates */
  float_tt ax,by,cz;
  unsigned natoms;
  atom *atoms; /* contains all the atoms within the super cell */
} superCellBox;

typedef boost::shared_ptr<superCellBox> superCellBoxPtr;

typedef struct atomBoxStruct {
  bool used;   /* indicate here whether this atom is used in the
		 particular problem */
  unsigned nx,ny,nz;
  float_tt dx,dy,dz;
  float_tt B;
  complex_tt ***potential;   /* 3D array containg 1st quadrant of real space potential */
  float_tt ***rpotential;   /* 3D array containg 1st quadrant of real space potential */
} atomBox;

typedef boost::shared_ptr<atomBox> atomBoxPtr;


// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
  #define QSTEM_HELPER_DLL_IMPORT __declspec(dllimport)
  #define QSTEM_HELPER_DLL_EXPORT __declspec(dllexport)
  #define QSTEM_HELPER_DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define QSTEM_HELPER_DLL_IMPORT __attribute__ ((visibility ("default")))
    #define QSTEM_HELPER_DLL_EXPORT __attribute__ ((visibility ("default")))
    #define QSTEM_HELPER_DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define QSTEM_HELPER_DLL_IMPORT
    #define QSTEM_HELPER_DLL_EXPORT
    #define QSTEM_HELPER_DLL_LOCAL
  #endif
#endif

#endif // STEMTYPES_H
