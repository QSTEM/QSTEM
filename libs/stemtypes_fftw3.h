#ifndef STEMTYPES_H
#define STEMTYPES_H

#include "memory_fftw3.h"
#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/StdVector>

// #include "floatdef.h"
#include "fftw3.h"

//#include "boost/multi_array.h"

////////////////////////////////////////////////////////////////
#if FLOAT_PRECISION == 1
#define fftw_real float
#ifndef float_tt
#define float_tt  float
#endif
#else  // FLOAT_PRECISION
#define fftw_real double
#ifndef float_tt
#define float_tt  double
#endif
#endif  // FLOAT_PRECISION

using namespace Eigen;

// Only this section uses the Eigen namespace.
// By default, Eigen handles memory in Column-major format (the opposite of C/C++)
//  If any array accessing is done, we have to be very careful about that!
//  Translation was done in Notepad++ with this regex:
//  example: Mm[0,1] -> Mm(1,0)
// (\w+)\[([\w\d\+\-\*]+)\]\[([\w\d\+\-\*]+)\] -> \1\(\3,\2\)

// Forced float (32-bit) types (capital F or C)
//typedef Matrix< float, Dynamic, Dynamic> QSFMat;
//typedef Matrix< std::complex<float_tt>, Dynamic, Dynamic> QSCMat;
//typedef Matrix< float, Dynamic, 1> QSFVec; 

// Variable float (32 or 64-bit) types (lowercase f or c)
typedef Matrix< float_tt, Dynamic, Dynamic> QSfMat;
typedef Matrix< float_tt, 3, 3> QSf3Mat;
typedef Matrix< std::complex<float_tt>, Dynamic, Dynamic> QScMat;
typedef Matrix< int, Dynamic, Dynamic> QSiMat;
typedef Matrix< float_tt, Dynamic, 1> QSfVec;
typedef Matrix< float_tt, 3, 1> QSf3Vec;
typedef Array<float_tt, 3, 1> QSf3Arr;
typedef Matrix< int, Dynamic, 1> QSiVec;

typedef std::vector<QScMat, Eigen::aligned_allocator<QScMat>> QSVecOfcMat;
typedef std::vector<QSfMat, Eigen::aligned_allocator<QSfMat>> QSVecOffMat;

////////////////////////////////////////////////////////////////

typedef struct atomStruct {
  //float z,y,x;
  QSf3Arr pos; // x, y, z;
  // float dx,dy,dz;  // thermal displacements
  float dw;      // Debye-Waller factor
  float occ;     // occupancy
  float q;       // charge 
  int Znum;
} atom;

/* Planes will be defined by the standard equation for a plane, i.e.
 * a point (point) and 2 vectors (vect1, vect2)
 */
typedef struct planeStruct {
  QSf3Vec norm, vect1, vect2;
  QSf3Arr point;
} plane;

typedef struct grainBoxStruct {
  int amorphFlag;
  float_tt density,rmin, rFactor;  /* density, atomic distance, reduced atomic distance
				  * for amorphous material.  Red. r is for making a hex.
				  * closed packed structure, which will fill all space, 
				  * but will only be sparsely filled, and later relaxed.
				  */ 
  std::string name;
  // char *name;
  std::vector<atom> unitCell; /* definition of unit cell */
  int natoms;     /* number of atoms in unit cell */
  QSf3Arr cellDims;
  //float_tt ax,by,cz; /* unit cell parameters */
  float_tt alpha, beta, gamma; /* unit cell parameters */
  float_tt tiltx,tilty,tiltz;
  float_tt shiftx,shifty,shiftz;
  std::vector<plane> planes;
  // plane *planes;   /* pointer to array of bounding planes */
  float_tt sphereRadius, sphereX,sphereY,sphereZ; /* defines a sphere instead of a grain with straight edges */
  int nplanes; /* number of planes in array planes */
} grainBox;

typedef struct superCellBoxStruct {
  //float_tt cmx,cmy,cmz;  /* fractional center of mass coordinates */
  QSf3Arr cm; //x, y, z
  //float_tt ax,by,cz;
  QSf3Arr cellDims; //x, y, z
  int natoms;
  std::vector<atom> atoms;
  // atom *atoms; /* contains all the atoms within the super cell */
} superCellBox;

typedef struct atomBoxStruct {
  int used;   /* indicate here whether this atom is used in the
		 particular problem */
  int nx,ny,nz;
  float_tt dx,dy,dz;
  float_tt B_;
  
  QSVecOfcMat potential;
  QSVecOffMat rpotential;
} atomBox;

typedef struct detectorStruct {
  float_tt rInside,rOutside;
  float_tt k2Inside,k2Outside;
  std::string name;
  QSfMat image;
  QSfMat image2;
  float_tt error;  
  float_tt shiftX,shiftY;
  int Navg;
} DETECTOR;

#endif // STEMTYPES_H

