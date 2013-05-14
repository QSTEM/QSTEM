#ifndef STEMTYPES_H
#define STEMTYPES_H

#include "memory_fftw3.h"
#include <vector>
#include <string>

////////////////////////////////////////////////////////////////////////
// define whether to use single or double precision
///////////////////////////////////////////////////////////////////////
#define FLOAT_PRECISION 1


#define BW (2.0F/3.0F)	/* bandwidth limit */
#define DOYLE_TURNER 0
#define WEICK_KOHL 1
#define CUSTOM 2
#define STEM    1
#define CBED    2
#define TEM     3
#define REFINE  4
#define MSCBED  5
#define TOMO    6

////////////////////////////////////////////////////////////////////////
// Define physical constants
////////////////////////////////////////////////////////////////////////
#define ELECTRON_CHARGE (1.6021773e-19)
#define PICO_AMPERE (1e-12/ELECTRON_CHARGE)
#define MILLISEC_PICOAMP (1e-3*PICO_AMPERE)

#include <Eigen/Core>
#include <Eigen/Dense>

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
typedef Matrix< int, Dynamic, 1> QSiVec;

////////////////////////////////////////////////////////////////

typedef struct atomStruct {
  float z,y,x;
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
  QSf3Vec norm, vect1, vect2, point;
} plane;

typedef struct grainBoxStruct {
  int amorphFlag;
  double density,rmin, rFactor;  /* density, atomic distance, reduced atomic distance
				  * for amorphous material.  Red. r is for making a hex.
				  * closed packed structure, which will fill all space, 
				  * but will only be sparsely filled, and later relaxed.
				  */ 
  std::string name;
  // char *name;
  std::vector<atom> unitCell; /* definition of unit cell */
  int natoms;     /* number of atoms in unit cell */
  float_tt ax,by,cz; /* unit cell parameters */
  float_tt alpha, beta, gamma; /* unit cell parameters */
  float_tt tiltx,tilty,tiltz;
  float_tt shiftx,shifty,shiftz;
  std::vector<plane> planes;
  // plane *planes;   /* pointer to array of bounding planes */
  float_tt sphereRadius, sphereX,sphereY,sphereZ; /* defines a sphere instead of a grain with straight edges */
  int nplanes; /* number of planes in array planes */
} grainBox;

typedef struct superCellBoxStruct {
  float_tt cmx,cmy,cmz;  /* fractional center of mass coordinates */
  float_tt ax,by,cz;
  int natoms;
  std::vector<atom> atoms;
  // atom *atoms; /* contains all the atoms within the super cell */
} superCellBox;

typedef struct atomBoxStruct {
  int used;   /* indicate here whether this atom is used in the
		 particular problem */
  int nx,ny,nz;
  float_tt dx,dy,dz;
  double B_;
  
  std::vector<QScMat> potential;
  std::vector<QSfMat> rpotential;
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


using namespace std;
vector<string> getElTable()
{
	vector<string> elTable(104);
	elTable.push_back("H");
	elTable.push_back("He");
	elTable.push_back("Li");
	elTable.push_back("Be");
	elTable.push_back("B");
	elTable.push_back("N");
	elTable.push_back("C");
	elTable.push_back("O");
	elTable.push_back("F");
	elTable.push_back("Ne");
	elTable.push_back("Na");
	elTable.push_back("Mg");
	elTable.push_back("Al");
	elTable.push_back("Si");
	elTable.push_back("P");
	elTable.push_back("S");
	elTable.push_back("Cl");
	elTable.push_back("Ar");
	elTable.push_back("K");
	elTable.push_back("Ar");
	elTable.push_back("K");
	elTable.push_back("Ca");
	elTable.push_back("Sc");
	elTable.push_back("Ti");
	elTable.push_back("V");
	elTable.push_back("Cr");
	elTable.push_back("Mn");
	elTable.push_back("Fe");
	elTable.push_back("Co");
	elTable.push_back("Ni");
	elTable.push_back("Cu");
	elTable.push_back("Zn");
	elTable.push_back("Ga");
	elTable.push_back("Ge");
	elTable.push_back("As");
	elTable.push_back("Se");
	elTable.push_back("Br");
	elTable.push_back("Kr");
	elTable.push_back("Rb");
	elTable.push_back("Sr");
	elTable.push_back("Y");
	elTable.push_back("Zr");
	elTable.push_back("Nb");
	elTable.push_back("Mo");
	elTable.push_back("Tc");
	elTable.push_back("Ru");
	elTable.push_back("Rh");
	elTable.push_back("Pd");
	elTable.push_back("Ag");
	elTable.push_back("Cd");
	elTable.push_back("In");
	elTable.push_back("Sn");
	elTable.push_back("Sb");
	elTable.push_back("Te");
	elTable.push_back("I");
	elTable.push_back("Xe");
	elTable.push_back("Cs");
	elTable.push_back("Ba");
	elTable.push_back("La");
	elTable.push_back("Ce");
	elTable.push_back("Pr");
	elTable.push_back("Nd");
	elTable.push_back("Pm");
	elTable.push_back("Sm");
	elTable.push_back("Eu");
	elTable.push_back("Gd");
	elTable.push_back("Tb");
	elTable.push_back("Dy");
	elTable.push_back("Ho");
	elTable.push_back("Er");
	elTable.push_back("Tm");
	elTable.push_back("Yb");
	elTable.push_back("Lu");
	elTable.push_back("Hf");
	elTable.push_back("Ta");
	elTable.push_back("W");
	elTable.push_back("Re");
	elTable.push_back("Os");
	elTable.push_back("Ir");
	elTable.push_back("Pt");
	elTable.push_back("Au");
	elTable.push_back("Hg");
	elTable.push_back("Tl");
	elTable.push_back("Pb");
	elTable.push_back("Bi");
	elTable.push_back("Po");
	elTable.push_back("At");
	elTable.push_back("Rn");
	elTable.push_back("Fr");
	elTable.push_back("Ra");
	elTable.push_back("Ac");
	elTable.push_back("Th");
	elTable.push_back("Pa");
	elTable.push_back("U");
	elTable.push_back("Np");
	elTable.push_back("Pu");
	elTable.push_back("Am");
	elTable.push_back("Cm");
	elTable.push_back("Bk");
	elTable.push_back("Cf");
	elTable.push_back("Es");
	elTable.push_back("Fm");
	elTable.push_back("Md");
	elTable.push_back("No");
	elTable.push_back("Lr");
	return elTable;
}

#endif // STEMTYPES_H

