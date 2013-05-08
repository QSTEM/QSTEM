#ifndef STEMTYPES_H
#define STEMTYPES_H

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

// #include "floatdef.h"
#include "fftw3.h"

#include "boost/shared_ptr.hpp"
#include "boost/multi_array.hpp"

////////////////////////////////////////////////////////////////////////
// define whether to use single or double precision
///////////////////////////////////////////////////////////////////////
#define FLOAT_PRECISION 1

////////////////////////////////////////////////////////////////
#if FLOAT_PRECISION == 1
#define fftw_real float
#ifndef float_tt
#define float_tt  float
#endif
#define real      float
#else  // FLOAT_PRECISION
#define fftw_real double
#ifndef float_tt
#define float_tt  float
#endif
#define real      float
#endif  // FLOAT_PRECISION
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
  double normX,normY,normZ;
  double vect1X,vect1Y,vect1Z;
  double vect2X,vect2Y,vect2Z;
  double pointX,pointY,pointZ;
} plane;

typedef struct grainBoxStruct {
  int amorphFlag;
  double density,rmin, rFactor;  /* density, atomic distance, reduced atomic distance
				  * for amorphous material.  Red. r is for making a hex.
				  * closed packed structure, which will fill all space, 
				  * but will only be sparsely filled, and later relaxed.
				  */ 
  char *name;
  atom *unitCell; /* definition of unit cell */
  int natoms;     /* number of atoms in unit cell */
  double ax,by,cz; /* unit cell parameters */
  double alpha, beta, gamma; /* unit cell parameters */
  double tiltx,tilty,tiltz;
  double shiftx,shifty,shiftz;
  plane *planes;   /* pointer to array of bounding planes */
  double sphereRadius, sphereX,sphereY,sphereZ; /* defines a sphere instead of a grain with straight edges */
  int nplanes; /* number of planes in array planes */
} grainBox;

typedef struct superCellBoxStruct {
  double cmx,cmy,cmz;  /* fractional center of mass coordinates */
  double ax,by,cz;
  int natoms;
  atom *atoms; /* contains all the atoms within the super cell */
} superCellBox;

typedef struct atomBoxStruct {
  int used;   /* indicate here whether this atom is used in the
		 particular problem */
  int nx,ny,nz;
  float_tt dx,dy,dz;
  double B;
#if FLOAT_PRECISION == 1
  fftwf_complex ***potential;   /* 3D array containg 1st quadrant of real space potential */
  float_tt ***rpotential;   /* 3D array containg 1st quadrant of real space potential */
#else
  fftw_complex ***potential;   /* 3D array containg 1st quadrant of real space potential */
  float_tt ***rpotential;   /* 3D array containg 1st quadrant of real space potential */
#endif
} atomBox;

typedef struct detectorStruct {
  float_tt rInside,rOutside;
  float_tt k2Inside,k2Outside;
  char name[32];
  float_tt **image;        // place for storing avg image = sum(data)/Navg
  float_tt **image2;        // we will store sum(data.^2)/Navg 
  float_tt error;  
  float_tt shiftX,shiftY;
  int Navg;
} DETECTOR;

/**************** storage/memory allocation structures ***************/
template <class T> class fftw_allocator
{
public:
    typedef T                 value_type;
    typedef value_type*       pointer;
    typedef const value_type* const_pointer;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
    typedef std::size_t       size_type;
    typedef std::ptrdiff_t    difference_type;

  template <class U> 
    struct rebind { typedef fftw_allocator<U> other; };

    fftw_allocator() {}
    fftw_allocator(const fftw_allocator&) {}
    
  template <class U> 
    fftw_allocator(const fftw_allocator<U>&) {}
    
    ~fftw_allocator() {}

    pointer address(reference x) const 
        { return &x; }
        
    const_pointer address(const_reference x) const 
        { return x; }

    pointer allocate(size_type n, const_pointer = 0) 
    {
        void* p = fftw_malloc(n * sizeof(T));
        if (!p)
            throw std::bad_alloc();
        return static_cast<pointer>(p);
    }

    void deallocate(pointer p, size_type) 
        { fftw_free(p); }

    size_type max_size() const 
        { return static_cast<size_type>(-1) / sizeof(T); }

    void construct(pointer p, const value_type& x) 
        { new(p) value_type(x); }
    
    void destroy(pointer p) 
        { p->~value_type(); }

    void operator=(const fftw_allocator&x) 
        { }
};


/******************* Vector and array type definitions *****************/
// vectors aligned for SIMD instructions
typedef boost::multi_array<float_tt, 1, fftw_allocator<float_tt>> float1D_type;
typedef boost::multi_array<double, 1, fftw_allocator<double>> double1D_type;

// 2D arrays aligned for SIMD instructions
typedef boost::multi_array<float_tt, 2, fftw_allocator<float_tt>> float2D_type;
typedef boost::multi_array<double, 2, fftw_allocator<double>> double2D_type;
typedef boost::multi_array<int, 2, fftw_allocator<int>> int2D_type;
typedef boost::multi_array<std::complex<float_tt>, 2, fftw_allocator<std::complex<float_tt>>> complex2D_type;

// 3D arrays aligned for SIMD instructions
typedef boost::multi_array<float_tt, 3, fftw_allocator<float_tt>> float3D_type;
typedef boost::multi_array<std::complex<float_tt>, 3, fftw_allocator<std::complex<float_tt>>> complex3D_type;

#endif // STEMTYPES_H
