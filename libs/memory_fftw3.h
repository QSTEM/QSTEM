#ifndef MEMORY_H
#define MEMORY_H

// #include <stdlib.h>
#include <stdio.h>
#include "fftw3.h"

// precision is defined in here
#include "stemtypes_fftw3.h"

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
// float vector aligned for SIMD instructions
typedef boost::multi_array<float_tt, 1, fftw_allocator<float_tt>> float1D_type;
typedef boost::multi_array<double_tt, 1, fftw_allocator<float_tt>> double1D_type;

// float 2D array aligned for SIMD instructions
typedef boost::multi_array<float_tt, 2, fftw_allocator<float_tt>> float2D_type;
typedef boost::multi_array<double, 2, fftw_allocator<double>> double2D_type;
typedef boost::multi_array<int, 2, fftw_allocator<int>> int2D_type;

// float 3D array aligned for SIMD instructions
typedef boost::multi_array<float_tt, 3, fftw_allocator<float_tt>> float3D_type;



/************** Functions to get vectors/arrays out ****************/
boost::shared_ptr<float_1D_type> float1D(int size, std::string message);
boost::shared_ptr<double_1D_type> double1D(int size, std::string message);
boost::shared_ptr<float2D_type> float2D(int nx, int ny, std::string message);
boost::shared_ptr<double2D_type> double2D(int nx, int ny, std::string message);
boost::shared_ptr<int2D_type> int2D(int nx, int ny, std::string message);


/*---------------------------- long2D() -------------------------------*/
/*
	2D array allocator for type long
	make space for m[0...(nx-1)][0..(ny-1)]
	
	message = char[] with error message
	
*/
long **long2D( int nx, int ny, const char *message );
int **int2D( int nx, int ny, const char *message );

/*---------------------------- float2D() -------------------------------*/
/*
	2D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)]

*/
float_tt **float2D( int nx, int ny, const char *message );
float **float32_2D( int nx, int ny, const char *message );
float ***float32_3D( int nx, int ny,int nz, const char *message );

/*---------------------------- double2D() -------------------------------*/
/*
	2D array allocator for type double
	make space for m[0...(nx-1)][0..(ny-1)]

*/
float_tt ***float3D( int nx, int ny,int nz, const char *message );
double **double2D( int nx, int ny, const char *message );

fftw_complex  **complex2D(int nx, int ny, const char *message);
fftwf_complex **complex2Df(int nx, int ny, const char *message);  // single precision
fftw_complex  ***complex3D(int nx, int ny,int nz, const char *message);
fftwf_complex ***complex3Df(int nx, int ny,int nz, const char *message); // single precision

void **any2D( int nx, int ny,int size, const char *message );
void ***any3D( int nx, int ny,int nz,int size, const char *message );

#endif

