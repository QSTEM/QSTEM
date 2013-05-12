#ifndef MEMORY_H
#define MEMORY_H

// #include <stdlib.h>
#include <stdio.h>

#include "fftw3.h"
#include "boost/shared_ptr.hpp"
#include "boost/multi_array.hpp"
#include "stemtypes_fftw3.h"

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
typedef boost::multi_array<float_tt, 1, fftw_allocator<float_tt> > _float1DArray;
typedef boost::multi_array<double, 1, fftw_allocator<double> > _double1DArray;
typedef boost::multi_array<int, 1, fftw_allocator<int> > int1DArray;

// 2D arrays aligned for SIMD instructions
// we subclass the floating point classes to add on rotation and inversion methods
typedef boost::multi_array<float_tt, 2, fftw_allocator<float_tt> > _float2DArray;
typedef boost::multi_array<double, 2, fftw_allocator<double> > _double2DArray;
typedef _float2DArray::array_view<1>::type float1DView;
typedef _double2DArray::array_view<1>::type double1DView;

typedef boost::multi_array<int, 2, fftw_allocator<int> > int2DArray;

// the extra dimension here is because FFTW treats complex numbers as 2-element arrays.
//    Those arrays confuse boost::multi_array.  This provides equivalent functionality.
typedef boost::multi_array<float_tt, 3, fftw_allocator<float_tt> > complex2DArray;
typedef boost::multi_array<double, 3, fftw_allocator<double> > complexDouble2DArray;
typedef complex2DArray::array_view<2>::type complex1DView;
typedef complexDouble2DArray::array_view<2>::type complexDouble1DView;

// 3D arrays aligned for SIMD instructions
typedef boost::multi_array<float_tt, 3, fftw_allocator<float_tt> > float3DArray;
typedef boost::multi_array<double, 3, fftw_allocator<double> > double3DArray;
typedef float3DArray::array_view<2>::type float2DView;
typedef double3DArray::array_view<2>::type complexDouble2DView;
// the extra dimension here is because FFTW treats complex numbers as 2-element arrays.
//    Those arrays confuse boost::multi_array.  This provides equivalent functionality.
typedef boost::multi_array<float_tt, 4, fftw_allocator<float_tt> > complex3DArray;
typedef boost::multi_array<double, 4, fftw_allocator<double> > complexDouble3DArray;
typedef complex3DArray::array_view<3>::type complex2DView;
typedef complexDouble3DArray::array_view<3>::type complexeDouble2DView;

class float1DArray : public _float1DArray
{
public:
	float1DArray(){};
	float1DArray(int size);
	void rotate(float1DArray &output, float_tt x, float_tt y, float_tt z);
};

// TODO: template this somehow to avoid duplication of the rotate and invert methods
class double1DArray : public _double1DArray
{
public:
	double1DArray(){};
	double1DArray(int size);
	void rotate(double1DArray &output, double x, double y, double z);
};

class float2DArray: public _float2DArray
{
public:
	float2DArray(){};
	float2DArray(int nx, int ny);
	void rotate(float2DArray &output, float_tt x, float_tt y, float_tt z);
	void invert(float2DArray &output);
};

// TODO: template this somehow to avoid duplication of the rotate and invert methods
class double2DArray : public _double2DArray
{
public:
	double2DArray(){};
	double2DArray(int nx, int ny);
	void rotate(double2DArray &output, float_tt x, float_tt y, float_tt z);
	void invert(double2DArray &output);
};

/************** Functions to get vectors/arrays out ****************/
float1DArray float1D(int size, std::string message);
double1DArray double1D(int size, std::string message);
int1DArray int1D(int size, std::string message);

float2DArray float2D(int nx, int ny, std::string message);
double2DArray double2D(int nx, int ny, std::string message);
int2DArray int2D(int nx, int ny, std::string message);
complex2DArray complex2D(int nx, int ny, std::string message);

float3DArray float3D(int nx, int ny, int nz, std::string message);
complex3DArray complex3D(int nx, int ny,int nz, std::string message);

void **any2D( int nx, int ny,int size, const char *message );
void ***any3D( int nx, int ny,int nz,int size, const char *message );

#endif

