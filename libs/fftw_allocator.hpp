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

// File copied from REDHAWK Basic components fftlib 20140125.  Added double allocator.

/*
* This file is protected by Copyright. Please refer to the COPYRIGHT file distributed with this
* source distribution.
*
* This file is part of REDHAWK Basic Components fftlib library.
*
* REDHAWK Basic Components fftlib library is free software: you can redistribute it and/or modify it under the terms of
* the GNU General Public License as published by the Free Software Foundation, either
* version 3 of the License, or (at your option) any later version.
*
* REDHAWK Basic Components fftlib library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE. See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this
* program. If not, see http://www.gnu.org/licenses/.
*/

#ifndef FFTW_ALLOCATOR_H_
#define FFTW_ALLOCATOR_H_

#include <iostream>
#include <ostream>
#include <memory>
#include "fftw3.h"
#include <complex>
#include <vector>

// FFTWComplex class copied from VIGRA: https://github.com/ukoethe/vigra
// their copyright notice:
/*

The VIGRA License
=================
(identical to the MIT X11 License)

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:
                                                               
The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the
Software.
                                                               
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE. 

 */

typedef double fftw_real;

template <class T>
struct FFTWReal;

template <>
struct FFTWReal<fftw_complex>
{
    typedef double type;
};

template <>
struct FFTWReal<fftwf_complex>
{
    typedef float type;
};

template <>
struct FFTWReal<fftwl_complex>
{
    typedef long double type;
};

template <class T>
struct FFTWReal2Complex;

template <>
struct FFTWReal2Complex<double>
{
    typedef fftw_complex type;
    typedef fftw_plan plan_type;
};

template <>
struct FFTWReal2Complex<float>
{
    typedef fftwf_complex type;
    typedef fftwf_plan plan_type;
};

template <>
struct FFTWReal2Complex<long double>
{
    typedef fftwl_complex type;
    typedef fftwl_plan plan_type;
};


template <class Real = double>
class FFTWComplex
{
  public:
        /** The wrapped complex type
*/
      typedef typename FFTWReal2Complex<Real>::type complex_type;

        /** The complex' component type, as defined in '<TT>fftw3.h</TT>'
*/
    typedef Real value_type;

        /** reference type (result of operator[])
*/
    typedef value_type & reference;

        /** const reference type (result of operator[] const)
*/
    typedef value_type const & const_reference;

        /** iterator type (result of begin() )
*/
    typedef value_type * iterator;

        /** const iterator type (result of begin() const)
*/
    typedef value_type const * const_iterator;

        /** The norm type (result of magnitude())
*/
    typedef value_type NormType;

        /** The squared norm type (result of squaredMagnitde())
*/
    typedef value_type SquaredNormType;

        /** Construct from real and imaginary part.
Default: 0.
*/
    FFTWComplex(value_type const & re = 0.0, value_type const & im = 0.0)
    {
        data_[0] = re;
        data_[1] = im;
    }

        /** Copy constructor.
*/
    FFTWComplex(FFTWComplex const & o)
    {
        data_[0] = o.data_[0];
        data_[1] = o.data_[1];
    }

        /** Copy constructor.
*/
    template <class U>
    FFTWComplex(FFTWComplex<U> const & o)
    {
        data_[0] = (Real)o.real();
        data_[1] = (Real)o.imag();
    }

        /** Construct from plain <TT>fftw_complex</TT>.
*/
    FFTWComplex(fftw_complex const & o)
    {
        data_[0] = (Real)o[0];
        data_[1] = (Real)o[1];
    }

        /** Construct from plain <TT>fftwf_complex</TT>.
*/
    FFTWComplex(fftwf_complex const & o)
    {
        data_[0] = (Real)o[0];
        data_[1] = (Real)o[1];
    }

        /** Construct from plain <TT>fftwl_complex</TT>.
*/
    FFTWComplex(fftwl_complex const & o)
    {
        data_[0] = (Real)o[0];
        data_[1] = (Real)o[1];
    }

        /** Construct from std::complex.
*/
    template <class T>
    FFTWComplex(std::complex<T> const & o)
    {
        data_[0] = (Real)o.real();
        data_[1] = (Real)o.imag();
    }

        /** Assignment.
*/
    FFTWComplex& operator=(FFTWComplex const & o)
    {
        data_[0] = o.data_[0];
        data_[1] = o.data_[1];
        return *this;
    }

        /** Assignment.
*/
    template <class U>
    FFTWComplex& operator=(FFTWComplex<U> const & o)
    {
        data_[0] = (Real)o.real();
        data_[1] = (Real)o.imag();
        return *this;
    }

        /** Assignment.
*/
    FFTWComplex& operator=(fftw_complex const & o)
    {
        data_[0] = (Real)o[0];
        data_[1] = (Real)o[1];
        return *this;
    }

        /** Assignment.
*/
    FFTWComplex& operator=(fftwf_complex const & o)
    {
        data_[0] = (Real)o[0];
        data_[1] = (Real)o[1];
        return *this;
    }

        /** Assignment.
*/
    FFTWComplex& operator=(fftwl_complex const & o)
    {
        data_[0] = (Real)o[0];
        data_[1] = (Real)o[1];
        return *this;
    }

        /** Assignment.
*/
    FFTWComplex& operator=(double o)
    {
        data_[0] = (Real)o;
        data_[1] = 0.0;
        return *this;
    }

        /** Assignment.
*/
    FFTWComplex& operator=(float o)
    {
        data_[0] = (Real)o;
        data_[1] = 0.0;
        return *this;
    }

        /** Assignment.
*/
    FFTWComplex& operator=(long double o)
    {
        data_[0] = (Real)o;
        data_[1] = 0.0;
        return *this;
    }

        /** Assignment.
*/
    template <class T>
    FFTWComplex& operator=(std::complex<T> const & o)
    {
        data_[0] = (Real)o.real();
        data_[1] = (Real)o.imag();
        return *this;
    }

    reference re()
        { return data_[0]; }

    const_reference re() const
        { return data_[0]; }

    reference real()
        { return data_[0]; }

    const_reference real() const
        { return data_[0]; }

    reference im()
        { return data_[1]; }

    const_reference im() const
        { return data_[1]; }

    reference imag()
        { return data_[1]; }

    const_reference imag() const
        { return data_[1]; }

        /** Unary negation.
*/
    FFTWComplex operator-() const
        { return FFTWComplex(-data_[0], -data_[1]); }

        /** Squared magnitude x*conj(x)
*/
    SquaredNormType squaredMagnitude() const
        { return data_[0]*data_[0]+data_[1]*data_[1]; }

        /** Magnitude (length of radius vector).
*/
    NormType magnitude() const
        { return sqrt(squaredMagnitude()); }

        /** Phase angle.
*/
    value_type phase() const
        { return atan2(data_[1], data_[0]); }

        /** Access components as if number were a vector.
*/
    reference operator[](int i)
        { return data_[i]; }

        /** Read components as if number were a vector.
*/
    const_reference operator[](int i) const
        { return data_[i]; }

        /** Length of complex number (always 2).
*/
    int size() const
        { return 2; }

    iterator begin()
        { return data_; }

    iterator end()
        { return data_ + 2; }

    const_iterator begin() const
        { return data_; }

    const_iterator end() const
        { return data_ + 2; }

  private:
    complex_type data_;
};

template <class R>
inline bool operator ==(FFTWComplex<R> const &a, const FFTWComplex<R> &b) {
    return a.re() == b.re() && a.im() == b.im();
}

template <class R>
inline bool operator ==(FFTWComplex<R> const &a, double b) {
    return a.re() == b && a.im() == 0.0;
}

template <class R>
inline bool operator ==(double a, const FFTWComplex<R> &b) {
    return a == b.re() && 0.0 == b.im();
}

    /// not equal
template <class R>
inline bool operator !=(FFTWComplex<R> const &a, const FFTWComplex<R> &b) {
    return a.re() != b.re() || a.im() != b.im();
}

    /// not equal
template <class R>
inline bool operator !=(FFTWComplex<R> const &a, double b) {
    return a.re() != b || a.im() != 0.0;
}

    /// not equal
template <class R>
inline bool operator !=(double a, const FFTWComplex<R> &b) {
    return a != b.re() || 0.0 != b.im();
}

    /// add-assignment
template <class R>
inline FFTWComplex<R> & operator +=(FFTWComplex<R> &a, const FFTWComplex<R> &b) {
    a.re() += b.re();
    a.im() += b.im();
    return a;
}

    /// subtract-assignment
template <class R>
inline FFTWComplex<R> & operator -=(FFTWComplex<R> &a, const FFTWComplex<R> &b) {
    a.re() -= b.re();
    a.im() -= b.im();
    return a;
}

    /// multiply-assignment
template <class R>
inline FFTWComplex<R> & operator *=(FFTWComplex<R> &a, const FFTWComplex<R> &b) {
    typename FFTWComplex<R>::value_type t = a.re()*b.re()-a.im()*b.im();
    a.im() = a.re()*b.im()+a.im()*b.re();
    a.re() = t;
    return a;
}

    /// divide-assignment
template <class R>
inline FFTWComplex<R> & operator /=(FFTWComplex<R> &a, const FFTWComplex<R> &b) {
    typename FFTWComplex<R>::value_type sm = b.squaredMagnitude();
    typename FFTWComplex<R>::value_type t = (a.re()*b.re()+a.im()*b.im())/sm;
    a.im() = (b.re()*a.im()-a.re()*b.im())/sm;
    a.re() = t;
    return a;
}

    /// add-assignment with scalar double
template <class R>
inline FFTWComplex<R> & operator +=(FFTWComplex<R> &a, double b) {
    a.re() += (R)b;
    return a;
}

    /// subtract-assignment with scalar double
template <class R>
inline FFTWComplex<R> & operator -=(FFTWComplex<R> &a, double b) {
    a.re() -= (R)b;
    return a;
}

    /// multiply-assignment with scalar double
template <class R>
inline FFTWComplex<R> & operator *=(FFTWComplex<R> &a, double b) {
    a.re() *= (R)b;
    a.im() *= (R)b;
    return a;
}

    /// divide-assignment with scalar double
template <class R>
inline FFTWComplex<R> & operator /=(FFTWComplex<R> &a, double b) {
    a.re() /= (R)b;
    a.im() /= (R)b;
    return a;
}

    /// addition
template <class R>
inline FFTWComplex<R> operator +(FFTWComplex<R> a, const FFTWComplex<R> &b) {
    a += b;
    return a;
}

    /// right addition with scalar double
template <class R>
inline FFTWComplex<R> operator +(FFTWComplex<R> a, double b) {
    a += b;
    return a;
}

    /// left addition with scalar double
template <class R>
inline FFTWComplex<R> operator +(double a, FFTWComplex<R> b) {
    b += a;
    return b;
}

    /// subtraction
template <class R>
inline FFTWComplex<R> operator -(FFTWComplex<R> a, const FFTWComplex<R> &b) {
    a -= b;
    return a;
}

    /// right subtraction with scalar double
template <class R>
inline FFTWComplex<R> operator -(FFTWComplex<R> a, double b) {
    a -= b;
    return a;
}

    /// left subtraction with scalar double
template <class R>
inline FFTWComplex<R> operator -(double a, FFTWComplex<R> const & b) {
    return (-b) += a;
}

    /// multiplication
template <class R>
inline FFTWComplex<R> operator *(FFTWComplex<R> a, const FFTWComplex<R> &b) {
    a *= b;
    return a;
}

    /// right multiplication with scalar double
template <class R>
inline FFTWComplex<R> operator *(FFTWComplex<R> a, double b) {
    a *= b;
    return a;
}

    /// left multiplication with scalar double
template <class R>
inline FFTWComplex<R> operator *(double a, FFTWComplex<R> b) {
    b *= a;
    return b;
}

    /// division
template <class R>
inline FFTWComplex<R> operator /(FFTWComplex<R> a, const FFTWComplex<R> &b) {
    a /= b;
    return a;
}

    /// right division with scalar double
template <class R>
inline FFTWComplex<R> operator /(FFTWComplex<R> a, double b) {
    a /= b;
    return a;
}

    /// absolute value (= magnitude)
template <class R>
inline typename FFTWComplex<R>::NormType abs(const FFTWComplex<R> &a)
{
    return a.magnitude();
}

    /// pahse
template <class R>
inline R arg(const FFTWComplex<R> &a)
{
    return a.phase();
}

    /// real part
template <class R>
inline R real(const FFTWComplex<R> &a)
{
    return a.real();
}

    /// imaginary part
template <class R>
inline R imag(const FFTWComplex<R> &a)
{
    return a.imag();
}

    /// complex conjugate
template <class R>
inline FFTWComplex<R> conj(const FFTWComplex<R> &a)
{
    return FFTWComplex<R>(a.re(), -a.im());
}

    /// norm (= magnitude)
template <class R>
inline typename FFTWComplex<R>::NormType norm(const FFTWComplex<R> &a)
{
    return a.magnitude();
}

    /// squared norm (= squared magnitude)
template <class R>
inline typename FFTWComplex<R>::SquaredNormType squaredNorm(const FFTWComplex<R> &a)
{
    return a.squaredMagnitude();
}


namespace std {
  template <class Real>
  ostream & operator<<(ostream & s, FFTWComplex<Real> const & v)
  {
    s << std::complex<Real>(v.re(), v.im());
    return s;
  }
}

template<class Ty>
class FFTWAllocator
{
  public:
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef Ty *pointer;
    typedef const Ty *const_pointer;
    typedef Ty& reference;
    typedef const Ty& const_reference;
    typedef Ty value_type;
    
    pointer address(reference val) const
        { return &val; }
        
    const_pointer address(const_reference val) const
        { return &val; }
        
    template<class Other>
    struct rebind
    {
        typedef FFTWAllocator<Other> other;
    };
    
    FFTWAllocator() throw()
    {}
    
    template<class Other>
    FFTWAllocator(const FFTWAllocator<Other>& right) throw()
    {}
    
    template<class Other>
    FFTWAllocator& operator=(const FFTWAllocator<Other>& right)
    {
        return *this;
    }
    
    pointer allocate(size_type count, void * = 0)
    {
        return (pointer)fftw_malloc(count * sizeof(Ty));
    }
    
    void deallocate(pointer ptr, size_type count)
    {
        fftw_free(ptr);
    }
    
    void construct(pointer ptr, const Ty& val)
    {
        new(ptr) Ty(val);
        
    }
    
    void destroy(pointer ptr)
    {
        ptr->~Ty();
    }
    
    size_type max_size() const throw()
    {
      return size_t(-1);
      //return NumericTraits<std::ptrdiff_t>::max() / sizeof(Ty);
    }
};

// determines if memory from another
// allocator can be deallocated from this one

template<typename T>
inline bool operator==(FFTWAllocator<T> const& lhs, FFTWAllocator<T> const& rhs) { 
    return operator==(static_cast<T&>(lhs), 
                       static_cast<T&>(rhs)); 
}


#endif /* FFTW_ALLOCATOR_H_ */
