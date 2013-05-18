#ifndef SPLINES_H
#define SPLINES_H

/* "THE BEER-WARE LICENSE" (Revision 42): Devin Lane wrote this file. As long as you retain 
 * this notice you can do whatever you want with this stuff. If we meet some day, and you
 * think this stuff is worth it, you can buy me a beer in return. 
 
 http://shiftedbits.org/2011/01/30/cubic-spline-interpolation/#more-191
 */

#include <vector>
#include <iostream>

/** Templated on type of X, Y. X and Y must have operator +, -, *, /. Y must have defined
 * a constructor that takes a scalar. 
 * This is a basic cubic spline.  */
template <typename X, typename Y>
class CubicSpline {
public:
    /** An empty, invalid spline */
    CubicSpline() {}
    
    /** A spline with x and y values */
    CubicSpline(const std::vector<X>& x, const std::vector<Y>& y) {
        if (x.size() != y.size()) {
            std::cerr << "X and Y must be the same size " << std::endl;
            return;
        }
        
        if (x.size() < 3) {
            std::cerr << "Must have at least three points for interpolation" << std::endl;
            return;
        }
        
        typedef typename std::vector<X>::difference_type size_type;
        
        size_type n = y.size() - 1;
        
        std::vector<Y> b(n), d(n), a(n), c(n+1), l(n+1), u(n+1), z(n+1);
        std::vector<X> h(n+1);

        l[0] = Y(1);
        u[0] = Y(0);
        z[0] = Y(0);
        h[0] = x[1] - x[0];
            
        for (size_type i = 1; i < n; i++) {
            h[i] = x[i+1] - x[i];
            l[i] = Y(2 * (x[i+1] - x[i-1])) - Y(h[i-1]) * u[i-1];
            u[i] = Y(h[i]) / l[i];
            a[i] = (Y(3) / Y(h[i])) * (y[i+1] - y[i]) - (Y(3) / Y(h[i-1])) * (y[i] - y[i-1]);
            z[i] = (a[i] - Y(h[i-1]) * z[i-1]) / l[i];
        }
            
        l[n] = Y(1);
        z[n] = c[n] = Y(0);
        
        for (size_type j = n-1; j >= 0; j--) {
            c[j] = z[j] - u[j] * c[j+1];
            b[j] = (y[j+1] - y[j]) / Y(h[j]) - (Y(h[j]) * (c[j+1] + Y(2) * c[j])) / Y(3);
            d[j] = (c[j+1] - c[j]) / Y(3 * h[j]);
        }
        
        for (size_type i = 0; i < n; i++) {
            mElements.push_back(Element(x[i], y[i], b[i], c[i], d[i]));
        }        
    }
    virtual ~CubicSpline() {}
    
    Y operator[](const X& x) const {
        return interpolate(x);
    }
    
    Y interpolate(const X&x) const {
        if (mElements.size() == 0) return Y();
        
        typename std::vector<element_type>::const_iterator it;
		// looks for existing element closest to given X coordinate to be interpolated from
        it = std::lower_bound(mElements.begin(), mElements.end(), element_type(x));
        if (it != mElements.begin()) {
            it--;
        }   
            
        return it->eval(x);
    }
    
    std::vector<Y> operator[](const std::vector<X>& xx) const {
        return interpolate(xx);
    }
    
    /* Evaluate at multiple locations, assuming xx is sorted ascending */
    std::vector<Y> interpolate(const std::vector<X>& xx) const {
        if (mElements.size() == 0) return std::vector<Y>(xx.size());
        
        typename std::vector<X>::const_iterator it;
        typename std::vector<element_type>::const_iterator it2;
        it2 = mElements.begin();
        std::vector<Y> ys;
        for (it = xx.begin(); it != xx.end(); it++) {
            it2 = std::lower_bound(it2, mElements.end(), element_type(*it));
            if (it2 != mElements.begin()) {
                it2--;
            }
                
            ys.push_back(it2->eval(*it));
        }

        return ys;
    }

protected:
    
    class Element {
    public:
        Element(X _x) : x(_x) {}
        Element(X _x, Y _a, Y _b, Y _c, Y _d)
        : x(_x), a(_a), b(_b), c(_c), d(_d) {}
        
        Y eval(const X& xx) const {
            X xix(xx - x);
            return a + b * xix + c * (xix * xix) + d * (xix * xix * xix);
        }
        
        bool operator<(const Element& e) const {
            return x < e.x;
        }
        bool operator<(const X& xx) const {
            return x < xx;
        }
        
        X x;
        Y a, b, c, d;
    };
            
    typedef Element element_type;
    std::vector<element_type> mElements;
};

/** Templated on type of X, Y. X and Y must have operator +, -, *, /. Y must have defined
 * a constructor that takes a scalar. 
 *  Akima splines are meant to not oscillate wildly around noisy data or outliers 
 
 [1] Spline fit as in H.Akima, J. ACM 17(1970)p.589-602
		'A New Method of Interpolation and Smooth
		Curve Fitting Based on Local Procedures'

 [2] H.Akima, Comm. ACM, 15(1972)p.914-918
 
 */
template <typename X_t, typename Y_t>
class AkimaSpline {
public:
    /** An empty, invalid spline */
    AkimaSpline() {}
    
    /** A spline with x and y values */
    AkimaSpline(const std::vector<X_t>& _X, const std::vector<Y_t>& _Y)	: X(_X), Y(_Y){

		int iv = (int) X.size();
		int i;
		Y_t a, b;
        if (X.size() != Y.size()) {
            std::cerr << "X and Y must be the same size " << std::endl;
            return;
        }
        
        if (X.size() < 3) {
            std::cerr << "Must have at least three points for interpolation" << std::endl;
            return;
        }
        
        typedef typename std::vector<X_t>::difference_type size_type;


		XM = std::vector<Y_t>(X.size()+4);
		Z = std::vector<Y_t>(X.size()+1);

		//X[0]=2.0*X[1]-X[2];
		//Calculate Akima coefficients, a and b
		for (i=0; i<iv-1; i++)
		{
			//Shift i to i+2
			XM[i+2]=(Y[i+1]-Y[i])/(X[i+1]-X[i]);
		}
		XM[iv+2]=2*XM[iv+1]-XM[iv];
		XM[iv+3]=2*XM[iv+2]-XM[iv+1];
		XM[2]=2*XM[3]-XM[4];
		XM[1]=2*XM[2]-XM[3];
		for (i=0; i<X.size(); i++)  
		{
			a=fabs(XM[i+3]-XM[i+2]);
			b=fabs(XM[i+1]-XM[i]);
			if (a+b!=0) 
			{
				Z[i]=(a*XM[i+1]+b*XM[i+2])/(a+b);
				continue;
			}
			Z[i]=(XM[i+2]+XM[i+1])/2.0f;
		}
		for ( i = 0; i < X.size()-1; i++) {
            mElements.push_back(Element(X[i], X[i+1], XM[i+2], Y[i], Z[i], Z[i+1]));
        }        
	}
    virtual ~AkimaSpline() {}
    
    Y_t operator[](const X_t& x) {
        return interpolate(x);
    }
    
    Y_t interpolate(const X_t &x) {
        if (mElements.size() == 0) return Y_t();
        size_t i = X.size()-1;
        typename std::vector<element_type>::iterator it;
		// looks for existing element closest to given X coordinate to be interpolated from
        it = std::lower_bound(mElements.begin(), mElements.end(), element_type(x));
        if (it != mElements.begin()) {
            it--;
			i--;
        }   
            
        return it->eval(x);
    }
    
    std::vector<Y_t> operator[](std::vector<X_t>& xx) {
        return interpolate(xx);
    }
    
    /* Evaluate at multiple locations, assuming xx is sorted ascending */
    std::vector<Y_t> interpolate(const std::vector<X_t>& xx)  {
        if (mElements.size() == 0) return std::vector<Y_t>(xx.size());
        
        typename std::vector<X_t>::iterator it;
        typename std::vector<element_type>::iterator it2;
        it2 = mElements.begin();
        std::vector<Y_t> ys;
        for (it = xx.begin(); it != xx.end(); it++) {
            it2 = std::lower_bound(it2, mElements.end(), element_type(*it));
            if (it2 != mElements.begin()) {
                it2--;
            }
                
            ys.push_back(it2->eval(*it));
        }

        return ys;
    }

protected:

	std::vector<X_t> X;
	std::vector<Y_t> Y;
	std::vector<Y_t> XM;
	std::vector<Y_t> Z;
    
    class Element {
    public:
		Element(X_t _x) : x(_x) {}
		Element(X_t _x, X_t _xP1, Y_t _xmP2, Y_t _y, Y_t _z, Y_t _zP1) : 
		  x(_x), xP1(_xP1), y(_y), xmP2(_xmP2), z(_z), zP1(_zP1) {}
        
        Y_t eval(const X_t& xx) {
			double yy;
			b=xP1-x;
			a=xx - x;
			yy=y + z * a + (3.0*xmP2 - 2.0 * z - zP1) * a * a / b;
            return static_cast<float_tt>(yy+(z+zP1-2.0*xmP2)*a*a*a/(b*b));
        }
        
        bool operator<(const Element& e) const {
            return x < e.x;
        }
        bool operator<(const X_t& xx) const {
            return x < xx;
        }
        
        X_t x, xP1;
        Y_t y, xmP2, z, zP1, a, b;
	private:
		int i;
    };
            
    typedef Element element_type;
    std::vector<element_type> mElements;
};

#endif