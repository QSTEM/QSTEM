/* file: STEMutil.c */

#include <stdio.h>	/*  ANSI-C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#ifndef WIN32
#include <sys/time.h>
#endif

#include "defines.h"
#include "stemlib.h"
#include "stemutil.h"
#include "memory_fftw3.h"	/* memory allocation routines */
// #include "tiffsubs.h"
#include "matrixlib.h"
#include "readparams.h"
#include "fileio_fftw3.h"
// #include "fileio.h"
// #include "floatdef.h"

// #define WRITEK

#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#ifdef WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif

int feTableRead=0;	/* flag to remember if the param file has been read */
QSfMat fparams;	/* to get fe(k) parameters */
const int nl=3, ng=3;	/* number of Lorenzians and Gaussians */


/******************************************************
 * memcopy replaces memcpy, because memcpy seems to have a bug!
 */
void *memcopy(void *dest, const void *src, size_t n)
{
  int i;

  for (i=0;i<n;i++)
    ((char *)dest)[i] = ((char *)src)[i]; 

  return dest;
  /* memcpy(dest,src,n) */
}

/*-------------------- bessi0() ---------------*/
/*
    modified Bessel function I0(x)
    see Abramowitz and Stegun page 379

    x = (float_tt) real arguments

    12-feb-1997 E. Kirkland
 */
 float_tt bessi0( float_tt x )
 {
 	int i;
 	float_tt ax, sum, t;
 	
 	float_tt i0a[] = { 1.0f, 3.5156229f, 3.0899424f, 1.2067492f,
		0.2659732f, 0.0360768f, 0.0045813f };

 	float_tt i0b[] = { 0.39894228f, 0.01328592f, 0.00225319f,
 		-0.00157565f, 0.00916281f, -0.02057706f, 0.02635537f,
 		-0.01647633f, 0.00392377f};

	ax = fabs( x );
	if( ax <= 3.75 ) {
		t = x / 3.75f;
		t = t * t;
		sum = i0a[6];
		for( i=5; i>=0; i--) sum = sum*t + i0a[i]; 
	} else {
		t = 3.75f / ax;
		sum = i0b[8];
		for( i=7; i>=0; i--) sum = sum*t + i0b[i];
		sum = exp( ax ) * sum / sqrt( ax );
	}
	return( sum );

}  /* end bessi0() */

/*-------------------- bessk0() ---------------*/
/*
    modified Bessel function K0(x)
    see Abramowitz and Stegun page 380
    
    Note: K0(0) is not define and this function
        returns 1E20
 
    x = (float_tt) real arguments
    
    this routine calls bessi0() = Bessel function I0(x)
    
    12-feb-1997 E. Kirkland
 */
 float_tt bessk0( float_tt x )
 {
 	float_tt bessi0(float_tt);
 
 	int i;
 	float_tt ax, x2, sum;
 	float_tt k0a[] = { -0.57721566f, 0.42278420f, 0.23069756f,
 		 0.03488590f, 0.00262698f, 0.00010750f, 0.00000740f};
 	        
 	float_tt k0b[] = { 1.25331414f, -0.07832358f, 0.02189568f,
 		 -0.01062446f, 0.00587872f, -0.00251540f, 0.00053208f};

	ax = fabs( x );
	if( (ax > 0.0)  && ( ax <=  2.0 ) ){
		x2 = ax/2.0f;
		x2 = x2 * x2;
		sum = k0a[6];
		for( i=5; i>=0; i--) sum = sum*x2 + k0a[i];
		sum = -log(ax/2.0f) * bessi0(x) + sum;
	} else if( ax > 2.0f ) {
		x2 = 2.0f/ax;
		sum = k0b[6];
		for( i=5; i>=0; i--) sum = sum*x2 + k0b[i];
		sum = exp( -ax ) * sum / sqrt( ax );
	} else sum = 1.0e20f;
	return ( sum );

}  /* end bessk0() */



/*--------------------- wavelength() -----------------------------------*/
/*
	return the electron wavelength (in Angstroms)
	keep this is one place so I don't have to keep typing in these
	constants (that I can never remember anyhow)

	ref: Physics Vade Mecum, 2nd edit, edit. H. L. Anderson
		(The American Institute of Physics, New York) 1989
		page 4.

	kev = electron energy in keV

*/

float_tt wavelength( float_tt kev )
{ 
  /* electron wavelength in Angstroms */
  float_tt w = static_cast<float_tt>(HC/sqrt( kev * ( 2*EMASS + kev ) ));
  
  return( w );
  
}  /* end wavelength() */


float_tt getTime() {
#ifdef WIN32
	return (float_tt)time(NULL);
#else
  timev tv;
  timez tz;

  gettimeofday(&tv, &tz); 
  // printf("Time: %ldsec, %ld usec\n",tv.tv_sec,tv.tv_usec);
  return (float_tt)tv.tv_sec+0.000001*(float_tt)tv.tv_usec;
#endif
}


/*--------------------- cputim() -----------------------------------*/
/*
  retrieve current CPU time in seconds
*/

float_tt cputim()
{
  return ( ( (float_tt)clock() ) / ( (float_tt)CLOCKS_PER_SEC) );
  
}  /* end cputim() */



/*-------------------- rangauss() -------------------------- */
/*
  Return a normally distributed random number with 
  zero mean and unit variance using Box-Muller method
  
  ranflat() is the source of uniform deviates
  
  ref.  Numerical Recipes, 2nd edit. page 289
  
  added log(0) test 10-jan-1998 E. Kirkland
*/
float_tt rangauss( unsigned long *iseed )
{
  float_tt ranflat( unsigned long* );
  float_tt x1, x2, y;
  
  /* be careful to avoid taking log(0) */
  do{
    x1 = ranflat( iseed );
    x2 = ranflat( iseed );
    
  } while ( (x1 < 1.0e-30) || (x2 < 1.0e-30) );
  
  y = sqrt( - 2.0f * log(x1) ) * (float_tt)(cos( TWOPI * x2 ));
  
  return( y );
  
}  /* end rangauss() */


/*--------------------- sigma() -----------------------------------*/
/*
	return the interaction parameter sigma in radians/(kv-Angstroms)
	keep this in one place so I don't have to keep typing in these
	constants (that I can never remember anyhow)

	ref: Physics Vade Mecum, 2nd edit, edit. H. L. Anderson
		(The American Institute of Physics, New York) 1989
		page 4.

	kev = electron energy in keV

*/

float_tt sigma( float_tt kev )
{
  float_tt s, wavl, x;
  const float_tt emass=510.99906f; /* electron rest mass in keV */
  float_tt wavelength( float_tt kev );  /*  get electron wavelength */
  
  x = ( emass + kev ) / ( 2.0f*emass + kev);  
  wavl = wavelength( kev );
  //pi = 4.0f * atan( 1.0f );
  
  s = static_cast<float_tt>(2.0 * PI * x / (wavl*kev));  // 2*pi*kz*(1+kev/emaxx)/(2*emass+kev)
  
  return( s );
  
}  /* end sigma() */

/*---------------------------- ranflat -------------------------------*/
/*
  return a random number in the range 0.0->1.0
  with uniform distribution
  
  the 'Magic Numbers' are from 
  Numerical Recipes 2nd edit pg. 285
*/
float_tt ranflat( unsigned long *iseed )
{
  static unsigned long a=1366, c=150889L, m=714025L;
  
  *iseed = ( a * (*iseed) + c ) % m;
  
  return( ((float_tt) *iseed)/m);
  
}  /* end ranflat() */

/*--------------------- parlay() -----------------------------------*/
/*
  subroutine to parse the atomic layer stacking sequence 
  for use with multislice programs.
  
  This converts layer structure definition of the form:
      2(abc)d
  into a sequence of numerical indices where a=1, b=2, c=3, ...
  The main attraction is the repeat operator n(...) where
  the structure inside the parenthesis (...) is repeated n times.
  for instance the above structure is translated into:
      1 2 3 1 2 3 4
  The parenthesis may be nested up to 100 levels (determined by
  nlmax). For instance  5( 2(abc) 3(cba) ) is also a valid structure
  definition. This is a compact way of specifying the layer structure
  for a multislice calculation. tab's may be present in the structure
  anywhere (they are ignored).

  This is done by pushing the position of each '(' and its repeat
  count onto a stack. Each ')' pops the last entry from the stack
  and invokes a duplication process.

  Layers refer to the distinquishable subset of different
  types of layers, whereas slices refer to the way these
  layers are sequentially stacked.

  fortran version started 22-dec-1989 earl j. kirkland
  added nested parenthesis by stacking entry points
     31-jan-1990 ejk
  added tab handling (in ascii and ebcdic) 1-feb-1990 ejk
  minor changes to error messages to make rs6000's happy
      7-july-1992 ejk
  converted to ANSI-C 26-jun-1995 ejk
  added include of stdlib.h 6-july-1995 ejk

   c             = input character string
   islice[nsmax] = integer array to get stacking sequence indicies
   nsmax         = (integer) size of layer
   lmax          = (integer) maximum allowed layer index
   nslice        = (integer) number of layers
   returned value= (integer) success/error code
                       0 : success
                      -1 : layer out of range
                      -2 : missing left parenthesis
                      -3 : parenthesis nested too deep
                      -4 : bad repeat code
                      -5 : unmatched right parenthesis
                      -6 : invalid character
                      -7 : too many layers
                      -8 : incomplete stacking sequence

   fperr       = (int) if this is not a NULL then write
                   error messages 

  NOTE: islice and nslice may be modified by this routine

  cname determines the mapping of characters into numbers.

*/

#define NSTKMAX 100	/* maximum stack depth (i.e. maximum level
				 of nesting the parenthesis) */
#define NCHARS 52	/* maximum number of character symbols */

int parlay( const char c[], int islice[], int nsmax, int lmax,
	    int *nslice, int fperr )
{
  /* define our own symbol sequence */
  const char cname[] =
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  
  int ic, lenc, name, i, j, istack, n,
    ipoint[NSTKMAX], repeat[NSTKMAX];
  
  /*  initialize misc constants  */
  
  *nslice = 0;
  istack = 0;
  lenc = (int) strlen( c );  /* length of c */
  
  for( ic=0; ic<lenc; ic++ ) {	/* loop over all characters in c */
    
    /*  skip over embedded spaces, tabs and wierd characters */
    while( isspace( c[ic] )  ) ic++;
    
    /*  if its a character do this  */		
    if( isalpha( c[ic] ) ) {
      for(name=0; name<NCHARS; name++) 
	if( c[ic] == cname[name] ) break;
      if( name > lmax ) {
	if( fperr != 0 ) 
	  printf("Layer /%c/ out of range in the stacking sequence"
		 " in subroutine parlay.\n", c[ic] );
	return( -1 );
      } else {
	if( *nslice >= nsmax ) {
	  if( fperr != 0 ) 
	    printf("Too many layers generated in the stacking"
		   " sequence in parlay.\n");
	  return( -7 );
	}
	islice[ (*nslice)++ ] = name;
      }
      
      /*  if its a number then extract repeat count up to the '('
	  and save it until a ')' appears (i.e. push on the stack)
      */
    } else if ( isdigit( c[ic] ) ) {
      if( istack >= NSTKMAX ) {
	if( fperr != 0 ) 
	  printf("Parenthesis nested too deep in "
		 "the stacking sequence in subroutine parlay.\n");
	return( -3 );
      }
      repeat[istack] = atoi( &c[ic] ) - 1; 
      while( isdigit(c[ic]) || isspace(c[ic]) ) ic++;
      if( c[ic] != '(' ) {
	if( fperr != 0 ) 
	  printf("Missing left parenthesis in "
		 "the stacking sequence in subroutine parlay.'\n");
	for(i=0; i<=ic; i++) printf("%c", c[i]);
	printf("?\n");
	return( -2 );
      }
      ipoint[istack++] = *nslice;
      
      /*  if its a ')' then repeat back to the last '('
	  (i.e. pop the stack) */
      
    } else if( c[ic] == ')' ) {
      if( istack <= 0 ) {
	if( fperr != 0 ) 
	  printf("Unmatched right parenthesis in "
		 "the stacking sequence in subroutine parlay.\n");
	return( -5 );
      }
      n = *nslice;
      istack--;
      for( j=0; j<repeat[istack]; j++)
	for(i=ipoint[istack]; i<n; i++){
	  if( *nslice >= nsmax ) {
	    if( fperr != 0 ) 
	      printf("Too many layers generated in the stacking"
		     " sequence in parlay.\n");
	    return( -7 );
	  }
	  islice[ (*nslice)++ ] = islice[ i ];
	}
    } else {
      if( fperr != 0 ) 
	printf("Invalid character /%c/ encountered in "
	       "the stacking sequence in subroutine parlay.\n",
	       c[ic]);
      return( -6 );
    }
    
  } /* end for( ic... */
  
  if( istack != 0 ) {
    if( fperr != 0 ) 
      printf("incomplete stacking sequence in parlay.\n");
    return( -8 ); 
  } else return( 0 );
  
#undef NSTKMAX 
#undef NCHARS
  
}  /* end parlay() */

/*------------------------- powerof2() ------------------------*/
/*
	return the nearest power of 2 greater than or equal to 
	the argument
*/
/*
  long powerof2( long n )
  {
  int ln;
  long n2;
  
  ln = 1;
  n2 = 2;
  while( (n2 < n) && (ln<31) ) { n2*=2; ln++; }
  
  return( n2 );
  }
*/

double cubicInterpolate (double p[4], double x) {
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}

double bicubicInterpolate (double p[4][4], double x, double y) {
	double arr[4];
	arr[0] = cubicInterpolate(p[0], y);
	arr[1] = cubicInterpolate(p[1], y);
	arr[2] = cubicInterpolate(p[2], y);
	arr[3] = cubicInterpolate(p[3], y);
	return cubicInterpolate(arr, x);
}

/* function bicubic from Matlab toolbox
 * zz = values of function F(z,x) F(ix,iz) = zz[iz*Nx+ix];
 * s = z-coordinate
 * t = x-coordinate
 */ 
#define FX (ptr[xi-1]*x0 + ptr[xi]*x1 + ptr[xi+1]*x2 + ptr[xi+2]*x)
float_tt bicubic(QSfMat ff,int Nz, int Nx,float_tt z,float_tt x) {
  static float_tt x0,x1,x2,f;
  static float_tt *ptr;
  static int xi,zi;

  // Now interpolate using computationally efficient algorithm.
  // s and t are the x- and y- coordinates of the points we are looking for
  // it is assumed that the function ff is defined at equally spaced intervals
  // 0,1,...N-1 in both directions.
 
  xi = (int)x;
  zi = (int)z;
  x   = x-(float_tt)xi;
  z   = z-(float_tt)zi;
  
  
  x0  = ((2.0f-x)*x-1.0f)*x;
  x1  = (3.0f*x-5.0f)*x*x+2.0f;
  x2  = ((4.0f-3.0f*x)*x+1.0f)*x;
  x   = (x-1.0f)*x*x;


  if (Nz > 2) {
	  // these were ff[0][0]
    if ((zi > Nz-3) || (xi > Nx-3)) return 0.0;
    if ((zi < 1) && (xi < 1)) return ff(0,0);  // was ff[0][0]
    //if (zi < 1) return ff[0][xi];
	if (zi < 1) return ff(xi,0);
    //if (xi < 1) return ff[zi][0];  // was ff[zi][0]
	if (xi < 1) return ff(0,zi);  // was ff[zi][0]

    ptr = ff.row(zi-1).data();
    f   = FX * (((2.0f-z)*z-1)*z);
    ptr = ff.row(zi).data();
    f   += FX * ((3.0f*z-5.0f)*z*z+2.0f);
    ptr = ff.row(zi+1).data();
    f   += FX * (((4.0f-3.0f*z)*z+1.0f)*z);
    ptr = ff.row(zi+2).data();
    f   += FX * ((z-1.0f)*z*z);
    f   *= 0.25f;
  }
  else {
    ptr = ff.data();
    f = 0.5f * FX;    
  }
  return f;
}
#undef FX


int atomCompare(const void *atom1,const void *atom2) {
  /*  return (((*(atom *)atom1).z == (*(atom *)atom2).z) ? 0 : 
	  (((*(atom *)atom1).z > (*(atom *)atom2).z) ? 1 : -1));
  */
  /* Use the fact that z is the first element in the atom struct */
  return ((*(float_tt *)atom1 == *(float_tt *)atom2) ? 0 : 
	  ((*(float_tt *)atom1 > *(float_tt *)atom2) ? 1 : -1)); 
}
