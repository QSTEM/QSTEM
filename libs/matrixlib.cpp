#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include "matrixlib.h"
#include "stemtypes_fftw3.h"
// #include "stemlib.h"
#include "memory_fftw3.h"
// #include "nrutil.h" 

#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#ifdef WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif

#define TINY 0.0 // 1.0e-20; //A small number. 

static float sqrarg = 0.0f; 
#define IMIN(a,b) ((a) < (b) ? (a) : (b))    // minimum of 2 int values
#define FMAX(a,b) ((a) > (b) ? (a) : (b))    // maximum of a and b
#define SIGN(a,b) ((b) >= 0 ? (a) : -(a))    // Magnitude of a times sign of b.
#define SQR(a) (((sqrarg=(a)) == 0.0) ? 0.0 : sqrarg*sqrarg)

void svdcmp1(double2DArray a, int m, int n, double w[], double2DArray v) {
  /* Given a matrix a[1..m][1..n], this routine computes its singular value decomposition,
     A = U   W   V T . 
     The matrix U replaces a on output. The diagonal matrix of singular values 
     W is out- put as a vector w[1..n]. The matrix V (not the transpose V T ) is output
     as v[1..n][1..n].
  */
  // uses:  double pythag(double a, double b); 
  int flag,i,its,j,jj,k,l=0,nm=0; 
  double anorm,c,f,g,h,s,scale,x,y,z;
  double1DArray rv1 = double1D(n+1,"rv1"); 

  g=scale=anorm=0.0; // Householder reduction to bidiagonal form. 
  for (i=1;i<=n;i++) { 
    l=i+1; 
    rv1[i]=scale*g; 
    g=s=scale=0.0; 
    if (i <= m) { 
      for (k=i;k<=m;k++) 
	scale += fabs(a[k][i]); 
      if (scale) { 
	for (k=i;k<=m;k++) { 
	  a[k][i] /= scale; 
	  s += a[k][i]*a[k][i]; 
	} 
	f=a[i][i]; 
	g = -SIGN(sqrt(s),f); 
	h=f*g-s; 
	a[i][i]=f-g; 
	for (j=l;j<=n;j++) { 
	  for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j]; 
	  f=s/h; 
	  for (k=i;k<=m;k++) a[k][j] += f*a[k][i]; 
	} 
	for (k=i;k<=m;k++) a[k][i] *= scale; 
      } 
    } 
    w[i]=scale *g; 
    g=s=scale=0.0; 
    if(i<= m && i != n) { 
      for (k=l;k<=n;k++) scale += fabs(a[i][k]); 
      if (scale) {

	for (k=l;k<=n;k++) { 
	  a[i][k] /= scale; 
	  s += a[i][k]*a[i][k]; 
	} 
	f=a[i][l]; 
	g = -SIGN(sqrt(s),f); 
	h=f*g-s; 
	a[i][l]=f-g; 
	for (k=l;k<=n;k++) rv1[k]=a[i][k]/h; 
	for (j=l;j<=m;j++) { 
	  for (s=0.0,k=l;k<=n;k++) 
	    s += a[j][k]*a[i][k]; 
	  for (k=l;k<=n;k++) 
	    a[j][k] += s*rv1[k]; 
	} 
	for (k=l;k<=n;k++) a[i][k] *= scale; 
      } 
    } 
    anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i]))); 
  } 
  for (i=n;i>=1;i--) { // Accumulation of right-hand transformations. 
    if(i < n) { 
      if (g) { 
	for (j=l;j<=n;j++) // Double division to avoid possible under ow. 
	  v[j][i]=(a[i][j]/a[i][l])/g; 
	for (j=l;j<=n;j++) { 
	  for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j]; 
	  for (k=l;k<=n;k++) v[k][j] += s*v[k][i]; 
	} 
      } 
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0; 
    } 
    v[i][i]=1.0; 
    g=rv1[i]; 
    l=i; 
  } 
  for (i=IMIN(m,n);i>=1;i--) { // Accumulation of left-hand transformations. 
    l=i+1; 
    g=w[i]; 
    for (j=l;j<=n;j++) a[i][j]=0.0; 
    if (g) { 
      g=1.0/g; 
      for (j=l;j<=n;j++) { 
	for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j]; 
	f=(s/a[i][i])*g; 
	for (k=i;k<=m;k++) a[k][j] += f*a[k][i]; 
      } 
      for (j=i;j<=m;j++) a[j][i] *= g; 
    } 
    else for (j=i;j<=m;j++) a[j][i]=0.0; ++a[i][i]; 
  } 
  for (k=n;k>=1;k--) { 
    // Diagonalization of the bidiagonal form: Loop over singular values,
    // and over allowed iterations. 
    for (its=1;its<=30;its++) { 
      flag=1; 
      for (l=k;l>=1;l--) { // Test for splitting. 
	nm=l-1; // Note that rv1[1] is always zero. 
	if ((fabs(rv1[l])+anorm) == anorm) { 
	  flag=0; 
	  break; 
	} 
	if ((fabs(w[nm])+anorm) == anorm) break; 
      } 
      if (flag) { 
	c=0.0; // Cancellation of rv1[l], if l>1. 
	s=1.0; 
	for (i=l;i<=k;i++) {
	  f=s*rv1[i]; 
	  rv1[i]=c*rv1[i]; 
	  if ((fabs(f)+anorm) == anorm) 
	    break; 
	  g=w[i]; 
	  h=pythag(f,g); 
	  w[i]=h; 
	  h=1.0/h; 
	  c=g*h; 
	  s = -f*h; 
	  for (j=1;j<=m;j++) { 
	    y=a[j][nm]; 
	    z=a[j][i]; 
	    a[j][nm]=y*c+z*s; 
	    a[j][i]=z*c-y*s; 
	  } 
	} 
      } 
      z=w[k]; 
      if (l == k) { //Convergence. 
	if (z < 0.0) { // Singular value is made nonnegative. 
	  w[k] = -z; 
	  for (j=1;j<=n;j++) v[j][k] = -v[j][k]; 
	} 
	break; 
      } 
      if (its == 30) { 
	printf("no convergence in 30 svdcmp iterations"); 
	exit(0);
      }
      x=w[l]; // Shiftfrom bottom 2-by-2minor. 
      nm=k-1; 
      y=w[nm]; 
      g=rv1[nm]; 
      h=rv1[k]; 
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y); 
      g=pythag(f,1.0); 
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x; 
      c=s=1.0; 
      // Next QR transformation: 
      for (j=l;j<=nm;j++) { 
	i=j+1; 
	g=rv1[i]; 
	y=w[i]; 
	h=s*g; 
	g=c*g; 
	z=pythag(f,h); 
	rv1[j]=z; 
	c=f/z; 
	s=h/z; 
	f=x*c+g*s; 
	g = g*c-x*s; 
	h=y*s; 
	y *= c; 
	for (jj=1;jj<=n;jj++) { 
	  x=v[jj][j]; 
	  z=v[jj][i]; 
	  v[jj][j]=x*c+z*s; 
	  v[jj][i]=z*c-x*s; 
	} 
	z=pythag(f,h); 
	w[j]=z; 
	// Rotation can be arbitrary if z = 0. 
	if (z) { 
	  z=1.0/z; 
	  c=f*z; 
	  s=h*z; 
	} 
	f=c*g+s*y; 
	x=c*y-s*g;
	for (jj=1;jj<=m;jj++) { 
	  y=a[jj][j]; 
	  z=a[jj][i]; 
	  a[jj][j]=y*c+z*s; 
	  a[jj][i]=z*c-y*s; 
	} 
      } 
      rv1[l]=0.0; 
      rv1[k]=f; 
      w[k]=x; 
    } 
  } 
}


double pythag(double a, double b) {
  /* Computes (a 2 + b 2 ) 1=2 without destructive under ow or over ow.
   */ 
  double absa,absb; 
  absa=fabs(a); 
  absb=fabs(b); 
  if (absa > absb) 
    return absa*sqrt(1.0+SQR(absb/absa)); 
  else 
    return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))); 
}


/* This function will calculate the 3-dimensional vector cross product 
 * c = [a x b]
 */
void crossProduct(const double *a, const double *b, double *c) {
  c[0] = a[1]*b[2]-a[2]*b[1];
  c[1] = a[2]*b[0]-a[0]*b[2];
  c[2] = a[0]*b[1]-a[1]*b[0];
  return;
}

double dotProduct(const double *a, const double *b) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

/* vector difference c = a - b 
 * If revFlag < 0 we will treat the first parameter's (a) coordinates
 * in reversed order, i.e. a[0]=z, a[1]=y, a[2]=x.
 */

void vectDiff_f(float *a, double *b, double *c,int revFlag) {
  if (revFlag > 0) {
    c[0] = a[0]-b[0];
    c[1] = a[1]-b[1];
    c[2] = a[2]-b[2];
  }
  else {
    c[0] = a[2]-b[0];
    c[1] = a[1]-b[1];
    c[2] = a[0]-b[2];
  }
}




void showMatrix(double **M,int Nx, int Ny,char *name) {
  int i,j;

  printf("%s:\n",name);
  for (i=0;i<Nx;i++) {
    for (j=0;j<Ny;j++) printf("%10f ",M[i][j]);
    printf("\n");
  }
}

/* This function will determine whether the point (point) lies
 * above or below the plane (p)
 * The plane defining vectors vect1, vect2, norm, and point will be used
 * revFlag indicates whether the parameter point is given in 
 * forward(1) (x,y,z), or reversed(-1) (z,y,x) order.
 * This is important for using the reversed order in the atom struct.
 */
double findLambda(plane *p, float *point, int revFlag) {
	double2DArray M;
	double2DArray Minv;
	double1DArray diff;
	double lambda; /* dummy variable */

	M = double2D(3,3,"M");
    Minv = double2D(3,3,"Minv");
    diff = double1D(3,"diff");

	// TODO: Eigen will be fantastic here.

	/*
	printf("hello x=(%g %g %g)\n",p->pointX,p->pointY,p->pointZ);
	printf("hello p->norm=(%g %g %g), (%d %d %d)\n",p->normX,p->normY,p->normZ,
		(int)M[0],(int)M[1],(int)M[2]);
	*/

	M[0][0] = -((*p).normX); 
	M[0][3] = -((*p).normY); 
	M[0][6] = -((*p).normZ);
	M[0][1] = p->vect1X; M[1][1] = p->vect1Y; M[2][1] = p->vect1Z;
	M[0][2] = p->vect2X; M[1][2] = p->vect2Y; M[2][2] = p->vect2Z;

	// Eigen
	vectDiff_f(point,&(p->pointX),diff,revFlag);

	// Eigen
	M.invert(Minv);
	//inverse_3x3 (Minv,M);

	// Eigen
	lambda = dotProduct(Minv[0],diff);
	// printf("lambda: %g\n",lambda);

	return lambda;
}




/* This function will perform a 3D rotation about the angles
 * phi_x,phi_y,phi_z (in that sequence and about the respective orthogonal axes)
 * of the vector vectIn, and store the result in vectOut.
 * The angles will be written in a matrix, which will be saved for the next rotation,
 * in case the angles don't change.
 * The same vector can be specified as input, as well as output vector.
 */
void rotateVect(double *vectIn,double *vectOut, double phi_x, double phi_y, double phi_z) {
	// Christoph statically allocated these to avoid reallocation for multiple rotations.  
		// Is there a major performance difference without it?
  //static double **Mrot = NULL;
  double2DArray Mrot;
  double1DArray vectOutTemp;
  //static double *vectOutTemp = NULL;
  static double sphi_x=0, sphi_y=0, sphi_z=0;
  //  static double *vectOut = NULL;
  // printf("angles: %g %g %g\n",phi_x,phi_y,phi_z);

  //if (Mrot == NULL) {
    Mrot = double2D(3,3,"Mrot");
    //memset(Mrot[0],0,9*sizeof(double));
    Mrot[0][0] = 1;Mrot[1][1] = 1;Mrot[2][2] = 1;
    vectOutTemp = double1D(3,"vectOutTemp");
  //}
  if ((phi_x!=sphi_x) || (phi_y!=sphi_y) || (phi_z!=sphi_z)) {
    Mrot[0][0] = cos(phi_z)*cos(phi_y);
    Mrot[0][1] = cos(phi_z)*sin(phi_y)*sin(phi_x)-sin(phi_z)*cos(phi_x);
    Mrot[0][2] = cos(phi_z)*sin(phi_y)*cos(phi_x)+sin(phi_z)*sin(phi_x);

    Mrot[1][0] = sin(phi_z)*cos(phi_y);
    Mrot[1][1] = sin(phi_z)*sin(phi_y)*sin(phi_x)+cos(phi_z)*cos(phi_x);
    Mrot[1][2] = sin(phi_z)*sin(phi_y)*cos(phi_x)-cos(phi_z)*sin(phi_x);

    Mrot[2][0] = -sin(phi_y);
    Mrot[2][1] = cos(phi_y)*sin(phi_x);
    Mrot[2][2] = cos(phi_y)*cos(phi_x);

    sphi_x = phi_x;
    sphi_y = phi_y;
    sphi_z = phi_z;
  }
  vectOutTemp[0] = Mrot[0][0]*vectIn[0]+Mrot[0][1]*vectIn[1]+Mrot[0][2]*vectIn[2];
  vectOutTemp[1] = Mrot[1][0]*vectIn[0]+Mrot[1][1]*vectIn[1]+Mrot[1][2]*vectIn[2];
  vectOutTemp[2] = Mrot[2][0]*vectIn[0]+Mrot[2][1]*vectIn[1]+Mrot[2][2]*vectIn[2];
  memcpy(vectOut,vectOutTemp.data(),3*sizeof(double));

  return;
}

void rotateMatrix(double *matrixIn,double2DArray matrixOut, double phi_x, double phi_y, double phi_z) 
{
	int i,j,k;
	static double2DArray Mrot;
	static double2DArray matrixOutTemp;
	static double sphi_x=0, sphi_y=0, sphi_z=0;
	// static double *vectOut = NULL;
	// printf("angles: %g %g %g\n",phi_x,phi_y,phi_z);

	Mrot = double2D(3,3,"Mrot");
	memset(Mrot.data(),0,9*sizeof(double));
	Mrot[0][0] = 1;Mrot[1][1] = 1;Mrot[2][2] = 1;
	matrixOutTemp = double1D(9,"vectOutTemp");
	
	if ((phi_x!=sphi_x) || (phi_y!=sphi_y) || (phi_z!=sphi_z)) 
	{
		Mrot[0][0] = cos(phi_z)*cos(phi_y);
		Mrot[0][1] = cos(phi_z)*sin(phi_y)*sin(phi_x)-sin(phi_z)*cos(phi_x);
		Mrot[0][2] = cos(phi_z)*sin(phi_y)*cos(phi_x)+sin(phi_z)*sin(phi_x);

		Mrot[1][0] = sin(phi_z)*cos(phi_y);
		Mrot[1][1] = sin(phi_z)*sin(phi_y)*sin(phi_x)+cos(phi_z)*cos(phi_x);
		Mrot[1][2] = sin(phi_z)*sin(phi_y)*cos(phi_x)-cos(phi_z)*sin(phi_x);

		Mrot[2][0] = -sin(phi_y);
		Mrot[2][1] = cos(phi_y)*sin(phi_x);
		Mrot[2][2] = cos(phi_y)*cos(phi_x);

		sphi_x = phi_x;
		sphi_y = phi_y;
		sphi_z = phi_z;
	}
	memset(matrixOutTemp.data(),0,9*sizeof(double));
	// TODO: replace this with a simple matrix multiply using Eigen.
	for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) 
	{
		matrixOutTemp[i*3+j] += Mrot[i][k]*matrixIn[k*3+j];
	}
	matrixOut = matrixOutTemp;
	//memcpy(matrixOut,matrixOutTemp,9*sizeof(double));

return;
}

/* This function will create a basis vector set for the axes of the unit cell
 * which are defined in lengths and angles.
 * The rotationall degree of freedom will be fixed by the fact that 
 * vax = (ax,0,0)
 * vby = (by*cos(gamma),by*sin(gamma),0)
 * vcz = (c*cos(beta),(c*cos(alpha)-cos(beta)*cos(gamma))/sin(gamma),
 * c*sqrt(1-cos(alpha)^2-cos(beta)^2-cos(gamma)^2+2*cos(alpha)*cos(beta)*cos(gamma))/sin(gamma))
 * I assume that 
 * alpha = angle between by and cz, 
 * beta = angle between ax and cz, and
 * gamma = angle between by and ax.
 */
// #define SQR(x) ((x)*(x))
void makeCellVect(grainBox *grain, double *vax, double *vby, double *vcz) {
  
  if ((grain->alpha == 90) && (grain->beta == 90) && (grain->gamma == 90)) {
    printf("Orthogonal unit cell\n");
    vax[0] = grain->ax; vax[1] = 0; vax[2] = 0;  
    vby[0] = 0; vby[1] = grain->by;  vby[2] = 0;
    vcz[0] = 0; vcz[1] = 0; vcz[2] = grain->cz; 
  }
 
  else {

    vax[0] = grain->ax; vax[1] = 0; vax[2] = 0;
    
    vby[0] = grain->by*cos(grain->gamma*PI180); 
    vby[1] = grain->by*sin(grain->gamma*PI180);
    vby[2] = 0;
    
    vcz[0] = grain->cz*cos(grain->beta*PI180);
    vcz[1] = grain->cz*(cos(grain->alpha*PI180)-cos(grain->beta*PI180)*
			cos(grain->gamma*PI180))/sin(grain->gamma*PI180);
    if ((fabs(grain->alpha-90.0) < 1e-4) && (fabs(grain->beta-90.0) < 1e-4))
      vcz[2] = grain->cz;
    else  {// the general function still seems to have a bug in it.
      // printf("alpha: %g, beta: %g\n",grain->alpha-90.0,grain->beta-90.0);
      vcz[2] = grain->cz*(sqrt(1-SQR(cos(grain->alpha*PI180))-SQR(cos(grain->beta*PI180))
			       +2*cos(grain->alpha*PI180)*cos(grain->beta*PI180)*
			       cos(grain->gamma*PI180))/sin(grain->gamma*PI180));
    }
  }
}
/* just like the function above, but with angles from the MULS struct */

void makeCellVectMuls(MULS *muls, double *vax, double *vby, double *vcz) {
  
  if ((muls->cAlpha == 90) && (muls->cBeta == 90) && (muls->cGamma == 90)) {
    printf("Orthogonal unit cell\n");
    vax[0] = muls->ax; vax[1] = 0; vax[2] = 0;  
    vby[0] = 0; vby[1] = muls->by;  vby[2] = 0;
    vcz[0] = 0; vcz[1] = 0; vcz[2] = muls->c; 
  }
 
  else {

    vax[0] = muls->ax; vax[1] = 0; vax[2] = 0;
    
    vby[0] = muls->by*cos(muls->cGamma*PI180); 
    vby[1] = muls->by*cos(muls->cGamma*PI180);
    vby[2] = 0;
    
    vcz[0] = muls->c*cos(muls->cBeta*PI180);
    vcz[1] = muls->c*(cos(muls->cAlpha*PI180)-cos(muls->cBeta*PI180)*
		      cos(muls->cGamma*PI180))/sin(muls->cGamma*PI180);
    vcz[2] = muls->c*(sqrt(1-cos(muls->cAlpha*PI180)*cos(muls->cAlpha*PI180)
			     -cos(muls->cBeta*PI180)*cos(muls->cBeta*PI180)
			     +2*cos(muls->cAlpha*PI180)*cos(muls->cBeta*PI180)*
			   cos(muls->cGamma*PI180))/sin(muls->cGamma*PI180));
  }
}

double vectLength(double *vect) {
  return sqrt(vect[0]*vect[0]+vect[1]*vect[1]+vect[2]*vect[2]);
}


/* c = a*b */
void matrixProductInt(double **a,int Nxa, int Nya, int **b,int Nxb, int Nyb, double **c) {
  int i,j,k;

  if (Nya != Nxb) {
    printf("multiplyMatrix: Inner Matrix dimensions do not agree!\n");
    return;
  }

  /*
  for (i=0;i<Nxa;i++)
    for (j=0;j<Nyb;j++) {
      c[i][j] = 0.0;
      for (k=0;k<Nya;k++)
	c[i][j] += a[i][k]*b[k][j];
    }
  */
  for (i=0;i<Nxa;i++) for (j=0;j<Nyb;j++) {
      c[0][i*Nyb+j] = 0.0;
      for (k=0;k<Nya;k++) c[0][i*Nyb+j] += a[0][i*Nya+k] * (double)(b[0][k*Nyb+j]);
    }
}

