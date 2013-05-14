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

#include <iostream>
#include <Eigen/Dense>

#define TINY 0.0 // 1.0e-20; //A small number. 

static float sqrarg = 0.0f; 
#define IMIN(a,b) ((a) < (b) ? (a) : (b))    // minimum of 2 int values
#define FMAX(a,b) ((a) > (b) ? (a) : (b))    // maximum of a and b
#define SIGN(a,b) ((b) >= 0 ? (a) : -(a))    // Magnitude of a times sign of b.
#define SQR(a) (((sqrarg=(a)) == 0.0) ? 0.0 : sqrarg*sqrarg)

float_tt pythag(float_tt a, float_tt b) {
  /* Computes (a 2 + b 2 ) 1=2 without destructive under ow or over ow.
   */ 
  float_tt absa,absb; 
  absa=fabs(a); 
  absb=fabs(b); 
  if (absa > absb) 
    return static_cast<float_tt>(absa*sqrt(1.0+SQR(absb/absa))); 
  else 
    return static_cast<float_tt>(absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))); 
}

/* vector difference c = a - b 
 * If revFlag < 0 we will treat the first parameter's (a) coordinates
 * in reversed order, i.e. a[0]=z, a[1]=y, a[2]=x.
 */

void vectDiff_f(QSf3Vec &a, QSf3Vec &b, QSf3Vec &c, int revFlag) {
  if (revFlag > 0) {
	c = a.array()-b.array();
  }
  else {
	c = a.reverse().array()-b.array();
  }
}

void showMatrix(QSfMat m, std::string name) {
  std::cout << "Here is the matrix " << name << ":\n" << m << std::endl;
}

/* This function will determine whether the point (point) lies
 * above or below the plane (p)
 * The plane defining vectors vect1, vect2, norm, and point will be used
 * revFlag indicates whether the parameter point is given in 
 * forward(1) (x,y,z), or reversed(-1) (z,y,x) order.
 * This is important for using the reversed order in the atom struct.
 */
double findLambda(plane *p, QSf3Vec point, int revFlag) {
  QSf3Mat M, Minv;
  // TODO: we need to define M's size!!
  QSf3Vec diff;
  //static double **M=NULL;
  //static double **Minv=NULL;
  //static double *diff=NULL;
  double lambda; /* dummy variable */

  /*
  printf("hello x=(%g %g %g)\n",p->pointX,p->pointY,p->pointZ);
  printf("hello p->norm=(%g %g %g), (%d %d %d)\n",p->normX,p->normY,p->normZ,
	 (int)M[0],(int)M[1],(int)M[2]);
  */

  /*
  M(0,0) = -((*p).normX); 
  M(3,0) = -((*p).normY); 
  M(6,0) = -((*p).normZ);
  M(1,0) = p->vect1X; M(1,1) = p->vect1Y; M(1,2) = p->vect1Z;
  M(2,0) = p->vect2X; M(2,1) = p->vect2Y; M(2,2) = p->vect2Z;
  vectDiff_f(point,&(p->pointX),diff,revFlag);
  */
  
  M(0,0) = -(p->norm[0]); M(0,1) = -(p->norm[1]); M(0,2) = -(p->norm[2]);
  M(1,0) = p->vect1[0]; M(1,1) = p->vect1[1]; M(1,2) = p->vect1[2];
  M(2,0) = p->vect2[0]; M(2,1) = p->vect2[1]; M(2,2) = p->vect2[2];
  vectDiff_f(point,p->point,diff,revFlag);

  Minv = M.inverse();
  //inverse_3x3 (Minv[0],M[0]);
  // showMatrix(M,3,3,"M");
  // showMatrix(Minv,3,3,"Minv");
  // showMatrix(&diff,1,3,"diff");
  // ludcmp(M,3,index,&d);
  // lubksb(M,3,index,point);

  lambda = Minv.dot(diff);
  //  lambda = dotProduct(Minv[0],diff);

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
void rotateVect(QSf3Vec vectIn, QSf3Vec vectOut, float_tt phi_x, float_tt phi_y, float_tt phi_z) {
  double sphi_x=0, sphi_y=0, sphi_z=0;
  //  static double *vectOut = NULL;
  // printf("angles: %g %g %g\n",phi_x,phi_y,phi_z);

  QSf3Mat Mrot;
  Mrot.setZero();
  Mrot(0,0) = 1; Mrot(1,1) = 1; Mrot(2,2) = 1;
  QSf3Vec vectOutTemp;

  if ((phi_x!=sphi_x) || (phi_y!=sphi_y) || (phi_z!=sphi_z)) {
    Mrot(0,0) = cos(phi_z)*cos(phi_y);
    Mrot(1,0) = cos(phi_z)*sin(phi_y)*sin(phi_x)-sin(phi_z)*cos(phi_x);
    Mrot(2,0) = cos(phi_z)*sin(phi_y)*cos(phi_x)+sin(phi_z)*sin(phi_x);

    Mrot(0,1) = sin(phi_z)*cos(phi_y);
    Mrot(1,1) = sin(phi_z)*sin(phi_y)*sin(phi_x)+cos(phi_z)*cos(phi_x);
    Mrot(2,1) = sin(phi_z)*sin(phi_y)*cos(phi_x)-cos(phi_z)*sin(phi_x);

    Mrot(0,2) = -sin(phi_y);
    Mrot(1,2) = cos(phi_y)*sin(phi_x);
    Mrot(2,2) = cos(phi_y)*cos(phi_x);

    sphi_x = phi_x;
    sphi_y = phi_y;
    sphi_z = phi_z;
  }
  vectOutTemp = vectIn*Mrot;
  //vectOutTemp[0] = Mrot(0,0)*vectIn[0]+Mrot(1,0)*vectIn[1]+Mrot(2,0)*vectIn[2];
  //vectOutTemp[1] = Mrot(0,1)*vectIn[0]+Mrot(1,1)*vectIn[1]+Mrot(2,1)*vectIn[2];
  //vectOutTemp[2] = Mrot(0,2)*vectIn[0]+Mrot(1,2)*vectIn[1]+Mrot(2,2)*vectIn[2];
  vectOut = vectOutTemp;
  //memcpy(vectOut,vectOutTemp,3*sizeof(double));

  return;
}

void rotateMatrix(QSfMat &matrixIn, QSfMat &matrixOut, double phi_x, double phi_y, double phi_z) {
QSf3Mat Mrot;
// This was static - so it may have saved some calculation.  Consider doing it back that way...
double sphi_x=0, sphi_y=0, sphi_z=0;
// static double *vectOut = NULL;
// printf("angles: %g %g %g\n",phi_x,phi_y,phi_z);

Mrot.setZero();
Mrot(0,0) = 1;Mrot(1,1) = 1;Mrot(2,2) = 1;

if ((phi_x!=sphi_x) || (phi_y!=sphi_y) || (phi_z!=sphi_z)) {
	Mrot(0,0) = static_cast<float_tt>(cos(phi_z)*cos(phi_y));
	Mrot(1,0) = static_cast<float_tt>(cos(phi_z)*sin(phi_y)*sin(phi_x)-sin(phi_z)*cos(phi_x));
	Mrot(2,0) = static_cast<float_tt>(cos(phi_z)*sin(phi_y)*cos(phi_x)+sin(phi_z)*sin(phi_x));

	Mrot(0,1) = static_cast<float_tt>(sin(phi_z)*cos(phi_y));
	Mrot(1,1) = static_cast<float_tt>(sin(phi_z)*sin(phi_y)*sin(phi_x)+cos(phi_z)*cos(phi_x));
	Mrot(2,1) = static_cast<float_tt>(sin(phi_z)*sin(phi_y)*cos(phi_x)-cos(phi_z)*sin(phi_x));

	Mrot(0,2) = static_cast<float_tt>(-sin(phi_y));
	Mrot(1,2) = static_cast<float_tt>(cos(phi_y)*sin(phi_x));
	Mrot(2,2) = static_cast<float_tt>(cos(phi_y)*cos(phi_x));

	sphi_x = phi_x;
	sphi_y = phi_y;
	sphi_z = phi_z;
}
matrixOut = matrixIn*Mrot;
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
void makeCellVect(grainBox &grain, QSf3Mat &Mm) 
{
  if ((grain.alpha == 90) && (grain.beta == 90) && (grain.gamma == 90)) {
    printf("Orthogonal unit cell\n");
    Mm(0,0) = grain.ax; Mm(1,0) = 0; Mm(2,0) = 0;  
    Mm(0,1) = 0; Mm(1,1) = grain.by;  Mm(2,1) = 0;
    Mm(0,2) = 0; Mm(1,2) = 0; Mm(2,2) = grain.cz; 
  }
 
  else {

    Mm(0,0) = grain.ax; Mm(1,0) = 0; Mm(2,0) = 0;
    
    Mm(0,1) = static_cast<float_tt>(grain.by*cos(grain.gamma*PI180)); 
    Mm(1,1) = static_cast<float_tt>(grain.by*sin(grain.gamma*PI180));
    Mm(2,1) = 0;
    
    Mm(0,2) = static_cast<float_tt>(grain.cz*cos(grain.beta*PI180));
    Mm(1,2) = static_cast<float_tt>(grain.cz*(cos(grain.alpha*PI180)-cos(grain.beta*PI180)*
			cos(grain.gamma*PI180))/sin(grain.gamma*PI180));
    if ((fabs(grain.alpha-90.0) < 1e-4) && (fabs(grain.beta-90.0) < 1e-4))
      Mm(2,2) = grain.cz;
    else  {// the general function still seems to have a bug in it.
      // printf("alpha: %g, beta: %g\n",grain.alpha-90.0,grain.beta-90.0);
      Mm(2,2) = static_cast<float_tt>(grain.cz*(sqrt(1-SQR(cos(grain.alpha*PI180))-SQR(cos(grain.beta*PI180))
			       +2*cos(grain.alpha*PI180)*cos(grain.beta*PI180)*
			       cos(grain.gamma*PI180))/sin(grain.gamma*PI180)));
    }
  }
}
/* just like the function above, but with angles from the MULS struct */

void makeCellVectMuls(MULS &muls, QSf3Mat &Mm) {
  
  if ((muls.cAlpha == 90) && (muls.cBeta == 90) && (muls.cGamma == 90)) {
    printf("Orthogonal unit cell\n");
	Mm(0,0) = muls.ax; Mm(1,0) = 0; Mm(2,0) = 0;  // ax[0], ax[1], ax[2]
    Mm(0,1) = 0; Mm(1,1) = muls.by;  Mm(2,1) = 0; // by[0], by[1], by[2]
    Mm(0,2) = 0; Mm(1,2) = 0; Mm(2,2) = muls.c;  // cz[0], cz[1], cz[2]
  }
 
  else {

	   
    Mm(0,0) = muls.ax; Mm(1,0) = 0; Mm(2,0) = 0;
    
    Mm(0,1) = static_cast<float_tt>(muls.by*cos(muls.cGamma*PI180)); 
    Mm(1,1) = static_cast<float_tt>(muls.by*sin(muls.cGamma*PI180));
    Mm(2,1) = 0;
    
	Mm(0,2) = static_cast<float_tt>(muls.c*cos(muls.cBeta*PI180));
    Mm(1,2) = static_cast<float_tt>(muls.c*(cos(muls.cAlpha*PI180)-cos(muls.cBeta*PI180)*
			cos(muls.cGamma*PI180))/sin(muls.cGamma*PI180));

    Mm(2,2) = static_cast<float_tt>(muls.c*(sqrt(1-cos(muls.cAlpha*PI180)*cos(muls.cAlpha*PI180)
			     -cos(muls.cBeta*PI180)*cos(muls.cBeta*PI180)
			     +2*cos(muls.cAlpha*PI180)*cos(muls.cBeta*PI180)*
			   cos(muls.cGamma*PI180))/sin(muls.cGamma*PI180)));
  }
}

double vectLength(double *vect) {
  return sqrt(vect[0]*vect[0]+vect[1]*vect[1]+vect[2]*vect[2]);
}


