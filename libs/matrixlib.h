#ifndef MATRIXLIB_H
#define MATRIXLIB_H

#include "stemtypes_fftw3.h"
#include "data_containers.h"

#define PI 3.1415926535898
#define PI180 1.7453292519943e-2

double pythag(double a, double b);

/* vector functions:
 */
void crossProduct(const double *a, const double *b, double *c);
double dotProduct(const double *a, const double *b);
double findLambda(plane *p, float *point, int revFlag);
void showMatrix(double **M,int Nx, int Ny,char *name);
void vectDiff_f(float *a, double *b, double *c,int revFlag);
double vectLength(double *vect);
void makeCellVect(grainBox *grain, double *vax, double *vby, double *vcz);
void makeCellVectMuls(MULS *muls, double *vax, double *vby, double *vcz);
void rotateVect(double *vectIn,double *vectOut, double phi_x, double phi_y, double phi_z);
void rotateMatrix(double *matrixIn,double *matrixOut, double phi_x, double phi_y, double phi_z);

/* |vect| */
double vectLength(double *vect);

void matrixProductInt(double **a,int Nxa, int Nya, int **b,int Nxb, int Nyb, double **c);

/* c = a*b */
template <class T> void matrixProduct(T a,int Nxa, int Nya, T b,int Nxb, int Nyb, T c) {
  int i,j,k;

  if (Nya != Nxb) {
    printf("multiplyMatrix: Inner Matrix dimensions do not agree!\n");
    return;
  }

  for (i=0;i<Nxa;i++) for (j=0;j<Nyb;j++) {
      c[0][i*Nyb+j] = 0.0;
      for (k=0;k<Nya;k++) c[0][i*Nyb+j] += a[0][i*Nya+k] * b[0][k*Nyb+j];
    }
}

/******************************************************************************
 Routine:   det_3x3
 Input:     m - matrix (3x3) address
 Output:    returns the determinant of 'm'
******************************************************************************/
template <class T> T det_3x3 (const T *mat)
{
    T det;

    det = mat[0] * (mat[4] * mat[8] - mat[7] * mat[5])
        - mat[1] * (mat[3] * mat[8] - mat[6] * mat[5])
        + mat[2] * (mat[3] * mat[7] - mat[6] * mat[4]);

    return det;
}

/******************************************************************************
 Routine:   trans_3x3
 Input:     Mstarget,Msource - matrix (3x3)
******************************************************************************/
template <class T> void trans_3x3 (T *Mt, const T *Ms)
{
  int i,j;
  for (i=0;i<2;i++) for (j=0;j<3;j++) Mt[i*3+j] = Ms[j*3+i];
  return;
}


/******************************************************************************
 Routine:   inverse_3x3
 Input:     m - matrix (3x3) address
 Output:    returns the inverse matrix of 'm'
******************************************************************************/
template <class T> void inverse_3x3 (T *res, const T *a)
{
    T det = det_3x3 (a);

    // printf("det: %g\n",det);

    if (fabs (det) < 0.0005f)
    {
        res[0] = 1.0f;
        res[1] = 0.0f;
        res[2] = 0.0f;
        res[3] = 0.0f;
        res[4] = 1.0f;
        res[5] = 0.0f;
        res[6] = 0.0f;
        res[7] = 0.0f;
        res[8] = 1.0f;
        return;
    }

    det = 1.0f / det;
    res[0] =  (a[4]*a[8] - a[5]*a[7]) * det;
    res[1] = -(a[1]*a[8] - a[7]*a[2]) * det;
    res[2] =  (a[1]*a[5] - a[4]*a[2]) * det;
      
    res[3] = -(a[3]*a[8] - a[5]*a[6]) * det;
    res[4] =  (a[0]*a[8] - a[6]*a[2]) * det;
    res[5] = -(a[0]*a[5] - a[3]*a[2]) * det;
      
    res[6] =  (a[3]*a[7] - a[6]*a[4]) * det;
    res[7] = -(a[0]*a[7] - a[6]*a[1]) * det;
    res[8] =  (a[0]*a[4] - a[1]*a[3]) * det;
}

#endif
