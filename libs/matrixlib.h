#ifndef MATRIXLIB_H
#define MATRIXLIB_H

#include "stemtypes_fftw3.h"

#define PI 3.1415926535898
#define PI180 1.7453292519943e-2

void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);
double det_3x3 (const double *mat);
void inverse_3x3 (double *res, const double *a);
void trans_3x3 (double *Mt, const double *Ms);

// svdcmp1 uses the NR unit-offset vectors :-(
void svdcmp1(double **a, int m, int n, double w[], double **v);
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

/* c = a*b */
void matrixProduct(double **a,int Nxa, int Nya, double **b,int Nxb, int Nyb, double **c);
void matrixProductInt(double **a,int Nxa, int Nya, int **b,int Nxb, int Nyb, double **c);

#endif
