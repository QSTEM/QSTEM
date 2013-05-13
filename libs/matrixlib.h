#ifndef MATRIXLIB_H
#define MATRIXLIB_H

#include "stemtypes_fftw3.h"
#include "data_containers.h"

#define PI 3.1415926535898
#define PI180 1.7453292519943e-2

double pythag(double a, double b);

/* vector functions:
 */
double findLambda(plane &p, float *point, int revFlag);
void showMatrix(double **M,int Nx, int Ny,char *name);
void vectDiff_f(float *a, double *b, double *c,int revFlag);
float_tt vectLength(QSfVec &vect);
void makeCellVect(grainBox *grain, float_tt *vax, float_tt *vby, float_tt *vcz);
void makeCellVectMuls(MULS *muls, float_tt *vax, float_tt *vby, float_tt *vcz);
void rotateVect(QSfVec &vectIn, QSfVec &vectOut, float_tt phi_x, float_tt phi_y, float_tt phi_z);
void rotateMatrix(QSfMat &matrixIn, QSfMat matrixOut, float_tt phi_x, float_tt phi_y, float_tt phi_z);

/* |vect| */
double vectLength(double *vect);

#endif
