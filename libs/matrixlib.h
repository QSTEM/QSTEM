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
void showMatrix(QSfMat &m, std::string name);
void vectDiff_f(QSfVec &a, QSfVec &b, QSfVec &c,int revFlag);
float_tt vectLength(QSfVec &vect);
void makeCellVect(grainBox &grain, QSf3Mat &Mm);
void makeCellVectMuls(MULS &muls, QSf3Mat &Mm);
void rotateVect(QSf3Vec &vectIn, QSf3Vec &vectOut, float_tt phi_x, float_tt phi_y, float_tt phi_z);
void rotateMatrix(QSf3Mat &matrixIn, QSf3Mat &matrixOut, float_tt phi_x, float_tt phi_y, float_tt phi_z);

#endif
