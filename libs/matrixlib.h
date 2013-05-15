#ifndef MATRIXLIB_H
#define MATRIXLIB_H

#include "stemtypes_fftw3.h"
#include "data_containers.h"

float_tt pythag(float_tt a, float_tt b);

/* vector functions:
 */
float_tt findLambda(plane &p, QSf3Arr &point, int revFlag);
void showMatrix(QSfMat &m, std::string name);
void vectDiff_f(QSf3Arr &a, QSf3Arr &b, QSf3Arr &c,int revFlag);
float_tt vectLength(QSfVec &vect);
void makeCellVect(grainBox &grain, QSf3Mat &Mm);
void makeCellVectMuls(MULS &muls, QSf3Mat &Mm);
void rotateVect(QSf3Vec &vectIn, QSf3Vec &vectOut, float_tt phi_x, float_tt phi_y, float_tt phi_z);
void rotateMatrix(QSf3Mat &matrixIn, QSf3Mat &matrixOut, float_tt phi_x, float_tt phi_y, float_tt phi_z);

#endif
