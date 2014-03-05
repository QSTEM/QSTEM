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

#ifndef MATRIXLIB_H
#define MATRIXLIB_H

#include "stemtypes_fftw3.hpp"

#define PI 3.1415926535898
#define PI180 1.7453292519943e-2

namespace QSTEM
{

void ludcmp(float_tt **a, int n, int *indx, float_tt *d);
void lubksb(float_tt **a, int n, int *indx, float_tt b[]);
float_tt det_3x3 (const float_tt *mat);
void inverse_3x3 (float_tt *res, const float_tt *a);
void trans_3x3 (float_tt *Mt, const float_tt *Ms);

// svdcmp1 uses the NR unit-offset vectors :-(
void svdcmp1(float_tt **a, int m, int n, float_tt w[], float_tt **v);
float_tt pythag(float_tt a, float_tt b);

/* vector functions:
 */
void crossProduct(const float_tt *a, const float_tt *b, float_tt *c);
float_tt dotProduct(const float_tt *a, const float_tt *b);
float_tt findLambda(plane *p, float *point, int revFlag);
void showMatrix(float_tt **M,int Nx, int Ny,char *name);
void vectDiff_f(float *a, float_tt *b, float_tt *c,int revFlag);
float_tt vectLength(float_tt *vect);
void makeCellVect(grainBox *grain, float_tt *vax, float_tt *vby, float_tt *vcz);
//void makeCellVectMuls(MULS *muls, float_tt *vax, float_tt *vby, float_tt *vcz);
void rotateVect(float_tt *vectIn,float_tt *vectOut, float_tt phi_x, float_tt phi_y, float_tt phi_z);
void rotateMatrix(float_tt *matrixIn,float_tt *matrixOut, float_tt phi_x, float_tt phi_y, float_tt phi_z);

/* |vect| */
float_tt vectLength(float_tt *vect);

/* c = a*b */
void matrixProduct(float_tt **a,int Nxa, int Nya, float_tt **b,int Nxb, int Nyb, float_tt **c);
void matrixProductInt(float_tt **a,int Nxa, int Nya, int **b,int Nxb, int Nyb, float_tt **c);

} // end namespace QSTEM

#endif
