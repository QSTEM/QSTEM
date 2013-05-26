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

#ifndef CUSTOMSLICE_H
#define CUSTOMSLICE_H

//#ifdef __cplusplus
//extern "C"
//{
//#endif /* __cplusplus */

void make3DSlicesFT(MULS *muls);
double **reduceAndExpand(fftw_complex **fc,int Nz,int Nx,int zOversample,int *fNz,int *fNx);

//#ifdef __cplusplus
//}
//#endif /* __cplusplus */

#endif // CUSTOMSLICE_H
CUSTOMSLICE_H
