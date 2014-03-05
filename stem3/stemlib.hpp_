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

#ifndef STEMLIB_H
#define STEMLIB_H

// #define WIN

#include "stemtypes_fftw3.hpp"
#include "potential.hpp"
#include "wavefunctions.hpp"
#include "data_containers.hpp"


/**********************************************
 * This function creates a incident STEM probe 
 * at position (dx,dy)
 * with parameters given in muls
 *********************************************/
// int probe(MULS *muls,double dx, double dy);
void probe(MULS *muls, WavePtr wave, double dx, double dy);
void probePlot(MULS *muls, WavePtr wave);

void initSTEMSlices(MULS *muls, int nlayer);
void interimWave(MULS *muls,WavePtr wave,int slice);
//void collectIntensity(MULS *muls, WavePtr wave, int slices);
//void detectorCollect(MULS *muls, WavePtr wave);
void saveSTEMImages(MULS *muls);

void createAtomBox(MULS *muls, int Znum, atomBox *aBox);
void transmit(void **wave,void **trans,int nx, int ny,int posx,int posy);
void propagate_slow(WavePtr wave,int nx, int ny,MULS *muls);

WAVEFUNC initWave(int nx, int ny);
void readStartWave(WavePtr wave);
/******************************************************************
 * runMulsSTEM() - do the multislice propagation in STEM mode
 *
 * waver, wavei are expected to contain incident wave function 
 * the will be updated at return
 *****************************************************************/
int runMulsSTEM_old(MULS *muls,int lstart);
int runMulsSTEM(MULS *muls, WavePtr wave, PotPtr pot);
void writePix(char *outFile,complex_tt **pict,MULS *muls,int iz);
void fft_normalize(void **array,int nx, int ny);
void showPotential(complex_tt ***pot,int nz,int nx,int ny,
		   double dx,double dy,double dz);
void atomBoxLookUp(complex_tt *vlu,MULS *muls,int Znum,double x,double y,
			   double z,double B);
void writeBeams(MULS *muls, WavePtr wave,int ilayer, int absolute_slice);

/***********************************************************************************
 * old image read/write functions, may soon be outdated
 */
void readRealImage_old(float_tt **pix, int nx, int ny, float_tt *t, char *fileName);
void readImage_old(complex_tt **pix, int nx, int ny, float_tt *t, char *fileName);
void writeFloat_TtImage_old(float_tt **pix, int nx, int ny, float_tt t,char *fileName);
void writeImage_old(complex_tt **pix, int nx, int ny, float_tt t,char *fileName);

#endif /* STEMLIB_H */
