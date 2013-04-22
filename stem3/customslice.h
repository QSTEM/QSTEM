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
