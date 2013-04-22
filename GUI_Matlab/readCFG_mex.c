#pragma comment(lib, "user32.lib", "fftw3.lib")
#include "C:\Users\koch\Christoph\programs\stem3_src\lib\fileio_fftw3.h"
#include "C:\Users\koch\Christoph\programs\stem3_src\lib\matrixlib.h"
#include "mex.h"

#define BUF_LEN 256

#if NAN_EQUALS_ZERO
#define IsNonZero(d) ((d) != 0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d) != 0.0)
#endif

// compile as mex readCFG_mex.c D:\\Christoph\programs\lib\fileio_fftw3.c D:\\Christoph\programs\lib\memory_fftw3.c D:\\Christoph\programs\lib\matrixlib.c D:\\Christoph\programs\lib\readparams.c
// compile as mex readCFG_mex.c C:\Users\koch\Christoph\programs\stem3_src\lib\fileio_fftw3.c C:\Users\koch\Christoph\programs\stem3_src\lib\memory.c C:\Users\koch\Christoph\programs\stem3_src\lib\matrixlib.c C:\Users\koch\Christoph\programs\stem3_src\lib\readparams.c

/* The gateway routine. */
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
    char fileName[BUF_LEN];
    MULS muls;
    int dims[3],ix,iy;
    float *ptr;
    double *dptr;
    int natom = 0;
    atom *atomPtr;
    
    // Check for the proper number of arguments.
    if (nrhs < 1)
        mexErrMsgTxt("Please provide a file name.");
    if (nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");
    
    // Check data type of input argument.
    if (!(mxIsChar(prhs[0]))) {
        mexErrMsgTxt("Input array must be of type string.");
    }
    mxGetString(prhs[0],fileName,BUF_LEN);
    
  


  
    ///////////////////////////////////////////////////////
    // initialize muls
    muls.tds       = 0;
    muls.natom     = 0;
    muls.atoms     = NULL;
    muls.Znums     = NULL;
    muls.atomKinds = 0;
    muls.u2        = NULL;
    muls.u2avg     = NULL;
    muls.cubex     = 0;
    muls.cubey     = 0;
    muls.cubez     = 0;
    muls.nCellX    = 1;
    muls.nCellY    = 1;
    muls.nCellZ    = 1;
    muls.ctiltx    = 0;
    muls.ctilty    = 0;
    muls.ctiltz    = 0;
    muls.xOffset   = 0;
    muls.yOffset   = 0;

    

    // read structure:
    muls.atoms = readUnitCell(&natom,fileName,&muls,0);
    if (muls.atoms == NULL) {
        // if an error occured:
        dims[0] = 1; dims[1] = 1; dims[2] = 1;
        plhs[0] = mxCreateNumericArray(1,dims,mxSINGLE_CLASS, mxREAL);
        mexErrMsgTxt("Could not read input file.");
    }
    else {
        // mexPrintf("Found %d atoms in %s\n",natom,fileName);
        // copy the array of atom positions
        if (natom > 0) {
            dims[0] = natom; dims[1] = 7; dims[2] = 1;
            plhs[0] = mxCreateNumericArray(2,dims,mxSINGLE_CLASS, mxREAL);
            
            ptr = (float *)mxGetPr(plhs[0]);            
            for (atomPtr= muls.atoms, ix=natom-1;ix >= 0;ix--,atomPtr++) {
                // printf("Atom %d: (%f %f %f)\n",ix,atomPtr->x,atomPtr->y,atomPtr->z);
                ptr[ix        ] = atomPtr->x;
                ptr[ix+1*natom] = atomPtr->y;
                ptr[ix+2*natom] = atomPtr->z;
                ptr[ix+3*natom] = (float)atomPtr->Znum;
                ptr[ix+4*natom] = atomPtr->dw;
                ptr[ix+5*natom] = atomPtr->occ;
                ptr[ix+6*natom] = atomPtr->q;
            }
        }       
        else {
            dims[0] = 1; dims[1] = 1; dims[2] = 1;
            plhs[0] = mxCreateNumericArray(1,dims,mxSINGLE_CLASS, mxREAL);
        }
        // copy the metric matrix
        if (nlhs > 1) {
            dims[0] = 3; dims[1] =3; dims[2] = 1;
            plhs[1] = mxCreateNumericArray(2,dims,mxDOUBLE_CLASS, mxREAL);
            dptr = mxGetPr(plhs[1]);
            for (ix=0;ix<3;ix++) for (iy=0;iy<3;iy++) dptr[ix+3*iy] = muls.Mm[ix][iy];
        }

    }
  
///////////////////////////////////////////////
  
  
  return;
}
