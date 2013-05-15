/* file: fileio_fftw3.c */

#include <stdio.h>	/*  ANSI-C libraries */
#include <time.h>

#include <vector>
#include <string>
#include <algorithm>

#include "defines.h"
#include "stemtypes_fftw3.h"
#include "data_containers.h"
#include "memory_fftw3.h"	/* memory allocation routines */
#include "matrixlib.h"
#include "readparams.h"
#include "fileio_fftw3.h"

#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#ifdef WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif

#define MAX_MASS_INDEX 59
int idArraySize=0;
int *idArray = NULL;
int idArrayPtr = 0;
float_tt massArray[MAX_MASS_INDEX]={1.008f,                                         4.0026f,
6.939f,9.012f,      10.81f,12.01f,14.01f,16.00f,19.00f,20.18f,
23.00f,24.31f,      26.98f,28.09f,30.97f,32.06f,35.45f,39.95f,
39.10f,40.08f,
44.96f,47.90f,50.94f,52.00f,54.94f,55.85f,58.93f,58.71f,63.55f,65.38f,
69.72f,72.59f,74.92f,78.96f,79.90f,83.80f,
85.47f,87.62f,88.91f,91.22f,92.91f,95.94f,98.91f,10.07f,102.9f,106.4f,107.9f,112.4f,
114.8f,118.7f,121.8f,127.6f,126.9f,131.3f,
132.9054f,137.33f,138.9055f,178.49f,180.9479f};
/* so far this list goes up to Xe (from Gerthsen) .. Ta (webelememts) */
float_tt chargeTable[MAX_MASS_INDEX];



/********************************************************
* writePDB(atoms,natoms,fileName)
* This function will write the atomic coordinates in the 
* atoms array to the file <fileName> in pdb format, which
* is readable by AtomicEye
* The return value is 1 for success, or 0 for failure.
********************************************************/

int writePDB(std::vector<atom> atoms,int natoms,char *fileName,MULS *muls) {
	std::vector<std::string> elTable = ElTable::Get();
  FILE *fp;
  size_t j,i;
  std::string elem;
  float_tt ax,by,cz;
  
  if (natoms < 1) {
    printf("Atom array empty - no file written\n");
    return 1;
  }
  
  fp = fopen(fileName, "w" );
  if( fp == NULL ) {
    printf("Cannot open file %s\n",fileName);
    return 0;
  }
  
  ax = muls->cellDims[0];
  by=muls->cellDims[1];
  cz=muls->cellDims[2];
  
  fprintf(fp,"HEADER    libAtoms:Config_save_as_pdb; %d atoms\n",natoms);
  fprintf(fp,"CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1\n",
          muls->cellDims[0],muls->cellDims[1],muls->cellDims[2]);

  for (j=0;j<natoms;j++) {
    elem = elTable[atoms[j].Znum];
    fprintf(fp,"ATOM   %4d %s",j+1,elem);
    for (i=strlen(elem.c_str());i<13;i++)
      fprintf(fp," ");
    fprintf(fp,"1   ");
    fprintf(fp," %8.3f%8.3f%8.3f\n",atoms[j].pos[0],atoms[j].pos[1],atoms[j].pos[2]);
    /*
      fprintf(fp," %8.3f%8.3f%8.3f\n",atoms[j].pos[0]/muls->cellDims[0],
      atoms[j].pos[1]/muls->cellDims[1],atoms[j].pos[2]/muls->cellDims[2]);
    */
  } 
  
  fclose(fp);
  return 1;
}


////////////////////////////////////////////////////
// This function siomply removes those atoms from the list
// whose x- y- or z-position is exactly 1.0:
int removeRedundantAtoms(std::vector<atom> atoms,int natoms) {
  return 0;
}


#define MIN_EDGE_LENGTH 5.18 /* minimal allowed edge length in A
* going below this limit will crash 
* AtomEye.
*/

int writeCFG(std::vector<atom> atoms,int natoms,std::string fileName,MULS *muls) {
  std::vector<std::string> elTable = ElTable::Get();
  FILE *fp;
  int j;

  std::string elem;
  float_tt ax,by,cz;
  
  if (natoms < 1) {
    printf("Atom array empty - no file written\n");
    return 1;
  }
  
  fp = fopen(fileName.c_str(), "w" );
  if( fp == NULL ) {
    printf("Cannot open file %s\n",fileName);
    return 0;
  }
  
  // ax = muls->cellDims[0] > MIN_EDGE_LENGTH ? muls->cellDims[0] : MIN_EDGE_LENGTH;
  // by = muls->cellDims[1] > MIN_EDGE_LENGTH ? muls->cellDims[1] : MIN_EDGE_LENGTH;
  // cz = muls->cellDims[2]  > MIN_EDGE_LENGTH ? muls->cellDims[2]  : MIN_EDGE_LENGTH;
  ax = muls->cellDims[0];
  by = muls->cellDims[1];
  cz = muls->cellDims[2];
  
  fprintf(fp,"Number of particles = %d\n",natoms);
  fprintf(fp,"A = 1.0 Angstrom (basic length-scale)\n");
  fprintf(fp,"H0(1,1) = %g A\nH0(1,2) = 0 A\nH0(1,3) = 0 A\n",ax);
  fprintf(fp,"H0(2,1) = 0 A\nH0(2,2) = %g A\nH0(2,3) = 0 A\n",by);
  fprintf(fp,"H0(3,1) = 0 A\nH0(3,2) = 0 A\nH0(3,3) = %g A\n",cz);
  fprintf(fp,".NO_VELOCITY.\nentry_count = 6\n");
  printf("ax: %g, by: %g, cz: %g n: %d\n",muls->cellDims[0],muls->cellDims[1],muls->cellDims[2],natoms);
  
  
  elem = elTable[atoms[0].Znum];
  // printf("ax: %g, by: %g, cz: %g n: %d\n",muls->cellDims[0],muls->cellDims[1],muls->cellDims[2],natoms);
  fprintf(fp,"%g\n%s\n",atoms[0].Znum,elem);
  fprintf(fp,"%g %g %g %g %g %g\n",atoms[0].pos[0]/ax,atoms[0].pos[1]/by,atoms[0].pos[2]/cz,
          atoms[0].dw,atoms[0].occ,atoms[0].q);
  


  for (j=1;j<natoms;j++) {
    if (atoms[j].Znum != atoms[j-1].Znum) {
      elem = elTable[atoms[j].Znum];
      fprintf(fp,"%g\n%s\n",atoms[j].Znum,elem);
      // printf("%d: %g\n%s\n",j,2.0*atoms[j].Znum,elem);
    }
    fprintf(fp,"%g %g %g %g %g %g\n",atoms[j].pos[0]/ax,atoms[j].pos[1]/by,atoms[j].pos[2]/cz,
            atoms[j].dw,atoms[j].occ,atoms[j].q);
    // if (atoms[j].occ != 1) printf("Atom %d: occ = %g\n",j,atoms[j].occ);
  } 
  fclose(fp);
  
  return 1;
}


// write CFG file using atomic positions stored in pos, Z's in Znum and DW-factors in dw
// the unit cell is assumed to be cubic
int writeCFGFractCubic(float_tt *pos,int *Znum,float_tt *dw,int natoms,char *fileName,
                       float_tt a,float_tt b,float_tt c) {
  std::vector<std::string> elTable = ElTable::Get();
  FILE *fp;
  int j;
  std::string elem;
  float_tt ax,by,cz;
  
  if (natoms < 1) {
    printf("Atom array empty - no file written\n");
    return 1;
  }
  
  fp = fopen(fileName, "w" );
  if( fp == NULL ) {
    printf("Cannot open file %s\n",fileName);
    return 0;
  }
  
  ax = static_cast<float_tt>(a > MIN_EDGE_LENGTH ? a : MIN_EDGE_LENGTH);
  by = static_cast<float_tt>(b > MIN_EDGE_LENGTH ? b : MIN_EDGE_LENGTH);
  cz = static_cast<float_tt>(c > MIN_EDGE_LENGTH ? c : MIN_EDGE_LENGTH);
  
  fprintf(fp,"Number of particles = %d\n",natoms);
  fprintf(fp,"A = 1.0 Angstrom (basic length-scale)\n");
  fprintf(fp,"H0(1,1) = %g A\nH0(1,2) = 0 A\nH0(1,3) = 0 A\n",ax);
  fprintf(fp,"H0(2,1) = 0 A\nH0(2,2) = %g A\nH0(2,3) = 0 A\n",by);
  fprintf(fp,"H0(3,1) = 0 A\nH0(3,2) = 0 A\nH0(3,3) = %g A\n",cz);
  fprintf(fp,".NO_VELOCITY.\nentry_count = 4\n");
  printf("ax: %g, by: %g, cz: %g n: %d\n",ax,by,c,natoms);
  
  
  elem = elTable[Znum[0]];
  // printf("ax: %g, by: %g, cz: %g n: %d\n",muls->cellDims[0],muls->cellDims[1],muls->cellDims[2],natoms);
  fprintf(fp,"%g\n%s\n",Znum[0],elem);
  fprintf(fp,"%g %g %g %g\n",pos[0]*a/ax,pos[1]*b/by,pos[2]*c/cz,dw[0]);
  
  for (j=1;j<natoms;j++) {
    if (Znum[j] != Znum[j-1]) {
      elem = elTable[Znum[j]];
      fprintf(fp,"%g\n%s\n",Znum[j],elem);
      // printf("%d: %g\n%s\n",j,2.0*atoms[j].Znum,elem);
    }
    fprintf(fp,"%g %g %g %g\n",pos[3*j+0]*a/ax,pos[3*j+1]*b/by,pos[3*j+2]*c/cz,dw[j]);
  } 
  fclose(fp);
  
  return 1;
}


/*******************************************************************************
* int phononDisplacement: 
* This function will calculate the phonon displacement for a given atom i of the
* unit cell, which has been replicated to the larger cell (icx,icy,icz)
* The phonon displacement is either defined by the phonon data file, or, 
* the Einstein model, if the appropriate flags in the muls struct are set
* The displacement will be given in fractional coordinates of a single unit cell.
*
* Input parameters:
* Einstein-mode:
* need only: dw, Znum, atomCount
* atomCount: give statistics report, if 0, important only for non-Einstein mode
* maxAtom: total number of atoms (will be called first, i.e. atomCount=maxAtoms-1:-1:0)
*
* Phonon-file mode:
* ...
*
********************************************************************************/ 
//  phononDisplacement(u,muls,jChoice,icx,icy,icz,j,atoms[jChoice].dw,natom,jz);
//  j == atomCount
int phononDisplacement(QSf3Vec u, MULS &muls, int id, int icx, int icy,
                       int icz, int atomCount, float_tt dw, int maxAtom, int ZnumIndex) {
  int ix,iy,idd; // iz;
  static FILE *fpPhonon = NULL;
  static int Nk, Ns;        // number of k-vectors and atoms per primitive unit cell
  QSfVec massPrim;
  //	static float *massPrim;   // masses for every atom in primitive basis
  QSfMat omega;
  //	static float **omega;     // array of eigenvalues for every k-vector 
  //std::vector <QSfMat> EigVecs;
  QSVecOfcMat eigVecs;
  //	static fftwf_complex ***eigVecs;  // array of eigenvectors for every k-vector
  QSfMat kVecs;
  //	static float **kVecs;     // array for Nk 3-dim k-vectors
  // TODO: do these need to be forced to float_tt precision?
  QSfMat q1, q2;
  //	static float_tt **q1=NULL, **q2=NULL;
  int ik,lambda,icoord; // Ncells, nkomega;
  float_tt kR,kRi,kRr,wobble;
  float_tt ux=0, uy=0, uz=0;
  QSfVec u2;
  //	static float_tt *u2=NULL,*u2T,ux=0,uy=0,uz=0; // u2Collect=0; // Ttotal=0;
  // static double uxCollect=0,uyCollect=0,uzCollect=0;
  std::vector<int> u2Count;
  int runCount = 1,u2Size = -1;
  long iseed=0;
  QSf3Mat Mm, MmInv;
  //	double **Mm=NULL,**MmInv=NULL;
  // static double **MmOrig=NULL,**MmOrigInv=NULL;
  QSf3Vec uf, b;
  // axCell and the rest are unused - were used at one point in makeCellVectMuls.
  //static double *axCell,*byCell,*czCell;//*uf,*b;
  static float_tt wobScale = 0,sq3,scale=0;
  
  if (muls.tds == 0) return 0;

  // TODO: u2 and u2Count
  
  /***************************************************************************
   * We will give a statistics report, everytime, when we find atomCount == 0
   ***************************************************************************/
  
  if (atomCount == 0) {
    for (ix=0;ix<muls.atomKinds;ix++) {
      // u2Collect += u2[ix]/u2Count[ix];
      // uxCollect += ux/maxAtom; uyCollect += uy/maxAtom; uzCollect += uz/maxAtom;
      /*
        printf("STATISTICS: sqrt(<u^2>): %g, CM: (%g %g %g) %d atoms, wob: %g\n"
        "                         %g, CM: (%g %g %g) run %d\n",
        sqrt(u2/u2Count),ux/u2Count,uy/u2Count,uz/u2Count,u2Count,scale*sqrt(dw*wobScale),
        sqrt(u2Collect/runCount),uxCollect/runCount,uyCollect/runCount,uzCollect/runCount,
        runCount);
      */
      // printf("Count: %d %g\n",u2Count[ix],u2(ix]);
      u2(ix) /= u2Count[ix];
      if (runCount > 0) 
        muls.u2avg(ix) = sqrt(((runCount-1)*(muls.u2avg(ix)*muls.u2avg(ix))+u2(ix))/runCount);
      else
        muls.u2avg(ix) = sqrt(u2(ix));
      
      muls.u2(ix)    = sqrt(u2(ix));
      
      u2(ix)=0; u2Count[ix]=0;
    }
    runCount++;
    ux=0; uy=0; uz=0; 
  }
  
  // memcpy(Mm[0],muls.Mm[0],3*3*sizeof(float_tt));
  // We need to copy the transpose of muls.Mm to Mm.
  // we therefore cannot use the following command:
  // memcpy(Mm[0],muls.Mm[0],3*3*sizeof(float_tt));
  // or Mm = muls.Mm;
  Mm = muls.Mm.transpose();
  //for (ix=0;ix<3;ix++) for (iy=0;iy<3;iy++) Mm(iy,ix)=muls.Mm(ix,iy);
  // makeCellVectMuls(muls, axCell, byCell, czCell);
  // memcpy(MmOrig[0],Mm[0],3*3*sizeof(float_tt));
  // TODO: do we still need MMInv?

  MmInv = Mm.inverse();
  //inverse_3x3(MmInv[0],Mm[0]);

  if (ZnumIndex >= u2Size) 
{
    u2 = QSfVec(ZnumIndex+1);
    u2Count = std::vector<int>(ZnumIndex+1);
    // printf("created phonon displacements %d!\n",ZnumIndex);								
    //    u2 = (float_tt *)realloc(u2,(ZnumIndex+1)*sizeof(float_tt));
    //    u2Count = (int *)realloc(u2Count,(ZnumIndex+1)*sizeof(int));
    
    // printf("%d .... %d\n",ZnumIndex,u2Size);
    if (u2Size < 1) {
      for (ix=0;ix<=ZnumIndex;ix++) {
        u2(ix) = 0;
        u2Count[ix] = 0;
      }
    }
    else {
      for (ix=u2Size;ix<=ZnumIndex;ix++) {
        u2(ix) = 0;
        u2Count[ix] = 0;
      }
    }						  
    // printf("%d ..... %d\n",ZnumIndex,u2Size);
    
    u2Size = ZnumIndex+1;
  }  
  
  
  /***************************************************************************
   * Thermal Diffuse Scattering according to accurate phonon-dispersion or 
   * just Einstein model
   *
   * Information in the phonon file will be stored in binary form as follows:
   * Nk (number of k-points: 32-bit integer)
   * Ns (number of atomic species 32-bit integer)
   * M_1 M_2 ... M_Ns  (32-bit floats)
   * kx(1) ky(1) kz(1) (32-bit floats)
   * w_1(1) q_11 q_21 ... q_(3*Ns)1    (32-bit floats)
   * w_2(1) q_12 q_22 ... q_(3*Ns)2
   * :
   * w_(3*Ns)(1) q_1Ns q_2Ns ... q_(3*Ns)Ns
   * kx(2) ky(2) kz(2) 
   * :
   * kx(Nk) ky(Nk) kz(Nk)
   * :
   * 
   *
   * Note: only k-vectors in half of the Brioullin zone must be given, since 
   * w(k) = w(-k)  
   * also: 2D arrays will be read slowly varying index = first index (i*m+j)
   **************************************************************************/
  
  if (wobScale == 0) {
    wobScale = static_cast<float_tt> (1.0/(8*PID*PID));   
    sq3 = static_cast<float_tt> (1.0/sqrt(3.0));  /* sq3 is an additional needed factor which stems from
                           * int_-infty^infty exp(-x^2/2) x^2 dx = sqrt(pi)
                           * introduced in order to match the wobble factor with <u^2>
                           */
    scale = static_cast<float_tt> (sqrt(muls.tds_temp/300.0)) ;
    iseed = -(long)(time(NULL));
  }
  
  
  if ((muls.Einstein == 0) && (fpPhonon == NULL)) {
    if ((fpPhonon = fopen(muls.phononFile.c_str(),"r")) == NULL) {
      printf("Cannot find phonon mode file, will use random displacements!");
      muls.Einstein = 1;
      //muls.phononFile = NULL;
    }
    else {
      
      if (2*sizeof(float) != sizeof(fftwf_complex)) {
        printf("phononDisplacement: data type mismatch: fftw_complex != 2*float!\n");
        exit(0);
      }
      fread(&Nk,sizeof(int),1,fpPhonon);
      fread(&Ns,sizeof(int),1,fpPhonon);
      // TODO: ok if precision if 1, will break if precision is 2.
      massPrim = QSfVec(Ns);
      //      massPrim =(float *)malloc(Ns*sizeof(float_tt));  // masses for every atom in primitive basis
	  // TODO: verify file read
      fread(massPrim.data(),sizeof(float_tt),Ns,fpPhonon);
      kVecs = QSfMat(Nk, 3);
      omega = QSfMat(Nk, 3*Ns);
      //      kVecs = float32_2D(Nk,3,"kVecs");
      //      omega = float32_2D(Nk,3*Ns,"omega");          
      /* array of eigenvalues for every k-vector 
       * omega is given in THz, but the 2pi-factor
       * is still there, i.e. f=omega/2pi
       */
      eigVecs = QSVecOfcMat(Nk);
	  for (int i=0; i<Nk; i++)
	  {
		  // TODO: comment to test
		  eigVecs[i] = QScMat(3*Ns, 3*Ns);
	  }
      //      eigVecs = complex3Df(Nk,3*Ns,3*Ns,"eigVecs"); // array of eigenvectors for every k-vector
      
	  // TODO: is the eigVecs getting filled properly without transpose?
      for (ix=0;ix<Nk;ix++) {
        fread(kVecs.data()+ix,sizeof(float),3,fpPhonon);  // k-vector
		ik=ix;
        for (iy=0;iy<3*Ns;iy++) {
			fread(omega.data()+iy+3*ix,sizeof(float),1,fpPhonon);
			//fread(&omega(iy,ix),sizeof(float_tt),1,fpPhonon);
          //fread(omega(iy,ix),sizeof(float_tt),1,fpPhonon);
          //fread(eigVecs(iy,ix),2*sizeof(float),3*Ns,fpPhonon);
			//fread(&eigVecs[ik](iy, ix),2*sizeof(float),3*Ns,fpPhonon);
			fread(eigVecs[ix].row(iy).data(),2*sizeof(float),3*Ns,fpPhonon);
        }	
      }
      /*
        printf("Masses: ");
        for (ix=0;ix<Ns;ix++) printf(" %g",massPrim[ix]);
        printf("\n");
        for (ix=0;ix<3;ix++) {
        printf("(%5f %5f %5f):  ",kVecs(0,ix),kVecs(1,ix),kVecs(2,ix));
        for (idd=0;idd<Ns;idd++) for (iy=0;iy<3;iy++)
        printf("%6g ",omega(iy+3*idd,ix));
        printf("\n");
        }
      */      
      /* convert omega into q scaling factors, since we need those, instead of true omega:    */
      /* The 1/sqrt(2) term is from the dimensionality ((q1,q2) -> d=2)of the random numbers */
      for (ix=0;ix<Nk;ix++) {
        for (idd=0;idd<Ns;idd++) for (iy=0;iy<3;iy++) {
            // quantize the energy distribution:
            // tanh and exp give different results will therefore use exp
            // nkomega = (int)(1.0/tanh(THZ_HBAR_2KB*omega(iy+3*id,ix)/muls.tds_temp));
            // wobble  =      (1.0/tanh(THZ_HBAR_2KB*omega(iy+3*id,ix)/muls.tds_temp)-0.5);
            // nkomega = (int)(1.0/(exp(THZ_HBAR_KB*omega(iy+3*id,ix)/muls.tds_temp)-1)+0.5);
			if (omega(iy+3*idd,ix) > 1e-4) {
            //if (omega(iy+3*idd,ix) > 1e-4) {
              wobble = static_cast<float_tt> (muls.tds_temp>0 ? (1.0/(exp(THZ_HBAR_KB*omega(iy+3*idd,ix)/muls.tds_temp)-1)):0);
			  //wobble = muls.tds_temp>0 ? (1.0/(exp(THZ_HBAR_KB*omega(iy+3*idd,ix)/muls.tds_temp)-1)):0;
              // if (ix == 0) printf("%g: %d %g\n",omega(iy+3*id,ix),nkomega,wobble);
              wobble = static_cast<float_tt> (sqrt((wobble+0.5)/(2*PID*Nk*2*massPrim[idd]*omega(iy+3*idd,ix)* THZ_AMU_HBAR)));  
			  //wobble = sqrt((wobble+0.5)/(2*PID*Nk*2*massPrim[idd]*omega(iy+3*idd,ix)* THZ_AMU_HBAR));  
            }
            else wobble = 0;
            /* Ttotal += 0.25*massPrim[id]*((wobble*wobble)/(2*Ns))*
               omega(iy+3*id,ix)*omega(iy+3*id,ix)*AMU_THZ2_A2_KB;
            */
			//
			omega(iy+3*idd,ix) = wobble;
            //omega(iy+3*idd,ix) = wobble;
          }  // idd
        // if (ix == 0) printf("\n");
      }
      // printf("Temperature: %g K\n",Ttotal);
      // printf("%d %d %d\n",(int)(0.4*(float_tt)Nk/11.0),(int)(0.6*(float_tt)Nk),Nk);
      q1 = QSfMat(3*Ns, Nk);
      q2 = QSfMat(3*Ns, Nk);
      //      q1 = float_tt2D(3*Ns,Nk,"q1");
      //      q2 = float_tt2D(3*Ns,Nk,"q2");
      
    }
    fclose(fpPhonon);    
  }  // end of if phononfile
  
  // 
  // in the previous bracket: the phonon file is only read once.
  /////////////////////////////////////////////////////////////////////////////////////
  if ((muls.Einstein == 0) && (atomCount == maxAtom-1)) {
    if (Nk > 800)
      printf("Will create phonon displacements for %d k-vectors - please wait ...\n",Nk);
	// TODO: optimize loop out with Eigen - it is a transpose along with multiplication by scalar.
    for (lambda=0;lambda<3*Ns;lambda++) for (ik=0;ik<Nk;ik++) {
		q1(ik, lambda) = (omega(lambda,ik) * gasdev( &iseed ));
        //q1(ik,lambda) = (omega(lambda,ik) * gasdev( &iseed ));
		// set q2 equal below.
        //q2(ik,lambda) = (omega(lambda,ik) * gasdev( &iseed ));
      }
    // printf("Q: %g %g %g\n",q1(0,0),q1(8,5),q1(3,0));
  }
  q2=q1;
  /********************************************************************************
   * Do the Einstein model independent vibrations !!!
   *******************************************************************************/
  if (muls.Einstein) {	    
    /* convert the Debye-Waller factor to sqrt(<u^2>) */
    wobble = scale*sqrt(dw*wobScale);
    u(0) = (wobble*sq3 * gasdev( &iseed ));
    u(1) = (wobble*sq3 * gasdev( &iseed ));
    u(2) = (wobble*sq3 * gasdev( &iseed ));
    ///////////////////////////////////////////////////////////////////////
    // Book keeping:
    u2[ZnumIndex] += u(0)*u(0)+u(1)*u(1)+u(2)*u(2);
    ux += u(0); uy += u(1); uz += u(2);
    u2Count[ZnumIndex]++;
    
    
    /* Finally we must convert the displacement for this atom back into its fractional
     * coordinates so that we can add it to the current position in vector a
     */
    // TODO: optimize
	//  20130514 correct Eigen format, I think.
	uf = MmInv*u;
    //matrixProduct(&u,1,3,MmInv,3,3,&uf);
    // test:
    /*
      matrixProduct(&uf,1,3,Mm,3,3,&b);
      if (atomCount % 5 == 0) {
      printf("Z = %d, DW = %g, u=[%g %g %g]\n       u'=[%g %g %g]\n",muls.Znums[ZnumIndex],dw,u(0),u(1),u(2),b(0),b(1),b(2));
      showMatrix(Mm,3,3,"Mm");
      showMatrix(MmInv,3,3,"MmInv");
      }
    */
    // end test
    // Optimize?  How well does Eigen handle copies?
    //memcpy(u,uf,3*sizeof(float_tt));
	u = uf;
  }
  else {
    // id seems to be the index of the correct atom, i.e. ranges from 0 .. Natom
    printf("created phonon displacements %d, %d, %d %d (eigVecs: %d %d %d)!\n",ZnumIndex,Ns,Nk,id,Nk,3*Ns,3*Ns);
    /* loop over k and lambda:  */
    //memset(u,0,3*sizeof(float_tt));
	u.setZero();
    for (lambda=0;lambda<3*Ns;lambda++) for (ik=0;ik<Nk;ik++) {
        // if (kVecs(2,ik) == 0){
        //kR = 2*PID*(icx*kVecs(0,ik)+icy*kVecs(1,ik)+icz*kVecs(2,ik));
		// Eigen translation:
		kR = static_cast<float_tt>(2*PID*(icx*kVecs(0,ik)+icy*kVecs(1,ik)+icz*kVecs(2,ik)));
        //  kR = 2*PID*(blat(0,0)*kVecs(0,ik)+blat(1,0)*kVecs(1,ik)+blat(2,0)*kVecs(2,ik));
        kRr = cos(kR); kRi = sin(kR);
        for (icoord=0;icoord<3;icoord++) {
			// Eigen translation:
			u[icoord] += q1(ik, lambda)*(eigVecs[ik](icoord+3*id,lambda).real()*kRr-
                                       eigVecs[ik](icoord+3*id,lambda).imag()*kRi)-
					   q2(ik, lambda)*(eigVecs[ik](icoord+3*id,lambda).real()*kRi+
								       eigVecs[ik](icoord+3*id,lambda).imag()*kRr);
          //u[icoord] += q1(ik,lambda)*(eigVecs(lambda,ik)[icoord+3*id].real()*kRr-
                                       //eigVecs(lambda,ik)[icoord+3*id].imag()*kRi)-
					   //q2(ik,lambda)*(eigVecs(lambda,ik)[icoord+3*id].real()*kRi+
								       //eigVecs(lambda,ik)[icoord+3*id].imag()*kRr);
        }
      }
    // printf("u: %g %g %g\n",u(0),u(1),u(2));
    /* Convert the cartesian displacements back to reduced coordinates
     */ 
    ///////////////////////////////////////////////////////////////////////
    // Book keeping:
    u2[ZnumIndex] += u(0)*u(0)+u(1)*u(1)+u(2)*u(2);
    ux += u(0); uy += u(1); uz += u(2);
	// TODO: define ax, by, c in muls as QSf3Vec
    u(0) /= muls.cellDims[0];
    u(1) /= muls.cellDims[1];
    u(2) /= muls.cellDims[2];
    u2Count[ZnumIndex]++;
    
  } /* end of if Einstein */
  
  // printf("%d ... %d\n",ZnumIndex,u2Size);
  // printf("atomCount: %d (%d) %d %d\n",atomCount,muls.Einstein,ZnumIndex,u2Size);
  
  
  /*
  // Used for Debugging the memory leak on Aug. 20, 2010:
  if (_CrtCheckMemory() == 0) {
  printf("Found bad memory check in phononDisplacement! %d %d\n",ZnumIndex,muls.atomKinds);
  }
  */
  return 0;  
}


/***********************************************************************
* The following function returns the number of atoms in the specified
* CFG file and updates the cell parameters in the muls struct
***********************************************************************/
int readDATCellParams(MULS &muls, QSf3Mat Mm, std::string fileName) {
  int ncoord;
  char buf[256];
  float_tt a,b,c,alpha=90.0,beta=90.0,gamma=90.0;
  
  // printf("Paramter File pointer (1): %d\n",(int)getFp());
  parFpPush(); /* push the current parameter file file pointer 
                * on the stack to make this function totaly 
                * transparent 
                */
  if (!parOpen(fileName.c_str())) {
    printf("Could not open CFG input file %s\n",fileName);
    parFpPull();  /* restore old parameter file pointer */
    return 0;
  }
  // printf("Paramter File pointer (2): %d\n",(int)getFp());
  
  resetParamFile();  
  setComment('#');  
  if (readparam("Number of atoms =",buf,1)) sscanf(buf,"%d",&ncoord);
  
  if (readparam("a =",buf,1)) sscanf(buf,"%lf",&a);
  if (readparam("b =",buf,1)) sscanf(buf,"%lf",&b);
  if (readparam("c =",buf,1)) sscanf(buf,"%lf",&c);
  
  if (readparam("alpha =",buf,1)) sscanf(buf,"%lf",&alpha);
  if (readparam("beta =",buf,1)) sscanf(buf,"%lf",&beta);
  if (readparam("gamma =",buf,1)) sscanf(buf,"%lf",&gamma);
  
  setComment('%');
  parClose();   
  parFpPull();  /* restore old parameter file pointer */
  
  
  muls.cellDims[0] = a;
  muls.cellDims[1] = b;
  muls.cellDims[2]  = c;
  muls.cGamma = alpha;
  muls.cBeta  = beta;
  muls.cAlpha = gamma;
  // construct the unit cell metric from the lattice parameters and angles:
  makeCellVectMuls(muls, Mm);   
  if (ncoord < 1) {
    printf("Number of atoms in CFG file not specified!\n");
    ncoord = 0;
  }
  return ncoord;
}



/***********************************************************************
* The following function returns the number of atoms in the specified
* CFG file and updates the cell parameters in the muls struct
***********************************************************************/
int readCFGCellParams(MULS &muls, QSf3Mat &Mm, std::string fileName) {
  int ncoord;
  char buf[256];
  float_tt lengthScale;
  
  // printf("Paramter File pointer (1): %d\n",(int)getFp());
  parFpPush(); /* push the current parameter file file pointer 
                * on the stack to make this function totaly 
                * transparent 
                */
  if (!parOpen(fileName.c_str())) {
    printf("Could not open CFG input file %s\n",fileName);
    parFpPull();  /* restore old parameter file pointer */
    return 0;
  }
	// printf("Paramter File pointer (2): %d\n",(int)getFp());
  
  resetParamFile();  
  setComment('#');  
  if (readparam("Number of particles =",buf,1)) sscanf(buf,"%d",&ncoord);
  if (readparam("A =",buf,1)) sscanf(buf,"%lf",&lengthScale);
  
  // TODO: verify that .data is behaving as a ptr
  if (readparam("H0(1,1) =",buf,1)) sscanf(buf,"%lf",Mm.data()+0);
  if (readparam("H0(1,2) =",buf,1)) sscanf(buf,"%lf",Mm.data()+1);
  if (readparam("H0(1,3) =",buf,1)) sscanf(buf,"%lf",Mm.data()+2);
  
  if (readparam("H0(2,1) =",buf,1)) sscanf(buf,"%lf",Mm.data()+3);
  if (readparam("H0(2,2) =",buf,1)) sscanf(buf,"%lf",Mm.data()+4);
  if (readparam("H0(2,3) =",buf,1)) sscanf(buf,"%lf",Mm.data()+5);
  
  if (readparam("H0(3,1) =",buf,1)) sscanf(buf,"%lf",Mm.data()+6);
  if (readparam("H0(3,2) =",buf,1)) sscanf(buf,"%lf",Mm.data()+7);
  if (readparam("H0(3,3) =",buf,1)) sscanf(buf,"%lf",Mm.data()+8);

  // Eigen is column-major, but our parameter reading in was just done as row-major.
  //  transpose to swap.
  // TODO: verify.
  Mm.transposeInPlace();
  /*
    if (readparam(".NO_VELOCITY.",buf,1)) noVelocityFlag = 1; 
    if (readparam("entry_count =",buf,1)) sscanf(buf,"%lf",&entryCount);
    if (!noVelocityFlag) entryCount+=3;
  */
  setComment('%');
  parClose();   
  parFpPull();  /* restore old parameter file pointer */
  
  // 
  Mm*=lengthScale;
  //for (i=0;i<9;i++) Mm(i,0) *= lengthScale;
  
  // todo: these are simple cartesian distances - should be easy to optimize.
  muls.cellDims[0] = sqrt(Mm(0,0)*Mm(0,0)+Mm(1,0)*Mm(1,0)+Mm(2,0)*Mm(2,0));
  muls.cellDims[1] = sqrt(Mm(0,1)*Mm(0,1)+Mm(1,1)*Mm(1,1)+Mm(2,1)*Mm(2,1));
  muls.cellDims[2]  = sqrt(Mm(0,2)*Mm(0,2)+Mm(1,2)*Mm(1,2)+Mm(2,2)*Mm(2,2));

  muls.cGamma = atan2(Mm(1,1),Mm(0,1));
  //muls.cGamma = atan2(Mm(1,1),Mm(0,1));
  muls.cBeta = acos(Mm(0,2)/muls.cellDims[2]);
  //muls.cBeta = acos(Mm(0,2)/muls.cellDims[2]);
  muls.cAlpha = acos(Mm(1,2)*sin(muls.cGamma)/muls.cellDims[2]+cos(muls.cBeta)*cos(muls.cGamma));
  //muls.cAlpha = acos(Mm(1,2)*sin(muls.cGamma)/muls.cellDims[2]+cos(muls.cBeta)*cos(muls.cGamma));
  muls.cGamma /= (float)PI180;
  muls.cBeta  /= (float)PI180;
  muls.cAlpha /= (float)PI180;
  if (ncoord < 1) {
    printf("Number of atoms in CFG file not specified!\n");
    ncoord = 0;
  }
  return ncoord;
}

/***********************************************************************
* The following function returns the number of atoms in the specified
* CFG file and updates the cell parameters in the muls struct
***********************************************************************/
int readCSSRCellParams(MULS &muls, QSf3Mat &Mm, std::string fileName) {
  FILE *fp;
  char s1[32],s2[32],s3[32],buf[NCMAX];
  int spaceGrp = 1,ncoord;
  
  fp = fopen(fileName.c_str(), "r" );
  if( fp == NULL ) {
    printf("Cannot open file %s\n",fileName);
    return 0;
  }
  /* read the cell parameters from first and secons line */
  ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
  sscanf( buf, " %s %s %s",s1,s2,s3);
  muls.cellDims[0] = static_cast<float_tt>(atof(s1)); 
  muls.cellDims[1] = static_cast<float_tt>(atof(s2)); 
  muls.cellDims[2] = static_cast<float_tt>(atof(s3));
  
  ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
  sscanf( buf, " %s %s %s",s1,s2,s3);
  muls.cAlpha = static_cast<float_tt>(atof(s1)); 
  muls.cBeta  = static_cast<float_tt>(atof(s2)); 
  muls.cGamma = static_cast<float_tt>(atof(s3));
  // What Eigen is providing here is the value 
  makeCellVectMuls(muls, Mm);   
  
  /* check the space group: */
  spaceGrp = atoi(strstr(buf,"SPGR =")+strlen("SPGR ="));
  if (spaceGrp != 1) {
    printf("cannot interpret space group %d\n",spaceGrp);
    exit(0);
  }
  ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
  ncoord = atoi(buf);
  fclose(fp);
  return ncoord;
}

/*******************************************************************************
* This function reads the atomic position and element data for a single atom
* from a .dat file.  The atomic positions are given in reduced coordinates.
* Calling this function several times will read one atom after another from the 
* file.
* This function will return -1, if the end of file is reached prematurely, 0 otherwise.
*******************************************************************************/

int readNextDATAtom(atom &newAtom, int flag, std::string fileName) {
  int printFlag = 0;
  static FILE *fp=NULL;
  int noVelocityFlag = 1,entryCount = 3,element = 1;
  char buf[NCMAX];
  float_tt *atomData = NULL;
  char *str,elementStr[3];
  int j;
  
  if (flag < 0) {
    if (fp != NULL) {
      parClose();   
      parFpPull();  /* restore old parameter file pointer */      
      fp = NULL;
      setComment('%');
    }
    return -1;
  }
  
  if (fp == NULL) {
    parFpPush();  /* save old parameter file pointer */      
    if (!parOpen(fileName.c_str())) {
      printf("Could not open DAT input file %s\n",fileName);
      parFpPull();  /* restore old parameter file pointer */
      return -1;
    }
    resetParamFile();  
    // advance file pointer to last (known) header line
    readparam("gamma =",buf,1);
    fp = getFp();  /* get the file pointer from the parameter file routines */
    atomData = (float_tt *)malloc(entryCount*sizeof(float_tt));
  }
  
  element = 0;
  do {
    if (fgets(buf,NCMAX,fp) == NULL) {
      if (printFlag) printf("Found end of file!\n");
      return -1;
		}
    /* check, if this is a new element name */	
    // str = strnext(buf," \t");
    if (buf != NULL) {
      memcpy(elementStr,buf,2);
	  element = getZNumber(std::string(elementStr));
    }
  } while (element == 0);
  if (printFlag) printf("Found Z=%d\n",element);
  if (element > 0) {
    str = buf+2;
    // skip leading spaces:
    while (strchr(" \t",*str) != NULL) str++; 
    for (j=0;j<entryCount;j++) {
      if (str==NULL) {
        printf("readNextCFGatom: Error: incomplete data line: >%s<\n",buf);
        return -1;
      }
      atomData[j] = static_cast<float_tt>(atof(str)); str=strnext(str," \t");
    }		
  }
  
  newAtom.Znum = element;
  newAtom.pos[0]    = static_cast<float_tt>(atomData[0]);
  newAtom.pos[1]    = static_cast<float_tt>(atomData[1]);
  newAtom.pos[2]    = static_cast<float_tt>(atomData[2]);
  newAtom.dw   = static_cast<float_tt>(0.45*28.0/(float_tt)(2.0*element));	
  newAtom.occ  = static_cast<float_tt>(1.0);
  newAtom.q    = static_cast<float_tt>(0.0);
  printf("Atom: %d (%g %g %g), occ=%g, q=%g\n",newAtom.Znum,newAtom.pos[0],newAtom.pos[1],newAtom.pos[2],newAtom.occ,newAtom.q);	  
  
  return 0;
}


/*******************************************************************************
* This function reads the atomic position and element data for a single atom
* from a .cfg file.  The atomic positions are given in reduced coordinates.
* Calling this function several times will read one atom after another from the 
* file.
* This function will return -1, if the end of file is reached prematurely, 0 otherwise.
*******************************************************************************/

int readNextCFGAtom(atom &newAtom, int flag, std::string fileName) {
	static FILE *fp=NULL;
	static int noVelocityFlag = 1,entryCount = 3,element = 1;
	char buf[NCMAX];
	float_tt *atomData = NULL;
	float_tt mass = 28;
	char *str = NULL;
	int j;

	if (flag < 0) {
		if (fp != NULL) {
			parClose();   
			parFpPull();  /* restore old parameter file pointer */      
			fp = NULL;
			setComment('%');
		}
		return -1;
	}

	if (fp == NULL) {
		parFpPush();  /* save old parameter file pointer */      
		if (!parOpen(fileName.c_str())) {
			printf("Could not open CFG input file %s\n",fileName);
			parFpPull();  /* restore old parameter file pointer */
			return -1;
		}
		resetParamFile();  
		if (readparam(".NO_VELOCITY.",buf,1)) noVelocityFlag = 1; 
		else noVelocityFlag = 0;
		if (readparam("entry_count =",buf,1)) sscanf(buf,"%d",&entryCount);
		if (!noVelocityFlag) entryCount+=3;
		fp = getFp();  /* get the file pointer from the parameter file routines */
		atomData = (float_tt *)malloc((entryCount+1)*sizeof(float_tt));
	}

	if (fgets(buf,NCMAX,fp) == NULL) return -1;
	/* check, if this is a new mass number */
	str = strnext(buf," \t");
	if ((atof(buf) >= 1.0) && ((str==NULL) || (*str == '#'))) {
		mass = static_cast<float_tt>(atof(buf));
		// printf("nV: %d, eC: %d (%g)\n",noVelocityFlag, entryCount,atof(buf));
		if (fgets(buf,NCMAX,fp) == NULL) return -1;    
		element = getZNumber(std::string(buf));
		// printf("*** found element %d (%s %d) ***\n",element,buf,strlen(buf));
		if (fgets(buf,NCMAX,fp) == NULL) return -1;
	}
	str = buf;
	// skip leading spaces:
	while (strchr(" \t",*str) != NULL) str++; 
	for (j=0;j<entryCount;j++) {
		if (str==NULL) {
			printf("readNextCFGatom: Error: incomplete data line: >%s<\n",buf);
			return -1;
		}
		atomData[j] = static_cast<float_tt>(atof(str)); str=strnext(str," \t");
	}


	newAtom.Znum = element;
	newAtom.pos[0]    = atomData[0];
	newAtom.pos[1]    = atomData[1];
	newAtom.pos[2]    = atomData[2];
	// newAtom->dw   = 0.45*28.0/((float_tt)(2*element));	
	// printf("Element: %d, mass=%g\n",element,mass);
	newAtom.dw   = static_cast<float_tt>(0.45*28.0/mass);	
	newAtom.occ  = 1.0f;
	newAtom.q    = 0.0f;
	// read the DW-factor
	if (entryCount > 3+3*(1-noVelocityFlag)) 
		newAtom.dw = atomData[3+3*(1-noVelocityFlag)];
	// read the atom's occupancy:
	if (entryCount > 4+3*(1-noVelocityFlag)) 
		newAtom.occ = atomData[4+3*(1-noVelocityFlag)];
	// read the atom's charge:
	if (entryCount > 5+3*(1-noVelocityFlag)) 
		newAtom.q = atomData[5+3*(1-noVelocityFlag)];
	// printf("Atom: %d (%g %g %g), occ=%g, q=%g\n",newAtom->Znum,newAtom->x,newAtom->y,newAtom->z,newAtom->occ,newAtom->q);	

	return 0;
}

/*******************************************************************************
* This function reads the atomic position and element data for a single atom
* from a .cssr file.  The atomic positions are given in reduced coordinates.
* Calling this function several times will read one atom after another from the 
* file.
* A ReadLine error will occur, if the end of file is reached prematurely.
*******************************************************************************/
int readNextCSSRAtom(atom &newAtom, int flag, std::string fileName) {
	static FILE *fp=NULL;
	char buf[NCMAX];
	int count;
	char element[32],s1[32],s2[32],s3[32];
	float_tt dw;

	if (flag < 0) {
		if (fp != NULL) fclose(fp);
		fp = NULL;
		return 0;
	}

	if (fp == NULL) {
		fp = fopen(fileName.c_str(), "r" );
		if( fp == NULL ) {
			printf("Cannot open file %s\n",fileName);
			exit( 0 );
		}
		// skip the header:
		ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
		ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
		ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
		ReadLine( fp, buf, NCMAX, "in ReadXYZcoord" );
	}

	ReadLine( fp, buf, NCMAX, "in ReadXYZcoord()" );
	/* for Si */
	/*    dw = 0.444; */
	sscanf(buf,"%d %s %s %s %s %d %d %d %d %d %d %d %d %lf",
		&count,element,s1,s2,s3,&count,&count,
		&count,&count,&count,&count,&count,&count,&dw);

	newAtom.pos[0] = static_cast<float_tt>(atof(s1));
	newAtom.pos[1] = static_cast<float_tt>(atof(s2));
	newAtom.pos[2] = static_cast<float_tt>(atof(s3));
	newAtom.occ = 1.0;
	newAtom.Znum = getZNumber(element); 
	newAtom.dw = dw;

	/*
	switch (newAtom->Znum) {
	case 14:  newAtom->dw = 0.45f;
	break;
	case 29:  newAtom->dw = 0.21f;
	break;
	default: newAtom->dw = 0.21f;
	}
	*/
	return 0;
}

// #define NCINMAX 500
// #define NPARAM	64    /* number of parameters */

// #undef NCINMAX 500
// #undef NPARAM	64    /* number of parameters */

/*
int readCubicCFG(float_tt **pos,float_tt **dw, int **Znums, float_tt *ax,float_tt *by,float_tt *cz, 
				 float_tt ctiltx, float_tt ctilty) {
					 atom *atoms;
					 int Natom;
					 int j;
					 boost::shared_ptr<MULS> mu = initMuls();

					 mu.atomKinds = 0;
					 mu.Znums = NULL;
					 mu.tds = 0;
					 mu.u2 = NULL;
					 mu.nCellX = 1;
					 mu.nCellY = 1;
					 mu.nCellZ = 1;
					 mu.ctiltx = ctiltx; mu.ctilty = ctilty; mu.ctiltz = 0;
					 mu.cubex = 0; mu.cubey = 0; mu.cubez = 0;

					 atoms = readUnitCell(&Natom,"demo.cfg",&mu,1);
					 *Znums = (int *)malloc(Natom*sizeof(int));
					 *pos   = (float_tt *)malloc(Natom*3*sizeof(float_tt));
					 *dw    = (float_tt *)malloc(Natom*sizeof(float_tt));
					 *ax = mu.ax;
					 *by = mu.by;
					 *cz = mu.c;
					 for (j=0;j<Natom;j++) {
						 (*Znums)[j] = atoms[j].Znum;
						 (*dw)[j] = 0; //  atoms[j].dw;
						 (*pos)[3*j+0] = atoms[j].y/mu.by;
						 (*pos)[3*j+1] = atoms[j].x/mu.ax;
						 (*pos)[3*j+2] = 1.0-atoms[j].z/mu.c;
						 if (j<10) {
							 printf("%2d: (%g %g %g) %d [%g]\n",j,(*pos)[3*j],(*pos)[3*j+1],(*pos)[3*j+2],(*Znums)[j],(*dw)[j]);
						 }
					 }
					 free(atoms);
					 return Natom;
}
*/

/*******************************************************
* This function reads in a .cssr or .pdb file and fills a 
* atom struct array with the files data
* .cssr: Cerius data format, also supported by Cerius2
* .cfg:  MD Simulations Extended Configuration format, 
* supported by AtomEye
* atomic positions have to be in FRACTIONAL coordinates!!!
******************************************************/
#define FORMAT_UNKNOWN 0
#define FORMAT_CSSR 1
#define FORMAT_CFG 2
#define FORMAT_PDB 3
#define FORMAT_XYZ 4
#define FORMAT_DAT 5

////////////////////////////////////////////////////////////////////////
// replicateUnitCell
// 
// Replicates the unit cell NcellX x NCellY x NCellZ times
// applies phonon displacement and removes vacancies and atoms appearing
// on same position:
// ncoord is the number of atom positions that has already been read.
// memory for the whole atom-array of size natom has already been allocated
// but the sites beyond natom are still empty.
// TODO: why is natom a ptr?
void replicateUnitCell(int ncoord,int &natom, MULS &muls,std::vector<atom> atoms,int handleVacancies) {
	int i,j,i2,jChoice,ncx,ncy,ncz,icx,icy,icz,jz,jCell,jequal,jVac;
	int 	atomKinds = 0;
	float_tt totOcc;
	float_tt choice,lastOcc;
	QSf3Vec u;
	//float_tt *u;
	// seed for random number generation
	static long idum = -1;

	ncx = muls.nCellX;
	ncy = muls.nCellY;
	ncz = muls.nCellZ;

	atomKinds = muls.atomKinds;
	//////////////////////////////////////////////////////////////////////////////
	// Look for atoms which share the same position:
	jVac = 0;  // no atoms have been removed yet
	for (i=ncoord-1;i>=0;) {

		////////////////
		if ((handleVacancies) && (atoms[i].Znum > 0)) {
			totOcc = atoms[i].occ;
			for (jequal=i-1;jequal>=0;jequal--) {
				// if there is anothe ratom that comes close to within 0.1*sqrt(3) A we will increase 
				// the total occupany and the counter jequal.
				if ((fabs(atoms[i].pos[0]-atoms[jequal].pos[0]) < 1e-6) && (fabs(atoms[i].pos[1]-atoms[jequal].pos[1]) < 1e-6) && (fabs(atoms[i].pos[2]-atoms[jequal].pos[2]) < 1e-6)) {
					totOcc += atoms[jequal].occ;
				}
				else break;
			} // jequal-loop
		}
		else {
			jequal = i-1;
			totOcc = 1;
			// Keep a record of the kinds of atoms we are reading
		}
		if (jequal == i-1) {
			for (jz=0;jz<atomKinds;jz++)	if (muls.Znums[jz] == atoms[i].Znum) break;
		}

		////////////////
		u.setZero();
		//memset(u,0,3*sizeof(float_tt));
		/* replicate unit cell ncx,y,z times: */
		/* We have to start with the last atoms first, because once we added the displacements 
		* to the original unit cell (icx=icy=icz=0), we cannot use those positions			
		* as unit cell coordinates for the other atoms anymore
		*/
		// printf("Will start phonon displacement (%f)\n",muls.tds,muls.temperature);
		// for (jz=0;jz<muls.atomKinds;jz++)	if (atoms[i].Znum == muls.Znums[jz]) break;

		for (icx=ncx-1;icx>=0;icx--) {
			for (icy=ncy-1;icy>=0;icy--) {
				for (icz=ncz-1;icz>=0;icz--) {
					jCell = (icz+icy*ncz+icx*ncy*ncz)*ncoord;
					j = jCell+i;
					/* We will also add the phonon displacement to the atomic positions now: */
					atoms[j].dw = atoms[i].dw;
					atoms[j].occ = atoms[i].occ;
					atoms[j].q = atoms[i].q;
					atoms[j].Znum = atoms[i].Znum; 

					// Now is the time to remove atoms that are on the same position or could be vacancies:
					// if we encountered atoms in the same position, or the occupancy of the current atom is not 1, then
					// do something about it:
					jChoice = i;
					if ((totOcc < 1) || (jequal < i-1)) { // found atoms at equal positions or an occupancy less than 1!
						// ran1 returns a uniform random deviate between 0.0 and 1.0 exclusive of the endpoint values. 
						// 
						// if the total occupancy is less than 1 -> make sure we keep this
						// if the total occupancy is greater than 1 (unphysical) -> rescale all partial occupancies!
						if (totOcc < 1.0) choice = ran1(&idum);   
						else choice = totOcc*ran1(&idum);
						// printf("Choice: %g %g %d, %d %d\n",totOcc,choice,j,i,jequal);
						lastOcc = 0;
						for (i2=i;i2>jequal;i2--) {
							atoms[jCell+i2].dw = atoms[i2].dw;
							atoms[jCell+i2].occ = atoms[i2].occ;
							atoms[jCell+i2].q = atoms[i2].q;
							atoms[jCell+i2].Znum = atoms[i2].Znum; 

							// if choice does not match the current atom:
							// choice will never be 0 or 1(*totOcc) 
							if ((choice <lastOcc) || (choice >=lastOcc+atoms[i2].occ)) {
								// printf("Removing atom %d, Z=%d\n",jCell+i2,atoms[jCell+i2].Znum);
								atoms[jCell+i2].Znum =  0;  // vacancy
								jVac++;
							}
							else {
								jChoice = i2;
							}
							lastOcc += atoms[i2].occ;
						}

						// Keep a record of the kinds of atoms we are reading
						for (jz=0;jz<atomKinds;jz++) {
							if (muls.Znums[jz] == atoms[jChoice].Znum) break;
						}
					}
					// printf("i2=%d, %d (%d) [%g %g %g]\n",i2,jequal,jz,atoms[jequal].pos[0],atoms[jequal].pos[1],atoms[jequal].pos[2]);

					// this function does nothing, if muls.tds == 0
					// if (j % 5 == 0) printf("atomKinds: %d (jz = %d, %d)\n",atomKinds,jz,atoms[jChoice].Znum);
					phononDisplacement(u,muls,jChoice,icx,icy,icz,j,atoms[jChoice].dw,natom,jz);
					// printf("atomKinds: %d (jz = %d, %d)\n",atomKinds,jz,atoms[jChoice].Znum);

					for (i2=i;i2>jequal;i2--) {
						atoms[jCell+i2].pos[0] = atoms[i2].pos[0]+icx+u(0);
						atoms[jCell+i2].pos[1] = atoms[i2].pos[1]+icy+u(1);
						atoms[jCell+i2].pos[2] = atoms[i2].pos[2]+icz+u(2);
					}
				}  // for (icz=ncz-1;icz>=0;icz--)
			} // for (icy=ncy-1;icy>=0;icy--) 
		} // for (icx=ncx-1;icx>=0;icx--)
		i=jequal;
	} // for (i=ncoord-1;i>=0;)
	if ((jVac > 0 ) &&(muls.printLevel)) printf("Removed %d atoms because of occupancies < 1 or multiple atoms in the same place\n",jVac);

}


// #define printf mexPrintf
//
// This function reads the atomic positions from fileName and also adds 
// Thermal displacements to their positions, if muls.tds is turned on.
// TODO: why is natom a ptr?
std::vector<atom> readUnitCell(int &natom, std::string fileName, MULS &muls, int handleVacancies) {
	char error_string[100];
	int printFlag = 1;
	// char buf[NCMAX], *str,element[16];
	// FILE *fp;
	// float_t alpha,beta,gamma;
	int ncoord=0,ncx,ncy,ncz,icx,icy,icz,jz;
	// float_t dw,occ,dx,dy,dz,r;
	int i,format=FORMAT_UNKNOWN,atomKinds=0;
	// char s1[16],s2[16],s3[16];
	float_tt boxXmin=0,boxXmax=0,boxYmin=0,boxYmax=0,boxZmin=0,boxZmax=0;
	float_tt boxCenterX,boxCenterY,boxCenterZ,bcX,bcY,bcZ;
	float_tt x,y,z;
	//float_tt *u = NULL;
	QSf3Vec u;
	QSf3Mat Mm;
	//float_tt **Mm = NULL;
	std::vector<atom> atoms;
	//static atom *atoms = NULL;
	static int ncoord_old = 0;
	printFlag = muls.printLevel;

	Mm.setZero();
	//memset(Mm[0],0,9*sizeof(float_tt));
	muls.Mm = Mm;
	//u = (float_tt *)malloc(3*sizeof(float_tt));

	std::string::size_type idx;
	idx = fileName.rfind('.');
	std::string extension="";

	if(idx != std::string::npos)
	{
		extension = fileName.substr(idx+1);
	}

	/* figure out, whether we have  cssr, pdb, or cfg */
	if (extension==".cssr") {
		format = FORMAT_CSSR;
		ncoord = readCSSRCellParams(muls,Mm,fileName);
	}
	if (extension==".cfg") {
		format = FORMAT_CFG;
		ncoord = readCFGCellParams(muls,Mm,fileName);
		// return readCFGUnitCell(natom,fileName,muls);
	}
	if (extension == ".dat") {
		format = FORMAT_DAT;
		ncoord = readDATCellParams(muls,Mm,fileName);
		// return readCFGUnitCell(natom,fileName,muls);
	}
	if (format == FORMAT_UNKNOWN) {
		sprintf(error_string, "readUnitCell: Cannot read anything else than .cssr, .cfg, or .dat files (%s)!\n", fileName);
		throw std::exception(error_string);
		
	}

	if (ncoord == 0) {
		sprintf(error_string, "readUnitCell: Error reading configuration file %s - ncoord =0\n", fileName);
		throw std::exception(error_string);
	}

	ncx = muls.nCellX;
	ncy = muls.nCellY;
	ncz = muls.nCellZ;

	if (printFlag) {
		printf("Lattice parameters: ax=%g by=%g cz=%g (%d atoms)\n",
			muls.cellDims[0],muls.cellDims[1],muls.cellDims[2],ncoord);

		if ((muls.cubex == 0) || (muls.cubey == 0) || (muls.cubez == 0))	 
			printf("Size of Super-lattice: ax=%g by=%g cz=%g (%d x %d x %d)\n",
			muls.cellDims[0]*ncx,muls.cellDims[1]*ncy,muls.cellDims[2]*ncz,ncx,ncy,ncz);
		else
			printf("Size of Cube: ax=%g by=%g cz=%g\n",
			muls.cubex,muls.cubey,muls.cubez);
	}
	/************************************************************
	* now that we know how many coordinates there are
	* allocate the arrays 
	************************************************************/
	/*
	printf("%d atoms (%d %d %d) (%d %d %d)\n",ncoord,ncx,ncy,ncz,
	muls.nCellX,muls.nCellY,muls.nCellZ);
	*/
	if (ncoord_old != ncoord) {
		// clear the vector of atoms
		atoms.clear();
		//if (atoms != NULL) free(atoms);
		atoms = std::vector<atom>(ncoord*ncx*ncy*ncz);
		//atoms = (atom *)malloc(ncoord*sizeof(atom)*ncx*ncy*ncz);
		ncoord_old = ncoord;
	}
	if (atoms.size()==0) {
		throw std::exception("readUnitCell: Could not allocate memory for atoms!\n");
	}
	natom = ncoord*ncx*ncy*ncz;

	/*
	if (muls.atomKinds < 1) {
	muls.atomKinds = 1;
	muls.Znums = (int *)malloc(muls.atomKinds*sizeof(int));
	}
	*/
	atomKinds = 0;

	/***********************************************************
	* Read actual Data
	***********************************************************/
	for(jz=0,i=ncoord-1; i>=0; i--) {
		switch (format) {
		case FORMAT_CFG: 
		if (readNextCFGAtom(atoms[i],0,fileName) < 0) {
			throw std::exception("readUnitCell: number of atoms does not agree with atoms in file!\n");
		}
		break;
		case FORMAT_DAT: 

		if (readNextDATAtom(atoms[i],0,fileName) < 0) {
			throw std::exception("readUnitCell: number of atoms does not agree with atoms in file!\n");
			
		}
		break;

		case FORMAT_CSSR:
		readNextCSSRAtom(atoms[i],0,fileName);
		break;
		default: throw std::exception("Unrecognized file format passed to readUnitCell function");
		}

		if((atoms[i].Znum < 1 ) || (atoms[i].Znum > NZMAX)) {
			/* for (j=ncoord-1;j>=i;j--)
			printf("%2d: %d (%g,%g,%g)\n",j,atoms[j].Znum,atoms[j].pos[0],
			atoms[j].pos[1],atoms[j].pos[2]);
			*/
			sprintf(error_string, "Error: bad atomic number %d in file %s (atom %d [%d: %g %g %g])\n",
				atoms[i].Znum,fileName,i,atoms[i].Znum,atoms[i].pos[0],atoms[i].pos[1],atoms[i].pos[2]);
			throw std::exception(error_string);
		}

		// Keep a record of the kinds of atoms we are reading
		for (jz=0;jz<atomKinds;jz++)	if (muls.Znums[jz] == atoms[i].Znum) break;
		// allocate more memory, if there is a new element
		if (jz == atomKinds) {
			atomKinds++;
			if (atomKinds > muls.atomKinds) {
				muls.atomKinds = atomKinds;
				// printf("%d kinds (%d)\n",atomKinds,atoms[i].Znum);
			}  
			muls.Znums[jz] = atoms[i].Znum;
		}


	} // for 1=ncoord-1:-1:0  - we've just read all the atoms.
	if (muls.tds) {
		muls.u2 = QSfVec(atomKinds);
		muls.u2.setZero();
		muls.u2avg = QSfVec(atomKinds);
		muls.u2avg.setZero();
	}

		////////////////////////////////////////////////////////////////
	// Close the file for further reading, and restore file pointer 
	switch (format) {
		case FORMAT_CFG: 
			readNextCFGAtom(atoms[0],-1,"");
			break;
		case FORMAT_DAT: 
			readNextDATAtom(atoms[0],-1,"");
			break;
		case FORMAT_CSSR:
			readNextCSSRAtom(atoms[0],-1,"");
			break;
	}

	// First, we will sort the atoms by position:
	if (handleVacancies) {
		std::sort (atoms.begin(), atoms.end(), atomCompareZYX());
		//qsort((void *)atoms,ncoord,sizeof(atom),);
	}


	/////////////////////////////////////////////////////////////////
	// Compute the phonon displacement and remove atoms which appear 
	// twice or have some probability to be vacancies:
	if ((muls.cubex > 0) && (muls.cubey > 0) && (muls.cubez > 0)) {
		/* at this point the atoms should have fractional coordinates */
		// printf("Entering tiltBoxed\n");
		atoms = tiltBoxed(ncoord,natom,muls,atoms,handleVacancies);
		// printf("ncoord: %d, natom: %d\n",ncoord,natom);
	}
	else {  // work in NCell mode
		// atoms are in fractional coordinates so far, we need to convert them to 
		// add the phonon displacement in this condition, because there we can 
		// actually do the correct Eigenmode treatment.
		// but we will probably just do Einstein vibrations anyway:
		replicateUnitCell(ncoord,natom,muls,atoms,handleVacancies);
		/**************************************************************
		* now, after we read all of the important coefficients, we
		* need to decide if this is workable
		**************************************************************/
		natom = ncoord*ncx*ncy*ncz;
		if (1) { // ((Mm(0,0)*Mm(1,1)*Mm(2,2) == 0) || (Mm(1,0)!=0)|| (Mm(2,0)!=0)|| (Mm(0,1)!=0)|| (Mm(2,1)!=0)|| (Mm(0,2)!=0)|| (Mm(1,2)!=0)) {
				// printf("Lattice is not orthogonal, or rotated\n");
				for(i=0;i<natom;i++) {
					/*
					x = Mm(0,0)*atoms[i].pos[0]+Mm(1,0)*atoms[i].pos[1]+Mm(2,0)*atoms[i].pos[2];
					y = Mm(0,1)*atoms[i].pos[0]+Mm(1,1)*atoms[i].pos[1]+Mm(2,1)*atoms[i].pos[2];
					z = Mm(0,2)*atoms[i].pos[0]+Mm(1,2)*atoms[i].pos[1]+Mm(2,2)*atoms[i].pos[2];
					*/

					// This converts also to cartesian coordinates
					x = Mm(0,0)*atoms[i].pos[0]+Mm(0,1)*atoms[i].pos[1]+Mm(0,2)*atoms[i].pos[2];
					y = Mm(1,0)*atoms[i].pos[0]+Mm(1,1)*atoms[i].pos[1]+Mm(1,2)*atoms[i].pos[2];
					z = Mm(2,0)*atoms[i].pos[0]+Mm(2,1)*atoms[i].pos[1]+Mm(2,2)*atoms[i].pos[2];

					atoms[i].pos[0] = x;
					atoms[i].pos[1] = y;
					atoms[i].pos[2] = z;

				}      
		}
		/**************************************************************
		* Converting to cartesian coordinates
		*************************************************************/
		else { 
			for(i=0;i<natom;i++) {
				atoms[i].pos[0] *= muls.cellDims[0]; 
				atoms[i].pos[1] *= muls.cellDims[1]; 
				atoms[i].pos[2] *= muls.cellDims[2];
			}		 
		}
		// Now we have all the cartesian coordinates of all the atoms!
		muls.cellDims[0] *= ncx;
		muls.cellDims[1] *= ncy;
		muls.cellDims[2]  *= ncz;

		/***************************************************************
		* Now let us tilt around the center of the full crystal
		*/   
			
		bcX = ncx/2.0f;
		bcY = ncy/2.0f;
		bcZ = ncz/2.0f;
		u(0) = Mm(0,0)*bcX+Mm(0,1)*bcY+Mm(0,2)*bcZ;
		u(1) = Mm(1,0)*bcX+Mm(1,1)*bcY+Mm(1,2)*bcZ;
		u(2) = Mm(2,0)*bcX+Mm(2,1)*bcY+Mm(2,2)*bcZ;
		boxCenterX = u(0);
		boxCenterY = u(1);
		boxCenterZ = u(2);
		
		// rotateVect(u,u,muls.ctiltx,muls.ctilty,muls.ctiltz);  // simply applies rotation matrix
		// boxCenterXrot = u(0); boxCenterYrot = u(1);	boxCenterZrot = u(2);
		
		// Determine the size of the (rotated) super cell
		for (icx=0;icx<=ncx;icx+=ncx) for (icy=0;icy<=ncy;icy+=ncy) for (icz=0;icz<=ncz;icz+=ncz) {
			u(0) = Mm(0,0)*(icx-bcX)+Mm(0,1)*(icy-bcY)+Mm(0,2)*(icz-bcZ);
			u(1) = Mm(1,0)*(icx-bcX)+Mm(1,1)*(icy-bcY)+Mm(1,2)*(icz-bcZ);
			u(2) = Mm(2,0)*(icx-bcX)+Mm(2,1)*(icy-bcY)+Mm(2,2)*(icz-bcZ);
			rotateVect(u,u,muls.ctiltx,muls.ctilty,muls.ctiltz);  // simply applies rotation matrix
			// x = u(0)+boxCenterXrot; y = u(1)+boxCenterYrot; z = u(2)+boxCenterZrot;
			x = u(0)+boxCenterX; y = u(1)+boxCenterY; z = u(2)+boxCenterZ;
			if ((icx == 0) && (icy == 0) && (icz == 0)) {
				boxXmin = boxXmax = x;
				boxYmin = boxYmax = y;
				boxZmin = boxZmax = z;
			}
			else {
				boxXmin = boxXmin>x ? x : boxXmin; boxXmax = boxXmax<x ? x : boxXmax; 
				boxYmin = boxYmin>y ? y : boxYmin; boxYmax = boxYmax<y ? y : boxYmax; 
				boxZmin = boxZmin>z ? z : boxZmin; boxZmax = boxZmax<z ? z : boxZmax; 
			}
		}

		// printf("(%f, %f, %f): %f .. %f, %f .. %f, %f .. %f\n",muls.cellDims[0],muls.cellDims[1],muls.cellDims[2],boxXmin,boxXmax,boxYmin,boxYmax,boxZmin,boxZmax);


		if ((muls.ctiltx != 0) || (muls.ctilty != 0) || (muls.ctiltz != 0)) {			
			for(i=0;i<(natom);i++) {

				u(0) = atoms[i].pos[0]-boxCenterX; 
				u(1) = atoms[i].pos[1]-boxCenterY; 
				u(2) = atoms[i].pos[2]-boxCenterZ; 
				rotateVect(u,u,muls.ctiltx,muls.ctilty,muls.ctiltz);  // simply applies rotation matrix
				u(0) += boxCenterX;
				u(1) += boxCenterY; 
				u(2) += boxCenterZ; 
				atoms[i].pos[0] = u(0);
				atoms[i].pos[1] = u(1); 
				atoms[i].pos[2] = u(2); 
				// boxXmin = boxXmin>u(0) ? u(0) : boxXmin; boxXmax = boxXmax<u(0) ? u(0) : boxXmax; 
				// boxYmin = boxYmin>u(1) ? u(1) : boxYmin; boxYmax = boxYmax<u(1) ? u(1) : boxYmax; 
				// boxZmin = boxZmin>u(2) ? u(2) : boxZmin; boxZmax = boxZmax<u(2) ? u(2) : boxZmax; 
			}
		} /* if tilts != 0 ... */

		for(i=0;i<(natom);i++) {
			atoms[i].pos[0]-=boxXmin; 
			atoms[i].pos[1]-=boxYmin; 
			atoms[i].pos[2]-=boxZmin; 
		}
		muls.cellDims[0] = boxXmax-boxXmin;
		muls.cellDims[1] = boxYmax-boxYmin;
		muls.cellDims[2]  = boxZmax-boxZmin;

		// printf("(%f, %f, %f): %f .. %f, %f .. %f, %f .. %f\n",muls.cellDims[0],muls.cellDims[1],muls.cellDims[2],boxXmin,boxXmax,boxYmin,boxYmax,boxZmin,boxZmax);
		/*******************************************************************
		* If one of the tilts was more than 30 degrees, we will re-assign 
		* the lattice constants ax, by, and c by boxing the sample with a box 
		******************************************************************/

	// Offset the atoms in x- and y-directions:
	// Do this after the rotation!
	if ((muls.xOffset != 0) || (muls.yOffset != 0)) {
		for(i=0;i<natom;i++) {
			atoms[i].pos[0] += muls.xOffset; 
			atoms[i].pos[1] += muls.yOffset; 
		}		 
	}
	} // end of Ncell mode conversion to cartesian coords and tilting.
	// printf("Offset: (%f, %f)\n",muls.xOffset,muls.yOffset);

	/*
	printf("Atom box: (%g %g %g) .. (%g %g %g)\n",boxXmin,boxYmin,boxZmin,boxXmax,boxYmax,boxZmax);
	printf("ax: %g, bx: %g, c: %g\n",muls.cellDims[0],muls.cellDims[1],muls.cellDims[2]);

	printf("%d atoms read from file <%s>, %d atoms in model (Tilt: x: %g mrad, y: %g mrad).\n",
	ncoord,fileName,natom,muls.ctiltx,muls.ctilty);
	*/

	// initialize vibration amplitude counters:
	//////////////////////////////////////////////////////////////////////////////////////////



	return atoms;
} // end of readUnitCell



std::vector<atom> tiltBoxed(int ncoord,int &natom, MULS &muls, std::vector<atom> atoms, int handleVacancies) {
	int atomKinds = 0;
	int iatom,jVac,jequal,jChoice,i2,ix,iy,iz,atomCount = 0, atomSize;
	QSf3Mat Mm, Mminv, MpRed, MpRedInv;
	QSf3Mat MbPrim, MbPrimInv, MmOrig, MmOrigInv;
	QSf3Vec a, aOrig, b, bfloor, blat;
	static int oldAtomSize = 0;
	//float_tt x,y,z,dx,dy,dz; 
	QSf3Vec pos, dpos;
	float_tt totOcc,lastOcc,choice;
	std::vector<atom> unitAtoms;
	atom newAtom;
	int nxmin,nxmax,nymin,nymax,nzmin,nzmax,jz;
	// static FILE *fpPhonon = NULL;
	// FILE *fpu2;
	// static int Nk, Ns;     // number of k-vectors and atoms per primitive unit cell
	// static float **bPrim;   // basis of primitive unit cell
	// static float **posPrim; // atomic positions per primitive basis
	// static float *massPrim;  // masses for every atom in primitive basis
	// static float **omega;  // array of eigenvalues for every k-vector 
	// static fftw_complex ***eigVecs;  // array of eigenvectors for every k-vector
	// static float **kVecs;    // array for Nk 3-dim k-vectors
	// static float_tt **q1=NULL, **q2=NULL;
	int Ncells;
	// float_tt kR,kRi,kRr;
	// FILE *fp;
	//static float_tt u2=0;
	//static int u2Count = 0;
	// static long iseed=0;
	// column vectors
	QSf3Vec u, uf;
	static long idum = -1;

	// if (iseed == 0) iseed = -(long) time( NULL );
	Ncells = muls.nCellX * muls.nCellY * muls.nCellZ;

	/* calculate maximum length in supercell box, which is naturally 
	* the room diagonal:

	maxLength = sqrt(muls.cubex*muls.cubex+
	muls.cubey*muls.cubey+
	muls.cubez*muls.cubez);
	*/

	//dx = 0; dy = 0; dz = 0;
	dpos.setZero();
	dpos[0] = muls.xOffset;
	dpos[1] = muls.yOffset;
	/* find the rotated unit cell vectors .. 
	* muls does still hold the single unit cell vectors in ax,by, and c
	*/
	// makeCellVectMuls(muls, axCell, byCell, czCell);
	// We need to copy the transpose of muls.Mm to Mm.
	// we therefore cannot use the following command:
	// memcpy(Mm[0],muls.Mm[0],3*3*sizeof(float_tt));
	Mm.transposeInPlace();
	//for (ix=0;ix<3;ix++) for (iy=0;iy<3;iy++) Mm(iy,ix)=muls.Mm(ix,iy);

	MmOrig = Mm.transpose();
	MmOrigInv = MmOrig.inverse();
	/* remember that the angles are in rad: */
	rotateMatrix(Mm,Mm,muls.ctiltx,muls.ctilty,muls.ctiltz);
	/*
	// This is wrong, because it implements Mrot*(Mm'):
	rotateVect(axCell,axCell,muls.ctiltx,muls.ctilty,muls.ctiltz);
	rotateVect(byCell,byCell,muls.ctiltx,muls.ctilty,muls.ctiltz);
	rotateVect(czCell,czCell,muls.ctiltx,muls.ctilty,muls.ctiltz);
	*/
	Mminv = Mm.inverse();
	//inverse_3x3(Mminv[0],Mm[0]);  // computes Mminv from Mm!
	/* find out how far we will have to go in unit of unit cell vectors.
	* when creating the supercell by checking the number of unit cell vectors 
	* necessary to reach every corner of the supercell box.
	*/
	// showMatrix(MmOrig,3,3,"Morig");
	// printf("%d %d\n",(int)Mm, (int)MmOrig);
	a.setZero();
	//memset(a[0],0,3*sizeof(float_tt));
	// matrixProduct(a,1,3,Mminv,3,3,b);
	//matrixProduct(Mminv,3,3,a,3,1,b);
	//  20130514 correct eigen format, I think.
	b = Mminv*a;
	// showMatrix(Mm,3,3,"M");
	// showMatrix(Mminv,3,3,"M");
	nxmin = nxmax = (int)floor(b(0,0)-dpos[0]); 
	nymin = nymax = (int)floor(b(1,0)-dpos[1]); 
	nzmin = nzmax = (int)floor(b(2,0)-dpos[2]);
	for (ix=0;ix<=1;ix++) for (iy=0;iy<=1;iy++)	for (iz=0;iz<=1;iz++) {
		a(0,0)=ix*muls.cubex; a(1,0)=iy*muls.cubey; a(2,0)=iz*muls.cubez;

		// matrixProduct(a,1,3,Mminv,3,3,b);
		//matrixProduct(Mminv,3,3,a,3,1,b);
		//  20130514 correct eigen format, I think.
		b = Mminv*a;

		// showMatrix(b,1,3,"b");
		// TODO: can I do comparisons elementwise with vectors?
		if (nxmin > (int)floor(b(0,0)-dpos[0])) nxmin=(int)floor(b(0,0)-dpos[0]);
		if (nxmax < (int)ceil( b(0,0)-dpos[0])) nxmax=(int)ceil( b(0,0)-dpos[0]);
		if (nymin > (int)floor(b(1,0)-dpos[1])) nymin=(int)floor(b(1,0)-dpos[1]);
		if (nymax < (int)ceil( b(1,0)-dpos[1])) nymax=(int)ceil( b(1,0)-dpos[1]);
		if (nzmin > (int)floor(b(2,0)-dpos[2])) nzmin=(int)floor(b(2,0)-dpos[2]);
		if (nzmax < (int)ceil( b(2,0)-dpos[2])) nzmax=(int)ceil( b(2,0)-dpos[2]);	  
	}

	// nxmin--;nxmax++;nymin--;nymax++;nzmin--;nzmax++;
	unitAtoms = std::vector<atom>(ncoord);
	//unitAtoms = (atom *)malloc(ncoord*sizeof(atom));

	// TODO: are we not copying the entire atoms vector (is ncoord the same
	//    as atoms.size()?)
	// std::vector<atom>
	//memcpy(unitAtoms,atoms,ncoord*sizeof(atom));

	atomSize = (1+(nxmax-nxmin)*(nymax-nymin)*(nzmax-nzmin)*ncoord);
	if (atomSize != oldAtomSize) {
		atoms = std::vector<atom>(atomSize);
		//atoms = (atom *)realloc(atoms,atomSize*sizeof(atom));
		oldAtomSize = atomSize;
	}
	// showMatrix(Mm,3,3,"Mm");
	// showMatrix(Mminv,3,3,"Mminv");
	// printf("Range: (%d..%d, %d..%d, %d..%d)\n",
	// nxmin,nxmax,nymin,nymax,nzmin,nzmax);

	atomCount = 0;  
	jVac = 0;
	for (iatom=0;iatom<ncoord;) {
		// printf("%d: (%g %g %g) %d\n",iatom,unitAtoms[iatom].pos[0],unitAtoms[iatom].pos[1],
		//   unitAtoms[iatom].pos[2],unitAtoms[iatom].Znum);
		//memcpy(newAtom,unitAtoms[iatom],sizeof(atom));
		newAtom = unitAtoms[iatom];
		for (jz=0;jz<muls.atomKinds;jz++)	if (muls.Znums[jz] == newAtom.Znum) break;
		// allocate more memory, if there is a new element
		/*
		if (jz == atomKinds) {
			atomKinds++;
			if (atomKinds > muls.atomKinds) {
				muls.Znums = (int *)realloc(muls.Znums,atomKinds*sizeof(int));
				muls.atomKinds = atomKinds;
				// printf("%d kinds (%d)\n",atomKinds,atoms[i].Znum);
			}  
			muls.Znums[jz] = newAtom.Znum;
		}
		*/
		/////////////////////////////////////////////////////
		// look for atoms at equal position
		if ((handleVacancies) && (newAtom.Znum > 0)) {
			totOcc = newAtom.occ;
			for (jequal=iatom+1;jequal<ncoord;jequal++) {
				// if there is anothe ratom that comes close to within 0.1*sqrt(3) A we will increase 
				// the total occupany and the counter jequal.
				if ((fabs(newAtom.pos[0]-unitAtoms[jequal].pos[0]) < 1e-6) && (fabs(newAtom.pos[1]-unitAtoms[jequal].pos[1]) < 1e-6) && (fabs(newAtom.pos[2]-unitAtoms[jequal].pos[2]) < 1e-6)) {
					totOcc += unitAtoms[jequal].occ;
				}
				else break;
			} // jequal-loop
		}
		else {
			jequal = iatom+1;
			totOcc = 1;
		}



		// printf("%d: %d\n",atomCount,jz);
		for (ix=nxmin;ix<=nxmax;ix++) {
			for (iy=nymin;iy<=nymax;iy++) {
				for (iz=nzmin;iz<=nzmax;iz++) {
					// atom position in cubic reduced coordinates: 
					aOrig(0,0) = ix+newAtom.pos[0]; aOrig(1,0) = iy+newAtom.pos[1]; aOrig(2,0) = iz+newAtom.pos[2];

					// Now is the time to remove atoms that are on the same position or could be vacancies:
					// if we encountered atoms in the same position, or the occupancy of the current atom is not 1, then
					// do something about it:
					// All we need to decide is whether to include the atom at all (if totOcc < 1
					// of which of the atoms at equal positions to include
					jChoice = iatom;  // This will be the atom we wil use.
					if ((totOcc < 1) || (jequal > iatom+1)) { // found atoms at equal positions or an occupancy less than 1!
						// ran1 returns a uniform random deviate between 0.0 and 1.0 exclusive of the endpoint values. 
						// 
						// if the total occupancy is less than 1 -> make sure we keep this
						// if the total occupancy is greater than 1 (unphysical) -> rescale all partial occupancies!
						if (totOcc < 1.0) choice = ran1(&idum);   
						else choice = totOcc*ran1(&idum);
						// printf("Choice: %g %g %d, %d %d\n",totOcc,choice,j,i,jequal);
						lastOcc = 0;
						for (i2=iatom;i2<jequal;i2++) {
							// atoms[atomCount].Znum = unitAtoms[i2].Znum; 
							// if choice does not match the current atom:
							// choice will never be 0 or 1(*totOcc) 
							if ((choice <lastOcc) || (choice >=lastOcc+unitAtoms[i2].occ)) {
								// printf("Removing atom %d, Z=%d\n",jCell+i2,atoms[jCell+i2].Znum);
								// atoms[atomCount].Znum =  0;  // vacancy
								jVac++;
							}
							else {
								jChoice = i2;
							}
							lastOcc += unitAtoms[i2].occ;
						}
						// printf("Keeping atom %d (%d), Z=%d\n",jChoice,iatom,unitAtoms[jChoice].Znum);
					}
					// if (jChoice != iatom) memcpy(&newAtom,unitAtoms+jChoice,sizeof(atom));
					if (jChoice != iatom) {
						for (jz=0;jz<muls.atomKinds;jz++)	if (muls.Znums[jz] == unitAtoms[jChoice].Znum) break;
					}




					// here we need to call phononDisplacement:
					// phononDisplacement(u,muls,iatom,ix,iy,iz,atomCount,atoms[i].dw,natom,atoms[i].Znum);
					if (muls.Einstein == 1) {
						// phononDisplacement(u,muls,iatom,ix,iy,iz,1,newAtom.dw,10,newAtom.Znum);
						if (muls.tds) {
							phononDisplacement(u,muls,jChoice,ix,iy,iz,1,unitAtoms[jChoice].dw,atomSize,jz);
							a(0) = aOrig(0)+u(0); a(1) = aOrig(1)+u(1); a(2) = aOrig(2)+u(2);
						}
						else {
							a(0) = aOrig(0); a(1) = aOrig(1); a(2) = aOrig(2);
						}
					}
					else {
						printf("Cannot handle phonon-distribution mode for boxed sample yet - sorry!!\n");
						exit(0);
					}
					// matrixProduct(aOrig,1,3,Mm,3,3,b);
					//  20130514 correct eigen format, I think.
					b = Mm*aOrig;
					//matrixProduct(Mm,3,3,aOrig,3,1,b);

					// if (atomCount < 2) {showMatrix(a,1,3,"a");showMatrix(b,1,3,"b");}
					// b now contains atom positions in cartesian coordinates */
					pos = b + dpos;
					//x  = b(0)+dpos[0]; 
					//y  = b(1)+dpos[1]; 
					//z  = b(2)+dpos[2]; 
					if ((pos[0] >= 0) && (pos[0] <= muls.cubex) &&
						(pos[1] >= 0) && (pos[1] <= muls.cubey) &&
						(pos[2] >= 0) && (pos[2] <= muls.cubez)) {
							// matrixProduct(a,1,3,Mm,3,3,b);
							//  20130514 correct eigen format, I think.
							b=Mm*a;
							//matrixProduct(Mm,3,3,a,3,1,b);
							atoms[atomCount].pos = b + dpos;
							//atoms[atomCount].pos[0]		= b(0)+dpos[0]; 
							//atoms[atomCount].pos[1]		= b(1)+dpos[1]; 
							//atoms[atomCount].pos[2]		= b(2)+dpos[2]; 
							atoms[atomCount].dw		= unitAtoms[jChoice].dw;
							atoms[atomCount].occ	= unitAtoms[jChoice].occ;
							atoms[atomCount].q		= unitAtoms[jChoice].q;
							atoms[atomCount].Znum	= unitAtoms[jChoice].Znum;
								
							atomCount++;	
							/*
							if (unitAtoms[jChoice].Znum > 22)
								printf("Atomcount: %d, Z = %d\n",atomCount,unitAtoms[jChoice].Znum);
							*/
					}

				} /* iz ... */
			} /* iy ... */
		} /* ix ... */
		iatom = jequal;
	} /* iatom ... */
	if (muls.printLevel > 2) printf("Removed %d atoms because of multiple occupancy or occupancy < 1\n",jVac);
	muls.cellDims[0] = muls.cubex;
	muls.cellDims[1] = muls.cubey;
	muls.cellDims[2]  = muls.cubez;
	natom = atomCount;
	// call phononDisplacement again to update displacement data:
	phononDisplacement(u,muls,iatom,ix,iy,iz,0,newAtom.dw,natom,jz);


	//free(unitAtoms);
	return atoms;
}  // end of 'tiltBoxed(...)'

/*--------------------- ReadLine() -----------------------*/
/*
read a full line from a file and 
return length of line read

to bad this looks like Pascal but its the easiest
way to read just whole line because fscanf() ignores
end of line characters

fpread = pointer to file
cMax = length of data buffer cRead
cRead = char[] buffer to read into
mesg = error message to print if not successful
*/
int ReadLine( FILE* fpRead, char* cRead, int cMax, const char *mesg )
{
	if( fgets( cRead, cMax, fpRead) == NULL ) {
		return 0;
		/*   printf("error reading input file: %s\n", mesg);
		exit( 0 );
		*/
	}
	return (int)( strlen( cRead ) );

}  /* end ReadLine */

int getZNumber(std::string element) {
	std::vector<std::string> elTable = ElTable::Get();
	std::vector<std::string>::iterator it;
   it = std::find(elTable.begin(), elTable.end(), element);
   if( it==elTable.end() )
		return 0;	
	else
		return (int) distance(elTable.begin(),it);
}

#define CHARGE 0.0

void writeFrameWork(FILE *fp,superCellBox superCell) {
	std::vector<std::string> elTable = ElTable::Get();
	int i,id = 0,newId = 1;
	float_tt charge=0.0;

	if (fp != NULL) fprintf(fp,"frame1 1 framework\n");
	for (i=0;i<superCell.natoms;i++) {
		if (i==0) {
			if (idArray == NULL) {
				idArraySize = 2;
				idArray=(int *)malloc(idArraySize*sizeof(int));
				id = 0;
				idArrayPtr = 0;
			}
			else {
				/* check, if this Znum is already somewhere in the list, and if not,
				* add it to the list, extending its size, if necessary
				*/
				for (id=0;id<=idArrayPtr;id++) if (superCell.atoms[i].Znum == idArray[id]) break;
				if (id>idArrayPtr) {
					if (idArrayPtr >= idArraySize-1) {
						idArraySize *=2;
						idArray = (int *)realloc(idArray,idArraySize*sizeof(int));
					}
					idArray[++idArrayPtr] = superCell.atoms[i].Znum;
					newId = 1;
				}           
			}
			idArray[id] = superCell.atoms[i].Znum;
		}
		else {
			newId = 0;
			if (superCell.atoms[i].Znum != superCell.atoms[i-1].Znum) {
				/* check, if this Znum is already somewhere in the list, and if not,
				* add it to the list, extending its size, if necessary
				*/
				for (id=0;id<=idArrayPtr;id++) if (superCell.atoms[i].Znum == idArray[id]) break;
				if (id>idArrayPtr) {
					if (idArrayPtr >= idArraySize-1) {
						idArraySize *=2;
						idArray = (int *)realloc(idArray,idArraySize*sizeof(int));
					}
					idArray[++idArrayPtr] = superCell.atoms[i].Znum;
					newId = 1;
				}
			}
		}
		if (newId) {
			switch (superCell.atoms[i].Znum) {
	  case 14: charge = 4.0; break;
	  case 7:  charge = -3.0; break;
	  case 8:  charge = -2.0; break;
	  case 39: charge = 3.0; break;
	  default: charge = 0.0; 
			}

			if (fp != NULL) fprintf(fp,"%d %7.3f %7.3f %7.3f %6.3f %6.3f %s\n",id+1,superCell.atoms[i].pos[0],
				superCell.atoms[i].pos[1],superCell.atoms[i].pos[2],
				// TODO: massArray should probably be a vector...
				massArray[superCell.atoms[i].Znum-1],charge,
				elTable[superCell.atoms[i].Znum]);
		}
		else
			if (fp != NULL) fprintf(fp,"%d %7.3f %7.3f %7.3f\n",id+1,superCell.atoms[i].pos[0],
			superCell.atoms[i].pos[1],superCell.atoms[i].pos[2]);
	}
}

#define A 0.0     // 975.857
#define B 43684.0   // 425379 
#define C 3.4483  // 5.04748

/* This function will write the amorphous data to the MD input file starting at 
* atoms nstart and ending just before atom nstop.   
*/
void writeAmorphous(FILE *fp,superCellBox superCell,int nstart,int nstop) {
	std::vector<std::string> elTable = ElTable::Get();
	int i,j,id;
	int *idCountArray = NULL;
	float_tt charge,b;//,x,y,z;
	QSf3Arr pos;
	// char elem[8];

	memset(chargeTable,0,MAX_MASS_INDEX*sizeof(float_tt));

	printf("amorph: %d ..%d-1\n",nstart,nstop);

	if (idArray == NULL) {
		idArraySize = 2;
		idArray=(int *)malloc(idArraySize*sizeof(int));
		id = 0;
		idArrayPtr = 0;
		idArray[0] = superCell.atoms[nstart].Znum;
	}
	else {
		/* check, if this Znum is already somewhere in the list, and if not,
		* add it to the list, extending its size, if necessary
		*/
		for (id=0;id<=idArrayPtr;id++) if (superCell.atoms[nstart].Znum == idArray[id]) break;
		if (id>idArrayPtr) {
			if (idArrayPtr >= idArraySize-1) {
				idArraySize *=2;
				idArray = (int *)realloc(idArray,idArraySize*sizeof(int));
			}
			idArray[++idArrayPtr] = superCell.atoms[nstart].Znum;
		}           
	}
	idCountArray = (int *)malloc(idArraySize*sizeof(int));
	memset(idCountArray,0,idArraySize*sizeof(int));
	idCountArray[id] = 1;

	for (i=nstart+1;i<nstop;i++) {
		if (superCell.atoms[i].Znum != superCell.atoms[i-1].Znum) {
			/* check, if this Znum is already somewhere in the list, and if not,
			* add it to the list, extending its size, if necessary
			*/

			for (id=0;id<=idArrayPtr;id++) {
				if (superCell.atoms[i].Znum == idArray[id]) break;
			}
			if (id>idArrayPtr) {
				if (idArrayPtr >= idArraySize-1) {
					idArraySize *=2;
					idArray = (int *)realloc(idArray,idArraySize*sizeof(int));
					idCountArray = (int *)realloc(idCountArray,idArraySize*sizeof(int));
				}
				idArray[++idArrayPtr] = superCell.atoms[i].Znum;
			}
		}
		idCountArray[id] += 1;
	}

	/**************************************************************************************
	* Now that we found all the species present in the amorphous phase we can go ahead and 
	* make a list of them, together with the number of atoms present of each kind
	*/
	for (id=0;id<=idArrayPtr;id++) {
		if (idCountArray[id] > 0) { 
			if (idArray[id] >= MAX_MASS_INDEX) {
				printf("mass exceeds array!!! - extend array!\n");
				exit(0);
			}
			switch (idArray[id]) {
	  case 14: charge = 4.0; break;
	  case 7:  charge = -3.0; break;
	  case 8:  charge = -2.0; break;
	  case 39: charge = 3.0; break;
	  default: charge = 0.0; 
			}
			if (fp != NULL) {
				fprintf(fp,"%c%c %d\n",elTable[2*idArray[id]-2],elTable[2*idArray[id]-1],idCountArray[id]);
				fprintf(fp,"%d %7.3f %7.3f %7.3f %6.3f %6.3f %c%c\n",id+1,0.0,0.0,0.0,
					massArray[idArray[id]-1],charge,
					elTable[2*idArray[id]-2],elTable[2*idArray[id]-1]);    
			}
		}
	}
	if (fp != NULL) fprintf(fp,"end\n");
	/***********************************************************************
	* The list of the potential parameters is next
	*/
	if (fp != NULL) {
	fprintf(fp,"buckingham\n");
	for (id=0;id<=idArrayPtr;id++)
		for (j=id;j<=idArrayPtr;j++) {
			// The following parameters are from S. Garofalini, J. Am. Cer. Soc. 67, 133 (1987)
			switch (10*(id+1) + (j+1)) {
	  case 11: b = 0.0f; break;  // N -N
	  case 12: b = 4.5f; break;  // N -Si
	  case 13: b = 4.5f; break;  // N -Y
	  case 14: b = 0.0f; break;  // N -O
	  case 22: b = 1.8770f; break;  // Si-Si
	  case 23: b = 4.4200f; break;  // Si-Y
	  case 24: b = 2.9620f; break;  // Si-O
	  case 33: b = 9.9706f; break;  // Y -Y
	  case 34: b = 6.0802f; break;  // Y -O
	  case 44: b = 0.7254f; break;  // O -O		
	  default : b=0.0f;
			}
			b = b*6.022e4f;  // *1e-16*6.022e23*1e-3
			fprintf(fp,"%d %d %g %g %g\n",id+1,j+1,A,b,C);

		}
		fprintf(fp,"end\n");

		/*************************************************************************
		* We can now write down the size of the MD cell:
		*/
		fprintf(fp,"%g %g %g 90 90 90 1 1 1\n",superCell.cellDims[0],
			superCell.cellDims[1],superCell.cellDims[2]);

		/****************************************************************************
		* now list the position of each one of the atoms in the amorphous phase
		* in FRACTIONAL COORDINATES !!!
		*/
		for(i=nstart;i<nstop;i++) {
			pos = superCell.atoms[i].pos/superCell.cellDims;
			if (fabs(pos[0]-1.0) < 1e-5) pos[0] = 0.0;
			if (fabs(pos[1]-1.0) < 1e-5) pos[1] = 0.0;
			if (fabs(pos[2]-1.0) < 1e-5) pos[2] = 0.0;
			fprintf(fp,"%s %7.5f %7.5f %7.5f\n",elTable[superCell.atoms[i].Znum],pos[0],pos[1],pos[2]);
		}
		if (nstart > 0)
			fprintf(fp,"frame1 %7.5f %7.5f %7.5f 1 0 0 0\n",superCell.cm[0],superCell.cm[1],superCell.cm[2]);
		fprintf(fp,"end\n");
	}  // end of if fp != NULL
}


#define IA 16807 
#define IM 2147483647 
#define AM (1.0/IM) 
#define IQ 127773 
#define IR 2836 
#define MASK 123459876

/********************************************************************
* Minimal  random number generator of Park and Miller. 
* Returns a uniform random deviate between 0.0 and 1.0.
* Set idum to any integer, except MASK, to initialize the sequence.
*******************************************************************/
float ran(long *idum) { 
	static long k; 
	static float ans; 

	*idum ^= MASK; // XORing with MASK allows use of zero and other . 
	k=(*idum)/IQ;  // simple bit patterns for idum.
	*idum=IA*(*idum-k*IQ)-IR*k; // Compute idum=(IA*idum) % IM without over-  ows by Schrage s method. 
	if (*idum < 0) *idum += IM; 
	ans=static_cast<float_tt>(AM*(*idum)); // Convert idum to a  floating result. 
	*idum ^= MASK; // Unmask before return. 
	return ans; 
}


#define IA 16807 
#define IM 2147483647 
#define AM (1.0/IM) 
#define IQ 127773 
#define IR 2836 
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB) 
#define EPS 1.2e-7 
#define RNMX (1.0-EPS)
/**********************************************************************
*  Minimal  random number generator of Park and Miller with Bays-Durham shuffle 
* and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 
* (exclusive of the endpoint values). Call with idum a negative integer to initialize; 
* thereafter, do not alter idum between successive deviates in a sequence. 
* RNMX should approximate the largest  floating value that is less than 1.
*/
float_tt ran1(long *idum) { 
	int j; 
	long k; 
	static long iy=0; 
	static long iv[NTAB]; 
	float_tt temp; 
	if (*idum <= 0 || !iy) { // Initialize. 
		if (-(*idum) < 1) *idum=1; // Be sure to prevent  idum = 0. 
		else *idum = -(*idum); 
		for (j=NTAB+7;j>=0;j--) { // Load the shu e table (after 8 warm-ups). 
			k=(*idum)/IQ; 
			*idum=IA*(*idum-k*IQ)-IR*k; 
			if (*idum < 0) *idum += IM; 
			if (j < NTAB) iv[j] = *idum; 
		} 
		iy=iv[0]; 
	} 
	k=(*idum)/IQ; // Start here when not initializing. 
	*idum=IA*(*idum-k*IQ)-IR*k; // Compute idum=(IA*idum) % IM without overflows by Schrage s method. 
	if (*idum < 0) *idum += IM; 
	j=iy/NDIV; // Will be in the range 0..NTAB-1. 
	iy=iv[j];  // Output previously stored value and re ll the shu e table. 
	iv[j] = *idum; 
	if ((temp=static_cast<float_tt>(AM*iy)) > RNMX) return static_cast<float_tt>(RNMX); // Because users don't expect endpoint values. 
	else return temp; 
}


/*****************************************************************
* Gaussian distribution with unit variance
* idum must be initailized to a negative integer 
//  TODO: we can use some other random number generation library.  Boost?
****************************************************************/
float_tt gasdev(long *idum) 
/* Returns a normally distributed deviate with zero mean and unit variance, 
* using ran1(idum) as the source of uniform deviates. */
{ 
	// float ran1(long *idum); 
	static int iset=0; 
	static float gset; 
	float_tt fac,rsq,v1,v2; 
	if (*idum < 0) {
		iset=0; // Reinitialize. 
		//    printf("reinit gasdev\n");
	}
	if (iset == 0) { 
		/* We don t have an extra deviate handy, so 
		* pick two uniform numbers in the square extending from -1 to +1 in each direction, */
		do { 
			v1=2.0f*ran1(idum)-1.0f;  
			v2=2.0f*ran1(idum)-1.0f; 
			rsq=v1*v1+v2*v2;  // see if they are in the unit circle,
		} while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again. 
		fac=static_cast<float_tt>(sqrt(-2.0*log(rsq)/rsq)); 
		/* Now make the Box-Muller transformation to get two normal deviates. 
		*  Return one and save the other for next time. 
		*/
		gset=v1*fac; 
		iset=1; // Set flag. 
		return v2*fac; 
	} 
	else { //We have an extra deviate handy,so unset the flag, and return it. 
		iset=0;  
		return gset; 
	} 
}

void writeSTEMinput(char* stemFile,char *cfgFile,MULS *muls) {
	FILE *fp;
	char folder[64];


	strcpy(folder,cfgFile);
	folder[strlen(folder)-4] = 0;  // cut off .cfg ending

	if ((fp=fopen(stemFile,"w")) == NULL) {
		printf("writeSTEMinput: Error opening file %s\n",stemFile);
		exit(0);
	}

	fprintf(fp,"%% STEM input file for tomography series, cerated automatically by STEM\n\n");

	fprintf(fp,"mode: CBED\n");
	fprintf(fp,"print level: 1\nsave level: %d\n",muls->saveLevel);
	fprintf(fp,"filename: %s\n",cfgFile);
	fprintf(fp,"NCELLX: %d\nNCELLY: %d\nNCELLZ: %d/%d\n",
		muls->nCellX,muls->nCellY,muls->nCellZ,muls->cellDiv);
	fprintf(fp,"v0: %g\n",muls->v0);
	fprintf(fp,"tds: %s\n",muls->tds ? "yes" : "no");
	fprintf(fp,"temperature: %g\n",muls->tds_temp);
	fprintf(fp,"slice-thickness: %g\n",muls->sliceThickness);
	fprintf(fp,"periodicXY: no\n");	
	fprintf(fp,"periodicZ: no\n");


	fprintf(fp,"nx: %d\n",muls->nx); 
	fprintf(fp,"ny: %d\n",muls->ny);      
	fprintf(fp,"Cs: %g\n",muls->Cs*1e-7);	
	fprintf(fp,"C5: %g\n",muls->C5);		
	fprintf(fp,"Cc: %g\n",muls->Cc*1e-7);		
	fprintf(fp,"dV/V: %g\n",muls->dE_E);	
	fprintf(fp,"alpha: %g\n",muls->alpha);		
	fprintf(fp,"AIS aperture: %g %% A\n",muls->aAIS); 
	fprintf(fp,"smooth: %s\n",muls->ismoth ? "yes" : "no");		
	fprintf(fp,"defocus: %g \n",muls->df0);
	fprintf(fp,"Source Size (diameter): %g \n",2*muls->sourceRadius);
	fprintf(fp,"gaussian: %s\n", muls->gaussFlag ? "yes" : "no");
	fprintf(fp,"potential3D: %s \n", muls->potential3D ? "yes" : "no");
	fprintf(fp,"atom radius: %g \n", muls->atomRadius);
	fprintf(fp,"plot V(r)*r: yes \n");	
	fprintf(fp,"bandlimit f_trans: yes\n");	
	fprintf(fp,"save potential: no	\n");
	fprintf(fp,"one time integration: yes \n");

	fprintf(fp,"tomo tilt: %lf\n",muls->tomoTilt);


	fprintf(fp,"Display Gamma: 0 \n");    
	fprintf(fp,"Folder: %s\n",folder);
	fprintf(fp,"Runs for averaging: %d\n",muls->avgRuns); 
	fprintf(fp,"Structure Factors: DT  \n");
	fprintf(fp,"show Probe: %s \n",muls->showProbe ? "yes" : "no");	
	fprintf(fp,"propagation progress interval: 10 \n");
	fprintf(fp,"potential progress interval: 1000 \n");
	fprintf(fp,"beams: n \n");       
	fprintf(fp,"sequence: 1 1\n");

	fclose(fp);
}
