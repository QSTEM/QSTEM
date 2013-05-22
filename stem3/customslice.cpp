#include <stdio.h>	/*  ANSI-C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "stemtypes_fftw3.h"
#include "memory_fftw3.h"
#include "imagelib_fftw3.h"
#include "stemutil.h"
#include "customslice.h"
#include "fileio_fftw3.h"
#include "comparators.h"

#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#ifdef WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif

#define FTBOX_AX 40
#define FTBOX_CZ 40
#define FTBOX_NX 1024
#define FTBOX_NZ 1024

#define SQR(x) ((x)*(x))

void make3DSlicesFT(MULS *muls) {

  /************************************************************************
   * temporary variables needed always:
   */
  float_tt scale,ffr,ffi,arg,r;  // ffr,i = form factor
  float_tt s,x,y,z,rxy2,zStart;
  int ix,iy,iz,j,i; // count = 0; // j,Ninteg;
  int atKind;                   // specifies kind of atom in list
  float_tt timer0,timer,t;
  std::vector<atom> atoms;
  int perX = 0, perY = 0, perZ = 0; // set 1 for periodic, 0 for non-periodic
  int pX,pY,pZ;                  // used for looping through neighboring cells


  /*************************************************************************
   * static variables which are used every time:
   */
  // static float ***potential = NULL;         // array holding real space potential
  int Nxp,Nyp,Nzp;                   // size of potential array
  int Nx=FTBOX_NX, Nz=FTBOX_NX, Nxm,Nzm; // size and center of single atom potential box
  int Nxl,Nzl;                       // size of lookup array
  float_tt axp, byp,czp,dXp,dYp,dZp;    // model dimensions and resolution
  //float_tt ***potLUT = NULL;           // potential lookup table for all atoms used
  QSVecOffMat potLUT;
  float_tt xOversample,zOversample;       // oversampling rate for x- and z-direction
  QSfVec rcutoff;            // radius after which set potential to 0 (one for each atom)
  int divCount = 0;
  char buf[128];
  float_tt dX,dZ;                      // real space resol. of FT box
  imageStruct *header = NULL;
  char fileName[64];


  /******************************************
   * only needed during initialization: 
   */
  fftw_plan plan;                   // fftw array
  int fftMeasureFlag = FFTW_ESTIMATE; // fftw plan needed for FFT
  QScMat pot;          // single atom potential box
  float_tt ax,cz;                        // real space size of FT box
  float_tt dsX,dsZ;                     // rec. space size of FT box
  float_tt sx,sz,sx2,sz2,sz2r;
  float_tt sx2max,sz2max;

  timer0 = getTime();

  /**********************************************************
   * Initialization of parameter used in this function
   *********************************************************/
  //if (potLUT == NULL) {
    Nxp = muls->potNx;
    Nyp = muls->potNy;
    Nzp = muls->slices;
    // potential = float3D(Nzp,Nyp,Nxp,"potential");
    // memset(potential(0,0),0,Nzp*Nyp*Nxp*sizeof(float));

    // read the atomic structure and the size of the box:
    atoms = muls->atoms;
    // axp = muls->ax; byp = muls->by; 
    czp = muls->cellDims[2]/muls->cellDiv;
    axp = muls->potSizeX; byp = muls->potSizeY;
    // Now we need to adjust ax, Nx, cz, and Nz to work with the model parameters
    dXp = axp/Nxp;
    dYp = byp/Nyp;
    dZp = czp/Nzp;
    
    // find the spacing with which we will create the projected 
    // potential lookup table.
    // we need to space it quite finely in z-direction in order to remain 
    // somewhat accurate.
    ax = FTBOX_AX;  cz = FTBOX_CZ;
    dZ = dZp;
    zOversample = 1;
    while (dZ > 0.4) zOversample += 2, dZ = dZp/(float_tt)(zOversample);
    Nz = (int)(cz/dZ+0.5);    if (Nz/2 < 0.5*(float_tt)Nz) Nz++;
    cz =  (float_tt)Nz*dZ;
    
    dX = (dXp < dYp ? dXp : dYp);
    xOversample = 1;
    while (dX > 0.1) dX = dXp/(float_tt)(++xOversample);    
    Nx = (int)(ax/dX); if (Nx/2 < 0.5*(float_tt)Nx) Nx++;
    ax = (float_tt)Nx*dX;
    // sAtom.x = 0.5*ax;  sAtom.y = 0.5*ax;  sAtom.z = 0.5*cz;
        
    printf("box: (%g, %g, %g), sampled: (%g, %g, %g)\n",axp,byp,czp,dXp,dYp,dZp);
    printf("Nr=%d, Nz=%d (%d,%d times oversampling)  (%g, %g) (%g, %g)\n",
	   Nx,Nz,zOversample,xOversample,ax,cz,dX,dZ);
    
    // create an array, so that index=iz*Nx*Ny+iy*Nx+ix = [iz][iy][ix]
    pot = QScMat(Nz, Nx);
	// TODO: is this right?
    potLUT = QSVecOffMat(muls->atomKinds);
		//(float_tt ***)malloc(muls->atomKinds*sizeof(float_tt **));
    rcutoff = QSfVec(muls->atomKinds);
	rcutoff.setZero();
    //memset(rcutoff,0,muls->atomKinds*sizeof(float_tt));
    Nxm = Nx/2;   // Nym = Ny/2;
    dsX = 1.0f/ax; // dsY = 1.0/by; 
    dsZ = 1.0f/cz;

    
    sx2max = 10.8f*(Nxm*dsX);  sx2max = 1.0f/SQR(sx2max);
    // sy2max = 10.0*(Nym*dsY);  sy2max = 1.0/SQR(sy2max);
    if (Nz > 1) Nzm = Nz/2, sz2max = 1.0f/SQR(10.8f*(Nzm*dsZ));  
    else Nzm = 1, sz2max = 1.0f;
    printf("ds=(%g,%g)/A, smax=(%g,%g)\n",dsX,dsZ,Nxm*dsX,Nzm*dsZ);

    // Now we need to make sure that the scattering factor tapers off to zero at the end:
    // sTable[sfCSize-1] = 0.5* (Nxm*dsX > Nzm*dsZ ? Nxm*dsX : Nzm*dsZ);
    // sfC[sfCSize-1] = 0.0;

    timer = getTime();
    /*************************************************************
     * We will now calculate the real space potential for every
     * kind of atom used in this model
     ************************************************************/
    for (atKind = 0; atKind<muls->atomKinds;atKind++) 
	{
      //memset(pot[0],0,sizeof(fftw_complex)*Nz*Nx);
		pot.setZero();

      for (iz=0;iz<Nz;iz++) 
	  {
		sz = (iz < Nzm ? iz : iz-Nz)*dsZ, sz2 = SQR(sz), sz2r = sz2max*sz2;
		if ((Nz < 2) || (sz2r <=1.0)) 
		{
			for (ix=0;ix<Nx;ix++) 
			{
			sx = (ix < Nxm ? ix : ix-Nx)*dsX, sx2 = SQR(sx);	
				if (sz2r+sx2max*sx2 <=1.0) 
				{  // enforce ellipse equation:
					s = sqrt(sz2+sx2);
					arg = static_cast<float_tt>(TWOPI * (sx*(0.5f*ax)+sz*(0.5f*cz))); // place single atom in center of box
					ffr = cos(arg); ffi = sin(arg);	  
					// all the s are actually q, therefore S = 0.5*q = 0.5*s:
					pot(ix,iz) = std::complex<float_tt>(sfLUT(0.5f*s,atKind,muls)*ffr, 
														sfLUT(0.5f*s,atKind,muls)*ffi);
				}      
			}
		}
      }

      // new fftw3 code:
	  // TODO: Define different plans for different data types.
      plan = fftw_plan_dft_2d(Nz,Nx,(fftw_complex *)pot.data(),(fftw_complex *)pot.data(),FFTW_BACKWARD,fftMeasureFlag);
      fftw_execute(plan);
      fftw_destroy_plan(plan);
      // old fftw2 code:
      // plan = fftw2d_create_plan(Nz,Nx,FFTW_BACKWARD,fftMeasureFlag | FFTW_IN_PLACE);  
      // fftwnd_one(plan, pot[0], NULL); 
      // fftwnd_destroy_plan(plan); 
    
      /* see L.M. Peng, Micron 30, p. 625 (1999) for details on the scale factor
       * so that pot is the true electrostatic potential. 
       */
      // scale = 47.87658/(sqrt((float_tt)(Nx*Ny))*(float_tt)Nz);
      scale = 2*47.87658f/(sqrt((float_tt)(Nx*Nz)*ax*cz));
      rcutoff[atKind] = 50;  // 50 A cutoff radius, initially
      for (iz=0;iz<Nz;iz++) {
	for (ix=0;ix<Nx;ix++) {
	  r = sqrt(SQR(dZ*(iz-Nzm))+SQR(dX*(ix-Nxm)));
	  if (r > rcutoff[atKind]) pot(ix,iz) = QScf(0.0f, 0.0f);
	  else {
	    if (pot(ix,iz).real() < 0) {
	      rcutoff[atKind] = r;
	      pot(ix,iz) = QScf(0.0f, 0.0f);
	    }
	    else
	      pot(ix,iz) *= scale;
	  }
	}
      }  
      for (iz=0;iz<Nz;iz++) {
	for (ix=0;ix<Nx;ix++) {
	  r = sqrt(SQR(dZ*(iz-Nzm))+SQR(dX*(ix-Nxm)));
	  if (r > rcutoff[atKind]) pot(ix,iz)=QScf(0.0f, 0.0f);
	}
      }

      /* Now we need to integrate over z, if we only have a single slice
       */
      if (Nzp == 1) {
	for (ix=0;ix<Nx;ix++) for (iz=1;iz<Nz;iz++)
	  pot(ix,0) += pot(ix,iz);	
	Nz = 1;
      }
      /* Now we will create a double array to hold the real valued potential as
       * a lookup table for later interpolation of potentials for all the atoms
       * and pre-integrate its values.
       * We are assuming that zOversample is always an odd number!
       */  
      potLUT[atKind] = reduceAndExpand(pot,Nz,Nx,zOversample,Nzl,Nxl);      
    } // end of for atKind ...
  //} // if potential == NULL
  /*******************************************************************/
  
  /*******************************************************
   * initializing  cz, and trans
   *************************************************************/

  // TODO: Do we need to loop through trans and zero it out?
  //memset(muls->trans(0,0),0,Nzp*Nxp*Nyp*sizeof(fftw_complex));
  //if (muls->cz == NULL) muls->cz = float1D(Nzp,"cz");
  muls->cz = QSfVec(Nzp);
  for (i=0;i<Nzp;i++) muls->cz[i] = muls->sliceThickness;  					
  
  /********************************************************************
   * Now that we have the lookup table and are able to find interpolations in it, 
   * we can start add the potentials of all the atoms in the structure
   */
  atoms = muls->atoms;
  if (divCount == 0) {
	std::sort(atoms.begin(), atoms.end(), atomCompareZnum());
    //qsort(atoms,muls->natom,sizeof(atom),atomCompare());
    if ((*muls).cfgFile != "") {
      sprintf(buf,"%s/%s",muls->folder,muls->cfgFile);
      writeCFG(atoms,muls->natom,buf,muls);	
    }
  }
  zStart = (*muls).czOffset+(float_tt)divCount*czp/((float_tt)(*muls).cellDiv);
  divCount = (divCount + 1) % (*muls).cellDiv;
  
  if (Nzp > 1) {
    timer = getTime();  
    for (j=0;j<muls->natom;j++) {
      for (atKind = 0;muls->Znums[atKind]!=atoms[j].Znum;atKind++); // find atKind 
      for (iz=0;iz<Nzp;iz++) for (pZ = 0;pZ<=perZ;pZ++) {
	z = fabs((iz+0.5f)*dZp-(atoms[j].pos[2]-zStart)+pZ*czp);
	// if (z <= rcutoff[atKind]+dZp) {
	if (z <= rcutoff[atKind]+dZp) {
	  for (ix=0;ix<Nxp;ix++)  for (pX = 0;pX<=perX;pX++) {
	    x =((ix+0.5f)*dXp-atoms[j].pos[0]+pX*axp)+ muls->potOffsetX;
	    for (iy=0;iy<Nyp;iy++) for (pY = 0;pY<=perY;pY++)  {
	      y = ((iy+0.5f)*dYp-atoms[j].pos[1]+pY*byp) + muls->potOffsetY;
	      rxy2 = SQR(x)+SQR(y);
	      r = sqrt(rxy2+SQR(z));
	      if (r <=rcutoff[atKind]+dZp+dXp)
			  // adding only to real value!
		muls->trans[iz](iy,ix) += bicubic(potLUT[atKind],Nzl,Nxl,z/dZ+1.0f,sqrt(rxy2)/dX+1.0f);
	      //potential(iy,iz)[ix] += bicubic(potLUT[atKind],Nzl,Nxl,z/dZ+1.0,sqrt(rxy2)/dX+1.0);
	    }      
	  }
	}
      }
      /*
      if (j % 40 == 0) 
	if (getTime()-timer >= 10) {
	  timer += 10.0;
	  printf("%2d%% (%d sec)\n",(int)(100.0*(float_tt)(j)/((float_tt)muls->natom)),
		 (int)(getTime()-timer0));
	}
      */      
      
      if (j % 100 == 0) {
	if ((t=getTime() - timer) >= 10) {
	  t = t*(float_tt)muls->natom/(float_tt)j;
	  timer  += t;
	  printf("Potential integration: %d%% done, %d min, %d sec left\n",(int)(100.0*(float_tt)(j)/((float_tt)muls->natom)),
	       (int)(t/60.0),(int)((int)(t) % 60));
	}
      }
      
    }  // end of for j=0 ... natom
  }
    
  // save the potential file:
  if (muls->savePotential) {
    for (iz=0;iz<Nzp;iz++) {
      sprintf(fileName,"%s/%s%d.img",muls->folder,muls->fileBase,iz);
      // printf("Saving potential layer %d to file %s\n",iz,filename); 
	  if (header == NULL) header = makeNewHeaderCompact(1,Nxp,Nyp,dZp,dXp,dYp,0,std::vector<float_tt>(),NULL);
	  sprintf(buf,"Projected Potential (%d slices)",muls->slices);
      header->comment = buf;
      header->commentSize = (int)header->comment.length();
      writeComplexImage(muls->trans[iz],header,fileName);      
    } 
  } /* end of if savePotential ... */
  
  printf("Calculation took %.1f sec, rc[0]: %gA\n",(getTime()-timer0),rcutoff[0]);
}  // end of function


/******************************************************************
 * This function takes part of the complex (single precision) 
 * potential and puts it into a real double precision array.
 *
 * possible expansion algorithm:
 *
 *  % Expand z so interpolation is valid at the boundaries.
 * zz = zeros(size(arg1)+2);
 * zz(1,2:ncols+1) = 3*arg1(1,:)-3*arg1(2,:)+arg1(3,:);
 * zz(2:nrows+1,2:ncols+1) = arg1;
 * zz(nrows+2,2:ncols+1) = 3*arg1(nrows,:)-3*arg1(nrows-1,:)+arg1(nrows-2,:);
 * zz(:,1) = 3*zz(:,2)-3*zz(:,3)+zz(:,4);
 * zz(:,ncols+2) = 3*zz(:,ncols+1)-3*zz(:,ncols)+zz(:,ncols-1);
 * nrows = nrows+2; ncols = ncols+2;
 ******************************************************************/
QSfMat reduceAndExpand(QScMat fc,int Nz,int Nx,int zOversample,int &fNz,int &fNx) 
{
  int ix,iz,j,nx,nz,Ninteg,Nxm,Nzm;

  Nxm = Nx/2;
  Nzm = Nz/2;
  nx = Nxm;
  nz = (Nz == 1 ? 1 : Nzm+2);
  QSfMat ff(nz, nx);
  ff.setZero();
  //ff = float_tt2D(nz,nx,"ff");
  //memset(ff[0],0,nz*nx*sizeof(float_tt));
  fNx = nx;
  fNz = nz;

  Ninteg = (zOversample-1)/2;
  if (Ninteg*2+1 != zOversample) {
    printf("Error: zOversample = %d (%d) is not an odd number!\n",zOversample,Ninteg);
    exit(0);
  }

  if (Nz == 1) {
    for (ix=0;ix<Nxm;ix++) ff(ix+1,0) = fc(ix+Nxm,0).real();
    
    ff(0,0) = 3.0f*ff(1,0)-3.0f*ff(2,0)+ff(3,0);
    // ff(0,0) = ff(1,0)+(ff(1,0)-ff(2,0));
    return ff;
  }
  
  for (iz=0;iz<Nzm;iz++) for (ix=0;ix<Nxm;ix++) { 
    for (j=-Ninteg ; j<=Ninteg;j++)
	if ((iz+j+Nzm>=0) && (iz+j<Nzm))
	  ff(ix+1,iz+1) += fc(ix+Nxm,iz+Nzm+j).real();
  }
  // zz(:,1) = 3*zz(:,2)-3*zz(:,3)+zz(:,4);
  // zz(1,2:ncols+1) = 3*arg1(1,:)-3*arg1(2,:)+arg1(3,:);
  for (iz=1;iz<=Nzm;iz++) 
    // ff(0,iz) = 3.0*ff(1,iz)-3.0*ff(2,iz)+ff(3,iz);
    ff(0,iz) = ff(1,iz)+(ff(1,iz)-ff(2,iz));
  for (ix=0;ix<=Nxm;ix++) 
    // ff(ix,0) = 3.0*ff(ix,1)-3.0*ff(ix,2)+ff(ix,3);
    ff(ix,0) = ff(ix,1)+(ff(ix,1)-ff(ix,2));
  
  return ff;
}

