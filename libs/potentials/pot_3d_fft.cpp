#include "pot_3d_fft.hpp"

C3DFFTPotential::C3DFFTPotential(ConfigReaderPtr &configReader) : C3DPotential(configReader)
{
}

void C3DFFTPotential::makeSlices(int nlayer, char *fileName, atom *center)
{
  for (i = 0;i<nlayer;i++) if ((*muls).cz[0] != (*muls).cz[i]) break;
  if (i<nlayer) printf("Warning: slice thickness not constant, will give wrong results (iz=%d)!\n",i);
}


void C3DFFTPotential::AddAtomNonPeriodic()
{
  iAtomZ = (int)floor(atomZ/muls->sliceThickness+0.5);
      
  iax0 = iAtomX-iRadX < 0 ? 0 : iAtomX-iRadX;
  iax1 = iAtomX+iRadX >= muls->potNx ? muls->potNx-1 : iAtomX+iRadX;
  iay0 = iAtomY-iRadY < 0 ? 0 : iAtomY-iRadY;
  iay1 = iAtomY+iRadY >= muls->potNy ? muls->potNy-1 : iAtomY+iRadY;
  // if within the potential map range:
  if ((iax0 < muls->potNx) && (iax1 >= 0) && (iay0 < muls->potNy) && (iay1 >= 0)) {
    // define range of sampling from atomZ-/+atomRadius
    iaz0 = iAtomZ-iRadZ < 0 ? -iAtomZ : -iRadZ;
    iaz1 = iAtomZ+iRadZ >= muls->slices ? muls->slices-iAtomZ-1 : iRadZ;
    if ((iAtomZ+iaz0 <        muls->slices) && (iAtomZ+iaz1 >= 0)) {
      // retrieve the pointer for this atom
      atPotPtr = getAtomPotential3D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw,&nzSub,&Nr,&Nz_lut);
#if USE_Q_POT_OFFSETS
      // retrieve the pointer to the array of charge-dependent potential offset
      // This function will return NULL; if the charge of this atom is zero:
      atPotOffsPtr = getAtomPotentialOffset3D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw,&nzSub,&Nr,&Nz_lut,atoms[iatom].q);
#endif // USE_Q_POT_OFFSETS
      iOffsLimHi = Nr*(Nz_lut-1);
      iOffsLimLo = -Nr*(Nz_lut-1);
      iOffsStep = nzSub*Nr;

      // Slices around the slice that this atom is located in must be affected by this atom:
      // iaz must be relative to the first slice of the atom potential box.
      for (iax=iax0; iax <= iax1; iax++) {
        potPtr = &(muls->trans[iAtomZ+iaz0][iax][iay0][0]);

        //////////////////////////////////////////////////////////////////////
        // Computation of Radius must be made faster by using pre-calculated ddx
        // and LUT for sqrt:
        // use of sqrt slows down from 120sec to 180 sec.
        x2 = iax*dx - atomX; x2 *= x2;
        for (iay=iay0; iay <= iay1; iay++) {
          // printf("iax=%d, iay=%d\n",iax,iay);
          y2 = iay*dy - atomY; y2 *= y2;
          r = sqrt(x2+y2);
          // r = (x2+y2);
          ddr = r/dr;
          ir  = (int)floor(ddr);
          // add in different slices, once r has been defined
          if (ir < Nr-1) {
            ddr = ddr-(double)ir;
            ptr = potPtr;
            
            dOffsZ = (iAtomZ+iaz0-atomZ/muls->sliceThickness)*nzSub;
#if Z_INTERPOLATION
            iOffsZ = (int)dOffsZ;
            ddz = fabs(dOffsZ - (double)iOffsZ);
#else
            iOffsZ = (int)(dOffsZ+0.5);
#endif
            iOffsZ *= Nr;
            
            for (iaz=iaz0; iaz <= iaz1; iaz++) {
              potVal = 0;
              if (iOffsZ < 0) {
                if (iOffsZ > iOffsLimLo) {
                  // do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
                  potVal = (1-ddz)*((1-ddr)*atPotPtr[ir-iOffsZ+Nr][0]+ddr*atPotPtr[ir+1-iOffsZ+Nr][0])+
                    ddz *((1-ddr)*atPotPtr[ir-iOffsZ ][0]+ddr*atPotPtr[ir+1-iOffsZ ][0]);
#if USE_Q_POT_OFFSETS
                  // add the charge-dependent potential offset
                  if (atPotOffsPtr != NULL) {
                    potVal += atoms[iatom].q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir-iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1-iOffsZ+Nr][0])+
                                              ddz *((1-ddr)*atPotOffsPtr[ir-iOffsZ ][0]+ddr*atPotOffsPtr[ir+1-iOffsZ ][0]));                
                  }        
#endif // USE_Q_POT_OFFSETS
#else // Z_INTERPOLATION
                  potVal = (1-ddr)*atPotPtr[ir-iOffsZ+Nr][0]+ddr*atPotPtr[ir+1-iOffsZ+Nr][0];
#if USE_Q_POT_OFFSETS
                  // add the charge-dependent potential offset
                  if (atPotOffsPtr != NULL) {
                        potVal += atoms[iatom].q*((1-ddr)*atPotOffsPtr[ir-iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1-iOffsZ+Nr][0]);
                      }        
#endif // USE_Q_POT_OFFSETS
#endif // Z_INTERPOLATION
                    }
                  } // if iOffZ < 0
                  else {
                    // select the pointer to the right layer in the lookup box
                    if (iOffsZ < iOffsLimHi) {
                      // do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
                      potVal = (1-ddz)*((1-ddr)*atPotPtr[ir+iOffsZ][0]+ddr*atPotPtr[ir+1+iOffsZ][0])+
                        ddz *((1-ddr)*atPotPtr[ir+iOffsZ+Nr][0]+ddr*atPotPtr[ir+1+iOffsZ+Nr][0]);
#if USE_Q_POT_OFFSETS
                      // add the charge-dependent potential offset
                      if (atPotOffsPtr != NULL) {
                        potVal += atoms[iatom].q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0])+
                                                  ddz *((1-ddr)*atPotOffsPtr[ir+iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1+iOffsZ+Nr][0]));
                      }
#endif // USE_Q_POT_OFFSETS
#else // Z_INTERPOLATION
                      potVal = (1-ddr)*atPotPtr[ir+iOffsZ][0]+ddr*atPotPtr[ir+1+iOffsZ][0];
#if USE_Q_POT_OFFSETS
                      // add the charge-dependent potential offset
                      if (atPotOffsPtr != NULL) {
                        potVal += atoms[iatom].q*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0]);                                                                                                
                      }
#endif // USE_Q_POT_OFFSETS
#endif // Z_INTERPOLATION
                                
                    }
                  } // if iOffsZ >=0
                  *ptr += potVal; // ptr = potPtr = muls->trans[...]
                            
                  ptr += sliceStep;        // advance to the next slice
                  // add the remaining potential to the next slice:
                  // if (iaz < iaz1)        *ptr += (1-ddz)*potVal;
                  iOffsZ += iOffsStep;
                } // end of iaz-loop
              } // if ir < Nr
              // advance pointer to next complex potential point:
              potPtr+=2;
                        
            } // iay=iay0 .. iay1        
          } // iax=iax0 .. iax1
        } // iaz0+iAtomZ < muls->slices
      } // if within bounds        
}


void C3DFFTPotential::AddAtomPeriodic()
{
  iAtomZ = (int)floor(atomZ/muls->sliceThickness+0.5);
  iax0 = iAtomX-iRadX;
  iax1 = iAtomX+iRadX;
  iay0 = iAtomY-iRadY;
  iay1 = iAtomY+iRadY;

  // define range of sampling from atomZ-/+atomRadius
  iaz0 = iAtomZ-iRadZ < 0 ? -iAtomZ : -iRadZ;
  iaz1 = iAtomZ+iRadZ >= muls->slices ? muls->slices-iAtomZ-1 : iRadZ;
  // if (iatom < 2) printf("iatomZ: %d, %d cz=%g, %g: %d, %d\n",iAtomZ,iaz0,muls->sliceThickness,atomZ,(int)(-2.5-atomZ),(int)(atomZ+2.5));
  // printf("%d: iatomZ: %d, %d cz=%g, %g\n",iatom,iAtomZ,iaz0,muls->sliceThickness,atomZ);
  if ((iAtomZ+iaz0 < muls->slices) && (iAtomZ+iaz1 >= 0)) {
    // retrieve the pointer for this atom
    atPotPtr = getAtomPotential3D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw,&nzSub,&Nr,&Nz_lut);
#if USE_Q_POT_OFFSETS
    // retrieve the pointer to the array of charge-dependent potential offset
    // This function will return NULL; if the charge of this atom is zero:
    atPotOffsPtr = getAtomPotentialOffset3D(atoms[iatom].Znum,muls,muls->tds ? 0 : atoms[iatom].dw,&nzSub,&Nr,&Nz_lut,atoms[iatom].q);
#endif // USE_Q_POT_OFFSETS
    iOffsLimHi =        Nr*(Nz_lut-1);
    iOffsLimLo = -Nr*(Nz_lut-1);
    iOffsStep = nzSub*Nr;

    // Slices around the slice that this atom is located in must be affected by this atom:
    // iaz must be relative to the first slice of the atom potential box.
    for (iax=iax0; iax < iax1; iax++) {
      potPtr = &(muls->trans[iAtomZ+iaz0][(iax+2*muls->potNx) % muls->potNx][(iay0+2*muls->potNy) % muls->potNy][0]);
      // potPtr = &(muls->trans[iAtomZ-iaz0+iaz][iax][iay0][0]);
      x2 = iax*dx - atomX;        x2 *= x2;
      for (iay=iay0; iay < iay1; ) {
        // printf("iax=%d, iay=%d\n",iax,iay);
        y2 = iay*dy - atomY;        y2 *= y2;
        r = sqrt(x2+y2);
        ddr = r/dr;
        ir = (int)floor(ddr);
        // add in different slices, once r has been defined
        if (ir < Nr-1) {
          ddr = ddr-(double)ir;
          ptr = potPtr;
          // Include interpolation in z-direction as well (may do it in a very smart way here !):
                        
          dOffsZ = (iAtomZ+iaz0-atomZ/muls->sliceThickness)*nzSub;
#if Z_INTERPOLATION
          iOffsZ = (int)dOffsZ;
          ddz         = fabs(dOffsZ - (double)iOffsZ);
#else // Z_INTERPOLATION
          iOffsZ = (int)(dOffsZ+0.5);
#endif // Z_INTERPOLATION
          iOffsZ *= Nr;
                        
          for (iaz=iaz0; iaz <= iaz1; iaz++) {
            potVal = 0;
            // iOffsZ = (int)(fabs(iAtomZ+iaz-atomZ/muls->sliceThickness)*nzSub+0.5);
            if (iOffsZ < 0) {
              if (iOffsZ > iOffsLimLo) {
                // do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
                potVal = (1-ddz)*((1-ddr)*atPotPtr[ir-iOffsZ+Nr][0]+ddr*atPotPtr[ir+1-iOffsZ+Nr][0])+
                  ddz *((1-ddr)*atPotPtr[ir-iOffsZ][0]+ddr*atPotPtr[ir+1-iOffsZ][0]);
#if USE_Q_POT_OFFSETS
                // add the charge-dependent potential offset
                if (atPotOffsPtr != NULL) {
                  potVal += atoms[iatom].q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir-iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1-iOffsZ+Nr][0])+
                                            ddz *((1-ddr)*atPotOffsPtr[ir-iOffsZ ][0]+ddr*atPotOffsPtr[ir+1-iOffsZ ][0]));                
                }        
#endif // USE_Q_POT_OFFSETS
#else // Z_INTERPOLATION
                potVal = (1-ddr)*atPotPtr[ir-iOffsZ+Nr][0]+ddr*atPotPtr[ir+1-iOffsZ+Nr][0];
#if USE_Q_POT_OFFSETS
                // add the charge-dependent potential offset
                if (atPotOffsPtr != NULL) {
                  potVal += atoms[iatom].q*((1-ddr)*atPotOffsPtr[ir-iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1-iOffsZ+Nr][0]);
                }        
#endif // USE_Q_POT_OFFSETS
#endif // Z_INTERPOLATION
              }
            }
            else {
              // select the pointer to the right layer in the lookup box
              // printf("%4d: iOffsZ: %d, iaz: %d (slice: %d, pos: %g [%d .. %d])\n",iatom,iOffsZ,iaz,iAtomZ+iaz,atomZ,iaz0,iaz1);
              if (iOffsZ < iOffsLimHi) {
                // do the real part by linear interpolation in r-dimension:
#if Z_INTERPOLATION
                potVal = (1-ddz)*((1-ddr)*atPotPtr[ir+iOffsZ][0]+ddr*atPotPtr[ir+1+iOffsZ][0])+
                  ddz *((1-ddr)*atPotPtr[ir+iOffsZ+Nr][0]+ddr*atPotPtr[ir+1+iOffsZ+Nr][0]);
#if USE_Q_POT_OFFSETS
                // add the charge-dependent potential offset
                if (atPotOffsPtr != NULL) {
                  potVal += atoms[iatom].q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0])+
                                            ddz *((1-ddr)*atPotOffsPtr[ir+iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1+iOffsZ+Nr][0]));
                }
#endif // USE_Q_POT_OFFSETS
#else // Z_INTERPOLATION
                potVal = (1-ddr)*atPotPtr[ir+iOffsZ][0]+ddr*atPotPtr[ir+1+iOffsZ][0];
#if USE_Q_POT_OFFSETS
                // add the charge-dependent potential offset
                if (atPotOffsPtr != NULL) {
                  potVal += atoms[iatom].q*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0]);                                                                                                
                }
#endif // USE_Q_POT_OFFSETS
#endif // Z_INTERPOLATION
              }
            }
            *ptr += potVal;
            
            ptr += sliceStep; // advance to the next slice
            // add the remaining potential to the next slice:
            // if (iaz < iaz1) *ptr += (1-ddz)*potVal;
            iOffsZ += iOffsStep;
          }
        } // if ir < Nr-1
        // advance pointer to next complex potential point:
        potPtr+=2;
        // make imaginary part zero for now
        // *potPtr = 0;
        // potPtr++;
                      
        // wrap around when end of y-line is reached:
        if (++iay % muls->potNy == 0) potPtr -= 2*muls->potNy;                
      }
    }
  } // iaz0+iAtomZ < muls->slices
}


/********************************************************************************
* Create Lookup table for 3D potential due to neutral atoms
********************************************************************************/
#define PHI_SCALE 47.87658
complex_tt *C3DFFTPotential::GetAtomPotential3D(int Znum, MULS *muls,double B,int *nzSub,int *Nr,int *Nz_lut) {
  int ix,iy,iz,iiz,ind3d,iKind,izOffset;
  double zScale,kzmax,zPos,xPos;
#if FLOAT_PRECISION == 1
  fftwf_plan plan;
#else
  fftw_plan plan;
#endif
  static double f,phase,s2,s3,kmax2,smax2,kx,kz,dkx,dky,dkz; // ,dx2,dy2,dz2;
  static int nx,ny,nz,nzPerSlice;
  static complex_tt **atPot = NULL;
  static complex_tt *temp = NULL;
#if SHOW_SINGLE_POTENTIAL == 1
  ImageIOPtr imageio = ImageIOPtr();
  static complex_tt *ptr = NULL;
  static char fileName[256];
#endif
  static double *splinb=NULL;
  static double *splinc=NULL;
  static double *splind=NULL;


  // scattering factors in:
  // float scatPar[4][30]
  if (atPot == NULL) {
    splinb = double1D(N_SF, "splinb" );
    splinc = double1D(N_SF, "splinc" );
    splind = double1D(N_SF, "splind" );

    
    nx = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionX);
    ny = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionY);
    // The FFT-resolution in the z-direction must be high enough to avoid
    // artifacts due to premature cutoff of the rec. space scattering factor
    // we will therefore make it roughly the same as the x-resolution
    // However, we will make sure that a single slice contains an integer number
    // of sampling points.
    nzPerSlice = (int)floor(OVERSAMP_X*muls->sliceThickness/muls->resolutionX);
    // make nzPerSlice odd:
    if (2.0*(nzPerSlice >> 1) == nzPerSlice) nzPerSlice += 1;
    // Total number of z-positions should be twice that of atomRadius/sliceThickness
    nz = (2*(int)ceil(muls->atomRadius/muls->sliceThickness))*nzPerSlice;
    if (muls->printLevel > 1) printf("Will use %d sampling points per slice, total nz=%d (%d)\n",nzPerSlice,nz,nzPerSlice >> 1);

    dkx = 0.5*OVERSAMP_X/(nx*muls->resolutionX); // nx*muls->resolutionX is roughly 2*muls->atomRadius
    dky = dkx;
    dkz = nzPerSlice/(double)(nz*muls->sliceThickness);
    kmax2 = 0.5*nx*dkx/(double)OVERSAMP_X; // largest k that we'll admit
    smax2 = kmax2;

    printf("dkx = %g, nx = %d, kmax2 = %g\n",dkx,nx,kmax2);
    if (muls->printLevel > 1) printf("Cutoff scattering angle: kmax=%g (1/A), dk=(%g,%g %g)\n",kmax2,dkx,dky,dkz);
    scatPar[0][N_SF-1] = 1.2*smax2;
    scatPar[0][N_SF-2] = 1.1*smax2;
    scatPar[0][N_SF-3] = smax2;
    // adjust the resolution of the lookup table if necessary
    if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3]) {

      if (1) {
        // set additional scattering parameters to zero:
        for (ix = 0;ix < N_SF-10;ix++) {
          if (scatPar[0][N_SF-4-ix] < scatPar[0][N_SF-3]-0.001*(ix+1)) break;
          scatPar[0][N_SF-4-ix] = scatPar[0][N_SF-3]-0.001*(ix+1);
          for (iy=1; iy<N_ELEM;iy++) scatPar[iy][N_SF-4-ix] = 0;
        }
      }
      else {
        for (ix = 0;ix < 20;ix++) {
          if (scatPar[0][N_SF-4-ix] < scatPar[0][N_SF-3]) break;
          scatPar[0][N_SF-4-ix] = scatPar[0][N_SF-3];
          for (iy=1; iy<N_ELEM;iy++) scatPar[iy][N_SF-4-ix] = scatPar[iy][N_SF-3];        
        }
      }                
      if (muls->printLevel > 1) printf("getAtomPotential3D: set resolution of scattering factor to %g/A!\n",
                                       scatPar[0][N_SF-4-ix]);
    }        // end of if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3])
    smax2 *= smax2;
    kmax2 *= kmax2;
    // allocate a list of pointers for the element-specific potential lookup table
    atPot = (fftwf_complex **)malloc((NZMAX+1)*sizeof(fftwf_complex *));
    for (ix=0;ix<=NZMAX;ix++) atPot[ix] = NULL;
    temp = (fftwf_complex*) fftwf_malloc(nx*nz*sizeof(fftwf_complex));
  }
  // initialize this atom, if it has not been done yet:
  if (atPot[Znum] == NULL) {
    iKind = Znum;
    
    // setup cubic spline interpolation:
    splinh(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF);
    
    // allocate a 3D array:
    atPot[Znum] = (complex_tt*) fftw_malloc(nx*nz/4*sizeof(complex_tt));
    memset(temp,0,nx*nz*sizeof(complex_tt));
    kzmax         = dkz*nz/2.0;
    // define x-and z-position of atom center:
    // The atom should sit in the top-left corner,
    // however (nzPerSlice+1)/2 above zero in z-direction
    xPos = -2.0*PI*0.0; // or muls->resolutionX*nx/(OVERSAMP_X), if in center
    izOffset = (nzPerSlice-1)/2;
    zPos = -2.0*PI*(muls->sliceThickness/nzPerSlice*(izOffset));

    // What this look-up procedure will do is to supply V(r,z) computed from fe(q).
    // Since V(r,z) is rotationally symmetric we might as well compute
    // V(x,y,z) at y=0, i.e. V(x,z).
    // In order to do the proper 3D inverse FT without having to do a complete 3D FFT
    // we will pre-compute the qy-integral for y=0.

    // kzborder = dkz*(nz/(2*OVERSAMP_Z) -1);
    for (iz=0;iz<nz;iz++) {
      kz = dkz*(iz<nz/2 ? iz : iz-nz);        
      // We also need to taper off the potential in z-direction
      // in order to avoid cutoff artifacts.
      // zScale = fabs(kz) <= kzborder ? 1.0 : 0.5+0.5*cos(M_PI*(fabs(kz)-kzborder)/(kzmax-kzborder));
      // printf("iz=%d, kz=%g, zScale=%g ",iz,kz,zScale);
      for (ix=0;ix<nx;ix++) {
        kx = dkx*(ix<nx/2 ? ix : ix-nx);        
        s2 = (kx*kx+kz*kz);
        // if this is within the allowed circle:
        if (s2<smax2) {
          ind3d = ix+iz*nx;
          // f = fe3D(Znum,k2,muls->tds,1.0,muls->scatFactor);
          // multiply scattering factor with Debye-Waller factor:
          // printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
          f = seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF,sqrt(s2))*exp(-s2*B*0.25);
          // perform the qy-integration for qy <> 0:
          for (iy=1;iy<nx;iy++) {
            s3 = dkx*iy;
            s3 = s3*s3+s2;
            if (s3<smax2) {
              f += 2*seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF,sqrt(s3))*exp(-s3*B*0.25);
            }
            else break;
          }
          f *= dkx;
          // note that the factor 2 is missing in the phase (2pi k*r)
          // this places the atoms in the center of the box.
          phase        = kx*xPos + kz*zPos;
          temp[ind3d][0] = f*cos(phase); // *zScale
          temp[ind3d][1] = f*sin(phase); // *zScale
          // if ((kx==0) && (ky==0)) printf(" f=%g (%g, [%g, %g])\n",f,f*zScale,atPot[Znum][ind3d][0],atPot[Znum][ind3d][1]);
        }
      }
    } // for iz ...


#if SHOW_SINGLE_POTENTIAL
    // 0 thickness
    imageio = ImageIOPtr(new CImageIO(nz, nx, 0, dkz, dkx));
    // This scattering factor agrees with Kirkland's scattering factor fe(q)
    sprintf(fileName,"pot_rec_%d.img",Znum);
    imageio->SetThickness(muls->sliceThickness);
    imageio->WriteComplexImage((void**)temp, fileName);
#endif        
    // This converts the 2D kx-kz map of the scattering factor to a 2D real space map.
#if FLOAT_PRECISION ==1
    plan = fftwf_plan_dft_2d(nz,nx,temp,temp,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
#else
    plan = fftw_plan_dft_2d(nz,nx,temp,temp,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
#endif
    // dx2 = muls->resolutionX*muls->resolutionX/(OVERSAMP_X*OVERSAMP_X);
    // dy2 = muls->resolutionY*muls->resolutionY/(OVERSAMP_X*OVERSAMP_X);
    // dz2 = muls->sliceThickness*muls->sliceThickness/(OVERSAMP_Z*OVERSAMP_Z);
    // We also make sure that the potential touches zero at least somewhere. This will avoid
    // sharp edges that could produce ringing artifacts.
    // It is certainly debatable whether this is a good apprach, or not.
    // printf("Setting up %d x %d potential for Znum=%d, (atom kind %d), Omega=%g\n",nx,nz,Znum,iKind,dkx*dky*dkz);
    // min = atPot[Znum][ny/2+nx/2*ny+nz/2*nx*ny][0];
    for (ix=0;ix<nx/2;ix++) for (iz=0;iz<nz/2;iz++) {
        ind3d = ix+iz*nx/2;
        // Integrate over nzPerSlice neighboring layers here:::::::::
        for (zScale=0,iiz=-izOffset;iiz<=izOffset;iiz++) {
          if (iz+izOffset+iiz < nz/2) zScale += temp[ix+(iz+izOffset+iiz)*nx][0];
        }
        if (zScale < 0) zScale = 0;
        // assign the iz-th slice the sum of the 3 other slices:
        // and divide by unit cell volume (since this is in 3D):
        // Where does the '1/2' come from??? OVERSAMP_X*OVERSAMP_Y/8 = 1/2
        // if nothing has changed, then OVERSAMP_X=2 OVERSAMP_Z=18.
        // remember, that s=0.5*k;         
        // This potential will later again be scaled by lambda*gamma (=0.025*1.39139)
        atPot[Znum][ind3d][0] = 47.8658*dkx*dkz/(nz)*zScale;

        // *8*14.4*0.529=4*a0*e (s. Kirkland's book, p. 207)
        // 2*pi*14.4*0.529 = 7.6176;
        // if (atPot[Znum][ind3d][0] < min) min = atPot[Znum][ind3d][0];        
        atPot[Znum][ind3d][1]= 0;
      }
    // make sure we don't produce negative potential:
    // if (min < 0) for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) atPot[Znum][iy+ix*ny][0] -= min;
#if SHOW_SINGLE_POTENTIAL
    imageio = ImageIOPtr(new CImageIO(nz/2, nx/2, 0, muls->sliceThickness/nzPerSlice,
                                      muls->resolutionX/OVERSAMP_X));
    // This scattering factor agrees with Kirkland's scattering factor fe(q)
    imageio->SetThickness(nz*muls->sliceThickness/nzPerSlice);
    sprintf(fileName,"potential_rz_%d.img",Znum);
    ptr = atPot[Znum];
    imageio->WriteComplexImage((void**)ptr, fileName);
#endif        
    if (muls->printLevel > 1) printf("Created 3D (r-z) %d x %d potential array for Z=%d (%d, B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)\n",
                                     nx/2,nz/2,Znum,iKind,B,dkx,dky,dkz,izOffset);
  }
  *Nz_lut = nz/2;
  *nzSub = nzPerSlice;
  *Nr         = nx/2;
  return atPot[Znum];
}


/********************************************************************************
* Lookup function for 3D potential offset due to charged atoms (ions)
********************************************************************************/
fftwf_complex *C3DFFTPotential::GetAtomPotentialOffset3D(int Znum, MULS *muls,double B,int *nzSub,int *Nr,int *Nz_lut,float q) {
  int ix,iy,iz,iiz,ind3d,iKind,izOffset;
  double zScale,kzmax,zPos,xPos;
  fftwf_plan plan;
  static double f,phase,s2,s3,kmax2,kx,kz,dkx,dky,dkz; // ,dx2,dy2,dz2;
  static int nx,ny,nz,nzPerSlice;
  static fftwf_complex **atPot = NULL;
  static fftwf_complex *temp = NULL;
#if SHOW_SINGLE_POTENTIAL == 1
  ImageIOPtr imageio = ImageIOPtr();
  static fftwf_complex *ptr = NULL;
  static char fileName[256];
#endif
  static double *splinb=NULL;
  static double *splinc=NULL;
  static double *splind=NULL;


  // if there is no charge to this atom, return NULL:
  if (q == 0) return NULL;


  // scattering factors in:
  // float scatPar[4][30]
  if (atPot == NULL) {
    splinb = double1D(N_SF, "splinb" );
    splinc = double1D(N_SF, "splinc" );
    splind = double1D(N_SF, "splind" );


    nx = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionX);
    ny = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionY);
    // The FFT-resolution in the z-direction must be high enough to avoid
    // artifacts due to premature cutoff of the rec. space scattering factor
    // we will therefore make it roughly the same as the x-resolution
    // However, we will make sure that a single slice contains an integer number
    // of sampling points.
    nzPerSlice = (int)floor(OVERSAMP_X*muls->sliceThickness/muls->resolutionX);
    // make nzPerSlice odd:
    if (2.0*floor((double)(nzPerSlice >> 1)) == nzPerSlice) nzPerSlice += 1;
    // Total number of z-positions should be twice that of atomRadius/sliceThickness
    nz = (2*(int)ceil(muls->atomRadius/muls->sliceThickness))*nzPerSlice;
    if (muls->printLevel > 1) printf("Potential offset: will use %d sampling points per slice, total nz=%d (%d)\n",nzPerSlice,nz,nzPerSlice >> 1);

    dkx = 0.5*OVERSAMP_X/(nx*muls->resolutionX);
    dky = 0.5*OVERSAMP_X/(ny*muls->resolutionY);
    dkz = nzPerSlice/(double)(nz*muls->sliceThickness);
    kmax2 = 0.5*nx*dkx/(double)OVERSAMP_X;        // largest k that we'll admit

    // printf("Cutoff scattering angle:kmax=%g, smax=%g (1/A), dk=(%g,%g %g)\n",kmax2,S_SCALE*kmax2,dkx,dky,dkz);
    scatParOffs[0][N_SF-1] = 1.2*kmax2;
    scatParOffs[0][N_SF-2] = 1.1*kmax2;
    scatParOffs[0][N_SF-3] = kmax2;
    if (scatParOffs[0][N_SF-4] > scatParOffs[0][N_SF-3]) {

      if (1) {
        // set additional scattering parameters to zero:
        for (ix = 0;ix < N_SF-10;ix++) {
          if (scatParOffs[0][N_SF-4-ix] < scatParOffs[0][N_SF-3]-0.001*(ix+1)) break;
          scatParOffs[0][N_SF-4-ix] = scatParOffs[0][N_SF-3]-0.001*(ix+1);
          for (iy=1; iy<N_ELEM;iy++) scatParOffs[iy][N_SF-4-ix] = 0;        
        }
      }
      else {
        for (ix = 0;ix < 20;ix++) {
          if (scatParOffs[0][N_SF-4-ix] < scatParOffs[0][N_SF-3]) break;
          scatParOffs[0][N_SF-4-ix] = scatParOffs[0][N_SF-3];
          for (iy=1; iy<N_ELEM;iy++) scatParOffs[iy][N_SF-4-ix] = scatParOffs[iy][N_SF-3];        
        }
      }        
      if (muls->printLevel > 1) printf("getAtomPotentialOffset3D: reduced angular range of scattering factor to %g/A!\n",
                                       scatParOffs[0][N_SF-4-ix]);
    } // end of if (scatParOffs[0][N_SF-4] > scatParOffs[0][N_SF-3])
    kmax2 *= kmax2;

    atPot = (fftwf_complex **)malloc((NZMAX+1)*sizeof(fftwf_complex *));
    for (ix=0;ix<=NZMAX;ix++) atPot[ix] = NULL;
    temp = (fftwf_complex*)fftwf_malloc(nx*nz*sizeof(fftwf_complex));
  }
  // initialize this atom, if it has not been done yet:
  if (atPot[Znum] == NULL) {
#if USE_REZ_SFACTS
    iKind = Znum;
#else
    printf("Using charged atoms only works with scattering factors by Rez et al!\n",Znum);
    exit(0);
#endif


    // setup cubic spline interpolation:
    splinh(scatParOffs[0],scatParOffs[iKind],splinb,splinc,splind,N_SF);

    atPot[Znum] = (fftwf_complex*)fftwf_malloc(nx*nz/4*sizeof(fftwf_complex));
    memset(temp,0,nx*nz*sizeof(fftwf_complex));
    kzmax         = dkz*nz/2.0;
    // define x-and z-position of atom center:
    // The atom should sit in the top-left corner,
    // however (nzPerSlice+1)/2 above zero in z-direction
    xPos = -2.0*PI*0.0; // or muls->resolutionX*nx/(OVERSAMP_X), if in center
    izOffset = (nzPerSlice-1)/2;
    zPos = -2.0*PI*(muls->sliceThickness/nzPerSlice*(izOffset));

    // kzborder = dkz*(nz/(2*OVERSAMP_Z) -1);
    for (iz=0;iz<nz;iz++) {
      kz = dkz*(iz<nz/2 ? iz : iz-nz);
      // We also need to taper off the potential in z-direction
      // in order to avoid cutoff artifacts.
      // zScale = fabs(kz) <= kzborder ? 1.0 : 0.5+0.5*cos(M_PI*(fabs(kz)-kzborder)/(kzmax-kzborder));
      // printf("iz=%d, kz=%g, zScale=%g ",iz,kz,zScale);
      for (ix=0;ix<nx;ix++) {
        kx = dkx*(ix<nx/2 ? ix : ix-nx);        
        s2 = (kx*kx+kz*kz);
        // if this is within the allowed circle:
        if (s2<kmax2) {
          ind3d = ix+iz*nx;
          // f = fe3D(Znum,k2,muls->tds,1.0,muls->scatFactor);
          // multiply scattering factor with Debye-Waller factor:
          // printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
          f = seval(scatParOffs[0],scatParOffs[iKind],splinb,splinc,splind,N_SF,sqrt(s2))*exp(-s2*B*0.25);
          // perform the qy-integration for qy <> 0:
          for (iy=1;iy<nx;iy++) {
            s3 = dky*iy;
            s3 = s3*s3+s2;
            if (s3<kmax2) {
              f += 2*seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF,sqrt(s3))*exp(-s3*B*0.25);
            }
            else break;
          }
          f *= dkx;
          // note that the factor 2 is missing in the phase (2pi k*r)
          // this places the atoms in the center of the box.
          phase = kx*xPos + kz*zPos;
          temp[ind3d][0] = f*cos(phase);        // *zScale
          temp[ind3d][1] = f*sin(phase);        // *zScale
          // if ((kx==0) && (ky==0)) printf(" f=%g (%g, [%g, %g])\n",f,f*zScale,atPot[Znum][ind3d][0],atPot[Znum][ind3d][1]);
        }
      }
    } // for iz ...




#if SHOW_SINGLE_POTENTIAL
    imageio = ImageIOPtr(new CImageIO(nz, nx, 0, dkx, dkz, std::vector<double>(),
                                      "rec. space potential"));
    // This scattering factor agrees with Kirkland's scattering factor fe(q)
    imageio->SetThickness(muls->sliceThickness);
    imageio->WriteComplexImage((void**)temp, fileName);
#endif        

    plan = fftwf_plan_dft_2d(nz,nx,temp,temp,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    // dx2 = muls->resolutionX*muls->resolutionX/(OVERSAMP_X*OVERSAMP_X);
    // dy2 = muls->resolutionY*muls->resolutionY/(OVERSAMP_X*OVERSAMP_X);
    // dz2 = muls->sliceThickness*muls->sliceThickness/(OVERSAMP_Z*OVERSAMP_Z);
    // We also make sure that the potential touches zero at least somewhere. This will avoid
    // sharp edges that could produce ringing artifacts.
    // It is certainly debatable whether this is a good apprach, or not.
    // printf("Setting up %d x %d potential for Znum=%d, (atom kind %d), Omega=%g\n",nx,nz,Znum,iKind,dkx*dky*dkz);
    // min = atPot[Znum][ny/2+nx/2*ny+nz/2*nx*ny][0];
    for (ix=0;ix<nx/2;ix++) for (iz=0;iz<nz/2;iz++) {
        ind3d = ix+iz*nx/2;
        // Integrate over nzPerSlice neighboring layers here:::::::::
        for (zScale=0,iiz=-izOffset;iiz<=izOffset;iiz++) {
          if (iz+izOffset+iiz < nz/2) zScale += temp[ix+(iz+izOffset+iiz)*nx][0];
        }
        if (zScale < 0) zScale = 0;
        // assign the iz-th slice the sum of the 3 other slices:
        // and divide by unit cell volume (since this is in 3D):
        // Where does the '1/2' come from??? OVERSAMP_X*OVERSAMP_Y/8 = 1/2
        // if nothing has changed, then OVERSAMP_X=2 OVERSAMP_Z=18.
        // remember, that s=0.5*k;        
        // This potential will later again be scaled by lambda*gamma (=0.025*1.39139)
        atPot[Znum][ind3d][0] = 47.8658*dkx*dkz/(nz)*zScale;
        atPot[Znum][ind3d][1] = 0;
      }
#if SHOW_SINGLE_POTENTIAL
    imageio = ImageIOPtr(new CImageIO(nz/2, nx/2, 0, muls->resolutionX/OVERSAMP_X,
                                      muls->sliceThickness/nzPerSlice));
    // This scattering factor agrees with Kirkland's scattering factor fe(q)
    imageio->SetThickness(nz*muls->sliceThickness/nzPerSlice);
    sprintf(fileName,"potentialOffs_rz_%d.img",Znum);
    ptr = atPot[Znum];
    imageio->WriteComplexImage((void**)ptr, fileName);
#endif        
    if (muls->printLevel > 1) printf("Created 3D (r-z) %d x %d potential offset array for Z=%d (%d, B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)\n",
                                     nx/2,nz/2,Znum,iKind,B,dkx,dky,dkz,izOffset);
  }
  *Nz_lut = nz/2;
  *nzSub = nzPerSlice;
  *Nr         = nx/2;
  return atPot[Znum];
}
// #undef SHOW_SINGLE_POTENTIAL
// end of fftwf_complex *getAtomPotential3D(...)


