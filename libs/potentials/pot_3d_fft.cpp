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


#include "pot_3d_fft.hpp"

// for memset
#include <cstring>

namespace QSTEM
{

C3DFFTPotential::C3DFFTPotential() : C3DPotential()
{
}

C3DFFTPotential::C3DFFTPotential(const ConfigReaderPtr &configReader) : C3DPotential(configReader)
{
  for (unsigned i = 0;i<m_nslices;i++) 
    {
      if (m_cz[0] != m_cz[i])
        {
          printf("Warning: slice thickness not constant, will give wrong results (iz=%d)!\n",i);
        }
    }
}

void C3DFFTPotential::DisplayParams()
{
  CPotential::DisplayParams();
  printf("* Potential calculation: 3D (FFT method)");
}

void C3DFFTPotential::AddAtomNonPeriodic(std::vector<atom>::iterator &atom, 
                         float_tt atomBoxX, unsigned int iAtomX, 
                         float_tt atomBoxY, unsigned int iAtomY, 
                         float_tt atomZ)
{
  unsigned iAtomZ = (int)floor(atomZ/m_sliceThickness+0.5);
  
  unsigned nzSub, Nr, Nz_lut;
      
  unsigned iax0 = iAtomX-m_iRadX < 0 ? 0 : iAtomX-m_iRadX;
  unsigned iax1 = iAtomX+m_iRadX >= m_nx ? m_nx-1 : iAtomX+m_iRadX;
  unsigned iay0 = iAtomY-m_iRadY < 0 ? 0 : iAtomY-m_iRadY;
  unsigned iay1 = iAtomY+m_iRadY >= m_ny ? m_ny-1 : iAtomY+m_iRadY;
  // if within the potential map range:
  if ((iax0 < m_nx) && (iax1 >= 0) && (iay0 < m_ny) && (iay1 >= 0)) {
    // define range of sampling from atomZ-/+atomRadius
    int iaz0 = iAtomZ-m_iRadZ < 0 ? -iAtomZ : -m_iRadZ;
    unsigned iaz1 = iAtomZ+m_iRadZ >= m_nslices ? m_nslices-iAtomZ-1 : m_iRadZ;
    if ((iAtomZ+iaz0 <        m_nslices) && (iAtomZ+iaz1 >= 0)) {
      // retrieve the pointer for this atom
      ComplexVector atPotPtr;
      GetAtomPotential3D(atom->Znum,m_tds ? 0 : atom->dw,nzSub,Nr,Nz_lut,atPotPtr);
#if USE_Q_POT_OFFSETS
      // retrieve the pointer to the array of charge-dependent potential offset
      // This function will return NULL; if the charge of this atom is zero:
      ComplexVector atPotOffsPtr;
      GetAtomPotentialOffset3D(atom->Znum,muls,m_tds ? 0 : atoms[iatom].dw,nzSub,Nr,Nz_lut,atom->q, atPotOffsPtr);
#endif // USE_Q_POT_OFFSETS
      unsigned iOffsLimHi = Nr*(Nz_lut-1);
      int iOffsLimLo = -Nr*(Nz_lut-1);
      unsigned iOffsStep = nzSub*Nr;

      // Slices around the slice that this atom is located in must be affected by this atom:
      // iaz must be relative to the first slice of the atom potential box.
      for (unsigned iax=iax0; iax <= iax1; iax++) {
        unsigned idx = iax*m_ny+iay0;
        float_tt *potPtr = &(m_trans[iAtomZ+iaz0][idx][0]);

        //////////////////////////////////////////////////////////////////////
        // Computation of Radius must be made faster by using pre-calculated ddx
        // and LUT for sqrt:
        // use of sqrt slows down from 120sec to 180 sec.
        float_tt x2 = iax*m_dx - atomBoxX; x2 *= x2;
        for (unsigned iay=iay0; iay <= iay1; iay++) {
          // printf("iax=%d, iay=%d\n",iax,iay);
          float_tt y2 = iay*m_dy - atomBoxY; y2 *= y2;
          float_tt r = sqrt(x2+y2);
          // r = (x2+y2);
          float_tt ddr = r/m_dr;
          unsigned ir  = (int)floor(ddr);
          // add in different slices, once r has been defined
          if (ir < Nr-1) {
            ddr = ddr-(double)ir;
            float_tt *ptr = potPtr;
            
            float_tt dOffsZ = (iAtomZ+iaz0-atomZ/m_sliceThickness)*nzSub;
#if Z_INTERPOLATION
            unsigned iOffsZ = (unsigned)dOffsZ;
            ddz = fabs(dOffsZ - (float_tt)iOffsZ);
#else
            unsigned iOffsZ = (int)(dOffsZ+0.5);
#endif
            iOffsZ *= Nr;
            
            for (int iaz=iaz0; iaz <= iaz1; iaz++) {
              float_tt potVal = 0;
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
                        potVal += atom->q*((1-ddz)*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0])+
                                                  ddz *((1-ddr)*atPotOffsPtr[ir+iOffsZ+Nr][0]+ddr*atPotOffsPtr[ir+1+iOffsZ+Nr][0]));
                      }
#endif // USE_Q_POT_OFFSETS
#else // Z_INTERPOLATION
                      potVal = (1-ddr)*atPotPtr[ir+iOffsZ][0]+ddr*atPotPtr[ir+1+iOffsZ][0];
#if USE_Q_POT_OFFSETS
                      // add the charge-dependent potential offset
                      if (atPotOffsPtr != NULL) {
                        potVal += atom->q*((1-ddr)*atPotOffsPtr[ir+iOffsZ][0]+ddr*atPotOffsPtr[ir+1+iOffsZ][0]);                                                 
                      }
#endif // USE_Q_POT_OFFSETS
#endif // Z_INTERPOLATION
                                
                    }
                  } // if iOffsZ >=0
                  *ptr += potVal; // ptr = potPtr = m_trans[...]
                            
                  ptr += m_sliceStep;        // advance to the next slice
                  // add the remaining potential to the next slice:
                  // if (iaz < iaz1)        *ptr += (1-ddz)*potVal;
                  iOffsZ += iOffsStep;
                } // end of iaz-loop
              } // if ir < Nr
              // advance pointer to next complex potential point:
              potPtr+=2;
                        
            } // iay=iay0 .. iay1        
          } // iax=iax0 .. iax1
        } // iaz0+iAtomZ < m_slices
      } // if within bounds        
}

void C3DFFTPotential::AddAtomToSlices(std::vector<atom>::iterator &atom, float_tt atomX, float_tt atomY, float_tt atomZ)
{
  unsigned iAtomX = (int)floor(atomX/m_dx);        
  unsigned iAtomY = (int)floor(atomY/m_dy);

  if (m_periodicXY)
    {
      AddAtomPeriodic(atom, atomX, iAtomX, atomY, iAtomY, atomZ);
    }
  else
    {
      AddAtomNonPeriodic(atom, atomX, iAtomX, atomY, iAtomY, atomZ);
    }
}

void C3DFFTPotential::AddAtomPeriodic(std::vector<atom>::iterator &atom, 
                         float_tt atomBoxX, unsigned int iAtomX, 
                         float_tt atomBoxY, unsigned int iAtomY, 
                         float_tt atomZ)
{
  unsigned iAtomZ = (int)floor(atomZ/m_sliceThickness+0.5);
  unsigned iax0 = iAtomX-m_iRadX;
  unsigned iax1 = iAtomX+m_iRadX;
  unsigned iay0 = iAtomY-m_iRadY;
  unsigned iay1 = iAtomY+m_iRadY;

  unsigned nzSub, Nr, Nz_lut;

  // define range of sampling from atomZ-/+atomRadius
  int iaz0 = iAtomZ-m_iRadZ < 0 ? -iAtomZ : -m_iRadZ;
  unsigned iaz1 = iAtomZ+m_iRadZ >= m_nslices ? m_nslices-iAtomZ-1 : m_iRadZ;
  // if (iatom < 2) printf("iatomZ: %d, %d cz=%g, %g: %d, %d\n",iAtomZ,iaz0,m_sliceThickness,atomZ,(int)(-2.5-atomZ),(int)(atomZ+2.5));
  // printf("%d: iatomZ: %d, %d cz=%g, %g\n",iatom,iAtomZ,iaz0,m_sliceThickness,atomZ);
  if ((iAtomZ+iaz0 < m_nslices) && (iAtomZ+iaz1 >= 0)) {
    // retrieve the pointer for this atom
    ComplexVector atPotPtr;
    GetAtomPotential3D(atom->Znum,m_tds ? 0 : atom->dw,nzSub,Nr,Nz_lut, atPotPtr);
#if USE_Q_POT_OFFSETS
    // retrieve the pointer to the array of charge-dependent potential offset
    // This function will return NULL; if the charge of this atom is zero:
    ComplexVector atPotOffsPtr;
    getAtomPotentialOffset3D(atoms[iatom].Znum,muls,m_tds ? 0 : atoms[iatom].dw,&nzSub,&Nr,&Nz_lut,atoms[iatom].q, atPotOffsPtr);
#endif // USE_Q_POT_OFFSETS
    unsigned iOffsLimHi =        Nr*(Nz_lut-1);
    int iOffsLimLo = -Nr*(Nz_lut-1);
    unsigned iOffsStep = nzSub*Nr;

    // Slices around the slice that this atom is located in must be affected by this atom:
    // iaz must be relative to the first slice of the atom potential box.
    for (unsigned iax=iax0; iax < iax1; iax++) {
      unsigned idx = ((iax+2*m_nx) % m_nx)*m_ny + (iay0+2*m_ny) % m_ny;
      float_tt *potPtr = &(m_trans[iAtomZ+iaz0][idx][0]);
      // potPtr = &(m_trans[iAtomZ-iaz0+iaz][iax][iay0][0]);
      float_tt x2 = iax*m_dx - atomBoxX;        x2 *= x2;
      for (unsigned iay=iay0; iay < iay1; ) {
        // printf("iax=%d, iay=%d\n",iax,iay);
        float_tt y2 = iay*m_dy - atomBoxY;        y2 *= y2;
        float_tt r = sqrt(x2+y2);
        float_tt ddr = r/m_dr;
        unsigned ir = (int)floor(ddr);
        // add in different slices, once r has been defined
        if (ir < Nr-1) {
          ddr = ddr-(double)ir;
          float_tt *ptr = potPtr;
          // Include interpolation in z-direction as well (may do it in a very smart way here !):
                        
          float_tt dOffsZ = (iAtomZ+iaz0-atomZ/m_sliceThickness)*nzSub;
#if Z_INTERPOLATION
          unsigned iOffsZ = (int)dOffsZ;
          float_tt ddz         = fabs(dOffsZ - (double)iOffsZ);
#else // Z_INTERPOLATION
          unsigned iOffsZ = (int)(dOffsZ+0.5);
#endif // Z_INTERPOLATION
          iOffsZ *= Nr;
                        
          for (int iaz=iaz0; iaz <= iaz1; iaz++) {
            float_tt potVal = 0;
            // iOffsZ = (int)(fabs(iAtomZ+iaz-atomZ/m_sliceThickness)*nzSub+0.5);
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
            
            ptr += m_sliceStep; // advance to the next slice
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
        if (++iay % m_ny == 0) potPtr -= 2*m_ny;                
      }
    }
  } // iaz0+iAtomZ < m_slices
}


/********************************************************************************
* Create Lookup table for 3D potential due to neutral atoms
********************************************************************************/
#define PHI_SCALE 47.87658
void C3DFFTPotential::GetAtomPotential3D(unsigned Znum, float_tt B,
                                                unsigned &nzSub,unsigned &Nr,unsigned &Nz_lut,
                                                ComplexVector &output) {
	/*
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
  */
  
  float_tt dkx, dkz;
  // largest k that we'll admit
  float_tt kmax2;

  ComplexVector temp;

  unsigned nx, ny, nz, nzPerSlice;
  std::vector<float_tt> splinb(N_SF), splinc(N_SF), splind(N_SF);

  // scattering factors in:
  // float scatPar[4][30]
  // TODO: this section meant to speed things up by skipping initialization.  Can we keep values
  //    as members instead, and move init to another function?
  if (m_atPot.size()==0) {
    nx = 2*OVERSAMPLING*(int)ceil(m_atomRadius/m_dx);
    ny = 2*OVERSAMPLING*(int)ceil(m_atomRadius/m_dy);
    temp.resize(nx*nz);
    // The FFT-resolution in the z-direction must be high enough to avoid
    // artifacts due to premature cutoff of the rec. space scattering factor
    // we will therefore make it roughly the same as the x-resolution
    // However, we will make sure that a single slice contains an integer number
    // of sampling points.
    nzPerSlice = (int)floor(OVERSAMPLING*m_sliceThickness/m_dx);
    // make nzPerSlice odd:
    if (2.0*(nzPerSlice >> 1) == nzPerSlice) nzPerSlice += 1;
    // Total number of z-positions should be twice that of atomRadius/sliceThickness
    nz = (2*(int)ceil(m_atomRadius/m_sliceThickness))*nzPerSlice;
    if (m_printLevel > 1) printf("Will use %d sampling points per slice, total nz=%d (%d)\n",nzPerSlice,nz,nzPerSlice >> 1);

    dkx = 0.5*OVERSAMPLING/(nx*m_dx); // nx*m_dx is roughly 2*m_atomRadius
    dkz = nzPerSlice/(nz*m_sliceThickness);
    kmax2 = 0.5*nx*dkx/OVERSAMPLING; 
    // Don't square kmax2 yet!

    printf("dkx = %g, nx = %d, kmax2 = %g\n",dkx,nx,kmax2);
    if (m_printLevel > 1) printf("Cutoff scattering angle: kmax=%g (1/A), dk=(%g,%g %g)\n",kmax2,dkx,dkx,dkz);
    scatPar[0][N_SF-1] = 1.2*kmax2;
    scatPar[0][N_SF-2] = 1.1*kmax2;
    scatPar[0][N_SF-3] = kmax2;
    // adjust the resolution of the lookup table if necessary
    if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3]) {
      unsigned ix=0;
      if (1) {
        // set additional scattering parameters to zero:
        for (ix;ix < N_SF-10;ix++) {
          if (scatPar[0][N_SF-4-ix] < scatPar[0][N_SF-3]-0.001*(ix+1)) break;
          scatPar[0][N_SF-4-ix] = scatPar[0][N_SF-3]-0.001*(ix+1);
          for (unsigned iy=1; iy<N_ELEM;iy++) scatPar[iy][N_SF-4-ix] = 0;
        }
      }
      else {
        for (ix;ix < 20;ix++) {
          if (scatPar[0][N_SF-4-ix] < scatPar[0][N_SF-3]) break;
          scatPar[0][N_SF-4-ix] = scatPar[0][N_SF-3];
          for (unsigned iy=1; iy<N_ELEM;iy++) scatPar[iy][N_SF-4-ix] = scatPar[iy][N_SF-3];        
        }
      }                
      if (m_printLevel > 1) printf("getAtomPotential3D: set resolution of scattering factor to %g/A!\n",
                                       scatPar[0][N_SF-4-ix]);
    }        // end of if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3])
    // allocate a list of pointers for the element-specific potential lookup table
    //atPot = (fftwf_complex **)malloc((N_ELEM+1)*sizeof(fftwf_complex *));
    //for (unsigned ix=0;ix<=N_ELEM;ix++) atPot[ix] = NULL;
    //temp = (complex_tt*) fftw_malloc(nx*nz*sizeof(complex_tt));
  }

  // Now kmax2 is actually kmax**2.
  kmax2 *= kmax2;

  // initialize this atom, if it has not been done yet:
  if (m_atPot.count(Znum) == 0) {
    // setup cubic spline interpolation:
    splinh(scatPar[0],scatPar[Znum],&splinb[0],&splinc[0],&splind[0],N_SF);
    
    // allocate a 3D array:
    m_atPot[Znum]=ComplexVector(nx*nz/4);
    //atPot[Znum] = (complex_tt*) fftw_malloc(nx*nz/4*sizeof(complex_tt));
    memset((void*)&temp[0],0,nx*nz*sizeof(complex_tt));
    float_tt kzmax         = dkz*nz/2.0;
    // define x-and z-position of atom center:
    // The atom should sit in the top-left corner,
    // however (nzPerSlice+1)/2 above zero in z-direction
    float_tt xPos = -2.0*PI*0.0; // or m_dx*nx/(OVERSAMPLING), if in center
    unsigned izOffset = (nzPerSlice-1)/2;
    float_tt zPos = -2.0*PI*(m_sliceThickness/nzPerSlice*(izOffset));

    // What this look-up procedure will do is to supply V(r,z) computed from fe(q).
    // Since V(r,z) is rotationally symmetric we might as well compute
    // V(x,y,z) at y=0, i.e. V(x,z).
    // In order to do the proper 3D inverse FT without having to do a complete 3D FFT
    // we will pre-compute the qy-integral for y=0.

    // kzborder = dkz*(nz/(2*OVERSAMP_Z) -1);
    for (unsigned iz=0;iz<nz;iz++) {
      float_tt kz = dkz*(iz<nz/2 ? iz : iz-nz);        
      // We also need to taper off the potential in z-direction
      // in order to avoid cutoff artifacts.
      // zScale = fabs(kz) <= kzborder ? 1.0 : 0.5+0.5*cos(M_PI*(fabs(kz)-kzborder)/(kzmax-kzborder));
      // printf("iz=%d, kz=%g, zScale=%g ",iz,kz,zScale);
      for (unsigned ix=0;ix<nx;ix++) {
        float_tt kx = dkx*(ix<nx/2 ? ix : ix-nx);        
        unsigned s2 = (kx*kx+kz*kz);
        // if this is within the allowed circle:
        if (s2<kmax2) {
          unsigned ind3d = ix+iz*nx;
          // f = fe3D(Znum,k2,m_tds,1.0,m_scatFactor);
          // multiply scattering factor with Debye-Waller factor:
          // printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
          float_tt f = seval(scatPar[0],scatPar[Znum],&splinb[0],&splinc[0],&splind[0],N_SF,sqrt(s2))*exp(-s2*B*0.25);
          // perform the qy-integration for qy <> 0:
          for (unsigned iy=1;iy<nx;iy++) {
            float_tt s3 = dkx*iy;
            s3 = s3*s3+s2;
            if (s3<kmax2) {
              f += 2*seval(scatPar[0],scatPar[Znum],&splinb[0],&splinc[0],&splind[0],N_SF,sqrt(s3))*exp(-s3*B*0.25);
            }
            else break;
          }
          f *= dkx;
          // note that the factor 2 is missing in the phase (2pi k*r)
          // this places the atoms in the center of the box.
          float_tt phase        = kx*xPos + kz*zPos;
          temp[ind3d]=std::complex<float_tt>(f*cos(phase), // *zScale
                                   f*sin(phase)); // *zScale
          // if ((kx==0) && (ky==0)) printf(" f=%g (%g, [%g, %g])\n",f,f*zScale,atPot[Znum][ind3d][0],atPot[Znum][ind3d][1]);
        }
      }
    } // for iz ...


#if SHOW_SINGLE_POTENTIAL
    // 0 thickness
    imageio = ImageIOPtr(new CImageIO(nz, nx, 0, dkz, dkx));
    // This scattering factor agrees with Kirkland's scattering factor fe(q)
    sprintf(fileName,"pot_rec_%d.img",Znum);
    imageio->SetThickness(m_sliceThickness);
    imageio->WriteImage((void**)temp, fileName);
#endif        
    // This converts the 2D kx-kz map of the scattering factor to a 2D real space map.
#if FLOAT_PRECISION ==1
    fftwf_complex *ptr=(fftwf_complex *)&temp[0];
    fftwf_plan plan = fftwf_plan_dft_2d(nz,nx,ptr,ptr,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
#else
    fftw_complex *ptr=(fftw_complex *)&temp[0];
    fftw_plan plan = fftw_plan_dft_2d(nz,nx,ptr,ptr,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
#endif
    // We also make sure that the potential touches zero at least somewhere. This will avoid
    // sharp edges that could produce ringing artifacts.
    // It is certainly debatable whether this is a good apprach, or not.
    // printf("Setting up %d x %d potential for Znum=%d, (atom kind %d), Omega=%g\n",nx,nz,Znum,iKind,dkx*dky*dkz);
    // min = atPot[Znum][ny/2+nx/2*ny+nz/2*nx*ny][0];
    for (unsigned ix=0;ix<nx/2;ix++) for (unsigned iz=0;iz<nz/2;iz++) {
        float_tt zScale=0;
        unsigned ind3d = ix+iz*nx/2;
        // Integrate over nzPerSlice neighboring layers here:::::::::
        for (int iiz=-izOffset;iiz<=izOffset;iiz++) {
          if (iz+izOffset+iiz < nz/2) zScale += temp[ix+(iz+izOffset+iiz)*nx][0];
        }
        if (zScale < 0) zScale = 0;
        // assign the iz-th slice the sum of the 3 other slices:
        // and divide by unit cell volume (since this is in 3D):
        // Where does the '1/2' come from??? OVERSAMPLING*OVERSAMP_Y/8 = 1/2
        // if nothing has changed, then OVERSAMPLING=2 OVERSAMP_Z=18.
        // remember, that s=0.5*k;         
        // This potential will later again be scaled by lambda*gamma (=0.025*1.39139)
        // Sets real value to this; imaginary value to 0
        m_atPot[Znum][ind3d] = 47.8658*dkx*dkz/(nz)*zScale;
      }
    // make sure we don't produce negative potential:
    if (m_printLevel > 1) printf("Created 3D (r-z) %d x %d potential array for Z=%d (B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)\n",
                                     nx/2,nz/2,Znum,B,dkx,dkx,dkz,izOffset);
  }
  Nz_lut = nz/2;
  nzSub  = nzPerSlice;
  Nr     = nx/2;
  output = m_atPot[Znum];
}


/********************************************************************************
* Lookup function for 3D potential offset due to charged atoms (ions)
********************************************************************************/
void C3DFFTPotential::GetAtomPotentialOffset3D(unsigned Znum, float_tt B,unsigned &nzSub,unsigned &Nr,unsigned &Nz_lut,float_tt q, ComplexVector &output) {
  /*
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
  */

  // if there is no charge to this atom, return NULL:
  if (q == 0) return;
#if !USE_REZ_SFACTS
    printf("Using charged atoms only works with scattering factors by Rez et al!\n",Znum);
    exit(0);
#endif


  unsigned nx, ny, nz, nzPerSlice;
  float_tt dkx, dky, dkz;
  std::vector<float_tt> splinb(N_SF), splinc(N_SF), splind(N_SF);
  float_tt kmax2;

  ComplexVector temp;

  // scattering factors in:
  // float scatPar[4][30]
  if (m_offsetPot.size() ==0) {

    nx = 2*OVERSAMPLING*(int)ceil(m_atomRadius/m_dx);
    ny = 2*OVERSAMPLING*(int)ceil(m_atomRadius/m_dy);
    temp.resize(nx*nz);

    // The FFT-resolution in the z-direction must be high enough to avoid
    // artifacts due to premature cutoff of the rec. space scattering factor
    // we will therefore make it roughly the same as the x-resolution
    // However, we will make sure that a single slice contains an integer number
    // of sampling points.
    nzPerSlice = (int)floor(OVERSAMPLING*m_sliceThickness/m_dx);
    // make nzPerSlice odd:
    if (2.0*floor((double)(nzPerSlice >> 1)) == nzPerSlice) nzPerSlice += 1;
    // Total number of z-positions should be twice that of atomRadius/sliceThickness
    nz = (2*(int)ceil(m_atomRadius/m_sliceThickness))*nzPerSlice;
    if (m_printLevel > 1) printf("Potential offset: will use %d sampling points per slice, total nz=%d (%d)\n",nzPerSlice,nz,nzPerSlice >> 1);

    dkx = 0.5*OVERSAMPLING/(nx*m_dx);
    dky = 0.5*OVERSAMPLING/(ny*m_dy);
    dkz = nzPerSlice/(double)(nz*m_sliceThickness);
    kmax2 = 0.5*nx*dkx/(double)OVERSAMPLING;        // largest k that we'll admit

    // printf("Cutoff scattering angle:kmax=%g, smax=%g (1/A), dk=(%g,%g %g)\n",kmax2,S_SCALE*kmax2,dkx,dky,dkz);
    scatParOffs[0][N_SF-1] = 1.2*kmax2;
    scatParOffs[0][N_SF-2] = 1.1*kmax2;
    scatParOffs[0][N_SF-3] = kmax2;
    if (scatParOffs[0][N_SF-4] > scatParOffs[0][N_SF-3]) {
      unsigned ix=0;
      if (1) {
        // set additional scattering parameters to zero:
        for (ix;ix < N_SF-10;ix++) {
          if (scatParOffs[0][N_SF-4-ix] < scatParOffs[0][N_SF-3]-0.001*(ix+1)) break;
          scatParOffs[0][N_SF-4-ix] = scatParOffs[0][N_SF-3]-0.001*(ix+1);
          for (unsigned iy=1; iy<N_ELEM;iy++) scatParOffs[iy][N_SF-4-ix] = 0;        
        }
      }
      else {
        for (ix;ix < 20;ix++) {
          if (scatParOffs[0][N_SF-4-ix] < scatParOffs[0][N_SF-3]) break;
          scatParOffs[0][N_SF-4-ix] = scatParOffs[0][N_SF-3];
          for (unsigned iy=1; iy<N_ELEM;iy++) scatParOffs[iy][N_SF-4-ix] = scatParOffs[iy][N_SF-3];        
        }
      }        
      if (m_printLevel > 1) printf("getAtomPotentialOffset3D: reduced angular range of scattering factor to %g/A!\n",
                                       scatParOffs[0][N_SF-4-ix]);
    } // end of if (scatParOffs[0][N_SF-4] > scatParOffs[0][N_SF-3])
    kmax2 *= kmax2;

    //atPot = (fftwf_complex **)malloc((N_ELEM+1)*sizeof(fftwf_complex *));
    //for (unsigned ix=0;ix<=N_ELEM;ix++) atPot[ix] = NULL;
  }
  // initialize this atom, if it has not been done yet:
  if (m_offsetPot.count(Znum)==0) {
    // setup cubic spline interpolation:
    splinh(scatParOffs[0],scatParOffs[Znum],&splinb[0],&splinc[0],&splind[0],N_SF);

    m_offsetPot[Znum] = ComplexVector(nx*nz/4);
    memset((void*)&temp[0],0,nx*nz*sizeof(complex_tt));

    //complex_tt *temp=(complex_tt *)fftw_malloc(nx*nz*sizeof(complex_tt));
    //memset(temp, 0, nx*nz*sizeof(complex_tt));

    float_tt kzmax         = dkz*nz/2.0;
    // define x-and z-position of atom center:
    // The atom should sit in the top-left corner,
    // however (nzPerSlice+1)/2 above zero in z-direction
    float_tt xPos = -2.0*PI*0.0; // or m_dx*nx/(OVERSAMPLING), if in center
    unsigned izOffset = (nzPerSlice-1)/2;
    float_tt zPos = -2.0*PI*(m_sliceThickness/nzPerSlice*(izOffset));

    // kzborder = dkz*(nz/(2*OVERSAMP_Z) -1);
    for (unsigned iz=0;iz<nz;iz++) {
      float_tt kz = dkz*(iz<nz/2 ? iz : iz-nz);
      // We also need to taper off the potential in z-direction
      // in order to avoid cutoff artifacts.
      // zScale = fabs(kz) <= kzborder ? 1.0 : 0.5+0.5*cos(M_PI*(fabs(kz)-kzborder)/(kzmax-kzborder));
      // printf("iz=%d, kz=%g, zScale=%g ",iz,kz,zScale);
      for (unsigned ix=0;ix<nx;ix++) {
        float_tt kx = dkx*(ix<nx/2 ? ix : ix-nx);        
        float_tt s2 = (kx*kx+kz*kz);
        // if this is within the allowed circle:
        if (s2<kmax2) {
          unsigned ind3d = ix+iz*nx;
          // f = fe3D(Znum,k2,m_tds,1.0,m_scatFactor);
          // multiply scattering factor with Debye-Waller factor:
          // printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
          float_tt f = seval(scatParOffs[0],scatParOffs[Znum],&splinb[0],&splinc[0],&splind[0],N_SF,
                             sqrt(s2))*exp(-s2*B*0.25);
          // perform the qy-integration for qy <> 0:
          for (unsigned iy=1;iy<nx;iy++) {
            float_tt s3 = dky*iy;
            s3 = s3*s3+s2;
            if (s3<kmax2) {
              f += 2*seval(scatPar[0],scatPar[Znum],&splinb[0],&splinc[0],&splind[0],N_SF,sqrt(s3))*exp(-s3*B*0.25);
            }
            else break;
          }
          f *= dkx;
          // note that the factor 2 is missing in the phase (2pi k*r)
          // this places the atoms in the center of the box.
          float_tt phase = kx*xPos + kz*zPos;
          temp[ind3d][0] = f*cos(phase);        // *zScale
          temp[ind3d][1] = f*sin(phase);        // *zScale
          // if ((kx==0) && (ky==0)) printf(" f=%g (%g, [%g, %g])\n",f,f*zScale,atPot[Znum][ind3d][0],atPot[Znum][ind3d][1]);
        }
      }
    } // for iz....


#if SHOW_SINGLE_POTENTIAL
    imageio = ImageIOPtr(new CImageIO(nz, nx, 0, dkx, dkz, std::vector<double>(),
                                      "rec. space potential"));
    // This scattering factor agrees with Kirkland's scattering factor fe(q)
    imageio->SetThickness(m_sliceThickness);
    imageio->WriteImage((void**)temp, fileName);
#endif        

#if FLOAT_PRECISION ==1
    fftwf_complex *ptr=(fftwf_complex *)&temp[0];
    fftwf_plan plan = fftwf_plan_dft_2d(nz,nx,ptr,ptr,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
#else
    fftw_complex *ptr=(fftw_complex *)&temp[0];
    fftw_plan plan = fftw_plan_dft_2d(nz,nx,ptr,ptr,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
#endif
    // We also make sure that the potential touches zero at least somewhere. This will avoid
    // sharp edges that could produce ringing artifacts.
    // It is certainly debatable whether this is a good apprach, or not.
    // printf("Setting up %d x %d potential for Znum=%d, (atom kind %d), Omega=%g\n",nx,nz,Znum,iKind,dkx*dky*dkz);
    // min = atPot[Znum][ny/2+nx/2*ny+nz/2*nx*ny][0];
    for (unsigned ix=0;ix<nx/2;ix++) for (unsigned iz=0;iz<nz/2;iz++) {
        unsigned ind3d = ix+iz*nx/2;
        float_tt zScale=0;
        // Integrate over nzPerSlice neighboring layers here:::::::::
        for (int iiz=-izOffset;iiz<=izOffset;iiz++) {
          if (iz+izOffset+iiz < nz/2) zScale += temp[ix+(iz+izOffset+iiz)*nx][0];
        }
        if (zScale < 0) zScale = 0;
        // assign the iz-th slice the sum of the 3 other slices:
        // and divide by unit cell volume (since this is in 3D):
        // Where does the '1/2' come from??? OVERSAMPLING*OVERSAMP_Y/8 = 1/2
        // if nothing has changed, then OVERSAMPLING=2 OVERSAMP_Z=18.
        // remember, that s=0.5*k;        
        // This potential will later again be scaled by lambda*gamma (=0.025*1.39139)
        //    implicitly sets imaginary part to 0.
        m_offsetPot[Znum][ind3d] = 47.8658*dkx*dkz/(nz)*zScale;
      }
#if SHOW_SINGLE_POTENTIAL
    imageio = ImageIOPtr(new CImageIO(nz/2, nx/2, 0, m_dx/OVERSAMPLING,
                                      m_sliceThickness/nzPerSlice));
    // This scattering factor agrees with Kirkland's scattering factor fe(q)
    imageio->SetThickness(nz*m_sliceThickness/nzPerSlice);
    sprintf(fileName,"potentialOffs_rz_%d.img",Znum);
    ptr = atPot[Znum];
    imageio->WriteImage((void**)ptr, fileName);
#endif        
    if (m_printLevel > 1) printf("Created 3D (r-z) %d x %d potential offset array for Z=%d (B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)\n",
                                     nx/2,nz/2,Znum,B,dkx,dky,dkz,izOffset);
  } // end of creating atom if not exist...
  Nz_lut = nz/2;
  nzSub  = nzPerSlice;
  Nr     = nx/2;
  output = m_offsetPot[Znum];
}
// #undef SHOW_SINGLE_POTENTIAL
// end of fftwf_complex *getAtomPotential3D(...)

} // end namespace QSTEM