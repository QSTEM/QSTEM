#include "pot_2d_fft.hpp"

C2DFFTPotential::C2DFFTPotential(ConfigReaderPtr &configReader) : C2DPotential(configReader)
{
}


void C2DFFTPotential::makeSlices(int nlayer, char *fileName, atom *center)
{
  /* check whether we have constant slice thickness */
  for (unsigned i = 0;i<nlayer;i++) 
	  {
		  if (m_cz[0] != m_cz[i])
		  {
			printf("Warning: slice thickness not constant, will give wrong results (iz=%d)!\n",i);
			break;
		  }
  }
}

void C2DFFTPotential::AddAtomToSlices(std::vector<atom>::iterator &atom, float_tt atomX, float_tt atomY, float_tt atomZ)
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

void C2DFFTPotential::AddAtomNonPeriodic(std::vector<atom>::iterator &atom, 
                         float_tt atomBoxX, unsigned int iAtomX, 
                         float_tt atomBoxY, unsigned int iAtomY, 
                         float_tt atomZ)
{
  unsigned iAtomZ = (int)floor(atomZ/m_sliceThickness);
  unsigned iax0 = iAtomX-m_iRadX < 0 ? 0 : iAtomX-m_iRadX;
  unsigned iax1 = iAtomX+m_iRadX >= m_nx ? m_nx-1 : iAtomX+m_iRadX;
  unsigned iay0 = iAtomY-m_iRadY < 0 ? 0 : iAtomY-m_iRadY;
  unsigned iay1 = iAtomY+m_iRadY >= m_ny ? m_ny-1 : iAtomY+m_iRadY;
  // if within the potential map range:
  if ((iax0 < m_nx) && (iax1 >= 0) && (iay0 < m_ny) && (iay1 >= 0)) 
    {
      float_tt ddx = (-(double)iax0+(atomBoxX/m_dx-(double)m_iRadX))*(double)OVERSAMPLING;
      float_tt ddy = (-(double)iay0+(atomBoxY/m_dy-(double)m_iRadY))*(double)OVERSAMPLING;
      unsigned iOffsX = (int)floor(ddx);
      unsigned iOffsY = (int)floor(ddy);
      ddx -= (double)iOffsX;
      ddy -= (double)iOffsY;
      float_tt s11 = (1-ddx)*(1-ddy);
      float_tt s12 = (1-ddx)*ddy;
      float_tt s21 = ddx*(1-ddy);
      float_tt s22 = ddx*ddy;
      complex_tt *atPotPtr = GetAtomPotential2D(atom->Znum,m_tds ? 0 : atom->dw);
                  
      for (unsigned iax=iax0; iax < iax1; iax++) 
        {
          // printf("(%d, %d): %d,%d\n",iax,nyAtBox,(iOffsX+OVERSAMPLING*(iax-iax0)),iOffsY+iay1-iay0);
          // potPtr and ptr are of type (float *)
          potPtr = &(m_trans[iAtomZ][iax][iay0][0]);
          ptr = &(atPotPtr[(iOffsX+OVERSAMPLING*(iax-iax0))*nyAtBox+iOffsY][0]);
          for (unsigned iay=iay0; iay < iay1; iay++) 
            {
              *potPtr += s11*(*ptr)+s12*(*(ptr+2))+s21*(*(ptr+nyAtBox2))+s22*(*(ptr+nyAtBox2+2));
              
              potPtr++;
              // *potPtr = 0;
              potPtr++;
              ptr += 2*OVERSAMPLING;
            }
        }
    }
}


void C2DFFTPotential::AddAtomPeriodic(std::vector<atom>::iterator &atom, 
                         float_tt atomBoxX, unsigned int iAtomX, 
                         float_tt atomBoxY, unsigned int iAtomY, 
                         float_tt atomZ)
{
    unsigned iAtomZ = (int)floor(atomZ/m_sliceThickness);
    unsigned iax0 = iAtomX-m_iRadX+2*m_nx;
    unsigned iax1 = iAtomX+m_iRadX+2*m_nx;
    unsigned iay0 = iAtomY-m_iRadY+2*m_ny;
    unsigned iay1 = iAtomY+m_iRadY+2*m_ny;

    float_tt ddx = (-(double)iax0+(atomBoxX/m_dx-(double)(m_iRadX-2*m_nx)))*(double)OVERSAMPLING;
    float_tt ddy = (-(double)iay0+(atomBoxY/m_dy-(double)(m_iRadY-2*m_ny)))*(double)OVERSAMPLING;
    unsigned iOffsX = (int)floor(ddx);
    unsigned iOffsY = (int)floor(ddy);
    ddx -= (double)iOffsX;
    ddy -= (double)iOffsY;
                  
    float_tt s22 = (1-ddx)*(1-ddy);
    float_tt s21 = (1-ddx)*ddy;
    float_tt s12 = ddx*(1-ddy);
    float_tt s11 = ddx*ddy;

    complex_tt *atPotPtr = GetAtomPotential2D(atom->Znum,m_tds ? 0 : atom->dw);

    for (unsigned iax=iax0; iax < iax1; iax++) { // TODO: should use ix += OVERSAMPLING
      // printf("(%d, %d): %d,%d\n",iax,nyAtBox,(iOffsX+OVERSAMPLING*(iax-iax0)),iOffsY+iay1-iay0);
      // potPtr and ptr are of type (float *)
      //////////////////
      // Only the exact slice that this atom is located in is affected by this atom:
      int atPosX = (OVERSAMPLING*(iax-iax0)-iOffsX);
      if ((atPosX >= 0) && (atPosX < nyAtBox-1)) {
        ptr = &(atPotPtr[atPosX*nyAtBox-iOffsY][0]);
        for (unsigned iay=iay0; iay < iay1; iay++) {
          // wrap around when end of y-line is reached:
          int atPosY = (iay-iay0)*OVERSAMPLING-iOffsY;
          if ((atPosY < nyAtBox-1) && (atPosY >=0)) {
            // do the real part
            m_trans[iAtomZ][iax % m_nx][iay % m_ny][0] +=
              s11*(*ptr)+s12*(*(ptr+2))+s21*(*(ptr+nyAtBox2))+s22*(*(ptr+nyAtBox2+2));
          }
          // make imaginary part zero for now
          // *potPtr = 0; potPtr++;
          ptr += 2*OVERSAMPLING;
        }
      } // if atPosX within limits
    }
}

#define PHI_SCALE 47.87658
// #define SHOW_SINGLE_POTENTIAL 0
////////////////////////////////////////////////////////////////////////////
// This function should be used yet, because it computes the projected
// potential wrongly, since it doe not yet perform the projection!!!
complex_tt *C2DFFTPotential::GetAtomPotential2D(int Znum, double B) {
        int ix,iy,iz,ind,iKind;
        double min;
        fftwf_plan plan;
        static double f,phase,s2,s3,kmax2,kx,ky,dkx,dky;
        static int nx,ny;
        static fftwf_complex **atPot = NULL;
#if SHOW_SINGLE_POTENTIAL == 1
        ImageIOPtr imageio = ImageIOPtr();
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

                nx = 2*OVERSAMPLING*(int)ceil(m_atomRadius/m_dx);
                ny = 2*OVERSAMPLING*(int)ceil(m_atomRadius/m_dy);
                dkx = 0.5*OVERSAMPLING/((nx)*m_dx);
                dky = 0.5*OVERSAMPLING/((ny)*m_dy);
                kmax2 = 0.5*nx*dkx/(double)OVERSAMPLING; // largest k that we'll admit

                printf("Cutoff scattering angle:kmax=%g (1/A)\n",kmax2);
                scatPar[0][N_SF-1] = 1.2*kmax2;
                scatPar[0][N_SF-2] = 1.1*kmax2;
                scatPar[0][N_SF-3] = kmax2;
                if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3]) {

                        if (1) {
                                // set additional scattering parameters to zero:
                                for (ix = 0;ix < N_SF-10;ix++) {
                                        if (scatPar[0][N_SF-4-ix] < scatPar[0][N_SF-3]-0.001*(ix+1)) break;
                                        scatPar[0][N_SF-4-ix] = scatPar[0][N_SF-3]-0.001*(ix+1);
                                        for (iy=1; iy<N_ELEM;iy++) scatPar[iy][N_SF-4-ix] = 0;        
                                }
                        }

                        if (m_printLevel > 1) printf("getAtomPotential2D: reduced angular range of scattering factor to %g/A!\n",
                                scatPar[0][N_SF-4-ix]);
                } // end of if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3])
                kmax2 *= kmax2;



                atPot = (fftwf_complex **)malloc((NZMAX+1)*sizeof(fftwf_complex *));
                for (ix=0;ix<=NZMAX;ix++) atPot[ix] = NULL;
        }
        // initialize this atom, if it has not been done yet:
        if (atPot[Znum] == NULL) {
                iKind = Znum;
                // setup cubic spline interpolation:
                splinh(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF);

                atPot[Znum] = (fftwf_complex*) fftwf_malloc(nx*ny*sizeof(fftwf_complex));
                // memset(temp,0,nx*nz*sizeof(fftwf_complex));
                memset(atPot[Znum],0,nx*ny*sizeof(fftwf_complex));
                for (ix=0;ix<nx;ix++) {
                        kx = dkx*(ix<nx/2 ? ix : nx-ix);
                        for (iy=0;iy<ny;iy++) {
                                ky = dky*(iy<ny/2 ? iy : ny-iy);
                                s2 = (kx*kx+ky*ky);
                                // if this is within the allowed circle:
                                if (s2<kmax2) {
                                        ind = iy+ix*ny;
                                        // f = fe3D(Znum,k2,m_tds,1.0,m_scatFactor);
                                        // multiply scattering factor with Debye-Waller factor:
                                        // printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
                                        f = seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF,sqrt(s2))*exp(-s2*B*0.25);
                                        phase = PI*(kx*m_dx*nx+ky*m_dy*ny);
                                        atPot[Znum][ind][0] = f*cos(phase);
                                        atPot[Znum][ind][1] = f*sin(phase);
                                }
                        }
                }
#if SHOW_SINGLE_POTENTIAL == 1
                imageio = ImageIOPtr(new CImageIO(ny, nx, 0, dkx, dky, std::vector<double>(),
                "potential"));
                // This scattering factor agrees with Kirkland's scattering factor fe(q)
                imageio->SetThickness(m_sliceThickness);
                imageio->WriteComplexImage((void**)atPot[Znum], fileName);
#endif
                plan = fftwf_plan_dft_2d(nx,ny,atPot[Znum],atPot[Znum],FFTW_BACKWARD,FFTW_ESTIMATE);
                fftwf_execute(plan);
                fftwf_destroy_plan(plan);
                for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
                                atPot[Znum][iy+ix*ny][0] *= dkx*dky*(OVERSAMPLING*OVERSAMPLING);
                }
                // make sure we don't produce negative potential:
                // if (min < 0) for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) atPot[Znum][iy+ix*ny][0] -= min;
#if SHOW_SINGLE_POTENTIAL == 1
                imageio = ImageIOPtr(new CImageIO(nx, ny, 0, m_dx/OVERSAMPLING,
                        m_dy/OVERSAMPLING, std::vector<double>(), "potential"));
                // This scattering factor agrees with Kirkland's scattering factor fe(q)
                //imageio->SetThickness(nz*m_sliceThickness/nzPerSlice);
                sprintf(fileName,"potential_%d.img",Znum);
                imageio->WriteComplexImage((void**)atPot[Znum], fileName);
#endif
                printf("Created 2D %d x %d potential array for Z=%d (%d, B=%g A^2)\n",nx,ny,Znum,iKind,B);
        }
        return atPot[Znum];
}
#undef PHI_SCALE
#undef SHOW_SINGLE_POTENTIAL
