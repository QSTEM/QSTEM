
from traits.api import HasTraits

class AtomPotential(HasTraits):
    

    def get3DPotential(self):
        pass

"""
Create Lookup table for 3D potential due to neutral atoms
"""
#define PHI_SCALE 47.87658
fftwf_complex *getAtomPotential3D(int Znum, MULS *muls,double B,int *nzSub,int *Nr,int *Nz_lut) {
	int ix,iy,iz,iiz,ind3d,iKind,izOffset;
	double zScale,kzmax,zPos,xPos;
	fftwf_plan plan;
	static double f,phase,s2,s3,kmax2,smax2,kx,kz,dkx,dky,dkz; // ,dx2,dy2,dz2;
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

		dkx = 0.5*OVERSAMP_X/(nx*muls->resolutionX);  // nx*muls->resolutionX is roughly 2*muls->atomRadius
		dky = dkx;                                    
		dkz = nzPerSlice/(double)(nz*muls->sliceThickness);
		kmax2 = 0.5*nx*dkx/(double)OVERSAMP_X;  # largest k that we'll admit
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
		}	// end of if (scatPar[0][N_SF-4] > scatPar[0][N_SF-3])
		smax2 *= smax2;
		kmax2 *= kmax2;
		// allocate a list of pointers for the element-specific potential lookup table
		atPot = (fftwf_complex **)malloc((NZMAX+1)*sizeof(fftwf_complex *));
		for (ix=0;ix<=NZMAX;ix++) atPot[ix] = NULL;
		temp  = (fftwf_complex*) fftwf_malloc(nx*nz*sizeof(fftwf_complex));
	}
	// initialize this atom, if it has not been done yet:
	if (atPot[Znum] == NULL) {
		iKind = Znum;

		// setup cubic spline interpolation:
		splinh(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF);

		// allocate a 3D array:
		atPot[Znum] = (fftwf_complex*) fftwf_malloc(nx*nz/4*sizeof(fftwf_complex));
		memset(temp,0,nx*nz*sizeof(fftwf_complex));
		kzmax	  = dkz*nz/2.0; 
		// define x-and z-position of atom center:
		// The atom should sit in the top-left corner, 
		// however (nzPerSlice+1)/2 above zero in z-direction
		xPos = -2.0*PI*0.0;  // or muls->resolutionX*nx/(OVERSAMP_X), if in center
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
					// multiply scattering factor with Debye-Waller factor:
					// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
                                        # evaluate the spline
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
					phase	= kx*xPos + kz*zPos;
					temp[ind3d][0] = f*cos(phase);  // *zScale
					temp[ind3d][1] = f*sin(phase);  // *zScale
					// if ((kx==0) && (ky==0)) printf(" f=%g (%g, [%g, %g])\n",f,f*zScale,atPot[Znum][ind3d][0],atPot[Znum][ind3d][1]);
				}
			}
		} // for iz ...

#if SHOW_SINGLE_POTENTIAL
		imageio = ImageIOPtr(new CImageIO(nx, nz, dkz, dkx));
		# This scattering factor agrees with Kirkland's scattering factor fe(q)
		sprintf(fileName,"pot_rec_%d.img",Znum);
		imageio->SetThickness(muls->sliceThickness);
		imageio->WriteComplexImage((void**)temp, fileName);
#endif	  
		# This converts the 2D kx-kz  map of the scattering factor to a 2D real space map.
                real_space_map = ifft2(temp)
		#plan = fftwf_plan_dft_2d(nz,nx,temp,temp,FFTW_BACKWARD,FFTW_ESTIMATE);

		// We also make sure that the potential touches zero at least somewhere.  This will avoid 
		// sharp edges that could produce ringing artifacts.  
		// It is certainly debatable whether this is a good apprach, or not. 
		// printf("Setting up %d x %d potential for Znum=%d, (atom kind %d), Omega=%g\n",nx,nz,Znum,iKind,dkx*dky*dkz);
		// min = atPot[Znum][ny/2+nx/2*ny+nz/2*nx*ny][0];
		for (ix=0;ix<nx/2;ix++)  for (iz=0;iz<nz/2;iz++) {
			ind3d = ix+iz*nx/2;
			// Integrate over nzPerSlice neighboring layers here:::::::::
			for (zScale=0,iiz=-izOffset;iiz<=izOffset;iiz++) {
				if (iz+izOffset+iiz < nz/2) zScale += temp[ix+(iz+izOffset+iiz)*nx][0];
			}
			if (zScale < 0) zScale = 0;
			// assign the iz-th slice the sum of the 3 other slices:
			// and divide by unit cell volume (since this is in 3D):
			// Where does the '1/2' come from???  OVERSAMP_X*OVERSAMP_Y/8 = 1/2
			// if nothing has changed, then OVERSAMP_X=2 OVERSAMP_Z=18.
			// remember, that s=0.5*k; 	
			// This potential will later again be scaled by lambda*gamma (=0.025*1.39139)
			atPot[Znum][ind3d][0] = 47.8658*dkx*dkz/(nz)*zScale; 

			atPot[Znum][ind3d][1]= 0;
		}
#if SHOW_SINGLE_POTENTIAL
		imageio = ImageIOPtr(new CImageIO(nz/2, nx/2, muls->sliceThickness/nzPerSlice, 
			muls->resolutionX/OVERSAMP_X));
		# This scattering factor agrees with Kirkland's scattering factor fe(q)
		imageio->SetThickness(nz*muls->sliceThickness/nzPerSlice);
		sprintf(fileName,"potential_rz_%d.img",Znum);
		ptr = atPot[Znum];
		imageio->WriteComplexImage((void**)ptr, fileName);
#endif	  
		if (muls->printLevel > 1) printf("Created 3D (r-z) %d x %d potential array for Z=%d (%d, B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)\n",
			nx/2,nz/2,Znum,iKind,B,dkx,dky,dkz,izOffset);
	}
	*Nz_lut = nz/2;
	*nzSub  = nzPerSlice;
	*Nr	  = nx/2;
	return atPot[Znum];
}
