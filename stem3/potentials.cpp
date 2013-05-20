#include "potentials.h"
#include "defines.h"

#define NRMAX	50	/* number of values in look-up-table in vzatomLUT */
#define RMIN	0.01	/* min r (in Ang) range of LUT for vzatomLUT() */
#define RMAX	5

//
std::vector<float_tt> Potential::m_splinr;

// initialization for 3D pots
void Potential::Potential3D(MULS *muls)
{
	int ix,iy,iz,iiz,ind3d,iKind,izOffset;
	double zScale,kzmax,zPos,xPos;
	fftwf_plan plan;
	static double f,phase,s2,s3,kmax2,smax2,kx,kz,dkx,dky,dkz; // ,dx2,dy2,dz2;
	static int nx,ny,nz,nzPerSlice;

	// scattering factors in:
	// float scatPar[4][30]

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

		atPot = QScMat(nx*nz, NZMAX);

		dkx = 0.5*OVERSAMP_X/(nx*muls->resolutionX);  // nx*muls->resolutionX is roughly 2*muls->atomRadius
		// dky = 0.5*OVERSAMP_X/(ny*muls->resolutionY);
		dky = dkx;                                    
		dkz = nzPerSlice/(double)(nz*muls->sliceThickness);
		kmax2 = 0.5*nx*dkx/(double)OVERSAMP_X;  // largest k that we'll admit
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
		//atPot = (fftwf_complex **)malloc((NZMAX+1)*sizeof(fftwf_complex *));
		for (ix=0;ix<=NZMAX;ix++) atPot[ix] = NULL;
		temp  = (fftwf_complex*) fftwf_malloc(nx*nz*sizeof(fftwf_complex));
}

//
void Potential::CreateAtPot(int Znum, bool useOffset)
{
	int ix,iy,iz,iiz,ind3d,iKind,izOffset;
	double zScale,kzmax,zPos,xPos;
	fftwf_plan plan;
	static double f,phase,s2,s3,kmax2,smax2,kx,kz,dkx,dky,dkz; // ,dx2,dy2,dz2;
	static int nx,ny,nz,nzPerSlice;

	iKind = Znum;

	// setup cubic spline interpolation:
	splinh(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF);

	// allocate a 3D array:
	//atPot[Znum] = QScMat(nx, nz); //(fftwf_complex*) fftwf_malloc(nx*nz/4*sizeof(fftwf_complex));
	//memset(temp,0,nx*nz*sizeof(fftwf_complex));
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
				phase	= kx*xPos + kz*zPos;
				temp[ind3d] = QScf(f*cos(phase), f*sin(phase)); // *zScale
				// if ((kx==0) && (ky==0)) printf(" f=%g (%g, [%g, %g])\n",f,f*zScale,atPot[Znum][ind3d].real(),atPot[Znum][ind3d][1]);
			}
		}
	} // for iz ...

#if SHOW_SINGLE_POTENTIAL
	// This scattering factor agrees with Kirkland's scattering factor fe(q)
	if (header == NULL) 
		header = makeNewHeaderCompact(1,nz,nx,0,dkx,dky,0,NULL,"rec. space potential");
	header->nx = nz;
	header->ny = nx;
	header->dx = dkz;
	header->dy = dkx;

	sprintf(fileName,"pot_rec_%d.img",Znum);
	header->t = muls->sliceThickness;
	writeImage((void **)&temp,header,fileName);
#endif	  
	// This converts the 2D kx-kz  map of the scattering factor to a 2D real space map.
	plan = fftwf_plan_dft_2d(nz,nx,(fftwf_complex*)temp.data(),(fftwf_complex*)temp.data(),FFTW_BACKWARD,FFTW_ESTIMATE);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	// dx2 = muls->resolutionX*muls->resolutionX/(OVERSAMP_X*OVERSAMP_X);
	// dy2 = muls->resolutionY*muls->resolutionY/(OVERSAMP_X*OVERSAMP_X);
	// dz2 = muls->sliceThickness*muls->sliceThickness/(OVERSAMP_Z*OVERSAMP_Z);
	// We also make sure that the potential touches zero at least somewhere.  This will avoid 
	// sharp edges that could produce ringing artifacts.  
	// It is certainly debatable whether this is a good apprach, or not. 
	// printf("Setting up %d x %d potential for Znum=%d, (atom kind %d), Omega=%g\n",nx,nz,Znum,iKind,dkx*dky*dkz);
	// min = atPot[Znum][ny/2+nx/2*ny+nz/2*nx*ny].real();
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
		atPot[Znum][ind3d] = QScf(47.8658*dkx*dkz/(nz)*zScale,0); 
			// *8*14.4*0.529=4*a0*e (s. Kirkland's book, p. 207)
		// 2*pi*14.4*0.529 = 7.6176;
		// if (atPot[Znum][ind3d].real() < min) min = atPot[Znum][ind3d].real();
		//atPot[Znum][ind3d][1]= 0;
	}
	// make sure we don't produce negative potential:
	// if (min < 0) for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) atPot[Znum][iy+ix*ny].real() -= min;
#if SHOW_SINGLE_POTENTIAL
		if (header == NULL) 
			header = makeNewHeaderCompact(1,nz/2,nx/2,0,muls->resolutionX/OVERSAMP_X,
			muls->sliceThickness/nzPerSlice,0,NULL,"potential");
		header->nx = nz/2;
		header->ny = nx/2;
		header->dx = muls->sliceThickness/nzPerSlice;
		header->dy = muls->resolutionX/OVERSAMP_X;
		header->t = nz*muls->sliceThickness/nzPerSlice;
		sprintf(fileName,"potential_rz_%d.img",Znum);
		ptr = atPot[Znum];
		writeImage((void **)&ptr,header,fileName);   
#endif	  
	if (muls->printLevel > 1) printf("Created 3D (r-z) %d x %d potential array for Z=%d (%d, B=%g, dkx=%g, dky=%g. dkz=%g,sps=%d)\n",
		nx/2,nz/2,Znum,iKind,B,dkx,dky,dkz,izOffset);
}

Potential::Potential() :
m_splines(std::vector<AkimaSpline<float_tt,float_tt>>(NZMAX+1))
{
	m_splinr=std::vector<float_tt>(NRMAX);
	m_splines=std::vector<AkimaSpline<float_tt,float_tt>>(NZMAX+1);
	float_tt dlnr = static_cast<float_tt>(log(RMAX/RMIN)/(NRMAX-1));
	for( int i=0; i<NRMAX; i++)
	{
		m_splinr[i] = static_cast<float_tt>(RMIN * exp( i * dlnr ));
	}
}

void Potential::GenerateSplineEntry(int Z)
{
	std::vector<float_tt> potential_values(m_splinr.size());

	// First, map the splinr std::vector to an Eigen Array (for vectorized operations in v3DAtom)
	float_tt* r_ptr = &m_splinr[0];
	Eigen::Map<QSfArr> e_splinr(r_ptr, NRMAX);

	// Map the potential_values to a separate Eigen Array to store the output
	float_tt* v_ptr = &potential_values[0];
	Eigen::Map<QSfArr> e_pots(v_ptr, NRMAX);
	// Do the computation.  NOTE: this is storing the data in potential_values, because the
	//    eigen array maps directly to its memory.
	//e_pots = v3Datom(Z, e_splinr, m_tdsFlag, m_scatFlag);

	// feed the data into the spline fitter
	m_splines[Z] = AkimaSpline<float_tt, float_tt>(m_splinr, potential_values);
	// record the fact that we know about this atom now.
	m_knownZvalues.push_back(Z);
}

float_tt Potential::LookUp(int Z, float_tt r)
{
	// check if spline has already been fit for this Z:
	if(std::find(m_knownZvalues.begin(), m_knownZvalues.end(), Z)==m_knownZvalues.end())
	{
		// If not, do so.
		GenerateSplineEntry(Z);
	}
	// evaluate the spline at radius r
	return m_splines[Z].interpolate(r);
}


/********************************************************************************
* Create Lookup table for 3D potential due to neutral atoms
********************************************************************************/
#define PHI_SCALE 47.87658
QScMat getAtomPotential3D(int Znum, MULS *muls,double B,int *nzSub,int *Nr,int *Nz_lut) {
	int ix,iy,iz,iiz,ind3d,iKind,izOffset;
	double zScale,kzmax,zPos,xPos;
	fftwf_plan plan;
	static double f,phase,s2,s3,kmax2,smax2,kx,kz,dkx,dky,dkz; // ,dx2,dy2,dz2;
	static int nx,ny,nz,nzPerSlice;
	static QScMat atPot;
	static fftwf_complex *temp = NULL;
#if SHOW_SINGLE_POTENTIAL == 1
	static imageStruct *header = NULL;
	static fftwf_complex *ptr = NULL;
	static char fileName[256];
#endif 

	// initialize this atom, if it has not been done yet:
	if (atPot[Znum] == NULL) {
		
	}
	*Nz_lut = nz/2;
	*nzSub  = nzPerSlice;
	*Nr	  = nx/2;
	return atPot[Znum];
}

/********************************************************************************
* Lookup function for 3D potential offset due to charged atoms (ions)
********************************************************************************/
QScMat getAtomPotentialOffset3D(int Znum, MULS *muls,double B,int *nzSub,int *Nr,int *Nz_lut,float q) {
	int ix,iy,iz,iiz,ind3d,iKind,izOffset;
	double zScale,kzmax,zPos,xPos;
	fftwf_plan plan;
	static double f,phase,s2,s3,kmax2,kx,kz,dkx,dky,dkz; // ,dx2,dy2,dz2;
	static int nx,ny,nz,nzPerSlice;
	static QScMat atPot;
	static fftwf_complex *temp = NULL;
#if SHOW_SINGLE_POTENTIAL == 1
	static imageStruct *header = NULL;
	static fftwf_complex *ptr = NULL;
	static char fileName[256];
#endif 

	// if there is no charge to this atom, return NULL:
	if (q == 0) return NULL;

	// initialize this atom, if it has not been done yet:
	if (atPot[Znum] == NULL) {
#if USE_REZ_SFACTS
		iKind = Znum;
#else
		printf("Using charged atoms only works with scattering factors by Rez et al!\n",Znum);
		exit(0);
#endif

	}
}
// #undef SHOW_SINGLE_POTENTIAL
// end of fftwf_complex *getAtomPotential3D(...)

#define PHI_SCALE 47.87658
// #define SHOW_SINGLE_POTENTIAL 0
////////////////////////////////////////////////////////////////////////////
// This function should be used yet, because it computes the projected
// potential wrongly, since it doe not yet perform the projection!!!
fftwf_complex *getAtomPotential2D(int Znum, MULS *muls,double B) {
	int ix,iy,iz,ind,iKind;
	double min;
	fftwf_plan plan;
	static double f,phase,s2,s3,kmax2,kx,ky,dkx,dky;
	static int nx,ny;
	static QScMat atPot;
#if SHOW_SINGLE_POTENTIAL == 1
	static imageStruct *header = NULL;
	static char fileName[256];
#endif 

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
					// f = fe3D(Znum,k2,muls->tds,1.0,muls->scatFactor);
					// multiply scattering factor with Debye-Waller factor:
					// printf("k2=%g,B=%g, exp(-k2B)=%g\n",k2,B,exp(-k2*B));
					f = seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF,sqrt(s2))*exp(-s2*B*0.25);
					/*
					for (iz=1;iz<nx;iz++) {
						s3 = 0.5*dkx*iz;
						s3 = s3*s3+s2;
						if (s3<kmax2) {
							f += 2*seval(scatPar[0],scatPar[iKind],splinb,splinc,splind,N_SF,sqrt(s3))*exp(-s3*B*0.25);
					} 
						else break;
					} 
					
					f *= 0.5*dkx;
					*/
					// phase	= kx*xPos + kz*zPos;
					// phase = PI*(kx*muls->resolutionX*nx/(OVERSAMP_X)+ky*muls->resolutionY*ny/(OVERSAMP_X));
					phase = PI*(kx*muls->resolutionX*nx+ky*muls->resolutionY*ny);
					atPot[Znum][ind][0] = f*cos(phase);
					atPot[Znum][ind][1] = f*sin(phase);
				}
			}
		}
#if SHOW_SINGLE_POTENTIAL == 1
		if (header == NULL) 
			header = makeNewHeaderCompact(1,nx,ny,0,dkx,dky,0,NULL,"potential");
		sprintf(fileName,"pot_rec_%d.img",Znum);
		header->dx = dkx;
		header->dy = dky;
		writeImage((void **)&(atPot[Znum]),header,fileName);
#endif    
		plan = fftwf_plan_dft_2d(nx,ny,atPot[Znum],atPot[Znum],FFTW_BACKWARD,FFTW_ESTIMATE);
		fftwf_execute(plan);
		fftwf_destroy_plan(plan);
		for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
				atPot[Znum][iy+ix*ny][0] *= dkx*dky*(OVERSAMP_X*OVERSAMP_X);  
		}
		/*
		for (min = atPot[Znum][0][0],ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
			if (sqrt((double)((ix-nx/2)*(ix-nx/2)+(iy-ny/2)*(iy-ny/2))) > nx/2)
				atPot[Znum][iy+ix*ny][0] = 0;
			else {

				// divide by unit cell area (volume, if in 3D):
				atPot[Znum][iy+ix*ny][0] *= dkx*dky; //  / (OVERSAMP_X*OVERSAMP_X);  
				if (atPot[Znum][iy+ix*ny][0] < 0) atPot[Znum][iy+ix*ny][0] = 0;
				if (atPot[Znum][iy+ix*ny][0] < min) min = atPot[Znum][iy+ix*ny][0];
			}
			atPot[Znum][iy+ix*ny][1]= 0;
		}
		*/
		/*
		printf("Found minimum potential value of %g ... subtracting it from 2D potential.\n",min);
		for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) {
			if (sqrt((double)((ix-nx/2)*(ix-nx/2)+(iy-ny/2)*(iy-ny/2))) <= nx/2)
				atPot[Znum][iy+ix*ny][0] -= min;
		}
		*/
		// make sure we don't produce negative potential:
		// if (min < 0) for (ix=0;ix<nx;ix++) for (iy=0;iy<ny;iy++) atPot[Znum][iy+ix*ny][0] -= min;
#if SHOW_SINGLE_POTENTIAL == 1
		if (header == NULL) 
			header = makeNewHeaderCompact(1,nx,ny,0,muls->resolutionX/OVERSAMP_X,
			muls->resolutionY/OVERSAMP_X,0,NULL,"potential");
		header->dx = muls->resolutionX/OVERSAMP_X;
		header->dy = muls->resolutionY/OVERSAMP_X;
		sprintf(fileName,"potential_%d.img",Znum);
		writeImage((void **)&(atPot[Znum]),header,fileName);
#endif    
		printf("Created 2D %d x %d potential array for Z=%d (%d, B=%g A^2)\n",nx,ny,Znum,iKind,B);
	}
	return atPot[Znum];
}
#undef PHI_SCALE
#undef SHOW_SINGLE_POTENTIAL
