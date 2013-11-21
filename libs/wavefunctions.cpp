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

#include "wavefunctions.hpp"

void CreateWaveFunctionDataSets(unsigned x, unsigned y, std::vector<unsigned> positions, std::string output_ext)
{
  CImageIO imageIO(x, y, "", output_ext);
  std::string potDataSetLabel = "Potential";
  std::string mulswavDataSetLabel = "mulswav";
  imageIO.CreateComplexDataSet(potDataSetLabel, positions);
  imageIO.CreateComplexDataSet(mulswavDataSetLabel, positions);
}

WAVEFUNC::WAVEFUNC(unsigned x, unsigned y, float_tt resX, float_tt resY, 
                   std::string input_ext, std::string output_ext) :
  //m_position(std::vector<unsigned>()),
  detPosX(0),
  detPosY(0),
  iPosX(0),
  iPosY(0),
  thickness(0.0),
  nx(x),
  ny(y),
  resolutionX(resX),
  resolutionY(resY),
  m_params(std::map<std::string, double>())
{
  diffpat = float2D(nx,ny,"diffpat");
  avgArray = float2D(nx,ny,"avgArray");

  // TODO: need to pass file extension through to this constructor
  m_imageIO=ImageIOPtr(new CImageIO(nx, ny, input_ext, output_ext));
	
  wave = complex2D(nx, ny, "wave");
#if FLOAT_PRECISION == 1
  fftPlanWaveForw = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD, FFTW_ESTIMATE);
  fftPlanWaveInv = fftwf_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD, FFTW_ESTIMATE);
#else
  fftPlanWaveForw = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_FORWARD,
                                     fftMeasureFlag);
  fftPlanWaveInv = fftw_plan_dft_2d(nx,ny,wave[0],wave[0],FFTW_BACKWARD,
                                    fftMeasureFlag);
#endif
}

WAVEFUNC::WAVEFUNC(ConfigReaderPtr &configReader)
  : detPosX(0)
  , detPosY(0)
  , iPosX(0)
  , iPosY(0)
  , thickness(0)
{
  configReader->ReadProbeArraySize(nx, ny);
  configReader->ReadResolution(resolutionX, resolutionY);
  configReader->ReadDoseParameters(beamCurrent, dwellTime);
  configReader->ReadVoltage(v0);
  electronScale = beamCurrent*dwellTime*MILLISEC_PICOAMP;

  // TODO: need to figure out how user is going to specify input/output formats
  WAVEFUNC(nx, ny, resolutionX, resolutionY, ".img", ".img");
}

void WAVEFUNC::_WriteWave(std::string &fileName, std::string comment,
                         std::map<std::string, double>params)
{
  params["dx"]=resolutionX;
  params["dy"]=resolutionY;
  params["Thickness"]=thickness;
  m_imageIO->WriteComplexImage((void **)wave, fileName, params, comment, m_position);
}

void WAVEFUNC::_WriteDiffPat(std::string &fileName, std::string comment,
                            std::map<std::string, double> params)
{
  params["dx"]=1.0/(nx*resolutionX);
  params["dy"]=1.0/(ny*resolutionY);
  params["Thickness"]=thickness;
  m_imageIO->WriteRealImage((void **)diffpat, fileName, params, comment, m_position);
}

  void WAVEFUNC::_WriteAvgArray(std::string &fileName, std::string comment, 
                               std::map<std::string, double> params)
{
  params["dx"]=1.0/(nx*resolutionX);
  params["dy"]=1.0/(ny*resolutionY);
  params["Thickness"]=thickness;
  m_imageIO->WriteRealImage((void **)avgArray, fileName, params, comment, m_position);
}

void WAVEFUNC::SetWavePosition(unsigned navg)
{
  m_position.resize(1);
  m_position[0]=navg;
}

void WAVEFUNC::SetWavePosition(unsigned posX, unsigned posY)
{
  detPosX=posX;
  detPosY=posY;
  m_position.resize(2);
  m_position[0]=posX;
  m_position[1]=posY;
}

void WAVEFUNC::ReadWave()
{
  m_position.clear();
  m_imageIO->ReadImage((void **)wave, waveFilePrefix, m_position);
}

void WAVEFUNC::ReadWave(unsigned navg)
{
  SetWavePosition(navg);
  m_imageIO->ReadImage((void **)wave, waveFilePrefix, m_position);

}

void WAVEFUNC::ReadWave(unsigned positionx, unsigned positiony)
{
  SetWavePosition(positionx, positiony);
  m_imageIO->ReadImage((void **)wave, waveFilePrefix, m_position);
}

void WAVEFUNC::ReadDiffPat()
{
  m_position.clear();
  m_imageIO->ReadImage((void **)diffpat, dpFilePrefix, m_position);
}

void WAVEFUNC::ReadDiffPat(unsigned navg)
{
  SetWavePosition(navg);
  m_imageIO->ReadImage((void **)diffpat, dpFilePrefix, m_position);
}

void WAVEFUNC::ReadDiffPat(unsigned positionx, unsigned positiony)
{
  SetWavePosition(positionx, positiony);
  m_imageIO->ReadImage((void **)diffpat, dpFilePrefix, m_position);
}

void WAVEFUNC::ReadAvgArray()
{
  m_position.clear();
  m_imageIO->ReadImage((void **)avgArray, avgFilePrefix, m_position);
}

void WAVEFUNC::ReadAvgArray(unsigned navg)
{
  SetWavePosition(navg);
  m_imageIO->ReadImage((void **)avgArray, avgFilePrefix, m_position);
}

void WAVEFUNC::ReadAvgArray(unsigned positionx, unsigned positiony)
{
  SetWavePosition(positionx, positiony);
  m_imageIO->ReadImage((void **)avgArray, avgFilePrefix, m_position);
}

/******************************************************************
* propagate_slow() 
* Propagates a wave
*****************************************************************/
WAVEFUNC::Propagate()
{
	int ixa, iya;
	float_tt wr, wi, tr, ti,ax,by;
	float_tt scale,t,dz; 
	static float_tt dzs=0;
    static std::vector<float_tt> propxr(nx), propxi(nx), propyr(ny), propyi(ny);
    float_tt wavlen;

	ax = m_resolutionX*m_nx;
	by = m_resolutionY*m_ny;
	dz = (*muls).cz[0];

	if (dz != dzs) {
          dzs = dz;
          scale = dz*PI;
          wavlen = wavelength(wave->v0);

          for( ixa=0; ixa<nx; ixa++) {
            wave->m_kx[ixa] = (ixa>nx/2) ? (float_tt)(ixa-nx)/ax : 
              (float_tt)ixa/ax;
            wave->m_kx2[ixa] = wave->m_kx[ixa]*wave->m_kx[ixa];
            t = scale * (wave->m_kx2[ixa]*wavlen);
            propxr[ixa] = (float_tt)  cos(t);
            propxi[ixa] = (float_tt) -sin(t);
          }
          for( iya=0; iya<ny; iya++) {
            wave->m_ky[iya] = (iya>ny/2) ? 
              (float_tt)(iya-ny)/by : 
              (float_tt)iya/by;
            wave->m_ky2[iya] = wave->m_ky[iya]*wave->m_ky[iya];
            t = scale * (wave->m_ky2[iya]*wavlen);
            propyr[iya] = (float_tt)  cos(t);
            propyi[iya] = (float_tt) -sin(t);
          }
          wave->m_k2max = nx/(2.0F*ax);
          if (ny/(2.0F*by) < wave->m_k2max ) wave->m_k2max = ny/(2.0F*by);
          wave->m_k2max = 2.0/3.0 * wave->m_k2max;
          // TODO: modifying shared value from multiple threads?
          wave->m_k2max = wave->m_k2max*wave->m_k2max;
	} 
	/* end of: if dz != dzs */
	/*************************************************************/

	/*************************************************************
	* Propagation
	************************************************************/
	for( ixa=0; ixa<nx; ixa++) {
          if( wave->m_kx2[ixa] < wave->m_k2max ) {
            for( iya=0; iya<ny; iya++) {
              if( (wave->m_kx2[ixa] + wave->m_ky2[iya]) < wave->m_k2max ) {
                
                wr = wave->wave[ixa][iya][0];
                wi = wave->wave[ixa][iya][1];
                tr = wr*propyr[iya] - wi*propyi[iya];
                ti = wr*propyi[iya] + wi*propyr[iya];
                wave->wave[ixa][iya][0] = tr*propxr[ixa] - ti*propxi[ixa];
                wave->wave[ixa][iya][1] = tr*propxi[ixa] + ti*propxr[ixa];

              } else
                wave->wave[ixa][iya][0] = wave->wave[ixa][iya][1] = 0.0F;
            } /* end for(iy..) */

          } else for( iya=0; iya<ny; iya++)
                   wave->wave[ixa][iya][0] = wave->wave[ixa][iya][1] = 0.0F;
	} /* end for(ix..) */
} /* end propagate */

/*------------------------ transmit() ------------------------*/
/*
transmit the wavefunction thru one layer 
(simply multiply wave by transmission function)

waver,i[ix][iy]  = real and imaginary parts of wavefunction
transr,i[ix][iy] = real and imag parts of transmission functions

nx, ny = size of array

on entrance waver,i and transr,i are in real space

only waver,i will be changed by this routine
*/
void WAVEFUNCTION::Transmit(PotPtr pot, unsigned sliceIdx) {
	double wr, wi, tr, ti;

	complex_tt **w,**t;
	w = (complex_tt **)m_wave;
	t = (complex_tt **)pot->trans[sliceIdx];

	/*  trans += posx; */
	for(unsigned ix=0; ix<m_nx; ix++) for(unsigned iy=0; iy<m_ny; iy++) {
		wr = w[ix][iy][0];
		wi = w[ix][iy][1];
		tr = t[ix+m_iPosX][iy+m_iPosY][0];
		ti = t[ix+m_iPosX][iy+m_iPosY][1];
		w[ix][iy][0] = wr*tr - wi*ti;
		w[ix][iy][1] = wr*ti + wi*tr;
	} /* end for(iy.. ix .) */
} /* end transmit() */



/**********************************************
* This function creates a incident STEM probe 
* at position (dx,dy)
* with parameters given in muls
*
* The following Abberation functions are being used:
* 1) ddf = Cc*dE/E + Cc2*(dE/E)^2,    
*    Cc, Cc2 = chrom. Abber. (1st, 2nd order) [1]
* 2) chi(qx,qy) = (2*pi/lambda)*{0.5*C1*(qx^2+qy^2)+
*                 0.5*C12a*(qx^2-qy^2)+
*                 C12b*qx*qy+
*                 C21a/3*qx*(qx^2+qy^2)+
*                 ... 
*                 +0.5*C3*(qx^2+qy^2)^2
*                 +0.125*C5*(qx^2+qy^2)^3
*                 ... (need to finish)
*
*
*    qx = acos(kx/K), qy = acos(ky/K) 
*
* References:
* [1] J. Zach, M. Haider, 
*    "Correction of spherical and Chromatic Abberation 
*     in a low Voltage SEM", Optik 98 (3), 112-118 (1995)
* [2] O.L. Krivanek, N. Delby, A.R. Lupini,
*    "Towards sub-Angstroem Electron Beams", 
*    Ultramicroscopy 78, 1-11 (1999)
*
*********************************************/
#define SMOOTH_EDGE 5 // make a smooth edge on AIS aperture over +/-SMOOTH_EDGE pixels
void WAVEFUNCTION::FormProbe()
{
	// static char *plotFile = "probePlot.dat",systStr[32];
	int ix, iy, nx, ny, ixmid, iymid;
	int CsDefAstOnly = 0;
	float rmin, rmax, aimin, aimax;
	// float **pixr, **pixi;
	double  kx, ky, ky2,k2, ktheta2, ktheta, k2max, v0, wavlen,ax,by,x,y,
		rx2, ry2,rx,ry, pi, scale, pixel,alpha,
		df, df_eff, chi1, chi2,chi3, sum, chi, time,r,phi;
	double gaussScale = 0.05;
	double envelope,delta,avgRes,edge;

	// FILE *fp=NULL;

	/* temporary fix, necessary, because fftw has rec. space zero 
	in center of image:
	*/
	nx = (*muls).nx;
	ny = (*muls).ny;
	ax = nx*(*muls).resolutionX; 
	by = ny*(*muls).resolutionY; 
	dx = ax-dx;
	dy = by-dy;
	gaussScale = (*muls).gaussScale;
	// average resolution:
	avgRes = sqrt(0.5*(muls->resolutionX*muls->resolutionX+muls->resolutionY*muls->resolutionY));
	edge = SMOOTH_EDGE*avgRes;

	/********************************************************
	* formulas from:
	* http://cimesg1.epfl.ch/CIOL/asu94/ICT_8.html
	*
	* dE_E = dE/E = energy spread of emitted electrons
	* dV_V = dV/V = acc. voltage fluctuations
	* dI_I = dI/I = lens current fluctuations
	* delta defocus in Angstroem (Cc in A)
	*******************************************************/
	delta = muls->Cc*muls->dE_E;
	if (muls->printLevel > 2) printf("defocus offset: %g nm (Cc = %g)\n",delta,muls->Cc);

	if (wave->wave == NULL) {
		printf("Error in probe(): Wave not allocated!\n");
		exit(0);
	}

	/**********************************************************
	*  Calculate misc constants  
	*********************************************************/  
	time = cputim( );
	pi = 4.0 * atan( 1.0 );

	rx = 1.0/ax;
	rx2 = rx * rx;
	ry = 1.0/by;
	ry2 = ry * ry;

	ixmid = nx/2;
	iymid = ny/2;

	// df = muls->df0;
	v0 = muls->v0;
	wavlen = 12.26/ sqrt( v0*1.e3 + v0*v0*0.9788 );

	/*  printf("Wavelength: %g A\n",wavlen);
	*/


	// chi2 = (*muls).Cs*0.5*wavlen*wavlen;
	// chi3 = (*muls).C5*0.25*wavlen*wavlen*wavlen*wavlen;
	/* delta *= 0.5*delta*pi*pi*wavlen*wavlen; */

	/* convert convergence angle from mrad to rad */
	alpha = 0.001*muls->alpha;
	k2max = sin(alpha)/wavlen;  /* = K0*sin(alpha) */
	k2max = k2max * k2max;

	/*   Calculate MTF 
	NOTE zero freg is in the bottom left corner and
	expandes into all other corners - not in the center
	this is required for FFT

	PIXEL = diagonal width of pixel squared
	if a pixel is on the apertur boundary give it a weight
	of 1/2 otherwise 1 or 0
	*/
	pixel = ( rx2 + ry2 );
	scale = 1.0/sqrt((double)nx*(double)ny);

	/*
	if ((muls.a33 == 0) && (muls.a31 == 0) && (muls.a44 == 0) && (muls.a42 == 0) &&
	(muls.a55 == 0) && (muls.a53 == 0) && (muls.a51 == 0) && 
	(muls.a66 == 0) && (muls.a64 == 0) && (muls.a62 == 0) && (muls.C5 == 0)) {
	CsDefAstOnly = 1;
	}
	*/

	for( iy=0; iy<ny; iy++) {
		ky = (double) iy;
		if( iy > iymid ) ky = (double) (iy-ny);
		ky2 = ky*ky*ry2;
		for( ix=0; ix<nx; ix++) {
			kx = (double) ix;
			if( ix > ixmid ) kx = (double) (ix-nx);
			k2 = kx*kx*rx2 + ky2;
			ktheta2 = k2*(wavlen*wavlen);
			ktheta = sqrt(ktheta2);
			phi = atan2(ry*ky,rx*kx);
			// compute the effective defocus from the actual defocus and the astigmatism: 
			// df_eff = df + muls->astigMag*cos(muls->astigAngle+phi);

			// chi = chi1*k2*(df_eff +chi2*k2)-2.0*pi*( (dx*kx/ax) + (dy*ky/by) );
			// defocus, astigmatism, and shift:
			chi = ktheta2*(muls->df0+delta + muls->astigMag*cos(2.0*(phi-muls->astigAngle)))/2.0;
			ktheta2 *= ktheta;  // ktheta^3 
			if ((muls->a33 > 0) || (muls->a31 > 0)) {
				chi += ktheta2*(muls->a33*cos(3.0*(phi-muls->phi33))+muls->a31*cos(phi-muls->phi31))/3.0;
			}	
			ktheta2 *= ktheta;   // ktheta^4
			if ((muls->a44 > 0) || (muls->a42 > 0) || (muls->Cs != 0)) {
				// chi += ktheta2*(muls->a33*cos(3*(phi-muls->phi33))+muls->a31*cos(phi-muls->phi31))/3.0;
				chi += ktheta2*(muls->a44*cos(4.0*(phi-muls->phi44))+muls->a42*cos(2.0*(phi-muls->phi42))+muls->Cs)/4.0;  
				//                     1/4*(a(4,4).*cos(4*(kphi-phi(4,4)))+a(4,2).*cos(2*(kphi-phi(4,2)))+c(4)).*ktheta.^4+...
			}
			ktheta2 *= ktheta;    // ktheta^5
			if ((muls->a55 > 0) || (muls->a53 > 0) || (muls->a51 > 0)) {
				chi += ktheta2*(muls->a55*cos(5.0*(phi-muls->phi55))+muls->a53*cos(3.0*(phi-muls->phi53))+muls->a51*cos(phi-muls->phi51))/5.0;
				//                     1/5*(a(5,5).*cos(5*(kphi-phi(5,5)))+a(5,3).*cos(3*(kphi-phi(5,3)))+a(5,1).*cos(1*(kphi-phi(5,1)))).*ktheta.^5+...
			}
			ktheta2 *= ktheta;    // ktheta^6
			if ((muls->a66 > 0) || (muls->a64 > 0) || (muls->a62 = 0) || (muls->C5 != 0)) {
				chi += ktheta2*(muls->a66*cos(6.0*(phi-muls->phi66))+muls->a64*cos(4.0*(phi-muls->phi64))+muls->a62*cos(2.0*(phi-muls->phi62))+muls->C5)/6.0;
				//                     1/6*(a(6,6).*cos(6*(kphi-phi(6,6)))+a(6,4).*cos(4*(kphi-phi(6,4)))+a(6,2).*cos(2*(kphi-phi(6,2)))+c(6)).*ktheta.^6);
			}

			chi *= 2*pi/wavlen;
			chi -= 2.0*pi*( (dx*kx/ax) + (dy*ky/by) );
			// include higher order aberrations


			if ( ( (*muls).ismoth != 0) && 
				( fabs(k2-k2max) <= pixel)) {
					wave->wave[ix][iy][0]= (float) ( 0.5*scale * cos(chi));
					wave->wave[ix][iy][1]= (float) (-0.5*scale* sin(chi));
			} 
			else if ( k2 <= k2max ) {
				wave->wave[ix][iy][0]= (float)  scale * cos(chi);
				wave->wave[ix][iy][1]= (float) -scale * sin(chi);
			} 
			else {
				wave->wave[ix][iy][0] = wave->wave[ix][iy][1] = 0.0f;
			}
		}
	}
	/* Fourier transform into real space */
	// fftwnd_one(muls->fftPlanInv, &(muls->wave[0][0]), NULL);
#if FLOAT_PRECISION == 1
	fftwf_execute(wave->fftPlanWaveInv);
#else
	fftw_execute(wave->fftPlanWaveInv);
#endif
	/**********************************************************
	* display cross section of probe intensity
	*/

	/* multiply with gaussian in Real Space in order to avoid artifacts */
	if (muls->gaussFlag) {
		for( ix=0; ix<nx; ix++) {
			for( iy=0; iy<ny; iy++) {
				r = exp(-((ix-nx/2)*(ix-nx/2)+(iy-ny/2)*(iy-ny/2))/(nx*nx*gaussScale));
				wave->wave[ix][iy][0] *= (float)r;
				wave->wave[ix][iy][1] *= (float)r;
			}
		}  
	}

	/* Apply AIS aperture in Real Space */
	// printf("center: %g,%g\n",dx,dy);
	if (muls->aAIS > 0) {
		for( ix=0; ix<nx; ix++) {
			for( iy=0; iy<ny; iy++) {
				x = ix*muls->resolutionX-dx;
				y = iy*muls->resolutionY-dy;
				r = sqrt(x*x+y*y);
				delta = r-0.5*muls->aAIS+edge;
				if (delta > 0) {
					wave->wave[ix][iy][0] = 0;
					wave->wave[ix][iy][1] = 0;
				}
				else if (delta >= -edge) {
					scale = 0.5*(1-cos(pi*delta/edge));
					wave->wave[ix][iy][0] = scale*wave->wave[ix][iy][0];
					wave->wave[ix][iy][1] = scale*wave->wave[ix][iy][1];
				}
			}
		}
	}

	/*  Normalize probe intensity to unity  */

	sum = 0.0;
	for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) 
		sum +=  wave->wave[ix][iy][0]*wave->wave[ix][iy][0]
	+ wave->wave[ix][iy][1]*wave->wave[ix][iy][1];

	scale = 1.0 / sum;
	scale = scale * ((double)nx) * ((double)ny);
	scale = (double) sqrt( scale );

	for( ix=0; ix<nx; ix++) 
		for( iy=0; iy<ny; iy++) {
			wave->wave[ix][iy][0] *= (float) scale;
			wave->wave[ix][iy][1] *= (float) scale;
		}

		/*  Output results and find min and max to echo
		remember that complex pix are stored in the file in FORTRAN
		order for compatability
		*/

		rmin = wave->wave[0][0][0];
		rmax = rmin;
		aimin = wave->wave[0][0][1];
		aimax = aimin;
		for( iy=0; iy<ny; iy++) {
			for( ix=0; ix<nx; ix++) {
				if( wave->wave[ix][iy][0] < rmin ) rmin = wave->wave[ix][iy][0];
				if( wave->wave[ix][iy][0] > rmax ) rmax = wave->wave[ix][iy][0];
				if( wave->wave[ix][iy][1] < aimin ) aimin = wave->wave[ix][iy][1];
				if( wave->wave[ix][iy][1] > aimax ) aimax = wave->wave[ix][iy][1];
			}
		}
		(*muls).rmin = rmin;
		(*muls).rmax = rmax;
		(*muls).aimin = aimin;
		(*muls).aimax = aimax;

		/**********************************************************/

}  /* end probe() */