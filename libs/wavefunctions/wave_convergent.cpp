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

#include "wave_convergent.hpp"

CConvergentWave::CConvergentWave(const ConfigReaderPtr &configReader) : CBaseWave(configReader)
{
  // TODO: where does beam current belong?
  //configReader->ReadDoseParameters(m_beamCurrent, m_dwellTime);
}

/** Copy constructor - used to copy wave just before dispatching multiple threads for STEM simulations */
CConvergentWave::CConvergentWave(const WavePtr& other) : CBaseWave(other)
{
  // TODO: need to copy arrays and anything pointed to - anything that needs to be thread-local
}

CConvergentWave::CConvergentWave() : CBaseWave()
{
}

void CConvergentWave::DisplayParams()
{
  CBaseWave::DisplayParams();

  printf("* Aperture half angle:  %g mrad\n",m_alpha);
  printf("* AIS aperture:         ");
  if (m_aAIS > 0) printf("%g A\n",m_aAIS);
  else printf("none\n");

  printf("* C_3 (C_s):            %g mm\n",m_Cs*1e-7);
  printf("* C_1 (Defocus):        %g nm%s\n",0.1*m_df0,
         (m_Scherzer == 1) ? " (Scherzer)" : (m_Scherzer==2) ? " (opt.)":"");
  printf("* Astigmatism:          %g nm, %g deg\n",0.1*m_astigMag,RAD2DEG*m_astigAngle);

  // more aberrations:
  if (m_a33 > 0)
    printf("* a_3,3:                %g nm, phi=%g deg\n",m_a33*1e-1,m_phi33*RAD2DEG);
  if (m_a31 > 0)
    printf("* a_3,1:                %g nm, phi=%g deg\n",m_a31*1e-1,m_phi31*RAD2DEG);
  
  if (m_a44 > 0)
    printf("* a_4,4:                %g um, phi=%g deg\n",m_a44*1e-4,m_phi44*RAD2DEG);
  if (m_a42 > 0)
    printf("* a_4,2:                %g um, phi=%g deg\n",m_a42*1e-4,m_phi42*RAD2DEG);

  if (m_a55 > 0)
    printf("* a_5,5:                %g um, phi=%g deg\n",m_a55*1e-4,m_phi55*RAD2DEG);
  if (m_a53 > 0)
    printf("* a_5,3:                %g um, phi=%g deg\n",m_a53*1e-4,m_phi53*RAD2DEG);
  if (m_a51 > 0)
    printf("* a_5,1:                %g um, phi=%g deg\n",m_a51*1e-4,m_phi51*RAD2DEG);

  if (m_a66 > 0)
    printf("* a_6,6:                %g um, phi=%g deg\n",m_a66*1e-7,m_phi66*RAD2DEG);
  if (m_a64 > 0)
    printf("* a_6,4:                %g um, phi=%g deg\n",m_a64*1e-7,m_phi64*RAD2DEG);
  if (m_a62 > 0)
    printf("* a_6,2:                %g um, phi=%g deg\n",m_a62*1e-7,m_phi62*RAD2DEG);
  if (m_C5 != 0)
    printf("* C_5:                  %g mm\n",m_C5*1e-7);
  
  printf("* C_c:                  %g mm\n",m_Cc*1e-7);

  printf("* Damping dE/E: %g / %g \n",sqrt(m_dE_E*m_dE_E+m_dV_V*m_dV_V+m_dI_I*m_dI_I)*m_v0*1e3,m_v0*1e3);

  /*
    // TODO: where does beam current belong?
  printf("* beam current:         %g pA\n",m_beamCurrent);
  printf("* dwell time:           %g msec (%g electrons)\n",
         m_dwellTime,m_electronScale);
  */
}

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
void CConvergentWave::FormProbe()
{
  // static char *plotFile = "probePlot.dat",systStr[32];
  int ix, iy, nx, ny, ixmid, iymid;
  int CsDefAstOnly = 0;
  float rmin, rmax, aimin, aimax;
  // float **pixr, **pixi;
  double dx, dy;
  double  kx, ky, ky2,k2, ktheta2, ktheta, k2max, v0, wavlen,ax,by,x,y,
    rx2, ry2,rx,ry, pi, scale, pixel,alpha,
    df, df_eff, chi1, chi2,chi3, sum, chi, time,r,phi;
  double envelope,delta,avgRes,edge;

  // FILE *fp=NULL;

  /* temporary fix, necessary, because fftw has rec. space zero 
     in center of image:
  */
  ax = m_nx*m_dx; 
  by = m_ny*m_dy; 
  dx = ax-dx;
  dy = by-dy;
  // average resolution:
  avgRes = sqrt(0.5*(m_dx*m_dx+m_dy*m_dy));
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
  delta = m_Cc*m_dE_E;
  if (m_printLevel > 2) printf("defocus offset: %g nm (Cc = %g)\n",delta,m_Cc);
  
  if (m_wave.size()==0) {
    printf("Error in probe(): Wave not allocated!\n");
    exit(0);
  }

  /**********************************************************
   *  Calculate misc constants  
   *********************************************************/  
  pi = 4.0 * atan( 1.0 );

  rx = 1.0/ax;
  rx2 = rx * rx;
  ry = 1.0/by;
  ry2 = ry * ry;

  ixmid = nx/2;
  iymid = ny/2;

  // df = m_df0;
  //v0 = m_v0;
  
  /*  printf("Wavelength: %g A\n",wavlen);
   */


  // chi2 = (*muls).Cs*0.5*wavlen*wavlen;
  // chi3 = (*muls).C5*0.25*wavlen*wavlen*wavlen*wavlen;
  /* delta *= 0.5*delta*pi*pi*wavlen*wavlen; */

  /* convert convergence angle from mrad to rad */
  alpha = 0.001*m_alpha;
  k2max = sin(alpha)/m_wavlen;  /* = K0*sin(alpha) */
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
    if ((m_a33 == 0) && (m_a31 == 0) && (m_a44 == 0) && (m_a42 == 0) &&
    (m_a55 == 0) && (m_a53 == 0) && (m_a51 == 0) && 
    (m_a66 == 0) && (m_a64 == 0) && (m_a62 == 0) && (m_C5 == 0)) {
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
      // df_eff = df + m_astigMag*cos(m_astigAngle+phi);

      // chi = chi1*k2*(df_eff +chi2*k2)-2.0*pi*( (dx*kx/ax) + (dy*ky/by) );
      // defocus, astigmatism, and shift:
      chi = ktheta2*(m_df0+delta + m_astigMag*cos(2.0*(phi-m_astigAngle)))/2.0;
      ktheta2 *= ktheta;  // ktheta^3 
      if ((m_a33 > 0) || (m_a31 > 0)) {
        chi += ktheta2*(m_a33*cos(3.0*(phi-m_phi33))+m_a31*cos(phi-m_phi31))/3.0;
      }	
      ktheta2 *= ktheta;   // ktheta^4
      if ((m_a44 > 0) || (m_a42 > 0) || (m_Cs != 0)) {
        // chi += ktheta2*(m_a33*cos(3*(phi-m_phi33))+m_a31*cos(phi-m_phi31))/3.0;
        chi += ktheta2*(m_a44*cos(4.0*(phi-m_phi44))+m_a42*cos(2.0*(phi-m_phi42))+m_Cs)/4.0;  
        //                     1/4*(a(4,4).*cos(4*(kphi-phi(4,4)))+a(4,2).*cos(2*(kphi-phi(4,2)))+c(4)).*ktheta.^4+...
      }
      ktheta2 *= ktheta;    // ktheta^5
      if ((m_a55 > 0) || (m_a53 > 0) || (m_a51 > 0)) {
        chi += ktheta2*(m_a55*cos(5.0*(phi-m_phi55))+m_a53*cos(3.0*(phi-m_phi53))+m_a51*cos(phi-m_phi51))/5.0;
        //                     1/5*(a(5,5).*cos(5*(kphi-phi(5,5)))+a(5,3).*cos(3*(kphi-phi(5,3)))+a(5,1).*cos(1*(kphi-phi(5,1)))).*ktheta.^5+...
      }
      ktheta2 *= ktheta;    // ktheta^6
      if ((m_a66 > 0) || (m_a64 > 0) || (m_a62 = 0) || (m_C5 != 0)) {
        chi += ktheta2*(m_a66*cos(6.0*(phi-m_phi66))+m_a64*cos(4.0*(phi-m_phi64))+m_a62*cos(2.0*(phi-m_phi62))+m_C5)/6.0;
        //                     1/6*(a(6,6).*cos(6*(kphi-phi(6,6)))+a(6,4).*cos(4*(kphi-phi(6,4)))+a(6,2).*cos(2*(kphi-phi(6,2)))+c(6)).*ktheta.^6);
      }

      chi *= 2*pi/wavlen;
      chi -= 2.0*pi*( (dx*kx/ax) + (dy*ky/by) );
      // include higher order aberrations


      if ( ( m_ismoth != 0) && 
           ( fabs(k2-k2max) <= pixel)) {
        m_wave[ix+m_nx*iy][0]= (float) ( 0.5*scale * cos(chi));
        m_wave[ix+m_nx*iy][1]= (float) (-0.5*scale* sin(chi));
      } 
      else if ( k2 <= k2max ) {
        m_wave[ix+m_nx*iy][0]= (float)  scale * cos(chi);
        m_wave[ix+m_nx*iy][1]= (float) -scale * sin(chi);
      } 
      else {
        m_wave[ix+m_nx*iy][0] = m_wave[ix+m_nx*iy][1] = 0.0f;
      }
    }
  }
  /* Fourier transform into real space */
  ToRealSpace();
  /**********************************************************
   * display cross section of probe intensity
   */
  
  /* multiply with gaussian in Real Space in order to avoid artifacts */
  if (m_gaussFlag) {
    for( ix=0; ix<nx; ix++) {
      for( iy=0; iy<ny; iy++) {
        r = exp(-((ix-m_nx/2)*(ix-m_nx/2)+(iy-m_ny/2)*(iy-m_ny/2))/(m_nx*m_nx*m_gaussScale));
        m_wave[ix+m_nx*iy][0] *= (float)r;
        m_wave[ix+m_nx*iy][1] *= (float)r;
      }
    }  
  }

  /* Apply AIS aperture in Real Space */
  // printf("center: %g,%g\n",dx,dy);
  if (m_aAIS > 0) {
    for( ix=0; ix<m_nx; ix++) {
      for( iy=0; iy<m_ny; iy++) {
        x = ix*m_dx-dx;
        y = iy*m_dy-dy;
        r = sqrt(x*x+y*y);
        delta = r-0.5*m_aAIS+edge;
        if (delta > 0) {
          m_wave[ix+m_nx*iy][0] = 0;
          m_wave[ix+m_nx*iy][1] = 0;
        }
        else if (delta >= -edge) {
          scale = 0.5*(1-cos(pi*delta/edge));
          m_wave[ix+m_nx*iy][0] = scale*m_wave[ix+m_nx*iy][0];
          m_wave[ix+m_nx*iy][1] = scale*m_wave[ix+m_nx*iy][1];
        }
      }
    }
  }

  /*  Normalize probe intensity to unity  */
  
  sum = 0.0;
  for( ix=0; ix<m_nx; ix++) 
    for( iy=0; iy<m_ny; iy++) 
      sum +=  m_wave[ix+m_nx*iy][0]*m_wave[ix+m_nx*iy][0]
        + m_wave[ix+m_nx*iy][1]*m_wave[ix+m_nx*iy][1];

  scale = 1.0 / sum;
  scale = scale * ((double)m_nx) * ((double)m_ny);
  scale = (double) sqrt( scale );

  for( ix=0; ix<m_nx; ix++) 
    for( iy=0; iy<m_ny; iy++) {
      m_wave[ix+m_nx*iy][0] *= (float) scale;
      m_wave[ix+m_nx*iy][1] *= (float) scale;
    }

  /*  Output results and find min and max to echo
      remember that complex pix are stored in the file in FORTRAN
      order for compatability
  */

  rmin = m_wave[0][0];
  rmax = rmin;
  aimin = m_wave[0][1];
  aimax = aimin;
  for( iy=0; iy<m_ny; iy++) {
    for( ix=0; ix<m_nx; ix++) {
      if( m_wave[ix+m_nx*iy][0] < rmin ) rmin = m_wave[ix+m_nx*iy][0];
      if( m_wave[ix+m_nx*iy][0] > rmax ) rmax = m_wave[ix+m_nx*iy][0];
      if( m_wave[ix+m_nx*iy][1] < aimin ) aimin = m_wave[ix+m_nx*iy][1];
      if( m_wave[ix+m_nx*iy][1] > aimax ) aimax = m_wave[ix+m_nx*iy][1];
    }
  }
  m_rmin = rmin;
  m_rmax = rmax;
  m_aimin = aimin;
  m_aimax = aimax;

  /**********************************************************/

}  /* end probe() */
