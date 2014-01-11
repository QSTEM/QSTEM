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

#include "wave_plane.hpp"

CPlaneWave::CPlaneWave(const ConfigReaderPtr &configReader) : WAVEFUNC(configReader)
{
}

void CPlaneWave::FormProbe()
{
  if ((m_btiltx == 0) && (m_btilty == 0)) {
      for (unsigned ix=0;ix<m_nx;ix++) for (unsigned iy=0;iy<m_ny;iy++) {
          m_wave[ix][iy][0] = 1; 
          m_wave[ix][iy][1] = 0;
        }
    }
    else {
      TiltBeam();
    }
}

void CPlaneWave::TiltBeam(bool tiltBack)
{
  if ((m_btiltx != 0) || (m_btilty != 0)) 
    {

      int direction = tiltBack ? -1 : 1;

      // produce a tilted wave function (btiltx,btilty):
      float_tt ktx = direction*2.0*M_PI*sin(m_btiltx)/GetWavelength();
      float_tt kty = direction*2.0*M_PI*sin(m_btilty)/GetWavelength();
      for (unsigned ix=0;ix<m_nx;ix++) {
        float_tt x = m_dx*(ix-m_nx/2);
        for (unsigned iy=0;iy<m_ny;iy++) {
          float_tt y = m_dy*(iy-m_ny/2);
          m_wave[ix][iy][0] = (float)cos(ktx*x+kty*y);	
          m_wave[ix][iy][1] = (float)sin(ktx*x+kty*y);
        }
      }
    }
}

inline void CPlaneWave::TiltBack()
{
  TiltBeam(true);
}
