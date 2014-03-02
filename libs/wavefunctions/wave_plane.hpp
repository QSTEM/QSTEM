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

#ifndef PLANE_WAVE_H
#define PLANE_WAVE_H

#include "wave_base.hpp"

namespace QSTEM
{

class QSTEM_HELPER_DLL_EXPORT CPlaneWave : public CBaseWave
{
public:
  CPlaneWave(const ConfigReaderPtr &configReader);
  CPlaneWave(const CPlaneWave& other);
  CPlaneWave();
  virtual void FormProbe();
  void TiltBeam(bool tiltBack=false);
  void TiltBack();
  virtual void DisplayParams();

  WavePtr Clone();

  // ReadImage is for TEM mode
  void ReadImage();
  void WriteImage();
protected:
  float_tt m_btiltx, m_btilty;     /* beam tilt, mrad */
  RealVector m_image;               /* Real-space image output */
private:
  friend class CWaveFactory;
  // Create an instance of this class, wrapped in a shared ptr
  //     This should not be inherited - any subclass needs its own implementation.
  static WavePtr Create(const ConfigReaderPtr &reader){
    return WavePtr(new CPlaneWave(reader));
  }  
};

typedef boost::shared_ptr<CPlaneWave> PlaneWavePtr;

}

#endif
