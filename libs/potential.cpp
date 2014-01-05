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


#include "potential.hpp"

#include "potentials/pot_2d.hpp"
#include "potentials/pot_2d_fft.hpp"
#include "potentials/pot_3d.hpp"
#include "potentials/pot_3d_fft.hpp"


PotPtr GetPotential(ConfigReaderPtr &configReader)
{
  bool _3D, fft;
  configReader->ReadPotentialCalculationParameters(fft, _3D);

  if (_3D)
    {
      if (fft)
        {
          return PotPtr(new C3DFFTPotential(configReader));
        }
      else
        {
          return PotPtr(new C3DPotential(configReader));
        }
    }
  else
    {
      if (fft)
        {
          return PotPtr(new C2DFFTPotential(configReader));
        }
      else
        {
          return PotPtr(new C2DPotential(configReader));
        }
    }
}
