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

#include "stdio.h"
#include <string.h>
#include "data_containers.hpp"
#include "memory_fftw3.hpp"

MULS::MULS():
  cellDiv(1),
  btiltx(0),btilty(0),tiltBack(true),
  centerSlices(false), sliceThickness(0),
  lpartl(0),
  atomRadius(5.0),
  saveFlag(0),
  sigmaf(0),
  dfdelt(0),
  acmax(0),
  acmin(0),
  aobj(0),
  aAIS(0),
  periodicZ(false),periodicXY(false),
  bandlimittrans(true),
  readPotential(false),
  savePotential(false),
  saveTotalPotential(false),
  plotPotential(false),
  fftpotential(true),
  avgRuns(10),
  storeSeries(false),
  potOffsetY(0),potOffsetX(0),
  tomoTilt(0),
  tomoStart(0),
  tomoStep(0),
  tomoCount(0),
  cz(NULL),
  scatFactor(DOYLE_TURNER),
  normHolog(0),
  gaussianProp(0),
  pendelloesung(NULL),
  onlyFresnel(NULL),
  showPhaseplate(NULL)
{

}
