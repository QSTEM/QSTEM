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

#ifndef EXPERIMENT_STEM_H
#define EXPERIMENT_STEM_H

#include "base.hpp"
#include "detectors.hpp"

class CExperimentSTEM : public CExperimentBase
{
public:
  CExperimentSTEM(const ConfigReaderPtr &configReader);
  void Run();
  void DisplayParams();
private:
  void SaveImages();
  virtual void CollectIntensity(unsigned absoluteSlice);
  unsigned m_stepsX, m_stepsY;
  DetectorMgrPtr m_detectors;
  float_tt m_scanXStart, m_scanXStop, m_scanYStart, m_scanYStop;
  // TODO: are these synonymous with m_stepsX?
  unsigned m_scanXN, m_scanYN;
  float_tt **m_avgArray;
};

#endif










