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

#ifndef EXPERIMENT_CBED_H
#define EXPERIMENT_CBED_H

#include "base.hpp"

class CExperimentCBED : public CExperimentBase
{
public:
  CExperimentCBED(const ConfigReaderPtr &configReader);
  void Run();
  void DisplayParams();
  virtual void WriteBeams(unsigned absoluteSlice);
protected:
  void CollectIntensity(unsigned absoluteSlice);
  void SaveImages();

  unsigned m_nbout;				/* number of recorded beams */
  float_tt **m_pendelloesung;
  bool m_lbeams;				/* flag indicating whether to record beams */	

  unsigned m_scanXStart, m_scanYStart;     /* The beam position on the sample */

  float_tt m_sourceRadius;                 /* The source radius in angstroms */

  bool m_showProbe;            /* if true, saves a plot of the probe */

  bool m_storeSeries;          

};

#endif
