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

#ifndef EXPERIMENT_BASE_H
#define EXPERIMENT_BASE_H

#include "experiment_interface.hpp"
#include "potential.hpp"
#include "wavefunctions.hpp"
#include "crystal.hpp"

class CExperimentBase : public IExperiment
{
public:
  CExperimentBase(const ConfigReaderPtr &configReader);
  virtual void DisplayProgress(int flag);
  virtual void DisplayParams();
  virtual void Run()=0;

protected:
  virtual void CollectIntensity(unsigned absoluteSlice)=0;
  virtual int RunMuls();
  virtual void InterimWave(int slice);
  virtual void SaveImages()=0;

  virtual void Transmit(unsigned sliceIdx);
  virtual void Propagate(float_tt dz);

  bool m_tds;
  unsigned m_avgRuns, m_avgCount;  // number of runs to average; runs currently averaged
  unsigned m_printLevel;

  boost::filesystem::path m_outputLocation;

  StructurePtr m_crystal;  // The structure of the sample (atom positions)
  WavePtr m_wave;		   // The electron wave (this may be copied for multiprocessing)
  PotPtr m_potential;      // The sample potential

  float_tt m_intIntensity;  // Integrated intensity from experiment - if too low, 
							// your wave array is too small, and the beam is being scattered beyond it.

  unsigned m_cellDiv;		// The number of sub-slabs that the supercell is divided into
  bool m_equalDivs;			// Whether or not all sub-slabs are the same size
  unsigned m_outputInterval;  // The number of slices between saving intermediate output files
  unsigned m_totalSliceCount; // The total number of slices that we've progressed through (all sub-slabs included)

  float_tt m_thickness;       // The total thickness of the sample at the current slice

  std::vector<float_tt> m_chisq;
  std::string m_mode;      // String representing the multislice mode (e.g. TEM, STEM, etc.)
};

#endif










