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
#include "potentials/pot_factory.hpp"
#include "wavefunctions.hpp"
#include "crystal.hpp"

class CExperimentBase : public IExperiment
{
public:
  CExperimentBase(const ConfigReaderPtr &configReader);
  virtual void DisplayProgress(int flag);
  virtual void DisplayParams();
  virtual void Run()=0;

  virtual void SaveImages()=0;

protected:
  virtual void PostSliceProcess(unsigned absoluteSlice){};  // Called in RunMuls after a slice is transmitted/propagated through.  Override as desired.
  
  virtual void CollectIntensity(unsigned absoluteSlice)=0;
  virtual int RunMuls(WavePtr wave);
  virtual void InterimWave(int slice);

  virtual void Transmit(WavePtr wave, unsigned sliceIdx);
  virtual void Propagate(WavePtr wave, float_tt dz);
  virtual void AddDPToAvgArray(const WavePtr &wave);

  void ReadAvgArray();
  void ReadAvgArray(unsigned navg);
  void ReadAvgArray(unsigned posX, unsigned posY);

  void _WriteAvgArray(std::string &fileName, std::string &comment, 
                      std::map<std::string, double> &params,
                      std::vector<unsigned> &position);

  inline void WriteAvgArray(std::string comment="Average Array", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    std::vector<unsigned> position;
    _WriteAvgArray(avgFilePrefix, comment, params, position);
  }
  inline void WriteAvgArray(unsigned navg, std::string comment="Average Array", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    std::vector<unsigned>position(1);
    position[0]=navg;
    _WriteAvgArray(avgFilePrefix, comment, params, position);
  }
  inline void WriteAvgArray(unsigned posX, unsigned posY, std::string comment="Average Array", 
                 std::map<std::string, double>params = std::map<std::string, double>())
  {
    std::vector<unsigned>position(2);
    position[0]=posX;
    position[1]=posY;
    _WriteAvgArray(avgFilePrefix, comment, params, position);
  }

  void fft_normalize(WavePtr wave);

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

  float_tt *m_avgArray;   // The averaged diffraction pattern (m_avgCount says how many are averaged currently)

  unsigned m_iPosX,m_iPosY;           /* integer offset for positioning probe within potential array */

  ImageIOPtr m_imageIO;

  std::vector<float_tt> m_propxr, m_propxi, m_propyr, m_propyi;
};

#endif










