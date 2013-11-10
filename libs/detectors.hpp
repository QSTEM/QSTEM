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

#ifndef DETECTORS_H
#define DETECTORS_H

#include "imagelib_fftw3.hpp"
#include "wavefunctions.hpp"

class Detector {
  ImageIOPtr m_imageIO;
public:
  int Navg;
  float_tt **m_image;        // place for storing avg image = sum(data)/Navg
  float_tt **m_image2;        // we will store sum(data.^2)/Navg 
  float_tt m_rInside,m_rOutside;
  float_tt m_k2Inside,m_k2Outside;
  char m_name[32];
  float_tt m_error;
  float_tt m_shiftX,m_shiftY;
  float_tt m_dx, m_dy;

public:
  Detector(int nx, int ny, float_tt resX, float_tt resY);
  boost::shared_ptr<Detector> Clone();
  void Initialize();
  void CollectIntensity(WavePtr &wave);
  void WriteImage(const char *fileName, const char *comment, std::map<std::string, double> &params,
                  std::vector<unsigned>position=std::vector<unsigned>());
};

typedef boost::shared_ptr<Detector> DetectorPtr;


// People will generally use this class - it will read config files to create
//   Detector objects as necessary.
class DetectorManager
{
  DetectorManager(ConfigReaderPtr &configReader);
  ~DetectorManager();
  void LoadDetectors(ConfigReaderPtr &configReader, std::vector<unsigned> &output_planes);
  // Use this for more generally setting output at set thicknesses (TODO: make sure intensity is 
  //    collected at these planes!!)
  void LoadDetectors(ConfigReaderPtr &configReader, std::vector<float_tt> &thicknesses);
  // Loads parameters for detectors from a config file
  void LoadDetectors(ConfigReaderPtr &configReader);
  // collects intensity from diffraction patterns to yield image intensity
  void CollectIntensity(int plane_idx, WavePtr &wave);
  // Saves detector images to files
  void SaveDetectors();
private:
  std::vector<std::vector<DetectorPtr> > m_detectors
};

typedef boost::shared_ptr<DetectorManager> DetectorMgrPtr;

#endif









