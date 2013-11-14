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

#include <map>
#include <string>
#include <vector>
#include <cmath>

#include "detectors.hpp"
#include "memory_fftw3.hpp"

Detector::Detector(int nx, int ny, float_tt resX, float_tt resY, float_tt wavelength)
  : m_error(0),
    m_shiftX(0),
    m_shiftY(0),
    m_Navg(0),
    m_nx(nx),
    m_ny(ny),
    m_dx(resX),
    m_dy(resY),
    m_wavelength(wavelength)
{
  m_image = float2D(nx,ny,"ADFimag");
  m_image2 = float2D(nx,ny,"ADFimag");
  
  // TODO: need way of passing file formats into constructor here
  m_imageIO=ImageIOPtr(new CImageIO(nx, ny));
}

// Creates a detector with the appropriately sized data arrays.  This is meant to then
// be populated with updated radius and shift values, and then initialized after that.
//  THE DETECTOR COMING OUT OF THIS IS NOT IMMEDIATELY VALID!  Make sure you are using
//    this in only your ReadDetectors method.
DetectorPtr Detector::Clone()
{
  return DetectorPtr(new Detector(m_nx, m_ny, m_dx, m_dy, m_wavelength));
}

void Detector::Initialize()
{
  /* determine v0 specific k^2 values corresponding to the angles */
  m_k2Inside = (float)(sin(m_rInside*0.001)/m_wavelength);
  m_k2Outside = (float)(sin(m_rOutside*0.001)/m_wavelength);
  /* calculate the squares of the ks */
  m_k2Inside *= m_k2Inside;
  m_k2Outside *= m_k2Outside;
}

void Detector::WriteImage(std::string &fileName, std::string &comment, std::map<std::string, double> &params,
                          std::vector<unsigned> &position)
{
  params["dx"]=m_dx;
  params["dy"]=m_dy;
  // Thickness is set externally and passed in on params
  m_imageIO->WriteRealImage((void **)m_image, fileName, params, comment, position);
}




DetectorManager::DetectorManager(ConfigReaderPtr &configReader)
{
  LoadDetectors(configReader);
}

void DetectorManager::LoadDetectors(ConfigReaderPtr &configReader, std::vector<float_tt> &thicknesses)
{
  int numDetectors;

  m_thicknesses=thicknesses;

  m_detectors.resize(m_thicknesses.size());
  configReader->ReadNumberOfDetectors(numDetectors);
  for (size_t plane_idx=0; plane_idx<m_thicknesses.size(); plane_idx++)
    {
      m_detectors[plane_idx].resize(numDetectors);
      for (size_t det_idx=0; det_idx<numDetectors; det_idx++)
        {
          DetectorPtr det = m_detectors[plane_idx][det_idx];
          configReader->ReadDetectorParameters(det_idx, det->m_rInside, det->m_rOutside, det->m_name, 
                                               det->m_shiftX, det->m_shiftY);
          det->Initialize();
        }
    }
}

void DetectorManager::LoadDetectors(ConfigReaderPtr &configReader, std::vector<unsigned> &output_planes)
{
  float_tt sliceThickness, dummyf;
  bool dummyb;
  unsigned dummyi;
  configReader->ReadSliceParameters(dummyb, sliceThickness, dummyi, dummyi, dummyf);
  // one extra plane for the final output
  std::vector<float_tt> thicknesses(output_planes.size(), float_tt());
  for (size_t plane_idx=0; plane_idx<thicknesses.size(); plane_idx++)
    {
      thicknesses[plane_idx]=output_planes[plane_idx]*sliceThickness;
    }
  return LoadDetectors(configReader, thicknesses);
}

void DetectorManager::LoadDetectors(ConfigReaderPtr &configReader)
{
  float_tt sliceThickness, dummyf;
  bool dummyb;
  unsigned nslices, outputInterval;
  int numDetectors;
  configReader->ReadSliceParameters(dummyb, sliceThickness, nslices, outputInterval, dummyf);
  // one extra plane for the final output
  unsigned nplanes = nslices/outputInterval;
  std::vector<float_tt> thicknesses(nplanes+1, float_tt());
  for (size_t plane_idx=0; plane_idx<nplanes; plane_idx++)
    {
      thicknesses[plane_idx]=(plane_idx+1)*sliceThickness;
    }
  thicknesses[nplanes]=nslices*sliceThickness;
  return LoadDetectors(configReader, thicknesses);
}

void DetectorManager::CollectIntensity(int plane_idx, WavePtr &wave)
{
  for (size_t det_idx=0; det_idx<m_detectors[plane_idx].size(); det_idx++)
    {
      m_detectors[plane_idx][det_idx]->CollectIntensity(wave);
    }
}

void DetectorManager::SaveDetectors(std::string &comment,
                                    std::map<std::string, double> &params)
{
  std::vector<unsigned>detector_id(2), final_id(1);
  for (size_t plane_idx=0; plane_idx<m_detectors.size(); plane_idx++)
    {
      for (size_t det_idx=0; det_idx<m_detectors[plane_idx].size(); det_idx++)
        {
          detector_id[0]=det_idx;
          detector_id[1]=plane_idx;
          final_id[0]=det_idx;
          DetectorPtr det = m_detectors[plane_idx][det_idx];
          params["Thickness"]=m_thicknesses[plane_idx];
          
          det->WriteImage(det->m_name, comment, params);
        }
    }
}




