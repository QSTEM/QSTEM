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

#include "detectors.hpp"
#include "memory_fftw3.hpp"

Detector::Detector(int nx, int ny, float_tt resX, float_tt resY)
  : m_error(0),
    m_shiftX(0),
    m_shiftY(0),
    m_Navg(0),
    m_nx(nx),
    m_ny(ny),
    m_dx(resX),
    m_dy(resY)
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
  return DetectorPtr(new Detector(m_nx, m_ny, m_dx, m_dy));
}

void Detector::Initialize()
{
  /* determine v0 specific k^2 values corresponding to the angles */
  k2Inside = (float)(sin(rInside*0.001)/(wavelength(v0)));
  k2Outside = (float)(sin(rOutside*0.001)/(wavelength(v0)));
  /* calculate the squares of the ks */
  k2Inside *= k2Inside;
  k2Outside *= k2Outside;
}

void Detector::WriteImage(const char *fileName, const char *comment, std::map<std::string, double> &params,
                          std::vector<unsigned> position)
{
  params["dx"]=dx;
  params["dy"]=dy;
  // Thickness is set externally and passed in on params
  m_imageIO->WriteRealImage((void **)image, fileName, params, std::string(comment), position);
}






DetectorManager::DetectorManager(ConfigReaderPtr &configReader)
{
  LoadDetectors(configReader);
}

DetectorManager::LoadDetectors(ConfigReaderPtr &configReader, std::vector<float_tt> &thicknesses)
{
  int numDetectors;

  m_detectors.resize(thicknesses.size());
  configReader->ReadNumberOfDetectors(numDetectors);
  for (size_t plane_idx=0; plane_idx<nplanes; plane_idx++)
    {
      m_detectors[plane_idx].resize(numDetectors);
      for (size_t det_idx=0; det_idx<numDetectors; det_idx++)
        {
          DetectorPtr det = m_detectors[plane_idx][det_idx];
          configReader->ReadDetectorParameters(det_idx, det->m_rInside, det->m_rOutside, det->m_name, 
                                               det->m_shiftX, det->m_shiftY);
          det->Initialize();
          det->SetThickness(thicknesses[det_idx]);
        }
    }
}

  DetectorManager::LoadDetectors(ConfigReaderPtr &configReader, std::vector<unsigned> &output_planes)
{
  float_tt sliceThickness, dummyf;
  bool dummyb;
  unsigned dummyi;
  configReader->ReadSliceParameters(dummyb, sliceThickness, dummyi, dummyi, dummyf);
  // one extra plane for the final output
  std::vector<float_tt> thicknesses(output_planes.size(), float_tt());
  for (size_t plane_idx=0; plane_idx<nplanes; plane_idx++)
    {
      thicknesses[plane_idx]=(plane_idx+1)*sliceThickness;
    }
  thicknesses[nplanes]=nslices*sliceThickness;
  return LoadDetectors(configReader, thicknesses);
}

DetectorManager::LoadDetectors(ConfigReaderPtr &configReader)
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

DetectorManager::CollectIntensity(int plane_idx, WavePtr &wave)
{
  for (size_t det_idx=0; det_idx<m_detectors[plane_idx].size(); det_idx++)
    {
      m_detectors[plane_idx][det_idx]->CollectIntensity(wave);
    }
}

DetectorManager::SaveDetectors()
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
          
          det->WriteImage(det->m_name, "", ,

}
















