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
    m_nx(nx),
    m_ny(ny),
    m_dx(resX),
    m_dy(resY),
    m_Navg(1),
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

void Detector::WriteImage(std::map<std::string, double> &params,
                          std::vector<unsigned> &position)
{
  params["dx"]=m_dx;
  params["dy"]=m_dy;
  std::string comment;
  // Thickness is set externally and passed in on params
  m_imageIO->WriteRealImage((void **)m_image, m_name, params, comment, position);
}

void Detector::CollectIntensity(const WavePtr &wave)
{
  /********************************************************************
   * collectIntensity(muls, wave, slice)
   * collect the STEM signal on the annular detector(s) defined in muls
   * and write the appropriate pixel in the image for each detector and thickness
   * The number of images is determined by the following formula:
   * muls->slices*muls->cellDiv/muls->outputInterval 
   * There are muls->detectorNum different detectors
   *******************************************************************/
  int ixs,iys,t;
  float_tt k2;
  float_tt intensity,scale,scaleCBED,scaleDiff,intensity_save;
  char fileName[256],avgName[256]; 
  float_tt **diffpatAvg = NULL;
  int tCount = 0;
  unsigned nx, ny, offsetX, offsetY;

  std::vector<std::vector<DetectorPtr> > detectors;

  wave->GetSizePixels(nx, ny);
  wave->GetPositionOffset(offsetX, offsetY);

  //scale=1;
  // TODO: come up with where electronScale belongs (might be STEM-specific?
  scale = m_electronScale/((double)(nx*ny)*(nx*ny));
  // scaleCBED = 1.0/(scale*sqrt((double)(muls->nx*muls->ny)));
  scaleDiff = 1.0/sqrt((double)(nx*ny));

  int position_offset = offsetY * m_nx + offsetX;

  // Multiply each image by its number of averages and divide by it later again:
  m_image[offsetX][offsetY]  *= m_Navg;	
  m_image2[offsetX][offsetY] *= m_Navg;	
  m_error = 0;

  /* add the intensities in the already 
     fourier transformed wave function */
  unsigned px = nx*ny;
  for (unsigned i = 0; i < px; i++) 
    {
        {
          unsigned ix = i%nx;
          unsigned iy = i/nx;
          k2 = wave->GetK2(ix, iy);
          intensity = wave->GetPixelIntensity(i);
          wave->SetDiffPatPixel((ix+nx/2)%nx,(iy+ny/2)%ny, intensity*scaleDiff);
          intensity *= scale;
          if ((k2 >= m_k2Inside) && (k2 <= m_k2Outside)) 
            {
              // detector in center of diffraction pattern:
              if ((m_shiftX == 0) && (m_shiftY == 0)) 
                {
                  m_image[offsetX][offsetY] += intensity;
                  // misuse the error number for collecting this pixels raw intensity
                  m_error += intensity;
                }
              /* special case for shifted detectors: */		
              else 
                {
                  intensity_save = intensity;
                  ixs = (ix+(int)m_shiftX+nx) % nx;
                  iys = (iy+(int)m_shiftY+ny) % ny;	    
                  intensity = scale * wave->GetPixelIntensity(ixs, iys);
                  m_image[offsetX][offsetY] += intensity;
                  // repurpose the error number for collecting this pixels raw intensity
                  m_error += intensity;
                  /* restore intensity, so that it will not be shifted for the other detectors */
                    intensity = intensity_save;
                  }
              } /* end of if k2 ... */
        } /* end of for iy=0... */
    } /* end of for ix = ... */

  
  
  // Divide each image by its number of averages again:
  // add intensity squared to image2 for this detector and pixel, then rescale:
  m_image2[offsetX][offsetY] += m_error*m_error;
  m_image2[offsetX][offsetY] /= m_Navg+1;	
    
  // do the rescaling for the average image:
  m_image[offsetX][offsetY] /= m_Navg+1;	
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

void DetectorManager::CollectIntensity(WavePtr &wave, int plane_idx)
{
  for (size_t det_idx=0; det_idx<m_detectors[plane_idx].size(); det_idx++)
    {
      m_detectors[plane_idx][det_idx]->CollectIntensity(wave);
    }
}

void DetectorManager::SaveDetectors(std::map<std::string, double> &params)
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
          
          det->WriteImage(params);
        }
    }
  /*
  int i, ix, islice;
  double intensity;
  static char fileName[256]; 
  //imageStruct *header = NULL;
  std::vector<DetectorPtr> detectors;
  float t;
  std::vector<unsigned> position(1,0);

  int tCount = (int)(ceil((double)((muls->slices * muls->cellDiv) / muls->outputInterval)));
  
  // Loop over slices (intermediates)
  for (islice=0; islice <= tCount; islice++)
    {
      if (islice<tCount)
        {
          t = ((islice+1) * muls->outputInterval ) * muls->sliceThickness;
        }
      else
        {
          t = muls->slices*muls->cellDiv*muls->sliceThickness;
        }
      detectors = muls->detectors[islice];
      // write the output STEM images:
      // This is done only after all pixels have completed, so that image is complete.
      for (i=0; i<muls->detectorNum; i++) 
        {
          // calculate the standard error for this image:
          detectors[i]->error = 0;
          intensity             = 0;
          for (ix=0; ix<muls->scanXN * muls->scanYN; ix++) 
            {
              detectors[i]->error += (detectors[i]->image2[0][ix]-
                                      detectors[i]->image[0][ix] * detectors[i]->image[0][ix]);
              intensity += detectors[i]->image[0][ix] * detectors[i]->image[0][ix];
            }
          detectors[i]->error /= intensity;
          position[0]=islice;
          sprintf(fileName,"%s/%s", muls->folder, detectors[i]->name);
          params["Thickness"]=t;
          params["Runs Averaged"]=(double)muls->avgCount+1;
          params["Error"]=(double)detectors[i]->error;

          // TODO: why is this here?  It's a second image.  Why aren't we just saving another image?
          //  for (ix=0; ix<muls->scanXN * muls->scanYN; ix++) 
          //  {
          //  detectors[i]->SetParameter(2+ix, (double)detectors[i]->image2[0][ix]);
          //  }

          // exclude the suffix if this is the last detector (at the final thickness)
          if (islice==tCount)
            detectors[i]->WriteImage(fileName, detectors[i]->name, params);
          else
            detectors[i]->WriteImage(fileName, detectors[i]->name, params, position);
        }
    }
*/
}

void DetectorManager::PrintDetectors()
{
  printf("* Number of detectors:  %d\n",m_detectors.size());
                
  for (size_t i=0;i<m_detectors.size();i++) {
    printf("* %d (\"%s\"):",i+1,m_detectors[0][i]->m_name.c_str());
    for (size_t j=0;j<14-m_detectors[0][i]->m_name.size();j++) printf(" ");
    printf(" %g .. %g mrad = (%.2g .. %.2g 1/A)\n",
           m_detectors[0][i]->m_rInside,
           m_detectors[0][i]->m_rOutside,
           m_detectors[0][i]->m_k2Inside,
           m_detectors[0][i]->m_k2Outside);
    if ((m_detectors[0][i]->m_shiftX != 0) ||(m_detectors[0][i]->m_shiftY != 0))
      printf("*   center shifted:     dkx=%g, dky=%g\n",
             m_detectors[0][i]->m_shiftX,m_detectors[0][i]->m_shiftY);
  }
}

