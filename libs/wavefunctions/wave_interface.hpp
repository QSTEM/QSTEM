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

#ifndef WAVE_INTERFACE_H
#define WAVE_INTERFACE_H

#include <map>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include "config_IO/config_reader_factory.hpp"
#include "stemtypes_fftw3.hpp"

class IWave;
typedef boost::shared_ptr<IWave> WavePtr;
typedef WavePtr (*CreateWaveFn)(const ConfigReaderPtr &reader);

// a structure for a probe/parallel beam wavefunction.
// Separate from mulsliceStruct for parallelization.
class IWave
{
public:
  virtual void CreateDataSets()=0;
  virtual void FormProbe()=0;

  virtual void DisplayParams()=0;

  virtual void ToRealSpace()=0;
  virtual void ToFourierSpace()=0;
  virtual bool IsRealSpace()=0;

  virtual WavePtr Clone()=0;

  virtual void GetSizePixels(unsigned &x, unsigned &y) const =0;
  virtual unsigned GetTotalPixels() const =0;
  virtual void GetResolution(float_tt &x, float_tt &y) const =0;
  virtual void GetPositionOffset(unsigned &x, unsigned &y) const =0;
  virtual float_tt GetK2(unsigned ix, unsigned iy) const =0;
  virtual float_tt GetKX2(unsigned ix) const =0;
  virtual float_tt GetKY2(unsigned iy) const =0;
  virtual float_tt GetK2Max() const =0;

  virtual void Resize(unsigned x, unsigned y) = 0;

  virtual float_tt GetVoltage()  const =0;
  virtual float_tt GetWavelength()  const =0;

  virtual float_tt GetPixelIntensity(unsigned i) const =0;
  virtual float_tt GetPixelIntensity(unsigned x, unsigned y) const =0;
  virtual float_tt GetDiffPatPixel(unsigned i)  const =0;
  virtual float_tt GetDiffPatPixel(unsigned x, unsigned y) const =0;
  virtual void SetDiffPatPixel(unsigned i, float_tt value) =0;
  virtual void SetDiffPatPixel(unsigned x, unsigned y, float_tt value) =0;

  virtual void ApplyTransferFunction(boost::shared_array<complex_tt> &wave)=0;

  void WriteBeams(unsigned absoluteSlice);

  virtual void WriteProbe()=0;

  // Methods for writing wavefunction
  // Method 1: no position reference
  virtual void WriteWave(std::string comment="Wavefunction")=0;
  // Method 2: pass an unsigned integer.  Use for either output at varying thickness, or navg
  virtual void WriteWave(unsigned navg, std::string comment="Wavefunction", 
                 std::map<std::string, double>params = std::map<std::string, double>())=0;
  // Method 3: pass two unsigned integers.  Use for output at varying STEM probe position
  virtual void WriteWave(unsigned posX, unsigned posY, std::string comment="Wavefunction", 
                 std::map<std::string, double>params = std::map<std::string, double>())=0;
  // Method 4: pass three unsigned integers.  Use for output at varying STEM probe position, and at navg/thickness
  virtual void WriteWave(unsigned posX, unsigned posY, unsigned posZ, std::string comment="Wavefunction", 
                 std::map<std::string, double>params = std::map<std::string, double>())=0;

  // Methods for writing diffraction pattern
  // Method 1: no position reference
  virtual void WriteDiffPat(std::string comment="Diffraction Pattern", 
                 std::map<std::string, double>params = std::map<std::string, double>())=0;
  // Method 2: pass an unsigned integer.  Use for either output at varying thickness, or navg
  virtual void WriteDiffPat(unsigned navg, std::string comment="Diffraction Pattern", 
                 std::map<std::string, double>params = std::map<std::string, double>())=0;
  // Method 3: pass two unsigned integers.  Use for output at varying STEM probe position
  virtual void WriteDiffPat(unsigned posX, unsigned posY, std::string comment="Diffraction Pattern", 
                 std::map<std::string, double>params = std::map<std::string, double>())=0;
    // Method 4: pass three unsigned integers.  Use for output at varying STEM probe position, and at navg/thickness
  virtual void WriteDiffPat(unsigned posX, unsigned posY, unsigned posZ, std::string comment="Diffraction Pattern", 
                 std::map<std::string, double>params = std::map<std::string, double>())=0;

  // People can change the wavefunction - for example, that's what we have to do when we
  //    transmit the wave through the sample's potential.
  virtual complex_tt *GetWavePointer()=0;
  // People should not directly change the diffraction pattern, since we'll re-calculate it when 
  //   the wavefunction changes.
  //   They can, however, access it.
  virtual const float_tt *GetDPPointer()=0;

  virtual float_tt GetIntegratedIntensity() const =0;

  virtual void ReadWave()=0;
  virtual void ReadWave(unsigned navg)=0;
  virtual void ReadWave(unsigned posX, unsigned posY)=0;
  virtual void ReadDiffPat()=0;
  virtual void ReadDiffPat(unsigned navg)=0;
  virtual void ReadDiffPat(unsigned posX, unsigned posY)=0;
};

#endif
