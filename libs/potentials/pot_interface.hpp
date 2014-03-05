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

#include "stemtypes_fftw3.hpp"
#include "config_IO/read_interface.hpp"

#ifndef POTENTIAL_INTERFACE_H
#define POTENTIAL_INTERFACE_H

namespace QSTEM
{

class IPotential;
typedef boost::shared_ptr<IPotential> PotPtr;
typedef PotPtr (*CreatePotentialFn)(void); 

class IPotential
{
public:
  // Any Potential class should define this as a private static member.  See pot_2d.hpp for example.
  // virtual PotPtr __stdcall Create()=0;

  // Resizes memory arrays, computes initial parameters based on start state
  virtual void Initialize()=0;
  // Same as Initialize, but before running Initialize(), it loads appropriate parameters from configReader
  virtual void Initialize(const ConfigReaderPtr &configReader){};


  virtual void DisplayParams(){};
  //virtual void initSTEMSlices();
  // encapsulates make slices and initSTEMslices - used to refresh the potential with a new structure (after a random
  //    shake)
  virtual void Refresh()=0;
  // You must provide a way of getting  slices somehow.  You must implement (at least) one of the following two
  //   functions in order for your Potential to be useful at all.
  // MakeSlices is called when you want to build the potential from atomic positions
  virtual void MakeSlices(int nlayer,char *fileName,atom *center){};
  // ReadPotential is used when you want to load the potential slices from some data files
  virtual void ReadPotential(std::string &fileName, unsigned subSlabIdx){};
  
  virtual unsigned GetNSlices() const =0;
  virtual float_tt GetSliceThickness() const =0;
  virtual float_tt GetSliceThickness(unsigned idx) const =0;
  virtual void GetSizePixels(unsigned &nx, unsigned &ny) const =0;

  virtual void WriteSlice(unsigned idx)=0;
  virtual void WriteProjectedPotential()=0;
  virtual void GetSlice(unsigned idx, ComplexVector &vec)=0;
  virtual complex_tt GetSlicePixel(unsigned iz, unsigned ix, unsigned iy)=0;
};

}

#endif
