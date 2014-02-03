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

#include "../stemtypes_fftw3.hpp"
#include "../imagelib_fftw3.hpp"
//#include "config_reader_factory.hpp"
#include "../memory_fftw3.hpp"
#include "../crystal.hpp"
#include "scatfactsRez.hpp"
#include <map>

#include "pot_interface.hpp"

#ifndef POTENTIAL_BASE_H
#define POTENTIAL_BASE_H

#define OVERSAMPLING 3
#define OVERSAMPLINGZ 3*OVERSAMPLING

class CPotential : public IPotential
{
public:
  CPotential();
  CPotential(unsigned nx, unsigned ny, unsigned nz, float_tt dx, float_tt dy, float_tt dz, float_tt atomRadius, float_tt v0);
  CPotential(const ConfigReaderPtr &configReader);
  //~CPotential();

  virtual void DisplayParams();

  void AtomBoxLookUp(complex_tt &val, int Znum, float_tt x, float_tt y, float_tt z, float_tt B);
  virtual void MakeSlices(int nlayer,char *fileName,atom *center);
  //virtual void initSTEMSlices();
  // encapsulates make slices and initSTEMslices - used to refresh the potential with a new structure (after a random
  //    shake)
  virtual void Refresh();
  // TODO: need abstracted structure reader
  //virtual void ReadAtoms();
  virtual void ReadPotential(std::string &fileName, unsigned subSlabIdx);
  virtual void CenterAtomZ(std::vector<atom>::iterator &atom, float_tt &z);
  virtual void AddAtomToSlices(std::vector<atom>::iterator &atom, float_tt atomX, float_tt atomY, float_tt atomZ)=0;
  void AddAtomRealSpace(std::vector<atom>::iterator &atom, float_tt atomX, float_tt atomY, float_tt atomZ);
  
  // *************************  Getters  ********************** 
  unsigned GetNSlices() const {return m_nslices;}
  float_tt GetSliceThickness() const {return m_sliceThickness;}
  float_tt GetSliceThickness(unsigned idx) const {return m_cz[idx];}
  void GetSizePixels(unsigned &nx, unsigned &ny) const;


  // *************************** Setters **********************
  //  For all of these, the potential needs to be recalculated after calling them.
  //    Each of them set a flag on this class indicating that, such that getting a slice will either recompute the potential for you,
  //       or raise an error (TODO!)
  // Set m_thickness, the uniform slice thickness.  Overrides m_cz.  Implicitly sets number of slices.
  void SetSliceThickness(float_tt thickness_Angstroms);
  // Set m_cz, the vector of thicknesses.  Implicitly also sets number of slices!
  void SetSliceThickness(std::vector<float_tt> thickness_Angstroms);
  // Set number of slices.  Implicitly sets m_thickness to the sub-slab height divided by nslices, and empties m_cz.
  void SetNSlices(unsigned slices);
  // Sets the structure being used to calculate the potential.
  void SetStructure(StructurePtr structure);


  void WriteSlice(unsigned idx);
  void WriteProjectedPotential();
  void GetSlice(unsigned idx, ComplexVector &vec){vec=m_trans[idx];}
  complex_tt GetSlicePixel(unsigned iz, unsigned ix, unsigned iy);

  atom GetAtom(unsigned idx){return m_atoms[idx];}

  // public members (should be moved to protected, and given getters/setters...
  bool m_potential3D;
  float_tt m_atomRadius;

protected:
  void Initialize();
  void SliceSetup();
  void ResizeSlices();
  void ReadSlice(const std::string &fileName, ComplexVector &slice, unsigned idx);
  virtual void _AddAtomRealSpace(std::vector<atom>::iterator &atom, 
                                 float_tt atomX, unsigned int ix, 
                                 float_tt atomY, unsigned int iy, 
                                 float_tt atomZ, unsigned int iatomZ)=0;


  ImageIOPtr m_imageIO;
  StructurePtr m_crystal;
  std::vector<ComplexVector> m_trans; //  The 3D specimen potential array as a vector of 1D vectors
  //complex_tt ***m_trans;    //  The 3D specimen potential array
  bool m_currentPotential;  // Indicates whether computed potential matches current parameters.  
                            //    Set to true after computing potential.  Reset to false when parameters change.
  unsigned m_nx, m_ny;    /* size of projected potential in pixels, possibly larger than wavefunc's nx/ny */
  float_tt m_dx, m_dy, m_dz;   // resolutions
  float_tt m_ddx, m_ddy, m_ddz;   // oversampled resolutions
  //
  int m_boxNx, m_boxNy, m_boxNz;
  float_tt m_v0;   // voltage
  std::map<unsigned, atomBoxPtr> m_atomBoxes;
  std::vector<atom> m_atoms;

  // *********** Slice parameters **********
  bool m_centerSlices;
  float_tt m_sliceThickness;
  float_tt m_zOffset; /* defines the offset for the first slice in fractional coordinates        */
  unsigned m_nslices;   // nslices is the number of slices PER SUB-SLAB!  Not the total.
  unsigned m_outputInterval;
  // vector of slice thicknesses.  Only really relevant when slice thickness is not uniform.
  std::vector<float_tt> m_cz;
  // vector of slice positions
  std::vector<float_tt> m_slicePos;

  bool m_tds;

  // m_c 
  float_tt m_c; // the thickness of the current sub-slab (in A)
  float_tt m_dr, m_atomRadius2;
  unsigned m_iRadX, m_iRadY, m_iRadZ, m_iRad2;
  bool m_periodicXY, m_periodicZ;

  float_tt m_offsetX, m_offsetY;

  // ********* multi-slab parameters ********
  unsigned m_cellDiv; // How many sub-slabs the model is divided into
  unsigned m_divCount; // How many sub-slabs we've already processed

  bool m_equalDivs;

  unsigned m_printLevel;
  unsigned m_displayPotCalcInterval;  /* show progress every .. atoms when computing potential */
  unsigned m_saveLevel;
  bool m_savePotential, m_saveProjectedPotential, m_plotPotential;

  int m_scatFactor;  // The scattering factor type.  One of: 0 (Doyle-Turner); 1 (Wieck-Kohl); 2 (Custom)
  
  std::string m_fileBase; // base filename for saving potential.  Will have slice index appended.
  std::string m_potFileBase; //base filename for loading potential from files.  Will have slice index appended.

  bool m_readPotential;

  float_tt sfLUT(float_tt s,int atKind);

  void splinh( float_tt x[], float_tt y[],
                   float_tt b[], float_tt c[], float_tt d[], int n);
  float_tt seval( float_tt *x, float_tt *y, float_tt *b, float_tt *c,
                  float_tt *d, int n, float_tt x0 );
};

#endif
