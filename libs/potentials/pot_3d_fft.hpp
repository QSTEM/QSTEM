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

#include "pot_3d.hpp"

namespace QSTEM
{

class C3DFFTPotential : public C3DPotential
{
public:
  C3DFFTPotential();
  C3DFFTPotential(const ConfigReaderPtr &configReader);
  virtual void DisplayParams();
  //virtual void makeSlices(int nlayer, char *fileName, atom *center);
  virtual void AddAtomToSlices(std::vector<atom>::iterator &atom, float_tt atomX, float_tt atomY, float_tt atomZ);
protected:
  virtual void AddAtomPeriodic(std::vector<atom>::iterator &atom, 
                         float_tt atomBoxX, unsigned int ix, 
                         float_tt atomBoxY, unsigned int iy, 
                         float_tt atomZ);
  virtual void AddAtomNonPeriodic(std::vector<atom>::iterator &atom, 
                         float_tt atomBoxX, unsigned int ix, 
                         float_tt atomBoxY, unsigned int iy, 
                         float_tt atomZ);
  void GetAtomPotential3D(unsigned Znum, float_tt B,unsigned &nzSub,unsigned &Nr,unsigned &Nz_lut,
                                 ComplexVector &output);
  void GetAtomPotentialOffset3D(unsigned Znum, float_tt B,unsigned &nzSub,unsigned &Nr,unsigned &Nz_lut,float_tt q,
                                ComplexVector &output);
private:
  // collection of atom potentials
  std::map<unsigned, ComplexVector> m_atPot;
  std::map<unsigned, ComplexVector> m_offsetPot;
private:	
  friend class CPotFactory;
  // Create an instance of this class, wrapped in a shared ptr
  //     This should not be inherited - any subclass needs its own implementation.
  static PotPtr Create() {return PotPtr(new C3DFFTPotential());}
  static PotPtr Create(const ConfigReaderPtr &configReader){return PotPtr(new C3DFFTPotential(configReader));}
};

}