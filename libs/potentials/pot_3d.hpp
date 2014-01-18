#ifndef POTENTIAL_3D_H
#define POTENTIAL_3D_H

#include "pot_base.hpp"

class C3DPotential : public CPotential
{
public:
	C3DPotential();
  C3DPotential(const ConfigReaderPtr &configReader);
  virtual void DisplayParams();
  virtual void atomBoxLookUp(complex_tt &val, int Znum, float_tt x, float_tt y, float_tt z, float_tt B);
  virtual void makeSlices(int nlayer, char *fileName, atom *center);
  void CenterAtomZ(std::vector<atom>::iterator &atom, float_tt &z);
  bool CheckAtomZInBounds(float_tt atomZ);
  virtual void AddAtomToSlices(std::vector<atom>::iterator &atom, 
                               float_tt atomX, float_tt atomY, float_tt atomZ);
  void _AddAtomRealSpace(std::vector<atom>::iterator &atom, 
                         float_tt atomBoxX, unsigned ix,
                         float_tt atomBoxY, unsigned iy,
                         float_tt atomZ, unsigned iAtomZ);
protected:
  unsigned m_sliceStep;  // number of float_tt to advance to next slice (2*m_nx*m_ny)

private:
  friend class CPotFactory;
  // Create an instance of this class, wrapped in a shared ptr
  //     This should not be inherited - any subclass needs its own implementation.
  static PotPtr Create() {return PotPtr(new C3DPotential());}
  static PotPtr Create(const ConfigReaderPtr &configReader){return PotPtr(new C3DPotential(configReader));}
};

#endif










