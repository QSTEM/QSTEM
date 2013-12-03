#ifndef POTENTIAL_3D_H
#define POTENTIAL_3D_H

#include "pot_base.hpp"

class C3DPotential : public CPotential
{
public:
  C3DPotential(ConfigReaderPtr &configReader);
  virtual void atomBoxLookUp(complex_tt &val, int Znum, float_tt x, float_tt y, float_tt z, float_tt B);
  virtual void makeSlices(int nlayer, char *fileName, atom *center);
  void CenterAtomZ(std::vector<atom>::iterator &atom, float_tt &z);
  bool CheckAtomZInBounds(float_tt atomZ);
  virtual void AddAtomToSlices(std::vector<atom>::iterator &atom, 
                               float_tt atomX, float_tt atomY, float_tt atomZ);
  void _AddAtomRealSpace(std::vector<atom>::iterator &atom, 
                         float_tt atomBoxX, unsigned int ix, 
                         float_tt atomBoxY, unsigned int iy, 
                         float_tt atomZ, unsigned int iatomZ);
};

#endif
