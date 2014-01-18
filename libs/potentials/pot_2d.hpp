#ifndef POTENTIAL_2D_H
#define POTENTIAL_2D_H

#include "pot_base.hpp"

class C2DPotential : public CPotential
{
public:
  C2DPotential();
  C2DPotential(const ConfigReaderPtr &configReader);
  virtual void Initialize();
  virtual void Initialize(const ConfigReaderPtr &configReader);
  virtual void DisplayParams();
  virtual void AtomBoxLookUp(complex_tt &val, int Znum, float_tt x, float_tt y, float_tt z, float_tt B);
  //virtual void MakeSlices(int nlayer, char *fileName, atom *center);
  bool CheckAtomZInBounds(float_tt atomZ);
  void CenterAtomZ(std::vector<atom>::iterator &atom, float_tt &z);
  virtual void AddAtomToSlices(std::vector<atom>::iterator &atom, 
                               float_tt atomX, float_tt atomY, float_tt atomZ);
  void _AddAtomRealSpace(std::vector<atom>::iterator &atom, 
                         float_tt atomBoxX, unsigned int ix, 
                         float_tt atomBoxY, unsigned int iy, 
                         float_tt atomZ, unsigned int iatomZ);
private:
  friend class CPotFactory;
  // Create an instance of this class, wrapped in a shared ptr
  //     This should not be inherited - any subclass needs its own implementation.
  static PotPtr Create(){return PotPtr(new C2DPotential());}
  static PotPtr Create(const ConfigReaderPtr &configReader){return PotPtr(new C2DPotential(configReader));}
};

#endif
