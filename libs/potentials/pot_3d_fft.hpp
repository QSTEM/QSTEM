#include "pot_3d.hpp"

class C3DFFTPotential : public C3DPotential
{
public:
  C3DFFTPotential(ConfigReaderPtr &configReader);
  virtual void DisplayParams();
  virtual void makeSlices(int nlayer, char *fileName, atom *center);
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
  complex_tt *GetAtomPotential3D(unsigned Znum, float_tt B,unsigned &nzSub,unsigned &Nr,unsigned &Nz_lut);
  complex_tt *GetAtomPotentialOffset3D(unsigned Znum, float_tt B,unsigned &nzSub,unsigned &Nr,unsigned &Nz_lut,float_tt q);

};
