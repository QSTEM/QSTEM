#include "pot_3d.hpp"

class C3DFFTPotential : public C3DPotential
{
public:
  C3DFFTPotential(ConfigReaderPtr &configReader);
  virtual void makeSlices(int nlayer, char *fileName, atom *center);
};
