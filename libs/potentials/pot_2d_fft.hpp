#include "pot_2d.hpp"

class C2DFFTPotential : public C2DPotential
{
public:
  C2DFFTPotential(ConfigReaderPtr &configReader);
  virtual void makeSlices(int nlayer, char *fileName, atom *center);
};
