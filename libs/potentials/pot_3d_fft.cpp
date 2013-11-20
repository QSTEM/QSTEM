#include "pot_3d_fft.hpp"

C3DFFTPotential::C3DFFTPotential(ConfigReaderPtr &configReader) : C3DPotential(configReader)
{
}

void C3DFFTPotential::makeSlices(int nlayer, char *fileName, atom *center)
{
  for (i = 0;i<nlayer;i++) if ((*muls).cz[0] != (*muls).cz[i]) break;
  if (i<nlayer) printf("Warning: slice thickness not constant, will give wrong results (iz=%d)!\n",i);
}










