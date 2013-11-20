#include "pot_2d_fft.hpp"

C2DFFTPotential::C2DFFTPotential(ConfigReaderPtr &configReader) : C2DPotential(configReader)
{
}


void C2DFFTPotential::makeSlices(int nlayer, char *fileName, atom *center)
{
  /* check whether we have constant slice thickness */
  for (i = 0;i<nlayer;i++) if ((*muls).cz[0] != (*muls).cz[i]) break;
  if (i<nlayer) printf("Warning: slice thickness not constant, will give wrong results (iz=%d)!\n",i);

  
}
