#include "potential.hpp"

PotPtr GetPotential(ConfigReaderPtr &configReader)
{
  bool _3D, fft;
  configReader->ReadPotentialCalculationParameters(fft, _3D);

  if (_3D)
    {
      if (fft)
        {
          return PotPtr(new C3DFFTPotential(configReader));
        }
      else
        {
          return PotPtr(new C3DPotential(configReader));
        }
    }
  else
    {
      if (fft)
        {
          return PotPtr(new C2DFFTPotential(configReader));
        }
      else
        {
          return PotPtr(new C2DPotential(configReader));
        }
    }
}
