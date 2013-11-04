#include "potential.hpp"

PotPtr GetPotential(bool _3D, bool fft, std::string cfg_file)
{
	if (_3D)
	{
		if (fft)
		{
			return PotPtr(new C3DFFTPotential(cfg_file));
		}
		else
		{
			return PotPtr(new C3DPotential(cfg_file));
		}
	}
	else
	{
		if (fft)
		{
			return PotPtr(new C2DFFTPotential(cfg_file));
		}
		else
		{
			return PotPtr(new C2DPotential(cfg_file));
		}
	}
}