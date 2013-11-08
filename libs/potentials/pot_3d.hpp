#ifndef POTENTIAL_3D_H
#define POTENTIAL_3D_H

#include "pot_base.hpp"

class C3DPotential : public CPotential
{
public:
	C3DPotential(ConfigReaderPtr &configReader);
	virtual void atomBoxLookUp(complex_tt &val, int Znum, float_tt x, float_tt y, float_tt z, float_tt B);
};

#endif
