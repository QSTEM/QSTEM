#ifndef POTENTIAL_2D_H
#define POTENTIAL_2D_H

#include "pot_base.hpp"

class C2DPotential : public CPotential
{
public:
	C2DPotential(std::string cfg_file);
	virtual void atomBoxLookUp(complex_tt &val, int Znum, float_tt x, float_tt y, float_tt z, float_tt B);
};

#endif