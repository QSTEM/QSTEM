#include "pot_2d.hpp"

C2DPotential::C2DPotential(ConfigReaderPtr &configReader) : CPotential(configReader)
{
	m_boxNz = 1;
}

void C2DPotential::atomBoxLookUp(complex_tt &sum, int Znum, float_tt x, float_tt y, float_tt z, float_tt B) 
{
	float_tt dx, dy;
	int ix, iy;

	// does the atom box lookup or calculation
	CPotential::atomBoxLookUp(sum, Znum, x, y, z, B);

	/***************************************************************
	* Do the trilinear interpolation
	*/
	sum[0] = 0.0;
	sum[1] = 0.0;
	if (x*x+y*y+z*z > m_radius2) {
		return;
	}
	x = fabs(x);
	y = fabs(y);
	ix = (int)(x/m_ddx);
	iy = (int)(y/m_ddy);
	dx = x-(float_tt)ix*m_ddx;
	dy = y-(float_tt)iy*m_ddy;

	if ((dx < 0) || (dy<0) ) {
		/* printf("Warning, dx(%g), dy(%g), dz(%g) < 0, (x=%g, y=%g, z=%g)\n",dx,dy,dz,x,y,z);
		*/
		if (dx < 0) dx = 0.0;
		if (dy < 0) dy = 0.0;
	}

	if (m_atomBoxes[Znum]->B > 0) {
		sum[0] = (1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->potential[0][ix][iy][0]+
			dx*m_atomBoxes[Znum]->potential[0][ix+1][iy][0])+
			dy*((1.0-dx)*m_atomBoxes[Znum]->potential[0][ix][iy+1][0]+
			dx*m_atomBoxes[Znum]->potential[0][ix+1][iy+1][0]);
		sum[1] = (1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->potential[0][ix][iy][1]+
			dx*m_atomBoxes[Znum]->potential[0][ix+1][iy][1])+
			dy*((1.0-dx)*m_atomBoxes[Znum]->potential[0][ix][iy+1][1]+
			dx*m_atomBoxes[Znum]->potential[0][ix+1][iy+1][1]);
	}
	else {
		sum[0] = (1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->rpotential[0][ix][iy]+
			dx*m_atomBoxes[Znum]->rpotential[0][ix+1][iy])+
			dy*((1.0-dx)*m_atomBoxes[Znum]->rpotential[0][ix][iy+1]+
			dx*m_atomBoxes[Znum]->rpotential[0][ix+1][iy+1]);
	}
}
