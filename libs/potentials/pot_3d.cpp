#include "pot_3d.hpp"

C3DPotential::C3DPotential(std::string cfg_file) : CPotential(cfg_file)
{
	m_boxNz = (int)(m_radius/m_ddz+2.0);
}

void C3DPotential::atomBoxLookUp(complex_tt &val, int Znum, float_tt x, float_tt y, float_tt z, float_tt B)
{
	float_tt dx, dy, dz;
	int ix, iy, iz;

	// does the atom box lookup or calculation
	CPotential::atomBoxLookUp(val, Znum, x, y, z, B);

	/***************************************************************
	* Do the trilinear interpolation
	*/
	val[0] = 0.0;
	val[1] = 0.0;
	if (x*x+y*y+z*z > m_radius2) {
		return;
	}
	x = fabs(x);
	y = fabs(y);
	z = fabs(z);
	ix = (int)(x/m_ddx);
	iy = (int)(y/m_ddy);
	iz = (int)(z/m_ddz);
	dx = x-(float_tt)ix*m_ddx;
	dy = y-(float_tt)iy*m_ddy;
	dz = z-(float_tt)iz*m_ddz;
	if ((dx < 0) || (dy<0) || (dz<0)) {
		/* printf("Warning, dx(%g), dy(%g), dz(%g) < 0, (x=%g, y=%g, z=%g)\n",dx,dy,dz,x,y,z);
		*/
		if (dx < 0) dx = 0.0;
		if (dy < 0) dy = 0.0;
		if (dz < 0) dz = 0.0;
	}
	
	if (m_atomBoxes[Znum]->B > 0) {
		val[0] = (1.0-dz)*((1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->potential[iz][ix][iy][0]+
			dx*m_atomBoxes[Znum]->potential[iz][ix+1][iy][0])+
			dy*((1.0-dx)*m_atomBoxes[Znum]->potential[iz][ix][iy+1][0]+
			dx*m_atomBoxes[Znum]->potential[iz][ix+1][iy+1][0]))+
			dz*((1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->potential[iz+1][ix][iy][0]+
			dx*m_atomBoxes[Znum]->potential[iz+1][ix+1][iy][0])+
			dy*((1.0-dx)*m_atomBoxes[Znum]->potential[iz+1][ix][iy+1][0]+
			dx*m_atomBoxes[Znum]->potential[iz+1][ix+1][iy+1][0]));
			val[1] = (1.0-dz)*((1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->potential[iz][ix][iy][1]+
			dx*m_atomBoxes[Znum]->potential[iz][ix+1][iy][1])+
			dy*((1.0-dx)*m_atomBoxes[Znum]->potential[iz][ix][iy+1][1]+
			dx*m_atomBoxes[Znum]->potential[iz][ix+1][iy+1][1]))+
			dz*((1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->potential[iz+1][ix][iy][1]+
			dx*m_atomBoxes[Znum]->potential[iz+1][ix+1][iy][1])+
			dy*((1.0-dx)*m_atomBoxes[Znum]->potential[iz+1][ix][iy+1][1]+
			dx*m_atomBoxes[Znum]->potential[iz+1][ix+1][iy+1][1]));
	}
	else {
		val[0] = (1.0-dz)*((1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->rpotential[iz][ix][iy]+
			dx*m_atomBoxes[Znum]->rpotential[iz][ix+1][iy])+
			dy*((1.0-dx)*m_atomBoxes[Znum]->rpotential[iz][ix][iy+1]+
			dx*m_atomBoxes[Znum]->rpotential[iz][ix+1][iy+1]))+
			dz*((1.0-dy)*((1.0-dx)*m_atomBoxes[Znum]->rpotential[iz+1][ix][iy]+
			dx*m_atomBoxes[Znum]->rpotential[iz+1][ix+1][iy])+
			dy*((1.0-dx)*m_atomBoxes[Znum]->rpotential[iz+1][ix][iy+1]+
			dx*m_atomBoxes[Znum]->rpotential[iz+1][ix+1][iy+1]));
	}
}
	