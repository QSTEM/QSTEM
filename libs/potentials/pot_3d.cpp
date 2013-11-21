#include "pot_3d.hpp"

C3DPotential::C3DPotential(ConfigReaderPtr &configReader) : CPotential(configReader)
{
	m_boxNz = (int)(m_atomRadius/m_ddz+2.0);
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
	if (x*x+y*y+z*z > m_atomRadius2) {
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
	
void C3DPotential::AddAtomToSlices(std::vector<atom>::iterator &atom, float_tt atomX, float_tt atomY, float_tt atomZ)
{
  if (!m_periodicZ && ) return;
  return AddAtomToSlicesRealSpaceLUT(atom, atomZ)
}


bool C3DPotential::CheckAtomZInBounds(float_tt atomZ)
{
  /*
   * c = the thickness of the current slab.
   *
   * if the z-position of this atom is outside the potential slab
   * we won't consider it and skip to the next
   */
  return ((atomZ - m_atomRadius > m_c) && (atomZ + m_atomRadius + m_sliceThickness >= 0));
}

void C3DPotential::AddAtomToSlicesRealSpaceLUT(std::vector<atom>::iterator &atom, float_tt atomZ)
{
  complex_tt dPot;
  /* calculate the range which we have left to cover with z-variation */
  /* iRadZ is the number of slices (rounded up) that this atom
   * will contribute to, given its current x,y-radius
   */
  iRadZ = (int)(sqrt(atomRadius2-r2sqr)/m_cz[0]+1.0);
  /* loop through the slices that this atoms contributes to */
  for (unsigned iaz=-iRadZ;iaz <=iRadZ;iaz++) {
    if (!m_periodicZ) {
      if (iaz+iAtomZ < 0) {
        if (-iAtomZ <= iRadZ) iaz = -iAtomZ;
        else break;
        if (abs(iaz)>nlayer) break;
      }
      if (iaz+iAtomZ >= nlayer)        break;
    }
    z = (double)(iAtomZ+iaz+0.5)*m_cz[0]-atomZ;
    /* shift into the positive range */
    iz = (iaz+iAtomZ+32*nlayer) % nlayer;        
    /* x,y,z is the true vector from the atom center
     * We can look up the proj potential at that spot
     * using trilinear extrapolation.
     */
    atomBoxLookUp(&dPot,atom->Znum,x,y,z,
                  m_tds ? 0 : atom->dw);
    // printf("access: %d %d %d\n",iz,ix,iy);
    m_trans[iz][ix][iy][0] += dPot[0];
    m_trans[iz][ix][iy][1] += dPot[1];         
  } /* end of for iaz=-iRadZ .. iRadZ */
}


void C3DPotential::CenterAtomZ(std::vector<atom>::iterator &atom, float_tt &z)
{
  CPotential::CenterAtomZ(atom, z);
  z -= m_sliceThickness;
}








