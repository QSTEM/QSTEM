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

bool C2DPotential::CheckAtomZInBounds(float_tt atomZ)
{
  /*
   * c = the thickness of the current slab.
   *
   * if the z-position of this atom is outside the potential slab
   * we won't consider it and skip to the next
   */
  return ((atomZ<m_c) && (atomZ>=0));
}

void C2DPotential::AddAtomToSlices(std::vector<atom>::iterator &atom, float_tt atomX, float_tt atomY,
                                               float_tt atomZ)
{
  // Note that if you override this method, you should do the following check to make sure the atom is in bounds.
  // skip atoms that are beyond the cell's boundaries
  if (!m_periodicZ)
    {
      if (atomZ > c) return;
      if ((atomZ >=0)
  

  AddAtomToSlicesRealSpaceLUT(atom, atomX, atomY, atomZ);
}

void C2DPotential::AddAtomToSlicesRealSpaceLUT(std::vector<atom>::iterator &atom, float_tt atomX, float_tt atomY,
                                               float_tt atomZ)
{
  complex_tt dPot;
  if (!periodicZ) {
    if (iAtomZ < 0) return;
    if (iAtomZ >= nlayer) return;        
  }                
  iz = (iAtomZ+32*nlayer) % nlayer;         /* shift into the positive range */
  atomBoxLookUp(&dPot,atom->Znum,x,y,0, m_tds ? 0 : atom->dw);
  z = (double)(iAtomZ+1)*m_cz[0]-atomZ;

  /* split the atom if it is close to the top edge of the slice */
  if ((z<0.15*(*muls).cz[0]) && (iz >0)) {
    m_trans[iz][ix][iy][0] += 0.5*dPot[0];
    m_trans[iz][ix][iy][1] += 0.5*dPot[1];
    m_trans[iz-1][ix][iy][0] += 0.5*dPot[0];
    m_trans[iz-1][ix][iy][1] += 0.5*dPot[1];                        
  }
  /* split the atom if it is close to the bottom edge of the slice */
  else {
    if ((z>0.85*m_cz[0]) && (iz < nlayer-1)) {
      m_trans[iz][ix][iy][0] += 0.5*dPot[0];
      m_trans[iz][ix][iy][1] += 0.5*dPot[1];        
      m_trans[iz+1][ix][iy][0] += 0.5*dPot[0];
      m_trans[iz+1][ix][iy][1] += 0.5*dPot[1];                 
    }
    else {
      m_trans[iz][ix][iy][0] += dPot[0];
      m_trans[iz][ix][iy][1] += dPot[1];        
    }
  }
}

void C2DPotential::CenterAtomZ(std::vector<atom>::iterator &atom, float_tt &z)
{
  CPotential::CenterAtomZ(atom, z);
  z += 0.5*muls->sliceThickness;
}
