/*
  QSTEM - image simulation for TEM/STEM/CBED
  Copyright (C) 2000-2010  Christoph Koch
  Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "pot_base.hpp"

CPotential::CPotential(ConfigReaderPtr &configReader)
{
  configReader->ReadProbeArraySize(m_nx, m_ny);
  configReader->ReadResolution(m_dx, m_dy);
  configReader->ReadVoltage(m_v0);
  configReader->ReadPotentialOutputParameters(m_savePotential, m_saveProjectedPotential, m_plotPotential);
}

void CPotential::Initialize()
{
	m_ddx = m_dx/(double)OVERSAMPLING;
	m_ddy = m_dy/(double)OVERSAMPLING;
	m_ddz = m_dz/(double)OVERSAMPLINGZ;

	/* For now we don't care, if the box has only small 
		* prime factors, because we will not fourier transform it
		* especially not very often.
		*/
	m_boxNx = (int)(m_radius/m_ddx+2.0);  
	m_boxNy = (int)(m_radius/m_ddy+2.0);  
}

/****************************************************************************
* function: atomBoxLookUp - looks up potential at position x, y, z, relative to atom center
*
* Znum = element
* x,y,z = real space position (in A)
* B = Debye-Waller factor, B=8 pi^2 <u^2>
***************************************************************************/
void CPotential::atomBoxLookUp(complex_tt &val, int Znum, float_tt x, float_tt y, float_tt z, float_tt B) 
{
  int boxNx,boxNy,boxNz;
  float_tt dx,dy,dz,ddx,ddy,ddz;
  int ix,iy,iz;
  float_tt maxRadius2;
  char fileName[256],systStr[256];
  int tZ, tnx, tny, tnz, tzOversample;  
  float_tt tdx, tdy, tdz, tv0, tB;
  FILE *fp;
  int numRead = 0,dummy;

  /* initialize all the atoms to non-used */
  if (!m_atomBoxes.count(Znum)) {
    m_atomBoxes[Znum]=atomBoxPtr(new atomBox());
    m_atomBoxes[Znum]->potential = NULL;
    m_atomBoxes[Znum]->rpotential = NULL;
    m_atomBoxes[Znum]->B = -1.0;

    m_radius2 = m_radius*m_radius;

		
    if (m_printLevel > 2)
      printf("Atombox has real space resolution of %g x %g x %gA (%d x %d x %d pixels)\n",
             ddx,ddy,ddz,boxNx,boxNy,boxNz);
  }
  // printf("Debugging: %d %g %g: %g\n",Znum,m_atomBoxes[Znum]->B,B,fabs(m_atomBoxes[Znum]->B - B));
  
  /* Creating/Reading a atombox for every new kind of atom, but only as needed */
  if (fabs(m_atomBoxes[Znum]->B - B) > 1e-6) {
    //  printf("Debugging 1 (%d: %.7g-%.7g= %.7g), %d\n",
    //	   Znum,m_atomBoxes[Znum]->B,B,fabs(m_atomBoxes[Znum]->B - B),fabs(m_atomBoxes[Znum]->B - B) > 1e-6);
    m_atomBoxes[Znum]->B = B;
    /* Open the file with the projected potential for this particular element
     */
    sprintf(fileName,"potential_%d_B%d.prj",Znum,(int)(100.0*B));
    if ((fp=fopen(fileName,"r")) == NULL) {
      sprintf(systStr,"scatpot %s %d %g %d %d %d %g %g %g %d %g",
              fileName,Znum,B,boxNx,boxNy,boxNz,ddx,ddy,ddz,OVERSAMPLINGZ,m_v0);
      if (m_printLevel > 2) {
        printf("Could not find precalculated potential for Z=%d,"
               " will calculate now.\n",Znum);
        printf("Calling: %s\n",systStr);
      }
      system(systStr);
      for (dummy=0;dummy < 10000;dummy++);
      if ((fp=fopen(fileName,"r")) == NULL) {
        printf("cannot calculate projected potential using scatpot - exit!\n");	
        exit(0);
      }  
      
    }
    fgets(systStr,250,fp);
    sscanf(systStr,"%d %le %d %d %d %le %le %le %d %le\n",
           &tZ, &tB, &tnx, &tny, &tnz, &tdx, &tdy, &tdz, &tzOversample, &tv0);
    /* If the parameters in the file don't match the current ones,
     * we need to create a new potential file
     */
    if ((tZ != Znum) || (fabs(tB-B)>1e-6) || (tnx != boxNx) || (tny != boxNy) || (tnz != boxNz) ||
        (fabs(tdx-ddx) > 1e-5) || (fabs(tdy-ddy) > 1e-5) || (fabs(tdz-ddz) > 1e-5) || 
        (tzOversample != OVERSAMPLINGZ) || (tv0 != m_v0)) {
      if (m_printLevel > 2) {
        printf("Potential input file %s has the wrong parameters\n",fileName);
        printf("Parameters:\n"
               "file:    Z=%d, B=%.3f A^2 (%d, %d, %d) (%.7f, %.7f %.7f) nsz=%d V=%g\n"
               "program: Z=%d, B=%.3f A^2 (%d, %d, %d) (%.7f, %.7f %.7f) nsz=%d V=%g\n"
               "will create new potential file, please wait ...\n",
               tZ,tB,tnx,tny,tnz,tdx,tdy,tdz,tzOversample,tv0,
               Znum,B,boxNx,boxNy,boxNz,ddx,ddy,ddz,OVERSAMPLINGZ,m_v0);
      }
      /* Close the old file, Create a new potential file now 
       */
      fclose(fp);
      sprintf(systStr,"scatpot %s %d %g %d %d %d %g %g %g %d %g",
              fileName,Znum,B,boxNx,boxNy,boxNz,ddx,ddy,ddz,OVERSAMPLINGZ,m_v0);
      system(systStr);
      if ((fp=fopen(fileName,"r")) == NULL) {
        printf("cannot calculate projected potential using scatpot - exit!\n");
        exit(0);
      }  
      fgets(systStr,250,fp);
    }

    /* Finally we can read in the projected potential
     */
    if (B == 0) {
      m_atomBoxes[Znum]->rpotential = float3D(boxNz,boxNx,boxNy,"atomBox");
      numRead = fread(m_atomBoxes[Znum]->rpotential[0][0],sizeof(float_tt),
                      (size_t)(boxNx*boxNy*boxNz),fp);	
    }
    else {
      m_atomBoxes[Znum]->potential = complex3D(boxNz,boxNx,boxNy,"atomBox");
      numRead = fread(m_atomBoxes[Znum]->potential[0][0],sizeof(complex_tt),
                      (size_t)(boxNx*boxNy*boxNz),fp);
    }
    
    /* writeImage_old(m_atomBoxes[Znum]->potential[0],boxNx,boxNy, 0.0,"potential.img");
       system("showimage potential.img");
    */
    fclose(fp);
    
    if (numRead == boxNx*boxNy*boxNz) {
      if (m_printLevel > 1)
        printf("Sucessfully read in the projected potential\n");
    }
    else {
      printf("error while reading potential file %s: read %d of %d values\n",
             fileName,numRead,boxNx*boxNy*boxNz);
      exit(0);
    }
  }	
}


void CPotential::ReadPotential(std::string &fileName)
{
  /*************************************************************************
   * read the potential that has been created externally!
   */
    for (i=(divCount+1)*muls->slices-1,j=0;i>=(divCount)*muls->slices;i--,j++) {
      ReadSlice(fileName, trans[j], i);
    }
    return;
}

void CPotential::ReadAtoms()
{
  if (muls->avgCount == 0) {
    // if this is the first run, the atoms have already been
    // read during initialization
    natom = (*muls).natom;
    atoms = (*muls).atoms;
  }
  else {
    /*
      the following function makes an array of natom atoms from
      the input file (with x,y,z,dw,occ);
    */
    // the last parameter is handleVacancies. If it is set to 1 vacancies
    
    // and multiple occupancies will be handled.
    atoms = readUnitCell(&natom,fileName,muls,1);
    if (muls->printLevel>=3)
      printf("Read %d atoms from %s, tds: %d\n",natom,fileName,muls->tds);
    muls->natom = natom;
    muls->atoms = atoms;
  }

minX = maxX = atoms[0].x;
    minY = maxY = atoms[0].y;
    minZ = maxZ = atoms[0].z;

    for (i=0;i<natom;i++) {
      if (atoms[i].x < minX) minX = atoms[i].x;
      if (atoms[i].x > maxX) maxX = atoms[i].x;
      if (atoms[i].y < minY) minY = atoms[i].y;
      if (atoms[i].y > maxY) maxY = atoms[i].y;
      if (atoms[i].z < minZ) minZ = atoms[i].z;
      if (atoms[i].z > maxZ) maxZ = atoms[i].z;
    }
    /*
      printf("Root of mean square TDS displacement: %f A (wobble=%g at %gK) %g %g %g\n",
      sqrt(u2/natom),wobble,(*muls).tds_temp,ux,uy,uz);
    */
    if (muls->printLevel >= 2) {
      printf("range of thermally displaced atoms (%d atoms): \n",natom);
      printf("X: %g .. %g\n",minX,maxX);
      printf("Y: %g .. %g\n",minY,maxY);
      printf("Z: %g .. %g\n",minZ,maxZ);
    }
    /*
      define the center of our unit cell by moving the atom specified
      by "center" at position (0.5,0.5,0.0)
    */
          
    if (center != NULL) {
      dx = (*muls).ax/2.0f - (*center).x;        
      dy = (*muls).by/2.0f - (*center).y;        
      dz = -(*center).z;
      for (i=0;i<natom;i++) {
        atoms[i].x += dx;
        if (atoms[i].x < 0.0f) atoms[i].x += (*muls).ax;
        else if (atoms[i].x > (*muls).ax) atoms[i].x -= (*muls).ax;
        atoms[i].y += dy;
        if (atoms[i].y < 0.0f) atoms[i].y += (*muls).by;
        else if (atoms[i].y > (*muls).by) atoms[i].y -= (*muls).by;
        atoms[i].z += dz;
        if (atoms[i].z < 0.0f) atoms[i].z += (*muls).c;
        else if (atoms[i].z > (*muls).c) atoms[i].z -= (*muls).c;
      }
    }

    /**********************************************************
     * Sort the atoms in z.
     *********************************************************/
    qsort(atoms,natom,sizeof(atom),atomCompare);
    if ((*muls).cfgFile != NULL) {
      sprintf(buf,"%s/%s",muls->folder,muls->cfgFile);
      // append the TDS run number
      if (strcmp(buf+strlen(buf)-4,".cfg") == 0) *(buf+strlen(buf)-4) = '\0';
      if (muls->tds) sprintf(buf+strlen(buf),"_%d.cfg",muls->avgCount);
      else sprintf(buf+strlen(buf),".cfg");
      
      // printf("Will write CFG file <%s> (%d)\n",buf,muls->tds)
      writeCFG(atoms,natom,buf,muls);
      
      if (muls->readPotential) {
        sprintf(buf,"nanopot %s/%s %d %d %d %s",muls->folder,muls->cfgFile,
                ny,nx,muls->slices*muls->cellDiv,muls->folder);
        system(buf);
      }
    }
}

/*****************************************************
* void make3DSlices()
*
* This function will create a 3D potential from whole
* unit cell, slice it, and make transr/i and propr/i
* Call this function with center = NULL, if you don't
* want the array to be shifted.
****************************************************/
void CPotential::makeSlices(MULS *muls,int nlayer,char *fileName, atom *center) 
{
  // FILE *fpu2;
  char filename[512];
  int natom,iatom,iz; /* number of atoms */
  atom *atoms;
  float_tt dx,dy,dz;
  float_tt c,atomX,atomY,atomZ;
  int i=0,j,nx,ny,ix,iy,iax,iay,iaz,sliceStep;
  int iAtomX,iAtomY,iAtomZ,iRadX,iRadY,iRadZ,iRad2;
  int iax0,iax1,iay0,iay1,iaz0,iaz1,nyAtBox,nyAtBox2,nxyAtBox,nxyAtBox2,iOffsX,iOffsY,iOffsZ;
  int nzSub,Nr,ir,Nz_lut;
  int iOffsLimHi,iOffsLimLo,iOffsStep;

  float_tt *slicePos;
  double z,x,y,r,ddx,ddy,ddr,dr,r2sqr,x2,y2,potVal,dOffsZ;
  char buf[BUF_LEN];
  FILE *sliceFp;
  float_tt minX,maxX,minY,maxY,minZ,maxZ;
  double atomRadius2;
  time_t time0,time1;
  float s11,s12,s21,s22;
  fftwf_complex        *atPotPtr;
  float *potPtr=NULL, *ptr;
  static int divCount = 0;
  ImageIOPtr imageIO = ImageIOPtr(new CImageIO(muls->potNx,muls->potNy));
  complex_tt dPot;
#if Z_INTERPOLATION
  double ddz;
#endif
#if USE_Q_POT_OFFSETS
  fftwf_complex        *atPotOffsPtr;
#endif

  // parameters to be saved to output files (e.g. thickness)
  std::map<std::string, double> params;

  /* return, if there is nothing to do */
  if (nlayer <1)
    return;

  nx = muls->potNx;
  ny = muls->potNy;

  /* we need to keep track of which subdivision of the unit cell we are in
   * If the cell is not subdivided, then muls.cellDiv-1 = 0.
   */
  if ((divCount == 0) || (muls->equalDivs))
    divCount = muls->cellDiv;
  divCount--;

  /* we only want to reread and shake the atoms, if we have finished the
   * current unit cell.
   */
  if (divCount == muls->cellDiv-1) {
    ReadAtoms();
    
  } /* end of if divCount==cellDiv-1 ... */
  else {
    natom = muls->natom;
    atoms = muls->atoms;
  }

  /**************************************************************
   *        setup the slices with their start and end positions
   *        then loop through all the atoms and add their potential to
   *        the slice that their potential reaches into (up to RMAX)
   *************************************************************/
  // c = (*muls).c/(float_tt)((*muls).cellDiv);
  c = muls->sliceThickness * muls->slices;
  dx = (*muls).resolutionX;
  dy = (*muls).resolutionY;
  dr = muls->resolutionX/OVERSAMP_X; // define step width in which radial V(r,z) is defined
  iRadX = (int)ceil((*muls).atomRadius/dx);
  iRadY = (int)ceil((*muls).atomRadius/dy);
  iRadZ = (int)ceil((*muls).atomRadius/muls->sliceThickness);
  iRad2 = iRadX*iRadX+iRadY*iRadY;
  atomRadius2 = (*muls).atomRadius * (*muls).atomRadius;
  
  if (muls->printLevel >= 3) {
    printf("Slab thickness: %gA z-offset: %gA (cellDiv=%d)\n",
           c,c*(float_tt)(muls->cellDiv-divCount-1),divCount);
  }        
  /*******************************************************
   * initializing slicPos, cz, and transr
   *************************************************************/
  if ((*muls).cz == NULL) {
    (*muls).cz = float1D(nlayer,"cz");
  }
  // sliceFp = fopen(sliceFile,"r");
  sliceFp = NULL;
  slicePos = float1D(nlayer,"slicePos");

  if (muls->sliceThickness == 0)
    (*muls).cz[0] = c/(float_tt)nlayer;
  else
    (*muls).cz[0] = muls->sliceThickness;
  slicePos[0] = (*muls).czOffset;
  /*
************************************************************/


  for (i=1;i<nlayer;i++) {
    if (sliceFp == NULL) (*muls).cz[i] = (*muls).cz[0];
    /* don't need to all be the same, yes they do for fast 3D-FFT method! */
    else {
      fgets(buf,BUF_LEN,sliceFp);
      (*muls).cz[i] = atof(buf);
    }
    slicePos[i] = slicePos[i-1]+(*muls).cz[i-1]/2.0+(*muls).cz[i]/2.0;
  }

  // reset the potential to zero:
  memset((void *)&(muls->trans[0][0][0][0]),0,
         muls->slices*muls->potNx*muls->potNy*sizeof(complex_tt));

  nyAtBox = 2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionY);
  nxyAtBox = nyAtBox*(2*OVERSAMP_X*(int)ceil(muls->atomRadius/muls->resolutionX));
  nyAtBox2 = 2*nyAtBox;
  nxyAtBox2 = 2*nxyAtBox;
  sliceStep = 2*muls->potNx*muls->potNy;

  /****************************************************************
   * Loop through all the atoms and add their potential
   * to the slices:                                                                        
   ***************************************************************/

  time(&time0);
  
  for (std::vector<atom>::iterator atom = m_atoms.begin();atom!=m_atoms.end();atom++) 
    {
    // make sure we skip vacancies:
    while (atom->Znum == 0) atom++;
    size_t iatom=atom-m_atoms.begin();
    if ((muls->printLevel >= 4) && (muls->displayPotCalcInterval > 0) && ((iatom+1) % (muls->displayPotCalcInterval)) == 0) 
      {
        printf("Adding potential for atom %d (Z=%d, pos=[%.1f, %.1f, %.1f])\n",iatom+1,atom->Znum,
               atom->x, atom->y, atom->z);
      }
    // printf("c=%g, slice thickness=%g, slices=%d, %d\n",c,muls->sliceThickness,muls->slices,muls->displayPotCalcInterval);
    /*
     * c = the thickness of the current slab.
     *
     * if the z-position of this atom is outside the potential slab
     * we won't consider it and skip to the next
     */
    /* cellDiv = number of times that the big super cell is divided into
     * less big ones, yet often still bigger than a single unit cell
     * (for saving memory)
     * divCount = counter of how many such semi-super cells we already
     * passed through.
     * c = height in A of one such semi-super cell.
     * Since cellDiv is always >=1, and divCount starts at 0, the actual position
     * of this atom with the super cell is given by:
     */
    /* c*(float_tt)((*muls).cellDiv-divCount-1) will pick the right super-cell
     * division in the big super-cell
     * The z-offset 0.5*cz[0] will position atoms at z=0 into the middle of the first
     * slice.
     */
    atomZ = atom->z-c*(float_tt)(muls->cellDiv-divCount-1) + muls->czOffset
      -(0.5*muls->sliceThickness*(1-muls->centerSlices));
    // make sure that slices are centered for 2D and 3D differently:
    if (muls->potential3D==0)        atomZ += 0.5*muls->sliceThickness;
    else atomZ -= muls->sliceThickness;

    /* Now we need to find the first atom that contributes to this slice */
    /* Make use of the fact that we sorted the atoms in z */
    if ((*muls).nonPeriodZ) {
      if (((*muls).potential3D) && (atomZ -(*muls).atomRadius > c)) break;        
      if (((*muls).potential3D==0) && (atomZ > c)) break;                
      do {
        // printf("z: %g c: %g\n",atomZ,c);
        if (((*muls).potential3D) && (atomZ+(*muls).atomRadius+muls->sliceThickness >=0)) break;
        if (((*muls).potential3D==0) && (atomZ >=0)) break;
        // atomZ is the z of the next atom minus ...
        atomZ = (++atom)->z-c*(float_tt)(muls->cellDiv-divCount-1) + muls->czOffset
                -(0.5*muls->sliceThickness*(1-muls->centerSlices));
        if (muls->potential3D==0)        atomZ += 0.5*muls->sliceThickness;
        else atomZ -= muls->sliceThickness;
      }
    }
    /* atom coordinates in cartesian coords
     * The x- and y-position will be offset by the starting point
     * of the actually needed array of projected potential
     */
    atomX = atom->x -(*muls).potOffsetX;
    atomY = atom->y -(*muls).potOffsetY;

    /* so far we need periodicity in z-direction.
     * This requirement can later be removed, if we
     * use some sort of residue slice which will contain the
     * proj. potential that we need to add to the first slice of
     * the next stack of slices
     */

    /*************************************************************
     * real space potential lookup table summation
     ************************************************************/
    if (muls->fftpotential) {
            
    }
          
    /**************************************************************************
     * Newer, even faster method based on FFT of tabulated scattering factors
     **************************************************************************/
    else {        /* fftpotential */
            
    } /* end of if (fftpotential) */
    } /* for iatom =0 ... */
  time(&time1);
  if (muls->printLevel) 
      printf("%g sec used for real space potential calculation (%g sec per atom)\n",difftime(time1,time0),difftime(time1,time0)/iatom);
    else
      if (muls->printLevel) printf("%g sec used for real space potential calculation\n",difftime(time1,time0));


  /*************************************************/
  /* Save the potential slices                                         */
  
  if (muls->savePotential) {
    for (iz = 0;iz<nlayer;iz++){
      // find the maximum value of each layer:
      potVal = muls->trans[iz][0][0][0];
      for (ddx=potVal,ddy = potVal,ix=0;ix<muls->potNy*muls->potNx;potVal = muls->trans[iz][0][++ix][0]) {
        if (ddy<potVal) ddy = potVal;
        if (ddx>potVal) ddx = potVal;
      }
      
      WriteSlice(iz);
      
      if (muls->printLevel >= 3)
        printf("Saving (complex) potential layer %d to file %s (r: %g..%g)\n",iz,filename,ddx,ddy);

      params["Thickness"]=muls->sliceThickness;
      sprintf(buf,"Projected Potential (slice %d)",iz);
      imageIO->WriteComplexImage((void **)muls->trans[iz], params, std::string(buf));
    } // loop through all slices
  } /* end of if savePotential ... */
  if (muls->saveTotalPotential) {
    if (tempPot == NULL) tempPot = float2D(muls->potNx,muls->potNy,"total projected potential");
    
    for (ix=0;ix<muls->potNx;ix++) for (iy=0;iy<muls->potNy;iy++) {
        tempPot[ix][iy] = 0;
        for (iz=0;iz<nlayer;iz++) tempPot[ix][iy] += muls->trans[iz][ix][iy][0];
      }

    for (ddx=tempPot[0][0],ddy = potVal,ix=0;ix<muls->potNy*muls->potNx;potVal = tempPot[0][++ix]) {
      if (ddy<potVal) ddy = potVal;
      if (ddx>potVal) ddx = potVal;
    }
    if (muls->printLevel >= 2)
      printf("Saving total projected potential to file %s (r: %g..%g)\n",filename,ddx,ddy);
    params["Thickness"]=muls->sliceThickness;
    sprintf(buf,"Projected Potential (sum of %d slices)",muls->slices);
    WriteProjectedPotential();
  }
} // end of make3DSlices


void AddAtomToSliceNonPeriodic()
{
  /* Warning: will assume constant slice thickness ! */
  /* do not round here: atomX=0..dx -> iAtomX=0 */
  /*
    iAtomX = (int)(atomX/dx);        
    if (atomX/dx < (float)iAtomX) iAtomX--; // in case iAtomX is negative
    iAtomY = (int)(atomY/dy);
    if (atomY/dy < (float)iAtomY) iAtomY--;
    iAtomZ = (int)(atomZ/(*muls).cz[0]);
    if (atomZ/(*muls).cz[0] < (float)iAtomZ) iAtomZ--;
  */
  iAtomX = (int)floor(atomX/dx);
  iAtomY = (int)floor(atomY/dy);
  iAtomZ = (int)floor(atomZ/muls->cz[0]);

  // printf("atomZ(%d)=%g(%d)\t",iatom,atomZ,iAtomZ);

  if (muls->displayPotCalcInterval > 0) {
    if ((muls->printLevel>=3) && ((iatom+1) % muls->displayPotCalcInterval == 0)) {
      printf("adding atom %d [%.3f %.3f %.3f (%.3f)], Z=%d\n",
             iatom+1,atomX+(*muls).potOffsetX,atomY+(*muls).potOffsetY,
             atoms[iatom].z,atomZ,atoms[iatom].Znum);
      /* (*muls).ax,(*muls).by,(*muls).potOffsetX,(*muls).potOffsetY); */
    }
  }

  for (iax = -iRadX;iax<=iRadX;iax++) {
    if ((*muls).nonPeriod) {
      if (iax+iAtomX < 0) {
        iax = -iAtomX;
        if (abs(iax)>iRadX) break;
      }
      if (iax+iAtomX >= nx)        break;
    }
    x = (double)(iAtomX+iax)*dx-atomX;
    ix = (iax+iAtomX+16*nx) % nx;        /* shift into the positive range */
    for (iay=-iRadY;iay<=iRadY;iay++) {
      if ((*muls).nonPeriod) {
        if (iay+iAtomY < 0) {
          iay = -iAtomY;
          if (abs(iay)>iRadY) break;
        }
        if (iay+iAtomY >= ny)        break;
      }
      y = (double)(iAtomY+iay)*dy-atomY;
      iy = (iay+iAtomY+16*ny) % ny;         /* shift into the positive range */
      r2sqr = x*x + y*y;
   if (r2sqr <= atomRadius2) 
     {
       AddAtomToPotential()
         }
}


/*------------------------ transmit() ------------------------*/
/*
transmit the wavefunction thru one layer 
(simply multiply wave by transmission function)

waver,i[ix][iy]  = real and imaginary parts of wavefunction
transr,i[ix][iy] = real and imag parts of transmission functions

nx, ny = size of array

on entrance waver,i and transr,i are in real space

only waver,i will be changed by this routine
*/
void CPotential::transmit(WavePtr wave, unsigned sliceIdx) {
	double wr, wi, tr, ti;

	complex_tt **w,**t;
	unsigned posx = wave->iPosX;
	unsigned posy = wave->iPosY;
	w = (complex_tt **)wave->wave;
	t = (complex_tt **)trans[sliceIdx];

	/*  trans += posx; */
	for(unsigned ix=0; ix<m_nx; ix++) for(unsigned iy=0; iy<m_ny; iy++) {
		wr = w[ix][iy][0];
		wi = w[ix][iy][1];
		tr = t[ix+posx][iy+posy][0];
		ti = t[ix+posx][iy+posy][1];
		w[ix][iy][0] = wr*tr - wi*ti;
		w[ix][iy][1] = wr*ti + wi*tr;
	} /* end for(iy.. ix .) */
} /* end transmit() */