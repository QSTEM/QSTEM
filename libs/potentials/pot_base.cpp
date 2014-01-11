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

const std::string kPotFileName = "potslice";
const int BUF_LEN = 256;

CPotential::CPotential(const ConfigReaderPtr &configReader)
{
  configReader->ReadProbeArraySize(m_nx, m_ny);
  configReader->ReadResolution(m_dx, m_dy);
  configReader->ReadVoltage(m_v0);
  configReader->ReadPotentialOutputParameters(m_savePotential, m_saveProjectedPotential, m_plotPotential);
  configReader->ReadAtomRadius(m_atomRadius);
  configReader->ReadSliceParameters(m_centerSlices, m_sliceThickness, m_nslices, m_outputInterval, m_zOffset);
  configReader->ReadSliceOffset(m_offsetX, m_offsetY);
  Initialize();
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
  m_boxNx = (int)(m_atomRadius/m_ddx+2.0);  
  m_boxNy = (int)(m_atomRadius/m_ddy+2.0);  

  m_c = m_sliceThickness * m_nslices;
  m_dr = m_dx/OVERSAMPLING; // define step width in which radial V(r,z) is defined
  m_iRadX = (int)ceil(m_atomRadius/m_dx);
  m_iRadY = (int)ceil(m_atomRadius/m_dy);
  m_iRadZ = (int)ceil(m_atomRadius/m_sliceThickness);
  m_iRad2 = m_iRadX*m_iRadX+m_iRadY*m_iRadY;
  m_atomRadius2 = m_atomRadius * m_atomRadius;

  m_cz.resize(m_nslices);
  m_slicePos.resize(m_nslices);

  if (m_sliceThickness == 0)
    m_cz[0] = m_c/(float_tt)m_nslices;
  else
    m_cz[0] = m_sliceThickness;

  m_slicePos[0] = m_zOffset;
}

void CPotential::DisplayParams()
{
  printf("*****************************************************\n");
  printf("********* Potential Parameters **********************\n");
  printf("*****************************************************\n");
  printf("* Print level:          %d\n",m_printLevel);
  printf("* Save level:           %d\n",m_saveLevel);
  if (m_savePotential)
    printf("* Potential file name:  %s\n",m_fileBase.c_str());
  printf("* Model Sampling:  %g x %g x %g A\n", m_dx,m_dy,m_sliceThickness);

  printf("* Pot. array offset:    (%g,%g,%g)A\n",m_offsetX,m_offsetY,m_zOffset);
  printf("* Potential periodic:   (x,y): %s, z: %s\n",
         (m_periodicXY) ? "yes" : "no",(m_periodicZ) ? "yes" : "no");
  if (k_fftMeasureFlag == FFTW_MEASURE)
    printf("* Potential array:      %d x %d (optimized)\n",m_nx,m_ny);
  else
    printf("* Potential array:      %d x %d (estimated)\n",m_nx,m_ny);
  printf("*                       %g x %gA\n",m_nx*m_dx,m_ny*m_dy);
  printf("* Scattering factors:   %d\n",m_scatFactor);

    printf("* Slices per division:  %d (%gA thick slices [%scentered])\n",
         m_nslices,m_sliceThickness,(m_centerSlices) ? "" : "not ");
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


void CPotential::ReadPotential(std::string &fileName, unsigned subSlabIdx)
{
  /*************************************************************************
   * read the potential that has been created externally!
   */
  unsigned slice_idx=0;
  for (unsigned i=(subSlabIdx+1)*m_nslices-1;i>=(subSlabIdx)*m_nslices;i--,slice_idx++) 
    {
      ReadSlice(fileName, m_trans[slice_idx], i);
    }
    return;
}

void CPotential::ReadSlice(std::string &fileName, complex_tt **, unsigned idx)
{
}

void CPotential::SliceSetup()
{
  FILE *sliceFp;
  char buf[BUF_LEN];
  /**************************************************************
   *        setup the slices with their start and end positions
   *        then loop through all the atoms and add their potential to
   *        the slice that their potential reaches into (up to RMAX)
   *************************************************************/
       
  for (unsigned i=1;i<m_nslices;i++) {
    if (sliceFp == NULL) m_cz[i] = m_cz[0];
    /* need to all be the same for fast 3D-FFT method, otherwise OK to be different */
    else {
      fgets(buf,BUF_LEN,sliceFp);
      m_cz[i] = atof(buf);
    }
    m_slicePos[i] = m_slicePos[i-1]+m_cz[i-1]/2.0+m_cz[i]/2.0;
  }
}

void CPotential::CenterAtomZ(std::vector<atom>::iterator &atom, float_tt &z)
{
  
  /*
   * Since cellDiv is always >=1, and divCount starts at 0, the actual position
   * of this atom with the super cell is given by:
   */
  /* c*(float_tt)((*muls).cellDiv-divCount-1) will pick the right super-cell
   * division in the big super-cell
   * The z-offset 0.5*cz[0] will position atoms at z=0 into the middle of the first
   * slice.
   */
    z = atom->z-m_c*(float_tt)(m_cellDiv-m_divCount-1) + m_zOffset
      -(0.5*m_sliceThickness*(1-(int)m_centerSlices));
}

/*****************************************************
* This function will create a 3D potential from whole
* unit cell, slice it, and make transr/i and propr/i
* Call this function with center = NULL, if you don't
* want the array to be shifted.
****************************************************/
void CPotential::MakeSlices(int nlayer,char *fileName, atom *center) 
{
  time_t time0, time1;
  /* return, if there is nothing to do */
  if (nlayer <1)
    return;

  /* we need to keep track of which subdivision of the unit cell we are in
   * If the cell is not subdivided, then m_cellDiv-1 = 0.
   */
  if ((m_divCount == 0) || (m_equalDivs))
    m_divCount = m_cellDiv;
  m_divCount--;

  /* we only want to reread and shake the atoms, if we have finished the
   * current unit cell.
   */
  if (m_divCount == m_cellDiv-1) {
    ReadAtoms();
  } /* end of if divCount==cellDiv-1 ... */

  SliceSetup();

  // reset the potential to zero:
  memset((void *)&(m_trans[0][0][0][0]),0,
         m_nslices*m_nx*m_ny*sizeof(complex_tt));

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
    if ((m_printLevel >= 4) && (m_displayPotCalcInterval > 0) && ((iatom+1) % (m_displayPotCalcInterval)) == 0) 
      {
        printf("Adding potential for atom %d (Z=%d, pos=[%.1f, %.1f, %.1f])\n",iatom+1,atom->Znum,
               atom->x, atom->y, atom->z);
      }
    /* atom coordinates in cartesian coords
     * The x- and y-position will be offset by the starting point
     * of the actually needed array of projected potential
     */
    float_tt atomX = atom->x - m_offsetX;
    float_tt atomY = atom->y - m_offsetY;
    float_tt atomZ;
    // make sure that slices are centered for 2D and 3D differently:
    CenterAtomZ(atom, atomZ);

    AddAtomToSlices(atom, atomX, atomY, atomZ);

    } /* for iatom =0 ... */
  time(&time1);
  if (m_printLevel) 
    printf("%g sec used for real space potential calculation (%g sec per atom)\n",difftime(time1,time0),difftime(time1,time0)/m_atoms.size());
    else
      if (m_printLevel) printf("%g sec used for real space potential calculation\n",difftime(time1,time0));


  /*************************************************/
  /* Save the potential slices                                         */
  
  if (m_savePotential) {
    for (unsigned iz = 0;iz<nlayer;iz++){
      // find the maximum value of each layer:
      float_tt potVal = m_trans[iz][0][0][0];
      float_tt ddx = potVal;
      float_tt ddy = potVal;
      for (unsigned ix=0;ix<m_ny*m_nx;potVal = m_trans[iz][0][++ix][0]) {
        if (ddy<potVal) ddy = potVal;
        if (ddx>potVal) ddx = potVal;
      }
            
      if (m_printLevel >= 3)
        printf("Saving (complex) potential layer %d to file (r: %g..%g)\n",iz,ddx,ddy);

      WriteSlice(iz);

      //imageIO->WriteComplexImage((void **)m_trans[iz], params, std::string(buf));
    } // loop through all slices
  } /* end of if savePotential ... */
  if (m_saveProjectedPotential) {
    WriteProjectedPotential();
  }
} // end of make3DSlices


void CPotential::AddAtomRealSpace(std::vector<atom>::iterator &atom, 
                                  float_tt atomX, float_tt atomY, float_tt atomZ)
{
  unsigned iatom = atom-m_atoms.begin();

  CenterAtomZ(atom, atomZ);
  
  /* Warning: will assume constant slice thickness ! */
  /* do not round here: atomX=0..dx -> iAtomX=0 */
  unsigned iAtomX = (int)floor(atomX/m_dx);
  unsigned iAtomY = (int)floor(atomY/m_dy);
  unsigned iAtomZ = (int)floor(atomZ/m_cz[0]);

  if (m_displayPotCalcInterval > 0) {
    if ((m_printLevel>=3) && ((iatom+1) % m_displayPotCalcInterval == 0)) {
      printf("adding atom %d [%.3f %.3f %.3f (%.3f)], Z=%d\n",
             iatom+1,atomX+m_offsetX,atomY+m_offsetY,
             atom->z,atomZ,atom->Znum);
    }
  }

  for (int iax = -m_iRadX;iax<=m_iRadX;iax++) {
    if (!m_periodicXY) {
      if (iax+iAtomX < 0) {
        iax = -iAtomX;
        if (abs(iax)>m_iRadX) break;
      }
      if (iax+iAtomX >= m_nx) break;
    }
    float_tt x = (iAtomX+iax)*m_dx-atomX;
    unsigned ix = (iax+iAtomX+16*m_nx) % m_nx;        /* shift into the positive range */
    for (int iay=-m_iRadY;iay<=m_iRadY;iay++) {
      if (!m_periodicXY) {
        if (iay+iAtomY < 0) {
          iay = -iAtomY;
          if (abs(iay)>m_iRadY) break;
        }
        if (iay+iAtomY >= m_ny)        break;
      }
      float_tt y = (iAtomY+iay)*m_dy-atomY;
      unsigned iy = (iay+iAtomY+16*m_ny) % m_ny;         /* shift into the positive range */
      float_tt r2sqr = x*x + y*y;
      if (r2sqr <= m_atomRadius2) 
        {
          // This (virtual) method is meant to be implemented by subclasses, 
          //    for specific functionality varying by dimensionality.
          _AddAtomRealSpace(atom, x, ix, y, iy, atomZ, iAtomZ);
        }
    }
  }
}

void CPotential::WriteSlice(unsigned idx)
{
    char buf[255];
    std::map<std::string, double> params;
    params["Thickness"]=m_sliceThickness;
    sprintf(buf,"Projected Potential (slice %d)",idx);
    std::string comment = buf;
    m_imageIO->WriteComplexImage((void **)m_trans[idx], kPotFileName, params, comment);
}

void CPotential::WriteProjectedPotential()
{
  std::map<std::string, double> params;
  char buf[255];
  float_tt **tempPot = float2D(m_nx,m_ny,"total projected potential");
  float_tt potVal=0;
    
  for (unsigned ix=0;ix<m_nx;ix++) for (unsigned iy=0;iy<m_ny;iy++) {
      tempPot[ix][iy] = 0;
      for (unsigned iz=0;iz<m_nslices;iz++) tempPot[ix][iy] += m_trans[iz][ix][iy][0];
    }

  float_tt ddx=tempPot[0][0], ddy=potVal;
  for (unsigned ix=0;ix<m_ny*m_nx;potVal = tempPot[0][++ix]) {
    if (ddy<potVal) ddy = potVal;
    if (ddx>potVal) ddx = potVal;
  }
  if (m_printLevel >= 2)
    printf("Saving total projected potential to file (r: %g..%g)\n",ddx,ddy);
  params["Thickness"]=m_sliceThickness;
  sprintf(buf,"Projected Potential (sum of %d slices)",m_nslices);
  std::string comment = buf;
  std::string fileName = "ProjectedPot";
  m_imageIO->WriteRealImage((void **)tempPot, fileName, params, comment);
}

/*
// TODO: this was taken from stemutils.  It seems to be used only in customslice, which then isn't used anywhere.
float_tt CPotential::sfLUT(float_tt s,int atKind)
{
   int i;
   double sf;
   static double *splinx=NULL;
   static double **spliny=NULL;
   static double **splinb=NULL;
   static double **splinc=NULL;
   static double **splind=NULL;
   static int sfSize = 0;
   static int atKinds = 0;
   static double maxK = 0;

   if(splinx == NULL) {
     // splinx = s;
     // spliny = sfC;
     sfSize = m_sfNk;
     splinx = m_sfkArray;
     spliny = m_sfTable;
     atKinds = m_atoms.size();
     splinb = double2D(atKinds,sfSize, "splinb" );
     splinc = double2D(atKinds,sfSize, "splinc" );
     splind = double2D(atKinds,sfSize, "splind" );
     maxK = splinx[sfSize-1];

     for (i=0;i<atKinds;i++)
       splinh(splinx,spliny[i],splinb[i],splinc[i],splind[i],sfSize);
   }
   
   
   // now that everything is set up find the
   //   scattering factor by interpolation in the table 
   
   if (s > maxK) return 0.0;     
   if (atKind < atKinds) {
     sf = seval(splinx,spliny[atKind],splinb[atKind],splinc[atKind],splind[atKind],sfSize,s);
     if (sf < 0) return 0.0;
     return(sf);
   }
   printf("sfLUT: invalid atom kind (%d) - exit!\n",atKind);
   exit(0);
}  // end sfLUT()
*/


/*------------------ splinh() -----------------------------*/
/*
        fit a quasi-Hermite cubic spline
        
        [1] Spline fit as in H.Akima, J. ACM 17(1970)p.589-602
                'A New Method of Interpolation and Smooth
                Curve Fitting Based on Local Procedures'

        [2] H.Akima, Comm. ACM, 15(1972)p.914-918

        E. Kirkland 4-JUL-85
        changed zero test to be a small nonzero number 8-jul-85 ejk
        converted to C 24-jun-1995 ejk

        The inputs are:
                x[n] = array of x values in ascending order, each X(I) must
                        be unique
                y[n] = array of y values corresponding to X(N)
                n = number of data points must be 2 or greater

        The outputs are (with z=x-x(i)):
                b[n] = array of spline coeficients for (x-x[i])
                c[n] = array of spline coeficients for (x-x[i])**2
                d[n] = array of spline coeficients for (x-x[i])**3
                ( x[i] <= x <= x[i+1] )
        To interpolate y(x) = yi + bi*z + c*z*z + d*z*z*z

        The coeficients b[i], c[i], d[i] refer to the x[i] to x[i+1]
        interval. NOTE that the last set of coefficients,
        b[n-1], c[n-1], d[n-1] are meaningless.
*/
void CPotential::splinh( float_tt x[], float_tt y[],
         float_tt b[], float_tt c[], float_tt d[], int n)
{
#define SMALL 1.0e-25

  int i, nm1, nm4;
  float_tt m1, m2, m3, m4, m5, t1, t2, m54, m43, m32, m21, x43;
  
  if( n < 4) return;
  
  /* Do the first end point (special case),
and get starting values */
  
  m5 = ( y[3] - y[2] ) / ( x[3] - x[2] );        /* mx = slope at pt x */
  m4 = ( y[2] - y[1] ) / ( x[2] - x[1] );
  m3 = ( y[1] - y[0] ) / ( x[1] - x[0] );
  
  m2 = m3 + m3 - m4;        /* eq. (9) of reference [1] */
  m1 = m2 + m2 - m3;
  
  m54 = fabs( m5 - m4);
  m43 = fabs( m4 - m3);
  m32 = fabs( m3 - m2);
  m21 = fabs( m2 - m1);
  
  if ( (m43+m21) > SMALL )
    t1 = ( m43*m2 + m21*m3 ) / ( m43 + m21 );
  else
    t1 = 0.5 * ( m2 + m3 );
  
  /* Do everything up to the last end points */
  
  nm1 = n-1;
  nm4 = n-4;
  
  for( i=0; i<nm1; i++) {
    
    if( (m54+m32) > SMALL )
      t2= (m54*m3 + m32*m4) / (m54 + m32);
    else
      t2 = 0.5* ( m3 + m4 );
    
    x43 = x[i+1] - x[i];
    b[i] = t1;
    c[i] = ( 3.0*m3 - t1 - t1 - t2 ) /x43;
    d[i] = ( t1 + t2 - m3 - m3 ) / ( x43*x43 );
    
    m1 = m2;
    m2 = m3;
    m3 = m4;
    m4 = m5;
    if( i < nm4 ) {
      m5 = ( y[i+4] - y[i+3] ) / ( x[i+4] - x[i+3] );
    } else {
      m5 = m4 + m4 - m3;
    }
    
    m21 = m32;
    m32 = m43;
    m43 = m54;
    m54 = fabs( m5 - m4 );
    t1 = t2;
  }
  
  return;
  
} /* end splinh() */

/*----------------------- seval() ----------------------*/
/*
        Interpolate from cubic spline coefficients

        E. Kirkland 4-JUL-85
        modified to do a binary search for efficiency 13-Oct-1994 ejk
        converted to C 26-jun-1995 ejk
        fixed problem on end-of-range 16-July-1995 ejk

        The inputs are:
                x[n] = array of x values in ascending order, each x[i] must
                        be unique
                y[n] = array of y values corresponding to x[n]
                b[n] = array of spline coeficients for (x-x[i])
                c[n] = array of spline coeficients for (x-x[i])**2
                d[n] = array of spline coeficients for (x-x[i])**3
                n = number of data points
                x0 = the x value to interpolate at
                (x[i] <= x <= x[i+1]) and all inputs remain unchanged

        The value returned is the interpolated y value.

        The coeficients b[i], c[i], d[i] refer to the x[i] to x[i+1]
        interval. NOTE that the last set of coefficients,
        b[n-1], c[n-1], d[n-1] are meaningless.
*/
float_tt CPotential::seval( float_tt *x, float_tt *y, float_tt *b, float_tt *c,
         float_tt *d, int n, float_tt x0 )
{
  int i, j, k;
  float_tt z, seval1;
  
  /* exit if x0 is outside the spline range */
  if( x0 <= x[0] ) i = 0;
  else if( x0 >= x[n-2] ) i = n-2;
  else {
    i = 0;
    j = n;
    do{ k = ( i + j ) / 2 ;
    if( x0 < x[k] ) j = k;
    else if( x0 >= x[k] ) i = k;
    } while ( (j-i) > 1 );
  }
  
  z = x0 - x[i];
  seval1 = y[i] + ( b[i] + ( c[i] + d[i] *z ) *z) *z;
  
  return( seval1 );
  
} /* end seval() */
