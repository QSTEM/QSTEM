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


