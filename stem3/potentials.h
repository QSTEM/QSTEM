#ifndef POTENTIALS_H
#define POTENTIALS_H

#include "stemtypes_fftw3.h"
#include "splines.h"
#include "data_containers.h"
#include <vector>

// Class for calculating and managing potentials
class Potential
{
	MULS *muls;
	// a set of logarithmic r values
	static std::vector<float_tt> m_splinr;
	std::vector<int> m_knownZvalues;
	std::vector<AkimaSpline<float_tt,float_tt> > m_potentialSplines;
	std::vector<AkimaSpline<float_tt,float_tt> > m_offsetSplines;
	QSVecOfcMat m_atPot;
	bool m_tdsFlag;
	int m_scatFlag;
	int m_Nr, m_Nz_lut, m_nzSub;

	void GenerateSplineEntry(int Z, float_tt charge=0);
	void CreateAtPot(int Znum, float_tt B, float_tt charge=0);

public:
	Potential(MULS *muls);
	QScMat GetPot(int Z, float_tt charge=0);
	float_tt SplineLookUp(int Z, float_tt R, float_tt charge=0);
	int GetNr() {return m_Nr;}
	int GetNzSub() {return m_nzSub;}
	int GetNzLUT() {return m_Nz_lut;}
};

#endif