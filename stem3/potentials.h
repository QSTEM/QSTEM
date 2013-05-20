// provide array variable names 
// - scatPar[N_ELEM][N_SF] and
// - scatParOffs[N_ELEM][N_SF]
// and also define N_SF and N_ELEM:
#include "scatfactsRez.h"
#include "stemtypes_fftw3.h"
#include "splines.h"
#include "data_containers.h"
#include <vector>

// Class for calculating and managing potentials
class Potential
{
	// a set of logarithmic r values
	static std::vector<float_tt> m_splinr;
	std::vector<int> m_knownZvalues;
	std::vector<AkimaSpline<float_tt,float_tt> > m_splines;
	QScMat atPot;
	bool m_tdsFlag;
	int m_scatFlag;

	void GenerateSplineEntry(int Z);

public:
	Potential();
	float_tt LookUp(int Z, float_tt R);
};