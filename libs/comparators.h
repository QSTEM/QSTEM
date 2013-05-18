#ifndef COMPARATORS_H
#define COMPARATORS_H

#include "stemtypes_fftw3.h"

struct atomCompareZnum {
	bool operator ()( atom const &atom1,atom const &atom2 ) const {
		return atom1.Znum<atom2.Znum;
	}
};

struct atomCompareZYX {
bool operator ()( atom const &atom1, const atom &atom2 ) const {
		if (atom1.pos[2]>atom2.pos[2]) return true;
		else if (atom1.pos[2]<atom2.pos[2]) return false;
		if (atom1.pos[1] > atom2.pos[1]) return true;
		else if (atom1.pos[1] < atom2.pos[1]) return false;
		if (atom1.pos[0] > atom2.pos[0]) return true;
		else if (atom1.pos[0] < atom2.pos[0]) return false;
		return false;
	}
};

#endif