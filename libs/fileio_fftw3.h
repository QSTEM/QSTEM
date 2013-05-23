#ifndef FILEIO_H
#define FILEIO_H

#include "data_containers.h"
#include "stemtypes_fftw3.h"

std::vector<atom> readUnitCell(int &natom,std::string fileName,MULS &muls,int handleVacancies);
int phononDisplacement(QSf3Vec &u, MULS &muls, int id, int icx, int icy,
                       int icz, int atomCount, float_tt dw, int maxAtom, int ZnumIndex);
void replicateUnitCell(int ncoord,int &natom,MULS &muls,std::vector<atom> atoms,int handleVacancies);
std::vector<atom> tiltBoxed(int ncoord,int &natom, MULS &muls,std::vector<atom> atoms,int handleVacancies);
int writePDB(std::vector<atom> atoms,int natoms,std::string fileName,MULS &muls);
int writeCFG(std::vector<atom> atoms,int natoms,std::string fileName,MULS *muls);
// write CFG file using atomic positions stored in pos, Z's in Znum and DW-factors in dw
// the unit cell is assumed to be cubic
int writeCFGFractCubic(float_tt *pos,int *Znum,float_tt *dw,int natoms,char *fileName,
		       float_tt a,float_tt b,float_tt c);
// Commented 2013-05-11 MCS - unused ATM.
//int readCubicCFG(float_tt **pos,float_tt **dw, int **Znum, float_tt *ax,float_tt *by,float_tt *cz,
		 //float_tt ctiltx, float_tt ctilty);

void writeSTEMinput(char* stemFile,char *cfgFile,MULS *muls);

/* Helper functions for above functions: */
int ReadLine( FILE* fpRead, char* cRead, int cMax, const char *mesg );
int readCFGCellParams(MULS &muls, QSf3Mat &Mm, std::string fileName);
int readCSSRCellParams(MULS &muls, QSf3Mat &Mm, std::string fileName);

void writeFrameWork(FILE *fp,superCellBox superCell);
void writeAmorphous(FILE *fp,superCellBox superCell,int nstart,int nstop);

float_tt gasdev(long *idum); 
float_tt ran(long *idum);
float_tt ran1(long *idum);

#endif /* FILEIO_H */
