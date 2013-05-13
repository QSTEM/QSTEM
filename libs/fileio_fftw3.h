#ifndef FILEIO_H
#define FILEIO_H

#include "data_containers.h"
#include "stemtypes_fftw3.h"

std::vector<atom> readUnitCell(int *natom,char *fileName,MULS *muls,int handleVacancies);
void replicateUnitCell(int ncoord,int *natom,MULS *muls,std::vector<atom> atoms,int handleVacancies);
std::vector<atom> tiltBoxed(int ncoord,int *natom, MULS *muls,std::vector<atom> atoms,int handleVacancies);
int writePDB(std::vector<atom> atoms,int natoms,char *fileName,MULS *muls);
int writeCFG(std::vector<atom> atoms,int natoms,char *fileName,MULS *muls);
// write CFG file using atomic positions stored in pos, Z's in Znum and DW-factors in dw
// the unit cell is assumed to be cubic
int writeCFGFractCubic(double *pos,int *Znum,double *dw,int natoms,char *fileName,
		       double a,double b,double c);
// Commented 2013-05-11 MCS - unused ATM.
//int readCubicCFG(double **pos,double **dw, int **Znum, double *ax,double *by,double *cz,
		 //double ctiltx, double ctilty);

void writeSTEMinput(char* stemFile,char *cfgFile,MULS *muls);

/* Helper functions for above functions: */
size_t ReadLine( FILE* fpRead, char* cRead, int cMax, const char *mesg );
int getZNumber(char *element);
int readCFGCellParams(MULS *muls, double **Mm, char *fileName);
int readCSSRCellParams(MULS *muls, double **Mm, char *fileName);

void writeFrameWork(FILE *fp,superCellBox superCell);
void writeAmorphous(FILE *fp,superCellBox superCell,int nstart,int nstop);

float_tt gasdev(long *idum); 
float_tt ran(long *idum);
float_tt ran1(long *idum);
int atomCompareZYX(const void *atPtr1,const void *atPtr2);
int atomCompareZnum(const void *atPtr1,const void *atPtr2);
#endif /* FILEIO_H */
