#include "stdio.h"
#include <string.h>
#include "data_containers.h"
#include "splines.h"

WAVEFUNC::WAVEFUNC(int x, int y) :
detPosX(0),
detPosY(0),
iPosX(0),
iPosY(0),
thickness(0.0),
intIntensity(0),
nx(0),
ny(0)
{

	std::string waveFile;
	std::string waveFileBase = "mulswav";
	nx = x;
	ny = y;
	diffpat = QSfMat::Zero(nx, ny);
	avgArray = QSfMat::Zero(nx, ny);
	wave = QScMat::Zero(nx, ny);
  
#if FLOAT_PRECISION == 1
	fftPlanWaveForw = fftwf_plan_dft_2d(nx,ny,(fftwf_complex*)wave.data(),(fftwf_complex*)wave.data(),FFTW_FORWARD, FFTW_ESTIMATE);
	fftPlanWaveInv = fftwf_plan_dft_2d(nx,ny,(fftwf_complex*)wave.data(),(fftwf_complex*)wave.data(),FFTW_BACKWARD, FFTW_ESTIMATE);
#else
	fftPlanWaveForw = fftw_plan_dft_2d(nx,ny,(fftwf_complex*)wave.data(),(fftwf_complex*)wave.data(),FFTW_FORWARD,
		fftMeasureFlag);
	fftPlanWaveInv = fftw_plan_dft_2d(nx,ny,(fftwf_complex*)wave.data(),(fftwf_complex*)wave.data(),FFTW_BACKWARD,
		fftMeasureFlag);
#endif
    waveFile = waveFileBase+".img";
    fileout = waveFile;
    fileStart = "mulswav.img";
}

WAVEFUNC::WAVEFUNC( WAVEFUNC& other ) :
iPosX(other.iPosX),
iPosY(other.iPosY),
detPosX(other.detPosX),
detPosY(other.detPosY),
nx(other.nx),
ny(other.ny),
fileStart(other.fileStart),
fileout(other.fileout),
avgName(other.avgName),
thickness(0.0),
intIntensity(0)
{
	diffpat = QSfMat::Zero(nx, ny);
	avgArray = QSfMat::Zero(nx, ny);

    wave = QScMat::Zero(nx,ny);

#if FLOAT_PRECISION == 1
	fftPlanWaveForw = fftwf_plan_dft_2d(nx,ny,(fftwf_complex*)wave.data(),(fftwf_complex*)wave.data(),FFTW_FORWARD, FFTW_ESTIMATE);
	fftPlanWaveInv = fftwf_plan_dft_2d(nx,ny,(fftwf_complex*)wave.data(),(fftwf_complex*)wave.data(),FFTW_BACKWARD, FFTW_ESTIMATE);
#else
	fftPlanWaveForw = fftw_plan_dft_2d(nx,ny,(fftwf_complex*)wave.data(),(fftwf_complex*)wave.data(),FFTW_FORWARD,
		fftMeasureFlag);
	fftPlanWaveInv = fftw_plan_dft_2d(nx,ny,(fftwf_complex*)wave.data(),(fftwf_complex*)wave.data(),FFTW_BACKWARD,
		fftMeasureFlag);
#endif
}

// reset the wave's thickness
void WAVEFUNC::ZeroWave(void)
{

}

MULS::MULS () 
{
	initMuls();
};

MULS::MULS (int slices)  
{
	// muls.areaAIS = 1.0;

	/* make multislice read the inout files and assign transr and transi: */
	//muls.trans = NULL;
	//muls.cz = NULL;  // (float_tt *)malloc(muls.slices*sizeof(float_tt));

	//muls.kx = NULL;
	//muls.kx2= NULL;
	//muls.ky = NULL;
	//muls.ky2= NULL;

	/****************************************************/
	/* copied from slicecell.c                          */
	//muls.pendelloesung = NULL;
}

void MULS::initMuls(int slices)
{
	initMuls();
	int sCount;
	/* general setup: */
	// cin2 should be a vector of strings - with each element being "a" followed by the slice index.  
	//     Is a a character (integer), and this is some kind of offset based on that?
	for (sCount =0;sCount<slices;sCount++)
		cin2[sCount] = 'a'+sCount;
	for (sCount = slices;sCount < NCINMAX;sCount++)
		cin2[sCount] = 0;

}

void MULS::initMuls()
{
	atomRadius= 5.0; /* radius in A for making the potential boxes */
	saveFlag=0;
	sigmaf=0; dfdelt=0; acmax=0; acmin=0; aobj=0; Cs=0; aAIS=0;
	// Tomography stuff - tomoCount = 0 indicates no Tomo simulation.
	tomoTilt=0; tomoStart=0; tomoStep=0; tomoCount=0;
	onlyFresnel=0;
	czOffset=0; /* defines the offset for the first slice in 
						fractional coordinates        */
	normHolog=0;
	gaussianProp=0;

	sparam = QSfVec::Zero(NPARAM);
}


//std::vector<std::string> ElTable::elements;
std::map<int, std::string> ElTable::ZToSymbol;
std::map<std::string, int> ElTable::SymbolToZ;
bool ElTable::isFilled(false);

ElTable::ElTable()
{
	Fill();
}

void ElTable::Fill()
{
	using namespace std;
	// offset by one so that elements can be addressed by their Z directly.
	ZToSymbol[0]="DUMMY"; SymbolToZ["DUMMY"]=0;
	ZToSymbol[1]="H";   SymbolToZ["H"]=1;
	ZToSymbol[2]="He";  SymbolToZ["He"]=2;
	ZToSymbol[3]="Li";  SymbolToZ["Li"]=3;
	ZToSymbol[4]="Be";  SymbolToZ["Be"]=4;
	ZToSymbol[5]="B";   SymbolToZ["B"]=5;
	ZToSymbol[6]="N";   SymbolToZ["N"]=6;
	ZToSymbol[7]="C";   SymbolToZ["C"]=7;
	ZToSymbol[8]="O";   SymbolToZ["O"]=8;
	ZToSymbol[9]="F";   SymbolToZ["F"]=9;
	ZToSymbol[10]="Ne"; SymbolToZ["Ne"]=10;
	ZToSymbol[11]="Na"; SymbolToZ["Na"]=11;
	ZToSymbol[12]="Mg"; SymbolToZ["Mg"]=12;
	ZToSymbol[13]="Al"; SymbolToZ["Al"]=13;
	ZToSymbol[14]="Si"; SymbolToZ["Si"]=14;
	ZToSymbol[15]="P";  SymbolToZ["P"]=15;
	ZToSymbol[16]="S";  SymbolToZ["S"]=16;
	ZToSymbol[17]="Cl"; SymbolToZ["Cl"]=17;
	ZToSymbol[18]="Ar"; SymbolToZ["Ar"]=18;
	ZToSymbol[19]="K";  SymbolToZ["K"]=19;
	ZToSymbol[20]="Ca"; SymbolToZ["Ca"]=20;
	ZToSymbol[21]="Sc"; SymbolToZ["Sc"]=21;
	ZToSymbol[22]="Ti"; SymbolToZ["Ti"]=22;
	ZToSymbol[23]="V";  SymbolToZ["V"]=23;
	ZToSymbol[24]="Cr"; SymbolToZ["Cr"]=24;
	ZToSymbol[25]="Mn"; SymbolToZ["Mn"]=25;
	ZToSymbol[26]="Fe"; SymbolToZ["Fe"]=26;
	ZToSymbol[27]="Co"; SymbolToZ["Co"]=27;
	ZToSymbol[28]="Ni"; SymbolToZ["Ni"]=28;
	ZToSymbol[29]="Cu"; SymbolToZ["Cu"]=29;
	ZToSymbol[30]="Zn"; SymbolToZ["Zn"]=30;
	ZToSymbol[31]="Ga"; SymbolToZ["Ga"]=31;
	ZToSymbol[32]="Ge"; SymbolToZ["Ge"]=32;
	ZToSymbol[33]="As"; SymbolToZ["As"]=33;
	ZToSymbol[34]="Se"; SymbolToZ["Se"]=34;
	ZToSymbol[35]="Br"; SymbolToZ["Br"]=35;
	ZToSymbol[36]="Kr"; SymbolToZ["Kr"]=36;
	ZToSymbol[37]="Rb"; SymbolToZ["Rb"]=37;
	ZToSymbol[38]="Sr"; SymbolToZ["Sr"]=38;
	ZToSymbol[39]="Y";  SymbolToZ["Y"]=39;
	ZToSymbol[40]="Zr"; SymbolToZ["Zr"]=40;
	ZToSymbol[41]="Nb"; SymbolToZ["Nb"]=41;
	ZToSymbol[42]="Mo"; SymbolToZ["Mo"]=42;
	ZToSymbol[43]="Tc"; SymbolToZ["Tc"]=43;
	ZToSymbol[44]="Ru"; SymbolToZ["Ru"]=44;
	ZToSymbol[45]="Rh"; SymbolToZ["Rh"]=45;
	ZToSymbol[46]="Pd"; SymbolToZ["Pd"]=46;
	ZToSymbol[47]="Ag"; SymbolToZ["Ag"]=47;
	ZToSymbol[48]="Cd"; SymbolToZ["Cd"]=48;
	ZToSymbol[49]="In"; SymbolToZ["In"]=49;
	ZToSymbol[50]="Sn"; SymbolToZ["Sn"]=50;
	ZToSymbol[51]="Sb"; SymbolToZ["Sb"]=51;
	ZToSymbol[52]="Te"; SymbolToZ["Te"]=52;
	ZToSymbol[53]="I";  SymbolToZ["I"]=53;
	ZToSymbol[54]="Xe"; SymbolToZ["Xe"]=54;
	ZToSymbol[55]="Cs"; SymbolToZ["Cs"]=55;
	ZToSymbol[56]="Ba"; SymbolToZ["Ba"]=56;
	ZToSymbol[57]="La"; SymbolToZ["La"]=57;
	ZToSymbol[58]="Ce"; SymbolToZ["Ce"]=58;
	ZToSymbol[59]="Pr"; SymbolToZ["Pr"]=59;
	ZToSymbol[60]="Nd"; SymbolToZ["Nd"]=60;
	ZToSymbol[61]="Pm"; SymbolToZ["Pm"]=61;
	ZToSymbol[62]="Sm"; SymbolToZ["Sm"]=62;
	ZToSymbol[63]="Eu"; SymbolToZ["Eu"]=63;
	ZToSymbol[64]="Gd"; SymbolToZ["Gd"]=64;
	ZToSymbol[65]="Tb"; SymbolToZ["Tb"]=65;
	ZToSymbol[66]="Dy"; SymbolToZ["Dy"]=66;
	ZToSymbol[67]="Ho"; SymbolToZ["Ho"]=67;
	ZToSymbol[68]="Er"; SymbolToZ["Er"]=68;
	ZToSymbol[69]="Tm"; SymbolToZ["Tm"]=69;
	ZToSymbol[70]="Yb"; SymbolToZ["Yb"]=70;
	ZToSymbol[71]="Lu"; SymbolToZ["Lu"]=71;
	ZToSymbol[72]="Hf"; SymbolToZ["Hf"]=72;
	ZToSymbol[73]="Ta"; SymbolToZ["Ta"]=73;
	ZToSymbol[74]="W";  SymbolToZ["W"]=74;
	ZToSymbol[75]="Re"; SymbolToZ["Re"]=75;
	ZToSymbol[76]="Os"; SymbolToZ["Os"]=76;
	ZToSymbol[77]="Ir"; SymbolToZ["Ir"]=77;
	ZToSymbol[78]="Pt"; SymbolToZ["Pt"]=78;
	ZToSymbol[79]="Au"; SymbolToZ["Au"]=79;
	ZToSymbol[80]="Hg"; SymbolToZ["Hg"]=80;
	ZToSymbol[81]="Tl"; SymbolToZ["Tl"]=81;
	ZToSymbol[82]="Pb"; SymbolToZ["Pb"]=82;
	ZToSymbol[83]="Bi"; SymbolToZ["Bi"]=83;
	ZToSymbol[84]="Po"; SymbolToZ["Po"]=84;
	ZToSymbol[85]="At"; SymbolToZ["At"]=85;
	ZToSymbol[86]="Rn"; SymbolToZ["Rn"]=86;
	ZToSymbol[87]="Fr"; SymbolToZ["Fr"]=87;
	ZToSymbol[88]="Ra"; SymbolToZ["Ra"]=88;
	ZToSymbol[89]="Ac"; SymbolToZ["Ac"]=89;
	ZToSymbol[90]="Th"; SymbolToZ["Th"]=90;
	ZToSymbol[91]="Pa"; SymbolToZ["Pa"]=91;
	ZToSymbol[92]="U";  SymbolToZ["U"]=92;
	ZToSymbol[93]="Np"; SymbolToZ["Np"]=93;
	ZToSymbol[94]="Pu"; SymbolToZ["Pu"]=94;
	ZToSymbol[95]="Am"; SymbolToZ["Am"]=95;
	ZToSymbol[96]="Cm"; SymbolToZ["Cm"]=96;
	ZToSymbol[97]="Bk"; SymbolToZ["Bk"]=97;
	ZToSymbol[98]="Cf"; SymbolToZ["Cf"]=98;
	ZToSymbol[99]="Es"; SymbolToZ["Es"]=99;
	ZToSymbol[100]="Fm"; SymbolToZ["Fm"]=100;
	ZToSymbol[101]="Md"; SymbolToZ["Md"]=101;
	ZToSymbol[102]="No"; SymbolToZ["No"]=102;
	ZToSymbol[103]="Lr"; SymbolToZ["Lr"]=103;

	isFilled=true;				   
}								   

std::string ElTable::GetSymbol(int element)
{
	std::map<int,std::string>::iterator it;
	if (!isFilled)
		Fill();
	it=ZToSymbol.find(element);
	return it->second;
}

int ElTable::GetZ(std::string element)
{
	std::map<std::string,int>::iterator it;
	if (!isFilled)
		Fill();
	std::string ext = element;
	ext.erase(std::remove_if(ext.begin(), ext.end(), isspace), ext.end());
	it=SymbolToZ.find(ext);
	return it->second;
}