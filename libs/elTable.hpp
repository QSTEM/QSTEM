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

#ifndef ELTABLE_H
#define ELTABLE_H

char *elTable = {
	"H HeLiBeB C N O F NeNaMgAlSiP S Cl"
	"ArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBr"
	"KrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTe"
	"I XeCsBaLaCePrNdPmSmEuGdTbDyHoErTm"
	"YbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn"
	"FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLr"
};

inline int getZNumber(char *element) {
	char *elem;
#ifndef WIN32
	char widows = '\n';
	widows = (char)((int)widows + 3);

	// fprintf(stderr,"getZnumber : element = *%s* \n",element); 
	element[2] = '\0'; 
	if ((atoi(element+1) != 0) || (element[1] == '\n')|| (element[1] == '\0') || element[1] == widows){
		element[1] = ' ';
		// fprintf(stderr,"getZnumber : element = *%s*  conversion \n",element); 
	} 
#else
	element[2] = '\0';
	if ((atoi(element+1) != 0) || (element[1] == '\n')|| (element[1] == '\0'))
		element[1] = ' ';
#endif
	if ((elem = strstr(elTable,element)) == NULL)
		return 0;	
	else
		return (int)(elem-elTable)/2+1;

}

#endif
