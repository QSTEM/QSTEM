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

#ifndef EXPERIMENT_TEM_H
#define EXPERIMENT_TEM_H

#include "base.hpp"

class CExperimentTEM : public CExperimentBase
{
public:
  CExperimentTEM(const ConfigReaderPtr &configReader);
  void Run();
  void DisplayParams();
private:
  void SaveImages();
  float_tt **m_pendelloesung;              /* pendelloesung plot */
  std::vector<int> m_hbeams,m_kbeams;	/* arrays to hold recorded 
					   beam indicies */
  bool m_lbeams;				/* flag indicating, whether 
					   to record beams */	
  char m_filebeam[512];		 	/* file, that beams get recorded in */
  unsigned m_nbout;				/* number of recorded beams */
};

#endif










