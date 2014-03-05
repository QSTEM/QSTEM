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

#ifndef STRUCTURE_INTERFACE_H
#define STRUCTURE_INTERFACE_H

#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>
#include "stemtypes_fftw3.hpp"

#include <cstring>

namespace QSTEM
{

class IStructureReader;
typedef boost::shared_ptr<IStructureReader> StructureReaderPtr;
typedef StructureReaderPtr (*CreateStructureReaderFn)(const boost::filesystem::path &filename);

class IStructureReader
{
public:
  virtual int ReadCellParams(float_tt **Mm)=0;
  virtual int ReadAtoms(std::vector<atom> &atoms)=0;
};


class IStructureWriter;
typedef boost::shared_ptr<IStructureWriter> StructureWriterPtr;
typedef StructureWriterPtr (*CreateStructureWriterFn)(const boost::filesystem::path &filename, 
                                                      float_tt ax, float_tt by, float_tt cz);

class IStructureWriter
{
public:
  virtual int Write(std::vector<atom> &atoms, std::string run_id)=0;
};

/*--------------------- ReadLine() -----------------------*/
/*
read a full line from a file and 
return length of line read

to bad this looks like Pascal but its the easiest
way to read just whole line because fscanf() ignores
end of line characters

fpread = pointer to file
cMax = length of data buffer cRead
cRead = char[] buffer to read into
mesg = error message to print if not successful
*/
inline int ReadLine( FILE* fpRead, char* cRead, int cMax, const char *mesg )
{
  if( fgets( cRead, cMax, fpRead) == NULL ) {
    return 0;
    /*   printf("error reading input file: %s\n", mesg);
         exit( 0 );
    */
  }
  return( strlen( cRead ) );
}  /* end ReadLine */

}
#endif
