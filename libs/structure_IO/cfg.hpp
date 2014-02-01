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

#ifndef CFG_STRUCTURE_H
#define CFG_STRUCTURE_H

#include "structureInterface.hpp"
#include <boost/filesystem.hpp>

class CCfgReader : public IStructureReader
{
public:
  CCfgReader(const boost::filesystem::path &structure_file);
  ~CCfgReader();
  
  int ReadCellParams(float_tt **Mm);
  int ReadAtoms(std::vector<atom> &atoms);
protected:
  int ReadNextAtom(atom *newAtomcom);
private:
  FILE *m_fp;
  std::vector<float_tt> m_atomData; // One line of atom data read from a file
  bool m_noVelocity;     // If true, the file does not provide velocity information
  unsigned m_entryCount; // The number of columns on any given line of data from the file
  bool m_isValid;  // If true, this object has a valid file to read atoms from.
  std::string m_buf;
private:
  friend class CStructureReaderFactory;
  // Create an instance of this class, wrapped in a shared ptr
  //     This should not be inherited - any subclass needs its own implementation.
  static StructureReaderPtr Create(const boost::filesystem::path &file){return StructureReaderPtr(new CCfgReader(file));}
};

class CCfgWriter : public IStructureWriter
{
public:
  CCfgWriter(const boost::filesystem::path &file, float_tt ax, float_tt by, float_tt cz);
  ~CCfgWriter();

  int Write(std::vector <atom> &atoms, unsigned run_number);
  int WriteFractCubic(double *pos,int *Znum,double *dw,int natoms,char *fileName,
                  double a,double b,double c);
private:
  FILE *m_fp;
  float_tt m_ax, m_by, m_cz;
private:
  friend class CStructureWriterFactory;
  // Create an instance of this class, wrapped in a shared ptr
  //     This should not be inherited - any subclass needs its own implementation.
  static StructureWriterPtr Create(const boost::filesystem::path &file, float_tt ax, float_tt by, float_tt cz){
    return StructureWriterPtr(new CCfgWriter(file, ax, by, cz));
  }  
};

#endif

















