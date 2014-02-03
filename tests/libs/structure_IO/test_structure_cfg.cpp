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
#define BOOST_TEST_MODULE TestCfgStructure
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <iostream>

#include "memory_fftw3.hpp"
#include "structure_IO/structure_factories.hpp"

struct cfgReaderFixture {
  cfgReaderFixture()
  {
    boost::filesystem::path cfgFile="SrTiO3.cfg";
    reader = CStructureReaderFactory::Get()->GetReader(cfgFile);
    //std::cout << "setup qsc config reader fixture" << std::endl; 
  }
  ~cfgReaderFixture()
  { //std::cout << "teardown qsc config reader fixture" << std::endl;
  }
  StructureReaderPtr reader;
};

BOOST_FIXTURE_TEST_SUITE(testCfgRead, cfgReaderFixture)

BOOST_AUTO_TEST_CASE(testReadMm)
{
  float_tt **Mm = float2D(3,3,"");
  reader->ReadCellParams(Mm);
  BOOST_CHECK_CLOSE(Mm[0][0], 3.905, 0.1);
  BOOST_CHECK_CLOSE(Mm[1][1], 3.905, 0.1);
  BOOST_CHECK_CLOSE(Mm[2][2], 3.905, 0.1);
}

BOOST_AUTO_TEST_CASE( testReadAtoms )
{
  std::vector<atom> atoms;
  reader->ReadAtoms(atoms);
  // Make sure that we read all the atoms in the file
  BOOST_CHECK_EQUAL(atoms.size(), 5);
  // Make sure that we read the first atom correctly
  BOOST_CHECK_CLOSE(atoms[0].x, 0, 0.05);
  BOOST_CHECK_CLOSE(atoms[0].y, 0, 0.05);
  BOOST_CHECK_CLOSE(atoms[0].z, 0, 0.05);
  BOOST_CHECK_CLOSE(atoms[0].dw, 0.6214, 0.05);  
  BOOST_CHECK_CLOSE(atoms[0].occ, 1.0, 0.05);
  BOOST_CHECK_CLOSE(atoms[0].q, 2.0, 0.05);
  BOOST_CHECK_CLOSE(atoms[0].mass, 76, 0.05);
  BOOST_CHECK_EQUAL(atoms[0].Znum, 38);
  // make sure second atom is read correctly - if it is, then the loop over atoms is behaving OK.
  BOOST_CHECK_CLOSE(atoms[1].x, 0.5, 0.05);
  BOOST_CHECK_CLOSE(atoms[1].y, 0.5, 0.05);
  BOOST_CHECK_CLOSE(atoms[1].z, 0.5, 0.05);
  BOOST_CHECK_CLOSE(atoms[1].dw, 0.4390, 0.05);  
  BOOST_CHECK_CLOSE(atoms[1].occ, 1.0, 0.05);
  BOOST_CHECK_CLOSE(atoms[1].q, 4.0, 0.05);
  BOOST_CHECK_CLOSE(atoms[1].mass, 44, 0.05);
  BOOST_CHECK_EQUAL(atoms[1].Znum, 22);
}

BOOST_AUTO_TEST_SUITE_END()


struct cfgWriterFixture {
  cfgWriterFixture()
  {
    boost::filesystem::path cfgFile="SrTiO3_test_write.cfg";
    writer = CStructureWriterFactory::Get()->GetWriter(cfgFile,3.905, 3.905, 3.905);
    //std::cout << "setup qsc config reader fixture" << std::endl; 
  }
  ~cfgWriterFixture()
  { //std::cout << "teardown qsc config reader fixture" << std::endl;
	  // delete the file now that we're done with it
    boost::filesystem::remove("SrTiO3_test_write.cfg");
  }
  StructureWriterPtr writer;
};


BOOST_FIXTURE_TEST_SUITE(testCfgWrite, cfgWriterFixture)

BOOST_AUTO_TEST_CASE( testWriteCfg )
{
  std::vector<atom> atoms(5);
  // input expects atom positions in true space; converts them into fractional coordinates
  //    before writing
  float_tt a=3.905, b=3.905, c=3.905;
  atoms[0] = atom(76, "Sr", 0*a, 0*b, 0*c, 0.6214, 1.0, 2.0);
  atoms[1] = atom(44, "Ti", 0.5*a, 0.5*b, 0.5*c, 0.4390, 1.0, 4.0);
  atoms[2] = atom(16, "O", 0*a, 0.5*b, 0.5*c, 0.7323, 1.0, -2.0);
  atoms[3] = atom(16, "O", 0.5*a, 0*b, 0.5*c, 0.7323, 1.0, -2.0);
  atoms[4] = atom(16, "O", 0.5*a, 0.5*b, 0*c, 0.7323, 1.0, -2.0);
  writer->Write(atoms,"1");
}

BOOST_AUTO_TEST_SUITE_END()













