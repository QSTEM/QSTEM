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
  BOOST_CHECK_EQUAL(atoms.size(), 5);
  BOOST_CHECK_CLOSE(atoms[0].x, 0, 0.05);
  BOOST_CHECK_CLOSE(atoms[0].y, 0, 0.05);
  BOOST_CHECK_CLOSE(atoms[0].z, 0, 0.05);
  BOOST_CHECK_CLOSE(atoms[0].dw, 0.6214, 0.05);  
  BOOST_CHECK_CLOSE(atoms[0].occ, 1.0, 0.05);
  BOOST_CHECK_CLOSE(atoms[0].q, 2.0, 0.05);
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
  }
  StructureWriterPtr writer;
};


BOOST_FIXTURE_TEST_SUITE(testCfgWrite, cfgWriterFixture)

BOOST_AUTO_TEST_CASE( testWriteCfg )
{
  std::vector<atom> atoms(5);
  atoms[0] = atom(76, "Sr", 0, 0, 0, 0.6214, 1.0, 2.0);
  atoms[1] = atom(44, "Ti", 0.5, 0.5, 0.5, 0.4390, 1.0, 4.0);
  atoms[2] = atom(16, "O", 0, 0.5, 0.5, 0.7323, 1.0, -2.0);
  atoms[3] = atom(16, "O", 0.5, 0, 0.5, 0.7323, 1.0, -2.0);
  atoms[4] = atom(16, "O", 0.5, 0.5, 0, 0.7323, 1.0, -2.0);
  writer->Write(atoms,1);
}

BOOST_AUTO_TEST_SUITE_END()













