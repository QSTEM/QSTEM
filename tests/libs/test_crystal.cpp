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


#define BOOST_TEST_MODULE TestCrystal
#include <boost/test/unit_test.hpp>

#include "crystal.hpp"

struct CrystalFixture {
  CrystalFixture()
  {
    configReader = CConfigReaderFactory::Get()->GetReader("stem_STO_4x4.qsc");
    cryst = StructurePtr(new CCrystal(configReader));
    //std::cout << "setup plane wave fixture" << std::endl; 
  }
  ~CrystalFixture()
  { 
    //std::cout << "teardown plane wave fixture" << std::endl; 
  }
  StructurePtr cryst;
  ConfigReaderPtr configReader;
};


BOOST_FIXTURE_TEST_SUITE(testCrystal, CrystalFixture)

BOOST_AUTO_TEST_CASE( testReadBaseAtoms )
{
  // we have already read the config file by default with the fixture.
  //  Make sure that it came out OK.
  float_tt alpha, beta, gamma;
  float_tt ax, by, cz;
  BOOST_CHECK_EQUAL(cryst->GetNumberOfCellAtoms(), 5);
  cryst->GetCellAngles(alpha, beta, gamma);
  BOOST_CHECK_CLOSE(alpha, 90, 0.05);
  BOOST_CHECK_CLOSE(beta, 90, 0.05);
  BOOST_CHECK_CLOSE(gamma, 90, 0.05);
  cryst->GetCellParameters(ax, by, cz);
  BOOST_CHECK_CLOSE(ax, 3.905, 0.05);
  BOOST_CHECK_CLOSE(by, 3.905, 0.05);
  BOOST_CHECK_CLOSE(cz, 3.905, 0.05);
}

BOOST_AUTO_TEST_CASE( testDuplicateAtoms )
{
  bool handleVacancies=false;
  cryst->ReplicateUnitCell(handleVacancies);
  BOOST_CHECK_EQUAL(cryst->GetNumberOfAtoms(),8840);
}

BOOST_AUTO_TEST_CASE( testPhononDisplacement )
{
  
}

BOOST_AUTO_TEST_CASE( testTiltBoxed )
{
  
}

BOOST_AUTO_TEST_CASE( testWriteStructure )
{
  
}

BOOST_AUTO_TEST_SUITE_END()
