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
#define BOOST_TEST_MODULE TestConfigReaders
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "config_readers.hpp"
#include "config_IO/read_interface.hpp"
#include "config_IO/read_qsc.hpp"

struct STEMQscFixture {
  STEMQscFixture()
  {
    std::string filename="stem_STO_4x4.qsc";
    configReader = GetConfigReader(filename);;
    //std::cout << "setup qsc config reader fixture" << std::endl; 
  }
  ~STEMQscFixture()
  { //std::cout << "teardown qsc config reader fixture" << std::endl;
  }
  ConfigReaderPtr configReader;
};

// Use the STEMQscFixture to run the base tests - could use any file, these traits are generic
BOOST_FIXTURE_TEST_SUITE (TestBaseQsc, STEMQscFixture)

BOOST_AUTO_TEST_CASE(testReadArraySize)
{
  unsigned nx, ny;
  configReader->ReadProbeArraySize(nx, ny);
  BOOST_CHECK_EQUAL(nx, 400);
  BOOST_CHECK_EQUAL(ny, 400);
}

BOOST_AUTO_TEST_CASE(testReadMode)
{
  std::string mode;
  configReader->ReadMode(mode);
  BOOST_CHECK_EQUAL(mode,"STEM");
}

BOOST_AUTO_TEST_CASE(testReadPrintLevel)
{
  unsigned printLevel;
  configReader->ReadPrintLevel(printLevel);
  BOOST_CHECK_EQUAL(printLevel,2);
}

BOOST_AUTO_TEST_CASE(testReadSaveLevel)
{
  unsigned saveLevel;
  configReader->ReadSaveLevel(saveLevel);
  BOOST_CHECK_EQUAL(saveLevel,1);
}

//BOOST_AUTO_TEST_CASE(

// Test image saving

// Test image reading

// Test probe calculation (parallel beam)

// Test probe calculation (convergent beam)

BOOST_AUTO_TEST_SUITE_END( )









