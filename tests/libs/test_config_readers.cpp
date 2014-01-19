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

struct TEMQscFixture {
  TEMQscFixture()
  {
    std::string filename="tem_STO_4x4.qsc";
    configReader = GetConfigReader(filename);;
    //std::cout << "setup qsc config reader fixture" << std::endl; 
  }
  ~TEMQscFixture()
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

BOOST_AUTO_TEST_CASE( testPotOutputInterval )
{
  unsigned displayPotCalcInterval;
  configReader->ReadPotentialOutputInterval(displayPotCalcInterval);
  BOOST_CHECK_EQUAL(displayPotCalcInterval, 1000);
}

BOOST_AUTO_TEST_CASE( testOutputName )
{
  std::string filename;
  configReader->ReadOutputName(filename);
  BOOST_CHECK_EQUAL(filename, "\"STO_4x4\"");
}

BOOST_AUTO_TEST_CASE( testNCells )
{
  unsigned nx, ny, nz;
  configReader->ReadNCells(nx, ny, nz);
  BOOST_CHECK_EQUAL(nx, 9);
  BOOST_CHECK_EQUAL(ny, 9);
  BOOST_CHECK_EQUAL(nz, 20);
}

BOOST_AUTO_TEST_CASE( testNSubSlabs )
{
  unsigned cellDiv;
  configReader->ReadNSubSlabs(cellDiv);
  BOOST_CHECK_EQUAL(cellDiv, 1);
}

BOOST_AUTO_TEST_CASE( testCrystalCubeAndTilt )
{
  float_tt tx, ty, tz, cx, cy, cz;
  bool adjustCubeSize=false;
  configReader->ReadCrystalCubeAndTilt(tx, ty, tz, cx, cy, cz, adjustCubeSize);
  BOOST_CHECK_CLOSE(tx, 0, 0.00005);
  BOOST_CHECK_CLOSE(ty, 0, 0.00005);
  BOOST_CHECK_CLOSE(tz, 0, 0.00005);
  BOOST_CHECK_CLOSE(cx, 35.145, 0.00005);
  BOOST_CHECK_CLOSE(cy, 35.145, 0.00005);
  BOOST_CHECK_CLOSE(cz, 78.100, 0.00005);
  BOOST_CHECK_EQUAL(adjustCubeSize, false);
}

BOOST_AUTO_TEST_CASE( testTemperatureData )
{
  bool doTDS, useEinstein;
  std::string phononFile;
  float_tt tds_temp;
  configReader->ReadTemperatureData(doTDS, tds_temp, phononFile, useEinstein);
  BOOST_CHECK_EQUAL(doTDS, true);
  BOOST_CHECK_CLOSE(tds_temp, 300, 0.00005);
  BOOST_CHECK_EQUAL(phononFile, "");
  BOOST_CHECK_EQUAL(useEinstein, true);

}

BOOST_AUTO_TEST_SUITE_END( )




BOOST_FIXTURE_TEST_SUITE(STEM_SPECIFIC_READER, STEMQscFixture)

BOOST_AUTO_TEST_CASE( testProbeProgressInterval )
{
  unsigned interval;
  configReader->ReadSTEMProgressInterval(interval);
  BOOST_CHECK_EQUAL(interval, 12);
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_FIXTURE_TEST_SUITE(TEM_SPECIFIC_READER, TEMQscFixture)

BOOST_AUTO_TEST_CASE( testBeamTilt )
{
  float_tt tx, ty;
  bool tiltBack;
  configReader->ReadBeamTilt(tx, ty, tiltBack);
  BOOST_CHECK_CLOSE(tx, 0, 0.00005);
  BOOST_CHECK_CLOSE(ty, 0, 0.00005);
  BOOST_CHECK_EQUAL(tiltBack, true);
}

BOOST_AUTO_TEST_SUITE_END()
