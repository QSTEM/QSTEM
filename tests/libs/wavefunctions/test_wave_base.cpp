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

#define BOOST_TEST_MODULE TestWaveBase
#include <boost/test/unit_test.hpp>

#include "wavefunctions/wave_interface.hpp"
#include "wavefunctions/wave_plane.hpp"

#include "config_IO/config_reader_factory.hpp"

#include "boost_collection_close.hpp"

// Although this uses the Plane wave fixture, we're testing fundamental methods
//   that both plane waves and convergent waves share (both inherit, and neither 
//   override the base class methods tested here.)

struct PlaneWaveFixture {
  PlaneWaveFixture()
  {
    wave = WavePtr(new CPlaneWave());
   configReader = CConfigReaderFactory::Get()->GetReader("tem_STO.qsc");
    //std::cout << "setup plane wave fixture" << std::endl; 
  }
  ~PlaneWaveFixture()
  { 
    //std::cout << "teardown plane wave fixture" << std::endl; 
  }
  WavePtr wave;
  ConfigReaderPtr configReader;
};

BOOST_FIXTURE_TEST_SUITE(TestBaseWave, PlaneWaveFixture)

BOOST_AUTO_TEST_CASE (testArrayAllocation)
{
  wave->Resize(512, 512);
  // Check array allocation
  BOOST_CHECK(wave->GetDPPointer() != NULL);
  BOOST_CHECK(wave->GetWavePointer() != NULL);
}

BOOST_AUTO_TEST_CASE(testReadCfgFromFile)
{
  BOOST_REQUIRE(configReader->IsValid());
  
  wave=WavePtr(new CPlaneWave(configReader));
  unsigned nx, ny;
  wave->GetSizePixels(nx, ny);
  BOOST_CHECK(nx == 400);
  BOOST_CHECK(ny == 400);
}

BOOST_AUTO_TEST_CASE(testSetPosition)
{
}

BOOST_AUTO_TEST_CASE( testReadWave )
{
  // our reference data is mulswav_16_2.img
  wave->ReadWave(16, 2);
}

BOOST_AUTO_TEST_CASE( testReadDiffPat )
{
  //wave->ReadDiffPat();
  
}

BOOST_AUTO_TEST_CASE( testWriteWave )
{
  
}

BOOST_AUTO_TEST_CASE( testWriteDiffPat )
{
  
}

BOOST_AUTO_TEST_SUITE_END()
