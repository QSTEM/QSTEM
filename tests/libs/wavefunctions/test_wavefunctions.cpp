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
#define BOOST_TEST_MODULE TestWavefunctions
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "wavefunctions/wave_interface.hpp"
#include "wavefunctions/wave_base.hpp"
#include "wavefunctions/wave_plane.hpp"
#include "wavefunctions/wave_convergent.hpp"

#include "config_readers.hpp"

struct PlaneWaveFixture {
  PlaneWaveFixture()
  {
    wave = WavePtr(new CPlaneWave());
    std::cout << "setup plane wave fixture" << std::endl; 
  }
  ~PlaneWaveFixture()
  { std::cout << "teardown plane wave fixture" << std::endl; }

  WavePtr wave;
};

BOOST_FIXTURE_TEST_SUITE (TestWave, PlaneWaveFixture)

BOOST_AUTO_TEST_CASE (testArrayAllocation)
{
  wave->Resize(512, 512);
  // Check array allocation
  BOOST_CHECK(wave->GetDPPointer() != NULL);
  BOOST_CHECK(wave->GetWavePointer() != NULL);
}

BOOST_AUTO_TEST_CASE(testReadCfgFromFile)
{
  ConfigReaderPtr configReader=GetConfigReader("stem_STO_4x4.qsc");
  wave=WavePtr(new CPlaneWave(configReader));
  unsigned nx, ny;
  wave->GetSizePixels(nx, ny);
  BOOST_CHECK(nx == 400);
  BOOST_CHECK(ny == 400);
}

// Test image saving

// Test image reading

// Test probe calculation (parallel beam)

// Test probe calculation (convergent beam)

BOOST_AUTO_TEST_SUITE_END( )


