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

#include <boost/test/unit_test.hpp>
#include <iostream>

#include "wavefunctions/wave_base.hpp"
#include "wavefunctions/wave_plane.hpp"
#include "wavefunctions/wave_convergent.hpp"


struct WaveFixture {
  WaveFixture():
    wave(WavePtr( new CWaveBase()))
  { 
    std::cout << "setup wave fixture" << std::endl; 
  }
  ~WaveFixture()
  { std::cout << "teardown wave fixture" << std::endl; }

  WavePtr wave;
};

BOOST_FIXTURE_TEST_SUITE (TestWave, WaveFixture)

BOOST_AUTO_TEST_CASE (testArrayAllocation)
{
  // Check array allocation
  BOOST_CHECK(wave->diffpat != NULL);
  BOOST_CHECK(wave->wave != NULL);
}

// Test image saving

// Test image reading

BOOST_AUTO_TEST_SUITE_END( )
