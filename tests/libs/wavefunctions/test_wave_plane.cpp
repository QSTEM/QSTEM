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
#define BOOST_TEST_MODULE TestPlaneWave
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "wavefunctions/wave_plane.hpp"

#include "config_IO/config_reader_factory.hpp"

using namespace QSTEM;

struct PlaneWaveFixture {
  PlaneWaveFixture()
  {
   configReader = CConfigReaderFactory::Get()->GetReader("tem_STO.qsc");
    wave = WavePtr(new CPlaneWave(configReader));
    //std::cout << "setup plane wave fixture" << std::endl; 
  }
  ~PlaneWaveFixture()
  { 
    //std::cout << "teardown plane wave fixture" << std::endl; 
  }
  WavePtr wave;
  ConfigReaderPtr configReader;
};

// Tests specific to the Plane Wave class
BOOST_FIXTURE_TEST_SUITE (TestPlaneWave, PlaneWaveFixture)

// Test probe calculation (parallel beam)




BOOST_AUTO_TEST_SUITE_END( )

