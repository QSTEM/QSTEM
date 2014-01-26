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
#define BOOST_TEST_MODULE TestConvergentWave
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <exception>

#include "wavefunctions/wave_convergent.hpp"

#include "config_IO/config_reader_factory.hpp"


struct ConvergentWaveFixture {
  ConvergentWaveFixture()
  {
   configReader = CConfigReaderFactory::Get()->GetReader("stem_STO_4x4.qsc");
    // just default to using the STEM config file here
    wave = WavePtr(new CConvergentWave(configReader));
    //std::cout << "setup convergent wave fixture" << std::endl; 
  }
  ~ConvergentWaveFixture()
  { 
    //std::cout << "teardown convergent wave fixture" << std::endl; 
  }

  WavePtr wave;
  ConfigReaderPtr configReader;
};

BOOST_FIXTURE_TEST_SUITE(TestConvergentWave, ConvergentWaveFixture)

// Test probe calculation (convergent beam)

// Test copy constructor
BOOST_AUTO_TEST_CASE( testCopy )
{
  
}

BOOST_AUTO_TEST_CASE( testFormProbe )
{
  
  wave->FormProbe();
  // TODO: how to best test?  
  // this will work if I turn the data into a vector
  //BOOST_CHECK_EQUAL_COLLECTIONS(values.begin(), values.end(), 
  //                            expected.begin(), expected.end());
}

BOOST_AUTO_TEST_SUITE_END()
