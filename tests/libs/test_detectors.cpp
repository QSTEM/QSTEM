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


#define BOOST_TEST_MODULE TestDetectors
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "detectors.hpp"

struct DetectorFixture {
  DetectorFixture() :
    det(DetectorPtr(new Detector(10, 10, 0.2f, 0.2f, 25.0f)))
  { 
    std::cout << "setup detector fixture" << std::endl; 
  }
  ~DetectorFixture()
  { std::cout << "teardown detector" << std::endl; }

  DetectorPtr det;
};

BOOST_FIXTURE_TEST_SUITE (TestDetector, DetectorFixture)

BOOST_AUTO_TEST_CASE (testArrayAllocation)
{
  // Check array allocation
  //BOOST_CHECK(det->image != NULL);
  //BOOST_CHECK(det->image2 != NULL);
}

BOOST_AUTO_TEST_CASE (testSetComment)
{
	
}

BOOST_AUTO_TEST_SUITE_END( )
