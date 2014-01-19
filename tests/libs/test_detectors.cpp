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
