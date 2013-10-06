#include <boost/test/unit_test.hpp>

#include "data_containers.hpp"
#include <iostream>

struct WaveFixture {
  WaveFixture():
    wave(WavePtr( new WAVEFUNC(10, 10, 1.0, 1.0)))
  { 
    std::cout << "setup wave fixture" << std::endl; 
  }
  ~WaveFixture()
  { std::cout << "teardown wave fixture" << std::endl; }

  WavePtr wave;
};

struct DetectorFixture {
  DetectorFixture() :
    det(DetectorPtr(new Detector(10, 10, 0.2f, 0.2f)))
  { 
    std::cout << "setup detector fixture" << std::endl; 
  }
  ~DetectorFixture()
  { std::cout << "teardown detector" << std::endl; }

  DetectorPtr det;
};

BOOST_FIXTURE_TEST_SUITE (TestWave, WaveFixture)

BOOST_AUTO_TEST_CASE (testArrayAllocation)
{
  // Check array allocation
  BOOST_CHECK(wave->diffpat != NULL);
  BOOST_CHECK(wave->avgArray != NULL);
  BOOST_CHECK(wave->wave != NULL);
}

// Test image saving

// Test image reading

BOOST_AUTO_TEST_SUITE_END( )


BOOST_FIXTURE_TEST_SUITE (TestDetector, DetectorFixture)

BOOST_AUTO_TEST_CASE (testArrayAllocation)
{
  // Check array allocation
  BOOST_CHECK(det->image != NULL);
  BOOST_CHECK(det->image2 != NULL);
}

BOOST_AUTO_TEST_CASE (testSetComment)
{
	
}

BOOST_AUTO_TEST_SUITE_END( )
