#include <boost/test/unit_test.hpp>

#include "data_containers.h"
#include <iostream>

int add(int i, int j)
{
    return i + j;
}
 
BOOST_AUTO_TEST_SUITE(Maths)
 
BOOST_AUTO_TEST_CASE(universeInOrder)
{
    BOOST_CHECK(add(2, 2) == 4);
}

BOOST_AUTO_TEST_SUITE_END( )


struct WaveFixture {
    WaveFixture() 
	{ 
		std::cout << "setup" << std::endl; 
		WavePtr wave = WavePtr( new WAVEFUNC(10, 10, 0.2f, 0.2f));
	}
    ~WaveFixture()
	{ std::cout << "teardown" << std::endl; }

    WavePtr wave;
};

struct DetectorFixture {
    DetectorFixture() 
	{ 
		std::cout << "setup" << std::endl; 
		DetectorPtr det = DetectorPtr( new Detector(10, 10, 0.2f, 0.2f));
	}
    ~DetectorFixture()
	{ std::cout << "teardown" << std::endl; }

    DetectorPtr det;
};
/*
BOOST_FIXTURE_TEST_SUITE (TestWave, WaveFixture) // name of the test suite is wavetest

BOOST_AUTO_TEST_CASE (testArrayAllocation)
{
  // Check array allocation
  //BOOST_CHECK(wave->diffpat != NULL);
  //BOOST_CHECK(wave->avgArray != NULL);
  //BOOST_CHECK(wave->wave != NULL);
	BOOST_CHECK(true);
}

// Test image saving

// Test image reading

BOOST_AUTO_TEST_SUITE_END( )


BOOST_FIXTURE_TEST_SUITE (TestDetector, DetectorFixture) // name of the test suite is wavetest

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
*/