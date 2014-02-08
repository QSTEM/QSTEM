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
#define BOOST_TEST_MODULE TestImgIO
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <iostream>

#include "data_IO/data_io_factories.hpp"

struct imgReaderFixture {
  imgReaderFixture()
  {
    dataReader = CDataReaderFactory::Get()->GetReader(".img");
    //std::cout << "setup qsc config reader fixture" << std::endl; 
  }
  ~imgReaderFixture()
  { //std::cout << "teardown qsc config reader fixture" << std::endl;
  }
  DataReaderPtr dataReader;
};

struct imgWriterFixture {
  imgWriterFixture()
  {
    dataWriter = CDataWriterFactory::Get()->GetWriter(".img");
    //std::cout << "setup qsc config reader fixture" << std::endl; 
  }
  ~imgWriterFixture()
  { //std::cout << "teardown qsc config reader fixture" << std::endl;
  }
  DataWriterPtr dataWriter;
};

BOOST_FIXTURE_TEST_SUITE(TestImgInput, imgReaderFixture)

BOOST_AUTO_TEST_CASE( testReadSize )
{
  std::string filename="diffAvg_0_16";
  unsigned nx, ny;
  dataReader->ReadSize(filename, nx, ny);
  BOOST_CHECK_EQUAL(nx, 400);
  BOOST_CHECK_EQUAL(ny, 400);
}

BOOST_AUTO_TEST_CASE( testReadComment )
{
  std::string filename="diffAvg_0_16";
  std::string comment;
  dataReader->ReadComment(filename, comment);
  BOOST_CHECK_EQUAL(comment, "Average Array");
}

BOOST_AUTO_TEST_CASE( testReadParameters )
{
  std::map<std::string, double> parameters;
  std::string filename="diffAvg_0_16";
  dataReader->ReadParameters(filename, parameters);
  // We don't have a meaningful way to know what parameters actually are, beyond these few:
  BOOST_CHECK_EQUAL(parameters.size(), 3);
  BOOST_CHECK_CLOSE(parameters["Thickness"],78.099,0.05);
  BOOST_CHECK_CLOSE(parameters["dx"],0.04, 0.05);
  BOOST_CHECK_CLOSE(parameters["dy"],0.04, 0.05);
}

BOOST_AUTO_TEST_CASE( testReadComplex )
{
  std::string filename="diffAvg_0_16";
  bool complex;
  dataReader->ReadComplex(filename, complex);
  BOOST_CHECK_EQUAL(complex, false);
}

BOOST_AUTO_TEST_CASE( testReadByteSize )
{
  std::string filename="diffAvg_0_16";
  unsigned byteSize;
  dataReader->ReadElementByteSize(filename, byteSize);
  BOOST_CHECK_EQUAL(byteSize, 4);
}

BOOST_AUTO_TEST_SUITE_END()




BOOST_FIXTURE_TEST_SUITE(TestImgOutput, imgWriterFixture)

BOOST_AUTO_TEST_CASE(testWriteHeader)
{
	std::vector<unsigned> size(2);
	std::map<std::string, double> pars;
}


BOOST_AUTO_TEST_CASE(testWriteFile)
{
	std::vector<unsigned> size(2,500);
	RealVector data(500*500,0);
	std::map<std::string, double> pars;
	pars["dx"]=0.5;
	pars["dy"]=0.5;
	pars["Thickness"]=10;
	std::string label = "test";
	std::string comment = "test_comment";
	dataWriter->WriteRealImage(data, size, label, comment, pars);
}

BOOST_AUTO_TEST_SUITE_END()
