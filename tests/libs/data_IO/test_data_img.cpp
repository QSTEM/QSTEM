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

struct imgFixture {
  imgFixture()
  {
    dataReader = CDataReaderFactory::Get()->GetReader(".img");
    dataWriter = CDataWriterFactory::Get()->GetWriter(".img");
    //std::cout << "setup qsc config reader fixture" << std::endl; 
  }
  ~imgFixture()
  { //std::cout << "teardown qsc config reader fixture" << std::endl;
  }
  DataReaderPtr dataReader;
  DataWriterPtr dataWriter;
};

BOOST_FIXTURE_TEST_SUITE(TestImgInput, imgFixture)

BOOST_AUTO_TEST_CASE( test_size )
{
  std::string filename="diffAvg_0_16";
  unsigned nx, ny;
  dataReader->ReadSize(filename, nx, ny);
  BOOST_CHECK_EQUAL(nx, 400);
  BOOST_CHECK_EQUAL(ny, 400);
}

BOOST_AUTO_TEST_CASE( testComment )
{
  std::string filename="diffAvg_0_16";
  std::string comment;
  dataReader->ReadComment(filename, comment);
  BOOST_CHECK_EQUAL(comment, "Average Array");
}

BOOST_AUTO_TEST_CASE( testParameters )
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

BOOST_AUTO_TEST_SUITE_END()









