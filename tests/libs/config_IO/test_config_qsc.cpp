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
#define BOOST_TEST_MODULE TestConfigReaders
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <iostream>

#include "config_IO/config_reader_factory.hpp"

struct STEMQscFixture {
  STEMQscFixture()
  {
   configReader = CConfigReaderFactory::Get()->GetReader("stem_STO_4x4.qsc");
    //std::cout << "setup qsc config reader fixture" << std::endl; 
  }
  ~STEMQscFixture()
  { //std::cout << "teardown qsc config reader fixture" << std::endl;
  }
  ConfigReaderPtr configReader;
};

struct TEMQscFixture {
  TEMQscFixture()
  {
   configReader = CConfigReaderFactory::Get()->GetReader("tem_STO.qsc");
   //std::cout << "setup qsc config reader fixture" << std::endl; 
  }
  ~TEMQscFixture()
  { //std::cout << "teardown qsc config reader fixture" << std::endl;
  }
  ConfigReaderPtr configReader;
};

/*
struct TomoQscFixture {
  TomoQscFixture()
  {
    std::string filename="tomo_STO.qsc";
    configReader = GetConfigReader(filename);;
    //std::cout << "setup qsc config reader fixture" << std::endl; 
  }
  ~TomoQscFixture()
  { //std::cout << "teardown qsc config reader fixture" << std::endl;
  }
  ConfigReaderPtr configReader;
};
*/

// Use the STEMQscFixture to run the base tests - could use any file, these traits are generic
BOOST_FIXTURE_TEST_SUITE (TestBaseQsc, STEMQscFixture)

BOOST_AUTO_TEST_CASE(testReadMode)
{
  std::string mode;
  configReader->ReadMode(mode);
  BOOST_CHECK_EQUAL(mode,"STEM");
}

BOOST_AUTO_TEST_CASE(testReadPrintLevel)
{
  unsigned printLevel;
  configReader->ReadPrintLevel(printLevel);
  BOOST_CHECK_EQUAL(printLevel,2);
}

BOOST_AUTO_TEST_CASE(testReadSaveLevel)
{
  unsigned saveLevel;
  configReader->ReadSaveLevel(saveLevel);
  BOOST_CHECK_EQUAL(saveLevel,1);
}

BOOST_AUTO_TEST_CASE( testPotOutputInterval )
{
  unsigned displayPotCalcInterval;
  configReader->ReadPotentialOutputInterval(displayPotCalcInterval);
  BOOST_CHECK_EQUAL(displayPotCalcInterval, 1000);
}

BOOST_AUTO_TEST_CASE( testOutputName )
{
  std::string filename;
  configReader->ReadOutputName(filename);
  BOOST_CHECK_EQUAL(filename, "\"STO_4x4\"");
}

BOOST_AUTO_TEST_CASE( testNCells )
{
  unsigned nx, ny, nz;
  configReader->ReadNCells(nx, ny, nz);
  BOOST_CHECK_EQUAL(nx, 9);
  BOOST_CHECK_EQUAL(ny, 9);
  BOOST_CHECK_EQUAL(nz, 20);
}

BOOST_AUTO_TEST_CASE( testNSubSlabs )
{
  unsigned cellDiv;
  configReader->ReadNSubSlabs(cellDiv);
  BOOST_CHECK_EQUAL(cellDiv, 1);
}

BOOST_AUTO_TEST_CASE( testCrystalCubeAndTilt )
{
  float_tt tx, ty, tz, cx, cy, cz;
  bool adjustCubeSize=false;
  configReader->ReadCrystalCubeAndTilt(tx, ty, tz, cx, cy, cz, adjustCubeSize);
  BOOST_CHECK_CLOSE(tx, 0, 0.00005);
  BOOST_CHECK_CLOSE(ty, 0, 0.00005);
  BOOST_CHECK_CLOSE(tz, 0, 0.00005);
  BOOST_CHECK_CLOSE(cx, 35.145, 0.00005);
  BOOST_CHECK_CLOSE(cy, 35.145, 0.00005);
  BOOST_CHECK_CLOSE(cz, 78.100, 0.00005);
  BOOST_CHECK_EQUAL(adjustCubeSize, false);
}

BOOST_AUTO_TEST_CASE( testTemperatureData )
{
  bool doTDS, useEinstein;
  std::string phononFile;
  float_tt tds_temp;
  configReader->ReadTemperatureData(doTDS, tds_temp, phononFile, useEinstein);
  BOOST_CHECK_EQUAL(doTDS, true);
  BOOST_CHECK_CLOSE(tds_temp, 300, 0.00005);
  BOOST_CHECK_EQUAL(phononFile, "");
  BOOST_CHECK_EQUAL(useEinstein, true);

}

BOOST_AUTO_TEST_CASE( testReadOffset )
{
  float_tt x=0, y=0;
  configReader->ReadSliceOffset(x, y);
  BOOST_CHECK_CLOSE(x, 0, 0.00005);
  BOOST_CHECK_CLOSE(y, 0, 0.00005);
}

BOOST_AUTO_TEST_CASE(testReadArraySize)
{
  unsigned nx, ny;
  configReader->ReadProbeArraySize(nx, ny);
  BOOST_CHECK_EQUAL(nx, 400);
  BOOST_CHECK_EQUAL(ny, 400);
}

BOOST_AUTO_TEST_CASE( testResolution )
{
  float_tt x=1, y=1;
  configReader->ReadResolution(x, y);
  BOOST_CHECK_CLOSE(x, 0.0625, 0.00005);
  BOOST_CHECK_CLOSE(y, 0.0625, 0.00005);
}

BOOST_AUTO_TEST_CASE( testVoltage )
{
  float_tt v0;
  configReader->ReadVoltage(v0);
  BOOST_CHECK_CLOSE(v0, 200, 0.00005);
}

BOOST_AUTO_TEST_CASE( testSliceParams )
{
  bool center=true; //if read properly, should return false.
  float_tt thickness=1, zOffset=0;
  unsigned nslices, outputInterval;
  configReader->ReadSliceParameters(center, thickness, nslices, outputInterval, zOffset);
  BOOST_CHECK_EQUAL(center, false);
  BOOST_CHECK_CLOSE(thickness, 1.925, 2);
  BOOST_CHECK_EQUAL(nslices, 40);
  BOOST_CHECK_EQUAL(outputInterval, 8);
  BOOST_CHECK_CLOSE(zOffset, 0.976, 0.00005);
}

BOOST_AUTO_TEST_CASE( testPeriodicParams )
{
  bool xy, z;
  configReader->ReadPeriodicParameters(xy, z);
  BOOST_CHECK_EQUAL(xy, false);
  BOOST_CHECK_EQUAL(z, false);
}

BOOST_AUTO_TEST_CASE( testBandLimitTrans )
{
  bool limit=true;
  // should flip it to false
  configReader->ReadBandLimitTrans(limit);
  BOOST_CHECK_EQUAL(limit, false);
}

BOOST_AUTO_TEST_CASE(testLoadPotential)
{
  bool load=false;
  configReader->ReadLoadPotential(load);
  BOOST_CHECK_EQUAL(load, false);
}

BOOST_AUTO_TEST_CASE( testPotOutputParams )
{
  bool savePot, saveProjPot, plotPot;
  configReader->ReadPotentialOutputParameters(savePot, saveProjPot, plotPot);
  BOOST_CHECK_EQUAL(savePot, false);
  BOOST_CHECK_EQUAL(saveProjPot, false);
  BOOST_CHECK_EQUAL(plotPot, false);
}

BOOST_AUTO_TEST_CASE( testPotCalcParams )
{
  bool _3D, fft;
  configReader->ReadPotentialCalculationParameters(fft, _3D);
  BOOST_CHECK_EQUAL(_3D, true);
  BOOST_CHECK_EQUAL(fft, true);
}

BOOST_AUTO_TEST_CASE( testAtomRadius )
{
  float_tt r;
  configReader->ReadAtomRadius(r);
  BOOST_CHECK_CLOSE(r, 5.0, 0.00005);
}

BOOST_AUTO_TEST_CASE( testStructureFactorType )
{
  std::string type;
  configReader->ReadStructureFactorType(type);
  BOOST_CHECK_EQUAL(type, "WK");
}

BOOST_AUTO_TEST_CASE( testAvgParameters )
{
  unsigned avgruns;
  bool storeSeries;
  configReader->ReadAverageParameters(avgruns, storeSeries);
  BOOST_CHECK_EQUAL(avgruns, 10);
  // TODO: need to check store series somewhere else
}

BOOST_AUTO_TEST_CASE( testStructureFilename )
{
  boost::filesystem::path filename;
  configReader->ReadStructureFileName(filename);
  BOOST_CHECK_EQUAL(filename, boost::filesystem::path("SrTiO3.cfg"));
}



BOOST_AUTO_TEST_SUITE_END( )




BOOST_FIXTURE_TEST_SUITE(STEM_SPECIFIC_READER, STEMQscFixture)

BOOST_AUTO_TEST_CASE( testProbeProgressInterval )
{
  unsigned interval;
  configReader->ReadSTEMProgressInterval(interval);
  BOOST_CHECK_EQUAL(interval, 12);
}

BOOST_AUTO_TEST_CASE( testScanParameters )
{
  float_tt xstart, xstop, ystart, ystop;
  unsigned nx, ny;
  configReader->ReadScanParameters(xstart, xstop, nx, ystart, ystop, ny);
  BOOST_CHECK_CLOSE(xstart, 11.67, 0.05);
  BOOST_CHECK_CLOSE(xstop, 19.71, 0.05);
  BOOST_CHECK_CLOSE(ystart, 11.45, 0.05);
  BOOST_CHECK_CLOSE(ystop, 19.32, 0.05);
  BOOST_CHECK_EQUAL(nx, 8);
  BOOST_CHECK_EQUAL(ny, 8);
}

BOOST_AUTO_TEST_CASE( testNumberOfDetectors )
{
  int num;
  configReader->ReadNumberOfDetectors(num);
  BOOST_CHECK_EQUAL(num, 4);
}

BOOST_AUTO_TEST_CASE( testDetectorParameters )
{
  float_tt rInside, rOutside, shiftx, shifty;
  std::string name;
  configReader->ReadDetectorParameters(1, rInside, rOutside, name, shiftx, shifty);
  BOOST_CHECK_CLOSE(rInside, 50, 0.05);
  BOOST_CHECK_CLOSE(rOutside, 70, 0.05);
  BOOST_CHECK_EQUAL(name, "detector2");
  BOOST_CHECK_CLOSE(shiftx, 0, 0.05);
  BOOST_CHECK_CLOSE(shifty, 0, 0.05);
}

BOOST_AUTO_TEST_CASE( testDoseParameters )
{
  float_tt beamCurrent, dwellTimeMs;
  configReader->ReadDoseParameters(beamCurrent, dwellTimeMs);
  BOOST_CHECK_CLOSE(beamCurrent, 1, .05);
  BOOST_CHECK_CLOSE(dwellTimeMs, 1.6021773e-4, .05);
}

BOOST_AUTO_TEST_CASE( testProbeParameters )
{
  float_tt dE_E=0, dI_I=0, dV_V=0, alpha, aAIS=0, sourceRadius=0;
  configReader->ReadProbeParameters(dE_E, dI_I, dV_V, alpha, aAIS, sourceRadius);
  BOOST_CHECK_CLOSE(dE_E, 0, 0.05);
  BOOST_CHECK_CLOSE(dI_I, 0, 0.05);
  BOOST_CHECK_CLOSE(dV_V, 0.000003, 0.05);
  BOOST_CHECK_CLOSE(alpha, 15, 0.05);
  BOOST_CHECK_CLOSE(aAIS, 0, 0.05);
  BOOST_CHECK_CLOSE(sourceRadius, 0, 0.05);
}

BOOST_AUTO_TEST_CASE( testProbeSmoothing )
{
  bool smooth, gaussFlag;
  float_tt gaussScale=0;
  configReader->ReadSmoothingParameters(smooth, gaussScale, gaussFlag);
  BOOST_CHECK_EQUAL(smooth, true);
  BOOST_CHECK_CLOSE(gaussScale, 0, 0.05);
  BOOST_CHECK_EQUAL(gaussFlag, false);
}

BOOST_AUTO_TEST_CASE( testAberrationAmplitudes )
{
  float_tt Cs=0, C5=0, Cc=0, df0=0, astig=0, \
    a33=0, a31=0, a42=0, a44=0, a55=0, \
    a53=0, a51=0, a66=0, a64=0, a62=0;
  std::string Scherzer;
  configReader->ReadAberrationAmplitudes(Cs, C5, Cc, df0, Scherzer, astig,
                                         a33, a31, a44, a42, a55, a53, a51,
                                         a66, a64, a62);
  BOOST_CHECK_CLOSE(Cs, 0.05, 0.00005);
  BOOST_CHECK_CLOSE(C5, 0, 0.05);
  BOOST_CHECK_CLOSE(Cc, 1.0, 0.00005);
  BOOST_CHECK_CLOSE(df0, -13.7, 0.00005);
  BOOST_CHECK_CLOSE(astig, 0, 0.05);
  BOOST_CHECK_CLOSE(a33, 0, 0.05);
  BOOST_CHECK_CLOSE(a31, 0, 0.05);
  BOOST_CHECK_CLOSE(a44, 0, 0.05);
  BOOST_CHECK_CLOSE(a42, 0, 0.05);
  BOOST_CHECK_CLOSE(a55, 0, 0.05);
  BOOST_CHECK_CLOSE(a53, 0, 0.05);
  BOOST_CHECK_CLOSE(a51, 0, 0.05);
  BOOST_CHECK_CLOSE(a66, 0, 0.05);
  BOOST_CHECK_CLOSE(a64, 0, 0.05);
  BOOST_CHECK_CLOSE(a62, 0, 0.05);
}


BOOST_AUTO_TEST_CASE( testAberrationAngles )
{
  float_tt astig=0, phi33=0, phi31=0, phi42=0, phi44=0, phi55=0, \
    phi53=0, phi51=0, phi66=0, phi64=0, phi62=0;
  configReader->ReadAberrationAngles(astig, phi33, phi31, phi44, phi42, phi55, phi53, phi51,
                                         phi66, phi64, phi62);
  BOOST_CHECK_CLOSE(astig, 0.0, 0.05);
  BOOST_CHECK_CLOSE(phi33, 0.0, 0.05);
  BOOST_CHECK_CLOSE(phi31, 0.0, 0.05);
  BOOST_CHECK_CLOSE(phi44, 0.0, 0.05);
  BOOST_CHECK_CLOSE(phi42, 0.0, 0.05);
  BOOST_CHECK_CLOSE(phi55, 0.0, 0.05);
  BOOST_CHECK_CLOSE(phi53, 0.0, 0.05);
  BOOST_CHECK_CLOSE(phi51, 0.0, 0.05);
  BOOST_CHECK_CLOSE(phi66, 0.0, 0.05);
  BOOST_CHECK_CLOSE(phi64, 0.0, 0.05);
  BOOST_CHECK_CLOSE(phi62, 0.0, 0.05);
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_FIXTURE_TEST_SUITE(TEM_SPECIFIC_READER, TEMQscFixture)

BOOST_AUTO_TEST_CASE( testBeamTilt )
{
  float_tt tx, ty;
  bool tiltBack;
  configReader->ReadBeamTilt(tx, ty, tiltBack);
  // 1 degree - output of function is in rad.
  BOOST_CHECK_CLOSE(tx, .01745, 0.05);
  // 2 degrees
  BOOST_CHECK_CLOSE(ty, 0.0349, 0.05);
  BOOST_CHECK_EQUAL(tiltBack, true);
}

BOOST_AUTO_TEST_CASE( testPendelloesung )
{
  std::vector<int> hbeams, kbeams;
  bool lbeams;
  unsigned nbout;
  configReader->ReadPendelloesungParameters(hbeams, kbeams, lbeams, nbout);
  // TODO: define tests for this
}

BOOST_AUTO_TEST_SUITE_END()



/*
BOOST_FIXTURE_TEST_SUITE(TOMO_READER, TomoQscFixture)

BOOST_AUTO_TEST_CASE( testTomoParameters )
{
  float_tt tilt, start, step;
  int count;
  float_tt zoomFactor;
  configReader->ReadTomoParameters(tilt, start, step, count, zoomFactor);
  //TODO: define tests for this
}

BOOST_AUTO_TEST_SUITE_END()
*/
