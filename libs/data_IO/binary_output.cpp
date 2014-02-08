#include "binary_output.hpp"

#include "boost/filesystem.hpp"   
using namespace boost::filesystem; 

CBinaryOutput::CBinaryOutput() : IDataWriter()
{
}

CBinaryOutput::~CBinaryOutput()
{
}

// TODO: decide what to do about run ids for non-qh5 datasets.
void CBinaryOutput::Initialize(std::string dirname, std::string run_id)
{
  char systStr[255];
  if (boost::filesystem::exists(dirname.c_str())) {
    printf(" (already exists)\n");
  }
  else {
    sprintf(systStr,"mkdir %s",dirname.c_str());
    system(systStr);
    printf(" (created)\n");
  }
}

void CBinaryOutput::DescribeFile(std::vector<unsigned> shape, unsigned element_size, 
                                 std::string label, 
                                 std::vector<unsigned> position, 
                                 std::string comment, std::map<std::string, double> parameters)
{
  // Write a text description of the file we wrote:
  std::stringstream filename;
  std::map<std::string, double>::iterator param;
  filename<<label;
  for (unsigned idx=0; idx<position.size(); idx++)
    {
      filename<<"_"<<idx;
    }
  filename<<".txt";
  std::ofstream outputFile;
  outputFile.open(filename.str().c_str(),std::ios::out);
  outputFile << "File comment:\n";
  outputFile << comment;
  outputFile << "File dimensions:\n";
  for (int dim=0; dim<shape.size(); dim++)
    {
      outputFile << dim << ", ";
    }
  outputFile << "\n";
  outputFile << "Data element size:\n";
  outputFile << element_size;
  outputFile << "\n";
  for (param=parameters.begin(); param!=parameters.end(); param++)
    {
      outputFile << param->first << ": " << param->second << "\n";
    }
  outputFile.close();
}

void CBinaryOutput::WriteRealVolume(float_tt *data, std::vector<unsigned> shape, std::string label,
                                    std::vector<unsigned> position, std::string comment, 
                                    std::map<std::string, double> parameters)
{
  WriteBlob(data, shape, label, position, parameters);
  DescribeFile(shape, sizeof(float_tt), label, position, comment, parameters);
}

void CBinaryOutput::WriteComplexVolume(complex_tt *data, std::vector<unsigned> shape, std::string label, 
                                       std::vector<unsigned> position, std::string comment, 
                                       std::map<std::string, double> parameters)
{
  WriteBlob(data, shape, label, position, parameters);
  DescribeFile(shape, sizeof(complex_tt), label, position, comment, parameters);
}

void CBinaryOutput::WriteRealImage(const RealVector &data, const std::vector<unsigned> &shape, const std::string &label, 
                                  const std::vector<unsigned> &position, const std::string &comment, 
								  const std::map<std::string, double> &parameters)
{
  WriteBlob(&data[0], shape, label, position, parameters);
  DescribeFile(shape, sizeof(float_tt), label, position, comment, parameters);
}

void CBinaryOutput::WriteComplexImage(const ComplexVector &data, const std::vector<unsigned> &shape, const std::string &label, 
                                  const std::vector<unsigned> &position, const std::string &comment, 
								  const std::map<std::string, double> &parameters)
{
  WriteBlob(&data[0], shape, label, position, parameters);
  DescribeFile(shape, sizeof(complex_tt), label, position, comment, parameters);
}

/*
void CBinaryOutput::WriteRealImage(QSfMat data, std::string label, std::vector<unsigned> position, 
                    std::map<std::string, double> parameters)
{
  std::vector<unsigned> shape(2);
  shape[0] = data.cols;
  shape[1] = data.rows;
  WriteBlob((float_tt *)data.data(), shape, label, position, parameters);
  DescribeFile(shape, label, position, parameters);
}


void CBinaryOutput::WriteComplexImage(QScMat data, std::string label, std::vector<unsigned> position, 
                       std::map<std::string, double> parameters)
{
  std::vector<unsigned> shape(2);
  shape[0] = data.cols;
  shape[1] = data.rows;
  WriteBlob((std::complex<float_tt> *)data.data(), shape, label, position, parameters);
  DescribeFile(shape, label, position, parameters);
}
*/
