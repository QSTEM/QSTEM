#include "binary_output.hpp"

CBinaryOutput::CBinaryOutput() : IDataWriter()
{
}

CBinaryOutput::~CBinaryOutput()
{
}

void CBinaryOutput::DescribeFile(std::vector<ulong> shape, ulong element_size, 
                                 std::string label, 
                                 std::vector<ulong> position, 
                                 std::string comment, std::map<std::string, double> parameters,
                                 std::vector<float_tt> resolution)
{
  // Write a text description of the file we wrote:
  std::stringstream filename;
  std::map<std::string, double>::iterator param;
  filename<<label;
  for (ulong idx=0; idx<position.size(); idx++)
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

void CBinaryOutput::WriteRealVolume(float_tt *data, std::vector<ulong> shape, std::string label,
                                    std::vector<ulong> position, std::string comment, 
                                    std::map<std::string, double> parameters,
                                    std::vector<float_tt> resolution)
{
  WriteBlob(data, shape, label, position, parameters);
  if (resolution == std::vector<float_tt>())
    {
      resolution = std::vector<float_tt>(shape.size(),1.0);
    }
  DescribeFile(shape, sizeof(float_tt), label, position, comment, parameters, resolution);
}

void CBinaryOutput::WriteComplexVolume(complex_tt *data, std::vector<ulong> shape, std::string label, 
                                       std::vector<ulong> position, std::string comment, 
                                       std::map<std::string, double> parameters,
                                       std::vector<float_tt> resolution)
{
  WriteBlob(data, shape, label, position, parameters);
  if (resolution == std::vector<float_tt>())
    {
      resolution = std::vector<float_tt>(shape.size(),1.0);
    }
  DescribeFile(shape, sizeof(complex_tt), label, position, comment, parameters, resolution);
}

void CBinaryOutput::WriteRealImage(float_tt **data, std::vector<ulong> shape, std::string label, 
                                   std::vector<ulong> position, std::string comment, 
                                   std::map<std::string, double> parameters,
                                   std::vector<float_tt> resolution)
{
  WriteBlob(data[0], shape, label, position, parameters);
  if (resolution == std::vector<float_tt>())
    {
      resolution = std::vector<float_tt>(shape.size(),1.0);
    }
  DescribeFile(shape, sizeof(float_tt), label, position, comment, parameters);
}

void CBinaryOutput::WriteComplexImage(complex_tt **data, std::vector<ulong> shape, std::string label, 
                                      std::vector<ulong> position, std::string comment,
                                      std::map<std::string, double> parameters,
                                      std::vector<float_tt> resolution)
{
  WriteBlob(data[0], shape, label, position, parameters);
  if (resolution == std::vector<float_tt>())
    {
      resolution = std::vector<float_tt>(shape.size(),1);
    }
  DescribeFile(shape, sizeof(complex_tt), label, position, comment, parameters);
}

/*
void CBinaryOutput::WriteRealImage(QSfMat data, std::string label, std::vector<ulong> position, 
                    std::map<std::string, double> parameters)
{
  std::vector<ulong> shape(2);
  shape[0] = data.cols;
  shape[1] = data.rows;
  WriteBlob((float_tt *)data.data(), shape, label, position, parameters);
  DescribeFile(shape, label, position, parameters);
}


void CBinaryOutput::WriteComplexImage(QScMat data, std::string label, std::vector<ulong> position, 
                       std::map<std::string, double> parameters)
{
  std::vector<ulong> shape(2);
  shape[0] = data.cols;
  shape[1] = data.rows;
  WriteBlob((std::complex<float_tt> *)data.data(), shape, label, position, parameters);
  DescribeFile(shape, label, position, parameters);
}
*/
