#ifndef BIN_OUTPUT_H
#define BIN_OUTPUT_H

#include "output_interface.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

#include <stdexcept>

/**
Binary output included as a fall-back for output writers that don't implement all features.
For example, the .IMG writer is not designed for 3D data, so it uses this class' WriteBlob method.

This is meant primarily as a base class, not for actual use...
 */
class CBinaryOutput : public IDataWriter
{
public:
  CBinaryOutput();
  ~CBinaryOutput();
  virtual void WriteRealVolume(float_tt *data, std::vector<ulong> shape, std::string label, 
                               std::vector<ulong> indices=std::vector<ulong>(), std::string comment=std::string(),
                               std::map<std::string, double> parameters=std::map<std::string, double>(),
                               std::vector<float_tt> resolution=std::vector<float_tt>());
  virtual void WriteComplexVolume(complex_tt *data, std::vector<ulong> shape, std::string label, 
                           std::vector<ulong> indices=std::vector<ulong>(), std::string comment=std::string(),
                           std::map<std::string, double> parameters=std::map<std::string, double>(),
                           std::vector<float_tt> resolution=std::vector<float_tt>());
  virtual void WriteRealImage(float_tt **data, std::vector<ulong> shape, std::string label, 
                              std::vector<ulong> indices=std::vector<ulong>(), std::string comment=std::string(),
                              std::map<std::string, double> parameters=std::map<std::string, double>(),
                              std::vector<float_tt> resolution=std::vector<float_tt>());
  virtual void WriteComplexImage(complex_tt **data, std::vector<ulong> shape, std::string label, 
                                 std::vector<ulong> indices=std::vector<ulong>(), std::string comment=std::string(),
                                 std::map<std::string, double> parameters=std::map<std::string, double>(),
                                 std::vector<float_tt> resolution=std::vector<float_tt>());

  /*
  virtual void WriteRealImage(QSfMat data, std::string label, std::vector<ulong> indices=std::vector<ulong>(), 
                      std::map<std::string, double> parameters=std::map<std::string, double>());
  virtual void WriteComplexImage(QScMat data, std::string label, std::vector<ulong> indices=std::vector<ulong>(), 
                         std::map<std::string, double> parameters=std::map<std::string, double>());
  */
protected:
  virtual void DescribeFile(std::vector<ulong> shape, ulong element_size, std::string label, 
                            std::vector<ulong> indices=std::vector<ulong>(), std::string comment="", 
                            std::map<std::string, double> parameters=std::map<std::string, double>(),
                            std::vector<float_tt> resolution=std::vector<float_tt>());

  template <typename dtype>
  void WriteBlob(dtype *data, std::vector<ulong> shape, std::string label, 
                   std::vector<ulong> indices=std::vector<ulong>(), 
                   std::map<std::string, double> parameters=std::map<std::string, double>())
    {
      std::stringstream filename;
      filename<<label;
      for (ulong idx=0; idx<indices.size(); idx++)
        {
          filename<<"_"<<idx;
        }
      filename<<".bin";
      std::fstream file(filename.str().c_str(), std::ios::out|std::ios::binary);
      if(!file.is_open()) 
        {
          throw std::runtime_error("WriteBlob: Could not open file for writing.");
        }

      ulong ndata=1;
      for (int dim=0; dim<shape.size(); dim++) ndata*=shape[dim];
      file.write(reinterpret_cast<const char*>(data),ndata*sizeof(dtype));
      file.close();
    }
};
#endif










