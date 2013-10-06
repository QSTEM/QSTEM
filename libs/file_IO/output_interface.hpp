#ifndef OUTPUT_INTERFACE_H
#define OUTPUT_INTERFACE_H

#include <vector>
#include <map>
#include "stemtypes_fftw3.hpp"
#include "boost/shared_ptr.hpp"

class IDataWriter
{
public:
  virtual void WriteRealVolume(float_tt *data, std::vector<ulong> shape, std::string label, 
                               std::vector<ulong> position=std::vector<ulong>(), std::string comment=std::string(),
                               std::map<std::string, double> parameters=std::map<std::string, double>(),
                               std::vector<float_tt> resolution=std::vector<float_tt>())=0;
  virtual void WriteComplexVolume(complex_tt *data, std::vector<ulong> shape, std::string label, 
                                  std::vector<ulong> position=std::vector<ulong>(), std::string comment=std::string(),
                                  std::map<std::string, double> parameters=std::map<std::string, double>(),
                                  std::vector<float_tt> resolution=std::vector<float_tt>())=0;
  virtual void WriteRealImage(float_tt **data, std::vector<ulong> shape, std::string label, 
                              std::vector<ulong> position=std::vector<ulong>(), std::string comment=std::string(),
                              std::map<std::string, double> parameters=std::map<std::string, double>(),
                              std::vector<float_tt> resolution=std::vector<float_tt>())=0;
  virtual void WriteComplexImage(complex_tt **data, std::vector<ulong> shape, std::string label, 
                                 std::vector<ulong> position=std::vector<ulong>(), std::string comment=std::string(),
                                 std::map<std::string, double> parameters=std::map<std::string, double>(),
                                 std::vector<float_tt> resolution=std::vector<float_tt>())=0;
  /*
  virtual void WriteRealFile(QSfMat &data, std::string label, std::vector<ulong> &position=std::vector<ulong>(), 
                             std::map<std::string, double> &parameters=std::map<std::string, double>())=0;
  virtual void WriteComplexFile(QScMat &data, std::string label, std::vector<ulong> &position=std::vector<ulong>(), 
                                std::map<std::string, double> &parameters=std::map<std::string, double>())=0;
  */
};

typedef boost::shared_ptr<IDataWriter> DataWriterPtr;

#endif
