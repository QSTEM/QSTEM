#ifndef IMG_OUTPUT_H
#define IMG_OUTPUT_H

#include "binary_output.hpp"

class CImgOutput : public CBinaryOutput
{
  int m_headerSize;
  int m_version;
public:
  CImgOutput();
  ~CImgOutput();
  // WriteVolume isn't implemented here - it uses the CBinaryOutput method instead.
  virtual void WriteRealImage(float_tt **data, std::vector<unsigned long> shape, std::string label, 
                              std::vector<unsigned long> position=std::vector<unsigned long>(), std::string comment=std::string(),
                              std::map<std::string, double> parameters=std::map<std::string, double>(),
                              std::vector<float_tt> resolution=std::vector<float_tt>());
  virtual void WriteComplexImage(complex_tt **data, std::vector<unsigned long> shape, std::string label, 
                                 std::vector<unsigned long> position=std::vector<unsigned long>(), std::string comment=std::string(),
                                 std::map<std::string, double> parameters=std::map<std::string, double>(),
                                 std::vector<float_tt> resolution=std::vector<float_tt>());
private:
  void WriteData(void **pix, bool is_complex, unsigned long dataSize, std::vector<unsigned long> shape, 
            std::string label, std::vector<unsigned long> position, std::string comment,
            std::map<std::string, double> parameters, std::vector<float_tt> resolution);
  /*
  void WriteRealImage(QSfMat &data, std::string label, std::vector<unsigned long> &position=std::vector<unsigned long>(), 
                      std::map<std::string, double> &parameters=std::map<std::string, double>());
  void WriteComplexImage(QScMat &data, std::string label, std::vector<unsigned long> &position=std::vector<unsigned long>(), 
                         std::map<std::string, double> &parameters=std::map<std::string, double>());
  */
};

#endif










