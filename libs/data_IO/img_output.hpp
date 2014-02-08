#ifndef IMG_OUTPUT_H
#define IMG_OUTPUT_H

#include "binary_output.hpp"

class CImgWriter : public CBinaryOutput
{
  int m_headerSize;
  int m_version;
public:
  CImgWriter();
  //~CImgWriter();
  virtual void Initialize(std::string dirOrFileName, std::string run_id="avg");
  // WriteVolume isn't implemented here - it uses the CBinaryWriter method instead.
  virtual void WriteRealImage(const RealVector &data, const std::vector<unsigned> &shape, const std::string &label, 
                                  const std::vector<unsigned> &position, const std::string &comment, 
								  const std::map<std::string, double> &parameters);
  virtual void WriteComplexImage(const ComplexVector &data, const std::vector<unsigned> &shape, const std::string &label, 
                                  const std::vector<unsigned> &position, const std::string &comment, 
								  const std::map<std::string, double> &parameters);
private:
  void WriteData(void *pix, bool is_complex, unsigned dataSize, const std::vector<unsigned> &shape, 
            const std::string &filebase, const std::vector<unsigned> &position, const std::string &comment,
            const std::map<std::string, double> &parameters);
  /*
  void WriteRealImage(QSfMat &data, std::string label, std::vector<unsigned> &position=std::vector<unsigned>(), 
                      std::map<std::string, double> &parameters=std::map<std::string, double>());
  void WriteComplexImage(QScMat &data, std::string label, std::vector<unsigned> &position=std::vector<unsigned>(), 
                         std::map<std::string, double> &parameters=std::map<std::string, double>());
  */
private:
  friend class CDataWriterFactory;
  // Create an instance of this class, wrapped in a shared ptr
  //     This should not be inherited - any subclass needs its own implementation.
  static DataWriterPtr Create(void){return DataWriterPtr(new CImgWriter());}
};

#endif










