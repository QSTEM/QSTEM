#include "img_output.hpp"

CImgOutput::CImgOutput() : CBinaryOutput(),
                           m_headerSize(56),
                           m_version(1)
{
}

CImgOutput::~CImgOutput()
{
}

void CImgOutput::WriteComplexImage(complex_tt **data, std::vector<unsigned long> shape, std::string label, 
                                     std::vector<unsigned long> position, std::string comment,
                                     std::map<std::string, double> parameters, std::vector<float_tt> resolution)
{
  unsigned long dataSize = 2*sizeof(float_tt);
  if (resolution == std::vector<float_tt>())
    {
      resolution = std::vector<float_tt>(shape.size(),1);
    }
  
  WriteData((void **)data, true, dataSize, shape, label, position, comment, parameters, resolution);
}

void CImgOutput::WriteRealImage(float_tt **data, std::vector<unsigned long> shape, std::string label, 
                                  std::vector<unsigned long> position, std::string comment,
                                  std::map<std::string, double> parameters, std::vector<float_tt> resolution)
{
  unsigned long dataSize = sizeof(float_tt);
  if (resolution == std::vector<float_tt>())
    {
      resolution = std::vector<float_tt>(shape.size(),1);
    }
  
  WriteData((void **)data, false, dataSize, shape, label, position, comment, parameters, resolution);
}

void CImgOutput::WriteData(void **pix, bool is_complex, unsigned long dataSize, std::vector<unsigned long> shape, 
                             std::string label, std::vector<unsigned long> position, std::string comment,
                             std::map<std::string, double> parameters, std::vector<float_tt> resolution)
{
  //FILE *fp;
  std::stringstream filename;
  char buf[200];
  filename<<label;
  for (unsigned long idx=0; idx<position.size(); idx++)
    {
      filename<<"_"<<idx;
    }
  filename<<".img";

  std::fstream file(filename.str().c_str(), std::ios::out|std::ios::binary);

  // Sychronize lengths of comments and parameters
  size_t paramSize = parameters.size();
  std::map<std::string, double>::iterator param;
  
  // TODO: need to parse parameters to pull out thickness, complexFlag, and dataSize.
  //    Hardcode version?  Things aren't going to change anymore...

  size_t commentSize = comment.size();
  double thickness = parameters["Thickness"];

  if(!file.is_open()) 
    {
      sprintf(buf,"WriteData: Could open file %s for writing\n",filename.str().c_str());
      throw std::runtime_error(buf);
    }

  file.write(reinterpret_cast <const char*> (&m_headerSize), 4);
  file.write(reinterpret_cast <const char*> (&paramSize), 4);
  file.write(reinterpret_cast <const char*> (&commentSize), 4);
  file.write(reinterpret_cast <const char*> (&shape[0]), 4);
  file.write(reinterpret_cast <const char*> (&shape[1]), 4);
  file.write(reinterpret_cast <const char*> ((int)is_complex), 4);
  file.write(reinterpret_cast <const char*> ((int)dataSize), 4);
  file.write(reinterpret_cast <const char*> (&m_version), 4);
  file.write(reinterpret_cast <const char*> (&thickness), 8);
  file.write(reinterpret_cast <const char*> (&resolution[0]), 8);
  file.write(reinterpret_cast <const char*> (&resolution[1]), 8);

  for (param=parameters.begin(); param!=parameters.end(); param++)
    {
      file.write(reinterpret_cast<const char*>(&(param->second)), sizeof(double));
    }
  
  file.write(comment.c_str(), commentSize);
  file.write(reinterpret_cast<char*>(pix[0]), shape[0]*shape[1]*dataSize);
  file.close();
}
