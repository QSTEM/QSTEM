#include "img_output.hpp"

CImgOutput::CImgOutput() : CBinaryOutput(),
                           m_headerSize(56),
                           m_version(1)
{
}

void CImgOutput::WriteComplexImage(complex_tt **data, std::vector<ulong> shape, std::string label, 
                                     std::vector<ulong> indices, std::string comment,
                                     std::map<std::string, double> parameters, std::vector<float_tt> resolution)
{
  ulong dataSize = 2*sizeof(float_tt);
  if (resolution == std::vector<float_tt>())
    {
      resolution = std::vector<float_tt>(shape.size(),1);
    }
  
  WriteData((void **)data, true, shape, dataSize, label, indices, parameters, resolution);
}

void CImgOutput::WriteRealImage(float_tt **data, std::vector<ulong> shape, std::string label, 
                                  std::vector<ulong> indices, std::string comment,
                                  std::map<std::string, double> parameters, std::vector<float_tt> resolution)
{
  dataSize = sizeof(float_tt);
  if (resolution == std::vector<float_tt>())
    {
      resolution = std::vector<float_tt>(shape.size(),1);
    }
  
  WriteData((void **)data, false, shape, dataSize, label, indices, parameters, resolution);
}

void CImgOutput::WriteData(float_tt **pix, bool is_complex, ulong dataSize, std::vector<ulong> shape, 
                             std::string label, std::vector<ulong> indices, std::string comment,
                             std::map<std::string, double> parameters, std::vector<float_tt> resolution)
{
  //FILE *fp;
  std::stringstream filename;
  filename<<label;
  for (ulong idx=0; idx<indices.size(); idx++)
    {
      filename<<"_"<<idx;
    }
  filename<<".img";

  std::fstream file(filename, std::ios::out|std::ios::binary);

  // Sychronize lengths of comments and parameters
  size_t paramSize = parameters.size();
  std::map<std::string, double>::iterator param;
  
  // TODO: need to parse parameters to pull out thickness, complexFlag, and dataSize.
  //    Hardcode version?  Things aren't going to change anymore...

  size_t commentSize = comment.size();

  if(!file.is_open()) 
    {
      sprintf(m_buf,"WriteData: Could open file %s for writing\n",fileName);
      throw std::runtime_error(m_buf);
    }

  file.write(reinterpret_cast <const char*> (&m_headerSize), 4);
  file.write(reinterpret_cast <const char*> (&paramSize), 4);
  file.write(reinterpret_cast <const char*> (&commentSize), 4);
  file.write(reinterpret_cast <const char*> (&shape[0]), 4);
  file.write(reinterpret_cast <const char*> (&shape[1]), 4);
  file.write(reinterpret_cast <const char*> (&(int)complexFlag), 4);
  file.write(reinterpret_cast <const char*> (&(int)dataSize), 4);
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
