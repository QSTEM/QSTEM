#include "img_output.hpp"

CImgOutput::CImgOutput() : CBinaryOutput(),
                           m_headerSize(56),
                           m_version(1)
{
}

CImgOutput::~CImgOutput()
{
}

void CImgOutput::WriteComplexImage(complex_tt **data, std::vector<unsigned> shape, std::string label, 
                                 std::vector<unsigned> position, std::string comment,
                                 std::map<std::string, double> parameters)
{
  unsigned dataSize = 2*sizeof(float_tt);
  WriteData((void **)data, true, dataSize, shape, label, position, comment, parameters);
}

void CImgOutput::WriteRealImage(float_tt **data, std::vector<unsigned> shape, std::string label, 
                                  std::vector<unsigned> position, std::string comment,
                                  std::map<std::string, double> parameters)
{
  unsigned dataSize = sizeof(float_tt);
    WriteData((void **)data, false, dataSize, shape, label, position, comment, parameters);
}

void CImgOutput::WriteData(void **pix, bool is_complex, unsigned dataSize, std::vector<unsigned> shape, 
                             std::string label, std::vector<unsigned> position, std::string comment,
                             std::map<std::string, double> parameters)
{
  //FILE *fp;
  std::stringstream filename;
  char buf[200];
  filename<<label;
  for (unsigned idx=0; idx<position.size(); idx++)
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
  double resX = parameters["dx"];
  double resY = parameters["dy"];

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
  file.write(reinterpret_cast <const char*> (&resX), 8);
  file.write(reinterpret_cast <const char*> (&resY), 8);

  for (param=parameters.begin(); param!=parameters.end(); param++)
    {
      file.write(reinterpret_cast<const char*>(&(param->second)), sizeof(double));
    }
  
  file.write(comment.c_str(), commentSize);
  file.write(reinterpret_cast<char*>(pix[0]), shape[0]*shape[1]*dataSize);
  file.close();
}
