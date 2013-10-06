#include "data_writers.hpp"

DataWriterPtr GetDataWriter(std::string extension)
{
  if (extension==".img")
    {
      return DataWriterPtr(new CImgOutput());
    }
  else
    {
      throw std::runtime_error("Unsupported file extension for data output");
    }
}
