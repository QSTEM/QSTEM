#include "data_io_factories.hpp"

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

DataReaderPtr GetDataReader(std::string extension)
{
  if (extension==".img")
    {
      return DataReaderPtr(new CImgInput());
    }
  else
    {
      throw std::runtime_error("Unsupported file extension for data input");
    }
}
