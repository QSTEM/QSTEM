#ifndef DATA_WRITERS_H
#define DATA_WRITERS_H

#include "output_interface.hpp"
#include "binary_output.hpp"
#include "img_output.hpp"

#include "input_interface.hpp"
#include "img_input.hpp"

#include <stdexcept>

DataWriterPtr GetDataWriter(std::string extension);
DataReaderPtr GetDataReader(std::string extension);

#endif
