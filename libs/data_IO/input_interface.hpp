/*
  QSTEM - image simulation for TEM/STEM/CBED
  Copyright (C) 2000-2010  Christoph Koch
  Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef INPUT_INTERFACE_H
#define INPUT_INTERFACE_H

#include <map>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

class IDataReader;
typedef boost::shared_ptr<IDataReader> DataReaderPtr;
typedef DataReaderPtr (*CreateDataReaderFn)();

// Inline definitions are here to allow specification of a filename without a position associated
class IDataReader
{
public:
  virtual void ReadImageData(const std::string &filename, const std::vector<unsigned> &position,
                             void *pix)=0;
  inline void ReadImageData(const std::string &filename, void *pix)
  {
    std::vector<unsigned> position;
    return ReadImageData(filename, position, pix);
  }
  inline void ReadImageData(const std::string &filename, unsigned position, void *pix)
  {
    std::vector<unsigned> posvec(1);
	posvec[0]=position;
    return ReadImageData(filename, posvec, pix);
  }


  virtual void ReadComment(const std::string &filename, const std::vector<unsigned> &position, 
                           std::string &comment)=0;
  inline void ReadComment(const std::string &filename, std::string &comment)
  {
    std::vector<unsigned> position;
    return ReadComment(filename, position, comment);
  }
inline void ReadComment(const std::string &filename, unsigned position, std::string &comment)
  {
    std::vector<unsigned> posvec(1);
	posvec[0]=position;
    return ReadComment(filename, posvec, comment);
  }


  virtual void ReadParameters(const std::string &filename,  const std::vector<unsigned> &position,
                              std::map<std::string, double> &parameters)=0;
  inline void ReadParameters(const std::string &filename, std::map<std::string, double> &parameters)
  {
    std::vector<unsigned> posvec;
    return ReadParameters(filename, posvec, parameters);
  }
inline void ReadParameters(const std::string &filename, unsigned position, std::map<std::string, double> &parameters)
  {
    std::vector<unsigned> posvec(1);
	posvec[0]=position;
    return ReadParameters(filename, posvec, parameters);
  }


  virtual void ReadSize(const std::string &filename, const std::vector<unsigned> &position,
                        unsigned &nx, unsigned &ny)=0;
  inline void ReadSize(const std::string &filename, unsigned &nx, unsigned &ny)
  {
    std::vector<unsigned> position;
    return ReadSize(filename, position, nx, ny);
  }
  inline void ReadSize(const std::string &filename, unsigned position, unsigned &nx, unsigned &ny)
  {
    std::vector<unsigned> posvec(1);
	posvec[0]=position;
    return ReadSize(filename, posvec, nx, ny);
  }

  virtual void ReadElementByteSize(const std::string &filename, const std::vector<unsigned> &position,
                                   unsigned &elementByteSize)=0;
  inline void ReadElementByteSize(const std::string &filename, unsigned &elementByteSize)
  {
    std::vector<unsigned> posvec;
    return ReadElementByteSize(filename, posvec, elementByteSize);
  }
  inline void ReadElementByteSize(const std::string &filename, unsigned position, unsigned &elementByteSize)
  {
    std::vector<unsigned> posvec(1);
	posvec[0]=position;
    return ReadElementByteSize(filename, posvec, elementByteSize);
  }


  virtual void ReadComplex(const std::string &filename, const std::vector<unsigned> &position, bool &complex)=0;
  inline void ReadComplex(const std::string &filename, bool &complex)
  {
    std::vector<unsigned> posvec;
    return ReadComplex(filename, posvec, complex);
  }
  inline void ReadComplex(const std::string &filename, unsigned position, bool &complex)
  {
    std::vector<unsigned> posvec(1);
	posvec[0]=position;
    return ReadComplex(filename, posvec, complex);
  }

  // One-shot catch-all
  virtual void ReadImage(const std::string &filename, const std::vector<unsigned> &position, void *pix, 
                         std::map<std::string, double> &params, std::string &comment)=0;
  inline void ReadImage(const std::string &filename, void *pix, std::map<std::string, double> &params,
                         std::string &comment)
  {
    std::vector<unsigned> posvec;
    return ReadImage(filename, posvec, pix, params, comment);
  }
  inline void ReadImage(const std::string &filename, unsigned position, void *pix, 
                        std::map<std::string, double> &params, std::string &comment)
  {
    std::vector<unsigned> posvec(1);
	posvec[0]=position;
    return ReadImage(filename, posvec, pix, params, comment);
  }
};


#endif













