#ifndef DATA_WRITERS_H
#define DATA_WRITERS_H

#include "output_interface.hpp"

#include "input_interface.hpp"

namespace QSTEM
{

// Factory for creating instances of IDataReader
class QSTEM_HELPER_DLL_EXPORT CDataReaderFactory
{
public:
    ~CDataReaderFactory() { m_FactoryMap.clear(); }

    static CDataReaderFactory *Get()
    {
        static CDataReaderFactory instance;
        return &instance;
    }

  void Register(const std::string &extension, CreateDataReaderFn pfnCreate);
  // Looks up which reader to get based on string mapping of registered readers
  DataReaderPtr GetReader(const std::string &animalName);

private:
  CDataReaderFactory();
  CDataReaderFactory(const CDataReaderFactory &) { }
  CDataReaderFactory &operator=(const CDataReaderFactory &) { return *this; }

  typedef std::map<std::string, CreateDataReaderFn> FactoryMap;
  FactoryMap m_FactoryMap;
};



class QSTEM_HELPER_DLL_EXPORT CDataWriterFactory
{
public:
    ~CDataWriterFactory() { m_FactoryMap.clear(); }

    static CDataWriterFactory *Get()
    {
        static CDataWriterFactory instance;
        return &instance;
    }

  void Register(const std::string &extension, CreateDataWriterFn pfnCreate);
  // Looks up which reader to get based on string mapping of registered readers
  DataWriterPtr GetWriter(const std::string &animalName);

private:
  CDataWriterFactory();
  CDataWriterFactory(const CDataWriterFactory &) { }
  CDataWriterFactory &operator=(const CDataWriterFactory &) { return *this; }

  typedef std::map<std::string, CreateDataWriterFn> FactoryMap;
  FactoryMap m_FactoryMap;
};

}

#endif
