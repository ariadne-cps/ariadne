#pragma once
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//! \ingroup InputOutputModule
//! Method to output data on file or other stream
//! Store an object using the field delimiter (data inside an object) and object
//! delimiter (delimiter between object); for example on csv file field
//! delimiter is ',', object delimiter is '\n'
<<<<<<< HEAD
class DataStream {
public:
  DataStream() = default;
  DataStream(const std::string &);

  //! Method to open a file
  //! @param oFile is the output file, default is /tmp/datastream.log
  //! @param append open file on-append mode (default is true)
  void open(const std::string &oFile = "/tmp/datastream.log",
            bool append = true);

  //! Parm field delimiter method
  void setFieldDelimiter(const std::string &);

  //! Parm object delimiter method
  void setObjDelimiter(const std::string &);

  //! Set the header on file (first line if csv file is considered)
  void setHeader(const std::vector<std::string> &);

  //! Close the file
  void close();

  //! Close the file but not write \n at end

  //! template method to append on a file the object. It must be defined field
  //! delimiter and object delimiter, otherwise an error is raised.
  template <class T> void append(const std::vector<T> &object);

  void write(const std::string &obj_str);

private:
  std::string outputFile;
  std::string fieldDelimiter;
  std::string objDelimiter;
  std::ofstream fout;
  std::vector<std::string> header;
=======
class DataStream
{
  public:
    DataStream() = default;
    DataStream(const std::string &);

    //! Method to open a file
    //! @param oFile is the output file, default is /tmp/datastream.log
    //! @param append open file on-append mode (default is true)
    void open(const std::string &oFile = "/tmp/datastream.log",
              bool append = true);

    //! Parm field delimiter method
    void setFieldDelimiter(const std::string &);

    //! Parm object delimiter method
    void setObjDelimiter(const std::string &);

    //! Set the header on file (first line if csv file is considered)
    void setHeader(const std::vector<std::string> &);

    //! Close the file
    void close();

    //! template method to append on a file the object. It must be defined field
    //! delimiter and object delimiter, otherwise an error is raised.
    template <class T> void append(const std::vector<T> &object);

    void write(const std::string &obj_str);

  private:
    std::string outputFile;
    std::string fieldDelimiter;
    std::string objDelimiter;
    std::ofstream fout;
    std::vector<std::string> header;
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
};

#include "impl/DataStream.i.hpp"
