#include "DataStream.hpp"

DataStream::DataStream(const std::string &oFile) { this->outputFile = oFile; }
void DataStream::open(const std::string &oFile, bool append) {
  // assert(header.size() > 0);
  assert(fieldDelimiter != "");
  assert(objDelimiter != "");
  outputFile = oFile;
  if (append)
    fout.open(outputFile, std::ios::out | std::ios::app);
  else
    fout.open(outputFile);
  if (header.size() == 0)
    return;
  unsigned i = 0;
  for (; i < header.size() - 1; ++i)
    fout << header[i] << fieldDelimiter;
  fout << header[i] << objDelimiter;
}

void DataStream::setFieldDelimiter(const std::string &fDelimiter) {
  fieldDelimiter = fDelimiter;
}
void DataStream::setObjDelimiter(const std::string &oDelimiter) {
  objDelimiter = oDelimiter;
}

void DataStream::setHeader(const std::vector<std::string> &h) { header = h; }

void DataStream::close() { fout.close(); }
