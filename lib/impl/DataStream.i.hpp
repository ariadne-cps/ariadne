<<<<<<< HEAD
template <class T> void DataStream::append(const std::vector<T> &object) {
  assert(fout);
  unsigned i = 0;
  for (; i < object.size() - 1; ++i)
    fout << object[i] << fieldDelimiter;
  fout << object[i] << objDelimiter;
}
namespace csv //< namespace used to create csv datastream automatically
{
static DataStream make_datastream() {
  DataStream ds;
  ds.setFieldDelimiter(",");
  ds.setObjDelimiter("\n");
  return ds;
=======
template <class T> void DataStream::append(const std::vector<T> &object)
{
    assert(fout);
    unsigned i = 0;
    for (; i < object.size() - 1; ++i)
        fout << object[i] << fieldDelimiter;
    fout << object[i] << objDelimiter;
}
namespace csv //< namespace used to create csv datastream automatically
{
static DataStream make_datastream()
{
    DataStream ds;
    ds.setFieldDelimiter(",");
    ds.setObjDelimiter("\n");
    return ds;
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
}
} //  namespace csv
