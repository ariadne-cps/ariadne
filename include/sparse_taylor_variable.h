#include <map>

class MultiIndex
{
 public:
  MultiIndex(uint as) : _as(as) { _d=0; for(uint i=0; i!=as; ++i) { _o[i]=0; } }
  ushort argument_size() const { return this->_as; }
  ushort degree() const { return this->_d; }
  ushort operator[](uint i) const { return this->_o[i]; }
  friend MultiIndex& operator+=(MultiIndex& a, const MultiIndex& b));
  friend bool operator==(const MultiIndex& a, const MultiIndex& b) const;
  friend bool operator<(const MultiIndex& a, const MultiIndex& b) const;
 private:
  uchar _as;
  uchar _d;
  uchar _o[6];
};

MultiIndex& operator+=(MultiIndex& a, const MultiIndex& b) {
  assert(a._as==b._as);
  a._d+=b._d;
  for(uint i=0; i!=a._as; ++i) {
    a._o[i]+=b._o[i];
  }
}

MultiIndex operator+(const MultiIndex& a, const MultiIndex& b) {
  MultiIndex c(a); c+=b; return c;
}

bool operator==(const MultiIndex& a, const MultiIndex& b) {
  if(a.argument_size()!=b.argument_size()) { return false; }
  for(uint i=0; i!=a.argument_size(); ++i) {
    if(a[i]!=b[i]) { return false; } }
  return true;
}
  
bool operator<(const MultiIndex& a, const MultiIndex& b) {
  if(a.degree()==b.degree()) {
    for(uint i=0; i!=a.argument_size(); ++i) {
      if(a[i]!=b[i]) { return a[i]>b[i]; } 
    }
    return false;
  } else {
    return a.degree()<b.degree();
  }
}




template<class X>
class SparseTaylorVariable 
{
  Differential<X>& operator+=(const Differential<X>& x);
  Differential<X>& operator*=(const X& c);
  friend Differential<X> operator*(const Differential<X>& x, const Differential<X>& y);
  friend Differential<X> compose(const Differential<X>& x, const Differential<X>& y);
 private:
  uint argument_size _as;
  std::map<MultiIndex,X> _data;
};

template<class X>
Differential<X>& Differential<X>::operator+=(const Differential<X>& x)
{
  for(const_iterator iter=x._data.begin(); iter!=x._data.end(); ++iter) {
    this->_data[iter->first]+=iter->second;
  }
}

template<class X>
Differential<X>& Differential<X>::operator*=(const X& c)
{
  for(const_iterator iter=x._data.begin(); iter!=x._data.end(); ++iter) {
    this->_data[iter->first]*=c;
  }
}

template<class X>
Differential<X> operator*(const Differential<X>& x, const Differential<X>& y)
{
  for(const_iterator xiter=x._data.begin(); xiter!=x._data.end(); ++xiter) {
    for(const_iterator yiter=y._data.begin(); yiter!=y._data.end(); ++yiter) {
      this->_data[xiter->first+yiter->first]+=xiter->second*xiter->second;
  }
}
