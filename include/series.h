#ifndef ARIADNE_SERIES_H
#define ARIADNE_SERIES_H

#include "array.h"

namespace Ariadne {


template<class X>
class Series
{
 public:
  explicit Series() : _data(1) { }
  explicit Series(uint d) : _data(d+1) { }
  uint degree() const { return this->_data.size()-1; }
  X& operator[](uint i) { return this->_data[i]; }
  const X& operator[](uint i) const { return this->_data[i]; }

  template<class R> Series<X>& operator+=(const Series<R>& y) { _data+=y._data; return *this; }
  template<class R> Series<X>& operator*=(const Series<R>& y) { *this=(*this)*y; return *this; }
  template<class R> Series<X>& operator+=(const R& c) { _data[0]+=c; return *this; }
  template<class R> Series<X>& operator*=(const R& c) { _data*=c; return *this; }

  static Series<X> rec(uint d, const X& x);
  static Series<X> pow(uint d, const X& x, int n);
  static Series<X> sqrt(uint d, const X& x);
  static Series<X> exp(uint d, const X& x);
  static Series<X> log(uint d, const X& x);
 private:
  array<X> _data;
};

template<class X>
Series<X> 
Series<X>::rec(uint d, const X& c) 
{
  Series<X> y(d);
  X mr = (-1)/c;
  for(uint i=0; i<=y.degree(); ++i) {
    y[i]=-Ariadne::pow(mr,i+1u);
  }
  return y;
}

template<class X>
Series<X> 
Series<X>::pow(uint d, const X& c, int k)
{

  uint n=k;
  Series<X> y(d);
  for(uint i=0; i<=std::min(uint(d),n); ++i) {
    uint j=n-i;
    y[i]=X(bin(n,j))*Ariadne::pow(c,j);
  }
  return y;
}

template<class X>
Series<X> 
Series<X>::sqrt(uint d, const X& c)
{
  Series<X> y(d);
  y[0]=Ariadne::sqrt(c);
  X mhr=-0.5/c;
  for(uint i=1; i<=y.degree(); ++i) {
    // Need to convert uint to int to prevent wraparound for 2*1u-3
    y[i]=((2*int(i)-3)*mhr)/i*y[i-1];
  }
  return y;
}

template<class X>
Series<X> 
Series<X>::exp(uint d, const X& c)
{
  Series<X> y(d);
  y[0]=Ariadne::exp(c);
  for(uint i=1; i<=y.degree(); ++i) {
    y[i]=y[i-1]/i;
  }
  return y;
}

template<class X>
Series<X> 
Series<X>::log(uint d, const X& c)
{
  Series<X> y(d);
  y[0]=Ariadne::log(c);
  X mr=(-1)/c;
  for(uint i=1; i<=y.degree();++i) {
    y[i]=-Ariadne::pow(mr,i)/i;
  }
  return y;
}


template<class X>
Series<X>
antiderivative(const Series<X>& x, const X& c)
{
  uint n=x.degree();
  Series<X> r(n+1);
  r[0]=c;
  for(uint i=0; i<=n; ++i) {
    r[i+1]=x[i]/(i+1);
  }
  return r;
}


template<class X>
Series<X> 
operator*(const Series<X>& x, const Series<X>& y) 
{
  //const Series<X>& x=*this;
  X zero=x[0]*y[0]*0.0;
  Series<X> r(std::min(x.degree(),y.degree()));
  for(uint k=0; k<=r.degree(); ++k) {
    r[k]=zero;
  }
  for(uint i=0; i<=r.degree(); ++i) {
    for(uint j=0; j<=r.degree()-i; ++j) {
      r[i+j]+=bin(i+j,i)*x[i]*y[j];
    }
  }
  return r;
}


template<class X>
std::ostream& operator<<(std::ostream& os, const Series<X>& x)
{
  os << "S";
  for(uint i=0; i<=x.degree(); ++i) {
    os << (i==0 ? '[' : ',') << x[i];
  }
  return os << "]";
}

} // namespace Ariadne


#endif /* ARIADNE_SERIES_H */

