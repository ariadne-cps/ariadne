

/* 
Compile with:
swig -python -c++ -I../include/ numeric.i
g++ -c -fpic numeric_wrap.cxx -I/usr/include/python2.4 -I../ -I../include/  
g++ -shared numeric_wrap.o  ../src/.libs/rational.o ../src/.libs/float64.o -lmpfr -lgmpxx -lgmp -o _numeric.so
*/

%module numeric 
%{
#include <sstream>
#include <string>
#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/float.h"
#include "numeric/interval.h"
#include "python/python_float.h"
using namespace Ariadne;
using Numeric::Integer;
using Numeric::Rational;
typedef Numeric::Interval<Float> Interval;
template<class T> const char* to_str(const T& t) {
  std::stringstream ss; ss << t; return ss.str().c_str(); }
%}

class Integer {
 public:
  Integer();
  Integer(int n);
  Integer(const Integer& n);
};
%extend Integer { 
  const char* __str__() { return to_str(*self); } 
  const char* __repr__() { return to_str(*self); } 
  Integer __add__(const Integer& n) { return *self + n; } 
  Integer __sub__(const Integer& n) { return *self - n; } 
  Integer __mul__(const Integer& n) { return *self * n; } 
  Integer __div__(const Integer& n) { return *self / n; } 
  Integer __add__(int n) { return *self + n; } 
  Integer __sub__(int n) { return *self - n; } 
  Integer __mul__(int n) { return *self * n; } 
  Integer __div__(int n) { return *self / n; } 
  Integer __radd__(int n) { return n + *self; } 
  Integer __rsub__(int n) { return n - *self; } 
  Integer __rmul__(int n) { return n * *self; } 
  Integer __rdiv__(int n) { return n / *self; } 
}

class Rational {
 public:
  Rational();
  Rational(int n);
  Rational(Integer n);
  Rational(int p,int q);
  Rational(Integer p,Integer q);
  Rational(const Rational& q);
};
%extend Rational { 
  const char* __str__() { return to_str(*self); } 
  const char* __repr__() { return to_str(*self); } 
  Rational __add__(const Rational& q) { return *self + q; } 
  Rational __sub__(const Rational& q) { return *self - q; } 
  Rational __mul__(const Rational& q) { return *self * q; } 
  Rational __div__(const Rational& q) { return *self / q; } 
  Rational __add__(const Integer& n) { return *self + n; } 
  Rational __sub__(const Integer& n) { return *self - n; } 
  Rational __mul__(const Integer& n) { return *self * n; } 
  Rational __div__(const Integer& n) { return *self / n; } 
  Rational __radd__(const Integer& n) { return n + *self; } 
  Rational __rsub__(const Integer& n) { return n - *self; } 
  Rational __rmul__(const Integer& n) { return n * *self; } 
  Rational __rdiv__(const Integer& n) { return n / *self; } 
  Rational __add__(int n) { return *self + n; } 
  Rational __sub__(int n) { return *self - n; } 
  Rational __mul__(int n) { return *self * n; } 
  Rational __div__(int n) { return *self / n; } 
  Rational __radd__(int n) { return n + *self; } 
  Rational __rsub__(int n) { return n - *self; } 
  Rational __rmul__(int n) { return n * *self; } 
  Rational __rdiv__(int n) { return n / *self; } 
  Rational __add__(double x) { return *self + x; } 
  Rational __sub__(double x) { return *self - x; } 
  Rational __mul__(double x) { return *self * x; } 
  Rational __div__(double x) { return *self / x; } 
  Rational __radd__(double x) { return x + *self; } 
  Rational __rsub__(double x) { return x - *self; } 
  Rational __rmul__(double x) { return x * *self; } 
  Rational __rdiv__(double x) { return x / *self; } 
}

class Float {
 public:
  Float();
  Float(int n);
  Float(double x);
  Float(const Float& q);
};


%extend Float { 
  const char* __str__() { return to_str(*self); } 
  const char* __repr__() { return to_str(*self); } 
  Interval __add__(const Float& x) { return *self + x; } 
  Interval __sub__(const Float& x) { return *self - x; } 
  Interval __mul__(const Float& x) { return *self * x; } 
  Interval __div__(const Float& x) { return *self / x; } 
  Interval __add__(double x) { return *self + x; } 
  Interval __sub__(double x) { return *self - x; } 
  Interval __mul__(double x) { return *self * x; } 
  Interval __div__(double x) { return *self / x; } 
  Interval __radd__(double x) { return x + *self; } 
  Interval __rsub__(double x) { return x - *self; } 
  Interval __rmul__(double x) { return x * *self; } 
  Interval __rdiv__(double x) { return x / *self; } 
  Interval __add__(int n) { return *self + n; } 
  Interval __sub__(int n) { return *self - n; } 
  Interval __mul__(int n) { return *self * n; } 
  Interval __div__(int n) { return *self / n; } 
  Interval __radd__(int n) { return n + *self; } 
  Interval __rsub__(int n) { return n - *self; } 
  Interval __rmul__(int n) { return n * *self; } 
  Interval __rdiv__(int n) { return n / *self; } 
  Rational __add__(const Rational& q) { return Rational(*self) + q; } 
  Rational __sub__(const Rational& q) { return Rational(*self) - q; } 
  Rational __mul__(const Rational& q) { return Rational(*self) * q; } 
  Rational __div__(const Rational& q) { return Rational(*self) / q; } 
  Rational __radd__(const Rational& q) { return q + Rational(*self); } 
  Rational __rsub__(const Rational& q) { return q - Rational(*self); } 
  Rational __rmul__(const Rational& q) { return q * Rational(*self); } 
  Rational __rdiv__(const Rational& q) { return q / Rational(*self); } 
}


class Interval {
 public:
  Interval();
  Interval(int n);
  Interval(double n);
  Interval(Float x);
  Interval(Rational q);
  Interval(double l,double u);
  Interval(Float l,Float u);
  Interval(const Interval& i);
};

%extend Interval { 
  const char* __str__() { return to_str(*self); } 
  const char* __repr__() { return to_str(*self); } 
  Interval __add__(const Interval& q) { return *self + q; } 
  Interval __sub__(const Interval& q) { return *self - q; } 
  Interval __mul__(const Interval& q) { return *self * q; } 
  Interval __div__(const Interval& q) { return *self / q; } 
  Interval __add__(const Float& x) { return *self + x; } 
  Interval __sub__(const Float& x) { return *self - x; } 
  Interval __mul__(const Float& x) { return *self * x; } 
  Interval __div__(const Float& x) { return *self / x; } 
  Interval __radd__(const Float& x) { return x + *self; } 
  Interval __rsub__(const Float& x) { return x - *self; } 
  Interval __rmul__(const Float& x) { return x * *self; } 
  Interval __rdiv__(const Float& x) { return x / *self; } 
  Interval __add__(int n) { return *self + n; } 
  Interval __sub__(int n) { return *self - n; } 
  Interval __mul__(int n) { return *self * n; } 
  Interval __div__(int n) { return *self / n; } 
  Interval __radd__(int n) { return Interval(n) + *self; } 
  Interval __rsub__(int n) { return Interval(n) - *self; } 
  Interval __rmul__(int n) { return Interval(n) * *self; } 
  Interval __rdiv__(int n) { return Interval(n) / *self; } 
  Interval __add__(double x) { return *self + x; } 
  Interval __sub__(double x) { return *self - x; } 
  Interval __mul__(double x) { return *self * x; } 
  Interval __div__(double x) { return *self / x; } 
  Interval __radd__(double x) { return Interval(x) + *self; } 
  Interval __rsub__(double x) { return Interval(x) - *self; } 
  Interval __rmul__(double x) { return Interval(x) * *self; } 
  Interval __rdiv__(double x) { return Interval(x) / *self; } 
}

