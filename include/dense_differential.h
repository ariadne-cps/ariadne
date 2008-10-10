#ifndef ARIADNE_DIFFERENTIAL_H
#define ARIADNE_DIFFERENTIAL_H

#include <gmpxx.h>
#include <cmath>
#include <limits>

#include "macros.h"
#include "array.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "series.h"
#include "multi_index.h"

namespace Ariadne {

uint fac(uint n);
uint bin(uint n, uint k);

class MultiIndex;
template<class X> class Array;
template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Series;
template<class X> class DenseDifferential;



/// The partial derivatives of a variable with respect to other variables.
template<class X>
class DenseDifferential
{
 public:
  typedef X ScalarType;

  /// Default constructor constructs a constant of degree zero.
  DenseDifferential();
  /// The constant zero of degree \a d in \a a arguments.
  DenseDifferential(uint a, uint d);
  /// The constant \a c of degree \a d in \a a arguments.
  DenseDifferential(uint a, uint d, X c);
  /// The \a i<sup>th</sup> variable \a x of degree \a d in \a a arguments.
  DenseDifferential(uint a, uint d, X x, uint i);
  /// A taylor variable of degree \a d in \a arguments, with values given by the array based at \a ptr.
  template<class XX> DenseDifferential(uint a, uint d, const XX* ptr);
  /// A taylor variable of degree \a d in \a arguments, with values given by the array based at \a ptr.
  template<class XX> DenseDifferential(const DenseDifferential<XX>& x);
  
  /// Assign from a constant.
  DenseDifferential<X>& operator=(const X& c);

  /// Equality operator.
  bool operator==(const DenseDifferential<X>& other) const;
  /// Inequality operator.
  bool operator!=(const DenseDifferential<X>& other) const;
  
  /// The number of variables of the argument.
  uint argument_size() const; 
  /// The degree (number of derivatives computed).
  uint degree() const; 
  /// The value of the quantity.
  const X& value() const;
  /// A reference to the value of the quantity.
  X& value();
  /// Set the value of the quantity.
  void set_value(const X&);
  /// The variation of the quantity with respect to the \a j<sup>th</sup> argument.
  const X& gradient(uint j) const;
  /// The array of derivative values.
  const Vector<X>& data() const;
  /// A reference to the array of derivative values.
  Vector<X>& data();
  /// A reference to the \a i<sup> th</sup> derivative \f$D^af=d^{|a|}f/dx_1^{a_1}\cdots dx_n^{a_n}\f$.
  X& operator[](const MultiIndex& a); 
  X& operator[](const uint& j);
  /// The \a i<sup> th</sup> derivative \f$D^af=d^{|a|}f/dx_1^{a_1}\cdots dx_n^{a_n}\f$.
  const X& operator[](const MultiIndex& a) const; 
  const X& operator[](const uint& j) const;
  
  /// Assign all elements of degree less than the degree of \a x to those of \a x .
  DenseDifferential<X>& assign(const DenseDifferential<X>& x);

  /// Add another variable.
  DenseDifferential<X>& operator+=(const DenseDifferential<X>& x);
  /// Subtract another variable.
  DenseDifferential<X>& operator-=(const DenseDifferential<X>& x);
  /// Add a constant.
  DenseDifferential<X>& operator+=(const X& c);
  /// Subtract a constant.
  DenseDifferential<X>& operator-=(const X& c);
  /// Multiply by a constant.
  DenseDifferential<X>& operator*=(const X& c);
  /// Divide by a constant.
  DenseDifferential<X>& operator/=(const X& c);

  /// The power series of the reciprocal function.
  static DenseDifferential<X> constant(uint as, uint d, const X& x);
  /// The power series of the reciprocal function.
  static DenseDifferential<X> variable(uint as, uint d, const X& x, uint i);

 public:
#ifdef DOXYGEN
  /// 
  friend bool operator<(const DenseDifferential<X>& x1, const DenseDifferential<X>& x2);

  ///
  friend DenseDifferential<X> operator+(const DenseDifferential<X>& x);
  ///
  friend DenseDifferential<X> operator-(const DenseDifferential<X>& x);
  ///
  friend DenseDifferential<X> operator+(const DenseDifferential<X>& x, const DenseDifferential<X>& y);
  ///
  friend DenseDifferential<X> operator-(const DenseDifferential<X>& x, const DenseDifferential<X>& y);
  ///
  friend DenseDifferential<X> operator*(const DenseDifferential<X>& x, const DenseDifferential<X>& y);
  ///
  friend DenseDifferential<X> operator/(const DenseDifferential<X>& x, const DenseDifferential<X>& y);
  ///
  friend DenseDifferential<X> min(const DenseDifferential<X>& x1, const DenseDifferential<X>& x2); 
  ///
  friend DenseDifferential<X> max(const DenseDifferential<X>& x1, const DenseDifferential<X>& x2); 
  ///
  friend DenseDifferential<X> abs(const DenseDifferential<X>& x);

  /// Reciprocal
  friend DenseDifferential<X> rec(const DenseDifferential<X>& x);
  /// Power
  friend DenseDifferential<X> pow(const DenseDifferential<X>& x, int k);
  /// Square root
  friend DenseDifferential<X> sqrt(const DenseDifferential<X>& x);
  /// Exponential
  friend DenseDifferential<X> exp(const DenseDifferential<X>& x); 
  /// Natural logarithm
  friend DenseDifferential<X> log(const DenseDifferential<X>& x); 
#endif
 public:
  uint _argument_size;
  uint _degree;
  Vector<X> _data;
};

template<class X> DenseDifferential<X> scalar_constant(uint as, uint d, const X& c);
template<class X> DenseDifferential<X> scalar_variable(uint as, uint d, const X& x, uint i);


/// The derivatives of a vector quantity with respect to independent variables.
template<class X> 
class DenseDifferentialVector
  : public Vector< DenseDifferential<X> >
{
  typedef Vector< DenseDifferential<X> > Data;
 public:
  DenseDifferentialVector()
    : Data(0,DenseDifferential<X>()) { } 
  DenseDifferentialVector(uint rs, uint as, uint d)
    : Data(rs,DenseDifferential<X>(as,d)) { } 
  template<class XX> DenseDifferentialVector(const Vector< DenseDifferential<XX> >& x) 
    : Data(x) { for(uint i=1; i!=x.size(); ++i) { assert(x[i].argument_size()==x[0].argument_size()); } }
  template<class XX> DenseDifferentialVector(uint rs, uint as, uint d, const XX* ptr)
    : Data(rs,DenseDifferential<X>(as,d)) 
  {
    for(uint i=0; i!=(*this).size(); ++i) {
      for(uint j=0; j!=(*this)[i]._data.size(); ++j) {
        (*this)[i]._data[j]=*ptr;
        ++ptr;
      }
    }
  } 

  DenseDifferentialVector(Vector<X> v, Matrix<X> A, uint d)
    : Data(A.row_size(),DenseDifferential<X>(A.column_size(),d))
  {
    for(uint i=0; i!=v.size(); ++i) {
      (*this)[i]._data[0]=v[i];
      for(uint j=0; j!=A.column_size(); ++j) {
        (*this)[i]._data[j+1]=A[i][j];
      } 
    }
  }

  uint result_size() const { return this->Data::size(); }
  uint argument_size() const { return (*this)[0].argument_size(); }
  uint degree() const { return (*this)[0].degree(); }

  Vector<X> value() const { 
    Vector<X> r(this->result_size()); 
    for(uint i=0; i!=r.size(); ++i) { 
      r[i]=(*this)[i].value(); } 
    return r; 
  }
  Matrix<X> jacobian() const { 
    Matrix<X> J(this->result_size(),this->argument_size()); 
    for(uint i=0; i!=J.row_size(); ++i) { 
      for(uint j=0; j!=J.column_size(); ++j) { 
        J[i][j]=(*this)[i]._data[j+1]; } }
    return J; 
  }

  void set_value(const Vector<X>& v) { 
    for(uint i=0; i!=v.size(); ++i) { 
      (*this)[i]._data[0]=v[i]; }
  }

  void set_jacobian(const Matrix<X>& J) { 
    for(uint i=0; i!=J.row_size(); ++i) { 
       for(uint j=0; j!=J.column_size(); ++j) { 
         (*this)[i]._data[j+1]=J[i][j]; } }
  }

};

template<class X, class Y> Y evaluate(const DenseDifferential<X>& y, const Vector<Y>& z);
template<class X> DenseDifferential<X> compose(const Series<X>& y, const DenseDifferential<X>& x);

template<class X> bool operator<(const DenseDifferential<X>& x1, const DenseDifferential<X>& x2);
template<class X> DenseDifferential<X> operator+(const DenseDifferential<X>& x);
template<class X> DenseDifferential<X> operator-(const DenseDifferential<X>& x);
template<class X> DenseDifferential<X> operator+(const DenseDifferential<X>& x, const DenseDifferential<X>& y);
template<class X> DenseDifferential<X> operator-(const DenseDifferential<X>& x, const DenseDifferential<X>& y);
template<class X> DenseDifferential<X> operator*(const DenseDifferential<X>& x, const DenseDifferential<X>& y);
template<class X> DenseDifferential<X> operator/(const DenseDifferential<X>& x, const DenseDifferential<X>& y);

template<class X> DenseDifferential<X> min(const DenseDifferential<X>& x1, const DenseDifferential<X>& x2); 
template<class X> DenseDifferential<X> max(const DenseDifferential<X>& x1, const DenseDifferential<X>& x2); 
template<class X> DenseDifferential<X> abs(const DenseDifferential<X>& x);

template<class X> DenseDifferential<X> rec(const DenseDifferential<X>& x);
template<class X> DenseDifferential<X> pow(const DenseDifferential<X>& x, int k);
template<class X> DenseDifferential<X> sqrt(const DenseDifferential<X>& x);
template<class X> DenseDifferential<X> exp(const DenseDifferential<X>& x); 
template<class X> DenseDifferential<X> log(const DenseDifferential<X>& x); 

template<class X, class XX> DenseDifferential<X> operator+(const XX& c, const DenseDifferential<X>& x) { 
  return DenseDifferential<X>(x)+=X(c); }
template<class X, class XX> DenseDifferential<X> operator+(const DenseDifferential<X>& x, const XX& c) { 
  return DenseDifferential<X>(x)+=X(c); }
template<class X, class XX> DenseDifferential<X> operator-(const XX& c, const DenseDifferential<X>& x) { 
  return DenseDifferential<X>(-x)+=X(c); }
template<class X, class XX> DenseDifferential<X> operator-(const DenseDifferential<X>& x, const XX& c) { 
  return DenseDifferential<X>(x)-=X(c); }
template<class X, class XX> DenseDifferential<X> operator*(const DenseDifferential<X>& x, const XX& c) { 
  return DenseDifferential<X>(x)*=c; }
template<class X, class XX> DenseDifferential<X> operator*(const XX& c, const DenseDifferential<X>& x) { 
  return DenseDifferential<X>(x)*=X(c); }
template<class X, class XX> DenseDifferential<X> operator/(const DenseDifferential<X>& x, const XX& c) { 
  return DenseDifferential<X>(x)*=(1.0/c); }

template<class X> std::ostream& operator<<(std::ostream& os, const DenseDifferential<X>& x);

//template<class X> DenseDifferential<X>& acc(DenseDifferential<X>& r, const DenseDifferential<X>& x, const DenseDifferential<X>& y);
//template<class X> DenseDifferential<X>& acc(DenseDifferential<X>& r, const X& c, const DenseDifferential<X>& x);




inline uint compute_polynomial_data_size(uint rs, uint as, uint d) { return rs*Ariadne::bin(d+as,as); }


template<class X> inline DenseDifferentialVector<X> differential_vector(uint rs, uint as, uint d, const X* ptr) {
  return DenseDifferentialVector<X>(rs,as,d,ptr); }

template<class X> DenseDifferentialVector<X> vector_variable(uint n, uint d, const Vector<X>& c);

template<class X> DenseDifferentialVector<X> project(const DenseDifferentialVector<X>& y, Slice rng);
template<class X> DenseDifferentialVector<X> restrict(const DenseDifferentialVector<X>& y, const Array<uint>& p);
template<class X> DenseDifferentialVector<X> expand(const DenseDifferentialVector<X>& y, const Array<uint>& p);
template<class X> DenseDifferentialVector<X> join(const DenseDifferentialVector<X>& x1, const DenseDifferentialVector<X>& x2);
template<class X> DenseDifferentialVector<X> join(const DenseDifferentialVector<X>& x1, const DenseDifferential<X>& x2);

template<class X> DenseDifferentialVector<X> operator-(const DenseDifferentialVector<X>& x, const DenseDifferentialVector<X>& y);
template<class X> DenseDifferentialVector<X> operator+(const DenseDifferentialVector<X>& v, const Vector<X>& c);
template<class X> DenseDifferentialVector<X> operator-(const DenseDifferentialVector<X>& v, const Vector<X>& c);

template<class X, class Y> Vector<Y> evaluate(const DenseDifferentialVector<X>& x, const Vector<Y>& y);

template<class X> DenseDifferentialVector<X> translate(const DenseDifferentialVector<X>& y, const Vector<X>& c);
template<class X> DenseDifferentialVector<X> compose(const DenseDifferentialVector<X>& y, const DenseDifferentialVector<X>& z);
template<class X> DenseDifferentialVector<X> inverse(const DenseDifferentialVector<X>& y);
template<class X> DenseDifferentialVector<X> implicit(const DenseDifferentialVector<X>& y);
template<class X> DenseDifferentialVector<X> derivative(const DenseDifferentialVector<X>& x, uint i);
template<class X> DenseDifferentialVector<X> antiderivative(const DenseDifferentialVector<X>& x, uint j);
//template<class X> DenseDifferentialVector<X> flow(const DenseDifferentialVector<X>& vf, const Vector<X>& x, uint ox);





template<class X> template<class XX> 
DenseDifferential<X>::DenseDifferential(uint a, uint d, const XX* ptr)
  : _argument_size(a), _degree(d), _data(compute_polynomial_data_size(1,a,d)) 
{
  for(uint i=0; i!=this->_data.size(); ++i) {
    this->_data[i]=ptr[i];
  }
}

template<class X> template<class XX> 
DenseDifferential<X>::DenseDifferential(const DenseDifferential<XX>& x)
  : _argument_size(x.argument_size()), _degree(x.degree()), _data(x.data()) 
{
}


template<class X>
DenseDifferential<X>::DenseDifferential(uint as, uint d, X x)
  : _argument_size(as), _degree(d), _data(compute_polynomial_data_size(1,as,d),0.0)
{
  _data[0]=x;
}


template<class X>
DenseDifferential<X>::DenseDifferential(uint as, uint d, X x, uint i)
  : _argument_size(as), _degree(d), _data(compute_polynomial_data_size(1,as,d),0.0)
{
  _data[0]=x;
  _data[i+1]=1.0;
}


template<class X>
DenseDifferential<X>& 
DenseDifferential<X>::operator=(const X& c) 
{
  this->_data[0]=c;
  for(uint i=1; i!=this->_data.size(); ++i) {
    this->_data[i]=0;
  }
  return *this;
}


template<class X>
bool 
DenseDifferential<X>::operator==(const DenseDifferential<X>& other) const
{
  return this->_argument_size==other._argument_size
    && this->_degree==other._degree
    && this->_data==other._data;
}



template<class X>
bool 
DenseDifferential<X>::operator!=(const DenseDifferential<X>& other) const
{
  return !(*this==other); 
}



template<class X>
X
evaluate(const DenseDifferential<X>& y, const Vector<X>& x) 
{
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  ARIADNE_ASSERT(y.argument_size()==x.size());
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;
  uint d=y.degree();
  uint ms=x.size();
  ARIADNE_ASSERT(d>=1);

  X zero = x[0]; zero*=0;
  X one = zero; one+=1;

  // Use inefficient brute-force approach with lots of storage...
  array< array< X > > val(ms, array< X >(d+1));
  for(uint j=0; j!=ms; ++j) {
    val[j][0]=one;
    val[j][1]=x[j];
    for(uint k=2; k<=d; ++k) {
      val[j][k]=val[j][k-1]*x[j];
    }
  }

  X r(zero);
  for(MultiIndex j(ms); j.degree()<=d; ++j) {
    X sf=fac(j);
    X t=one;
    for(uint k=0; k!=ms; ++k) {
      t=t*val[k][j[k]];
    }
    t*=X(y[j]/sf);
    r+=t;
  }
  return r;
}


template<class X, class Y>
Y
evaluate(const DenseDifferential<X>& y, const Vector<Y>& x) 
{
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  ARIADNE_ASSERT(y.argument_size()==x.size());
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;
  uint d=y.degree();
  uint ms=x.size();
  ARIADNE_ASSERT(d>=1);

  Y zero = x[0]; zero*=0;
  Y one = zero; one+=1;

  // Use inefficient brute-force approach with lots of storage...
  array< array< Y > > val(ms, array< Y >(d+1));
  for(uint j=0; j!=ms; ++j) {
    val[j][0]=one;
    val[j][1]=x[j];
    for(uint k=2; k<=d; ++k) {
      val[j][k]=val[j][k-1]*x[j];
    }
  }

  Y r(zero);
  for(MultiIndex j(ms); j.degree()<=d; ++j) {
    Y sf=fac(j);
    Y t=one;
    for(uint k=0; k!=ms; ++k) {
      t=t*val[k][j[k]];
    }
    t*=Y(y[j]/sf);
    r+=t;
  }
  return r;
}





template<class X>
DenseDifferential<X> 
min(const DenseDifferential<X>& x1, const DenseDifferential<X>& x2) 
{
  if(x1.value()==x2.value()) {
    ARIADNE_THROW(std::runtime_error,"min(Taylor<X> x1, Taylor<X> x2)","x1[0]==x2[0]");
  }
  return x1.value()<x2.value() ? x1 : x2;
}

  
template<class X>
DenseDifferential<X> 
max(const DenseDifferential<X>& x1,const DenseDifferential<X>& x2) 
{
  if(x1.value()==x2.value()) { 
    ARIADNE_THROW(std::runtime_error,"max(Taylor<X> x1, Taylor<X> x2)","x1[0]==x2[0]"); 
  }
  return x1.value()>x2.value() ? x1 : x2;
}

 
template<class X>
DenseDifferential<X> 
pos(const DenseDifferential<X>& x)
{
  return x;
}

 
template<class X>
DenseDifferential<X> 
neg(const DenseDifferential<X>& x)
{
  DenseDifferential<X> y(x.argument_size(),x.degree());
  for(uint n=0; n<y.data().size(); ++n) {
    y.data()[n] = -x.data()[n];
  }
  return y;
}

  
template<class X>
DenseDifferential<X> 
abs(const DenseDifferential<X>& x) 
{
  if(x.value()==0) { 
    ARIADNE_THROW(std::runtime_error,"abs(Taylor<X> x)","x[0]==0"); 
  }
  return x.value()>0 ? pos(x) : neg(x); 
}

 
template<class X>
DenseDifferential<X> 
add(const DenseDifferential<X>& x, const DenseDifferential<X>& y)
{
  assert(x.argument_size()==y.argument_size());
  DenseDifferential<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  for(uint n=0; n<z.data().size(); ++n) {
    z.data()[n] = x.data()[n]+y.data()[n];
  }
  return z;
}

 
template<class X>
DenseDifferential<X> 
sub(const DenseDifferential<X>& x, const DenseDifferential<X>& y)
{
  assert(x.argument_size()==y.argument_size());
  DenseDifferential<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  for(uint n=0; n<z.data().size(); ++n) {
    z.data()[n] = x.data()[n]-y.data()[n];
  }
  return z;
}

 
template<class X>
DenseDifferential<X> 
mul(const DenseDifferential<X>& x, const DenseDifferential<X>& y)
{
  assert(x.argument_size()==y.argument_size());
  DenseDifferential<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  acc(z,x,y);
  return z;
}

 
template<class X>
DenseDifferential<X> 
div(const DenseDifferential<X>& x, const DenseDifferential<X>& y)
{
  return mul(x,rec(y));
}

template<class X>
DenseDifferential<X> 
pow(const DenseDifferential<X>& x, int k)
{
  return compose(Series<X>::pow(x.degree(),x.value(),k),x);
}

 



template<class X>
DenseDifferential<X> 
rec(const DenseDifferential<X>& x)
{
  return compose(Series<X>::rec(x.degree(),x.value()),x);
}

  
template<class X>
DenseDifferential<X> 
sqrt(const DenseDifferential<X>& x) 
{
  return compose(Series<X>::sqrt(x.degree(),x.value()),x);
}

  
template<class X>
DenseDifferential<X> 
exp(const DenseDifferential<X>& x) 
{
  return compose(Series<X>::exp(x.degree(),x.value()),x);
}

  
template<class X>
DenseDifferential<X> 
log(const DenseDifferential<X>& x) 
{
  return compose(Series<X>::log(x.degree(),x.value()),x);
}

 
template<class X>
DenseDifferential<X> 
operator+(const DenseDifferential<X>& x)
{
  return pos(x);
}

 
template<class X>
DenseDifferential<X> 
operator-(const DenseDifferential<X>& x)
{
  return neg(x);
}

 
template<class X>
DenseDifferential<X> 
operator+(const DenseDifferential<X>& x, const DenseDifferential<X>& y)
{
  return add(x,y);
}

 
template<class X>
DenseDifferential<X> 
operator-(const DenseDifferential<X>& x, const DenseDifferential<X>& y)
{
  return sub(x,y);
}

 
template<class X>
DenseDifferential<X> 
operator*(const DenseDifferential<X>& x, const DenseDifferential<X>& y)
{
  return mul(x,y);
}

 
template<class X>
DenseDifferential<X> 
operator/(const DenseDifferential<X>& x, const DenseDifferential<X>& y)
{
  return div(x,y);
}







template<class X>
DenseDifferential<X>&
DenseDifferential<X>::operator+=(const DenseDifferential<X>& x)
{
  this->_data+=x._data;
  return *this;
}

template<class X>
DenseDifferential<X>&
DenseDifferential<X>::operator-=(const DenseDifferential<X>& x)
{
  this->_data-=x._data;
  return *this;
}


template<class X>
DenseDifferential<X>&
DenseDifferential<X>::operator+=(const X& c)
{
  this->_data[0]+=c;
  return *this;
}

template<class X>
DenseDifferential<X>&
DenseDifferential<X>::operator-=(const X& c)
{
  this->_data[0]-=c;
  return *this;
}

template<class X>
DenseDifferential<X>&
DenseDifferential<X>::operator*=(const X& c)
{
  this->_data*=c; 
  return *this;
}

template<class X>
DenseDifferential<X>&
DenseDifferential<X>::operator/=(const X& c)
{
  this->_data/=c; 
  return *this;
}

template<class X>
DenseDifferential<X>::DenseDifferential()
  : _argument_size(1), _degree(0), _data(1u,X(0)) 
{
}

 
template<class X>
DenseDifferential<X>::DenseDifferential(uint a, uint d)
  : _argument_size(a), _degree(d), _data(compute_polynomial_data_size(1,a,d),X(0))
{
}

 
template<class X>
DenseDifferential<X>
DenseDifferential<X>::constant(uint as, uint d, const X& c)  
{ 
  DenseDifferential<X> r(as,d); r._data[0]=c; return r;
}

 
template<class X>
DenseDifferential<X>
DenseDifferential<X>::variable(uint as, uint d, const X& x, uint i)  
{ 
  DenseDifferential<X> r(as,d); r._data[0]=x; r._data[1+i]=1; return r;
}

 
template<class X>
uint 
DenseDifferential<X>::argument_size() const 
{ 
  return this->_argument_size;
}

 
template<class X>
uint
DenseDifferential<X>::degree() const 
{ 
  return this->_degree;
}

 
template<class X>
const X&
DenseDifferential<X>::value() const 
{ 
  return this->_data[0];
}

 
template<class X>
X&
DenseDifferential<X>::value()  
{ 
  return this->_data[0];
}

template<class X>
void
DenseDifferential<X>::set_value(const X& x)  
{ 
  this->_data[0]=x;
}
 
template<class X>
const X&
DenseDifferential<X>::gradient(uint j) const 
{ 
  return this->_data[j+1u];
}

 
template<class X>
Vector<X>& 
DenseDifferential<X>::data()
{
  return this->_data; 
}

 
template<class X>
const Vector<X>& 
DenseDifferential<X>::data() const 
{
  return this->_data; 
}


 
template<class X>
X& 
DenseDifferential<X>::operator[](const MultiIndex& a) 
{ 
  return this->_data[a.position()]; 
}

 
template<class X>
const X& 
DenseDifferential<X>::operator[](const MultiIndex& a) const 
{ 
  return this->_data[a.position()]; 
}

template<class X>
X& 
DenseDifferential<X>::operator[](const uint& i) 
{ 
  ARIADNE_ASSERT(i<this->argument_size());
  return this->_data[i+1];
}

template<class X>
const X& 
DenseDifferential<X>::operator[](const uint& i) const 
{ 
  ARIADNE_ASSERT(i<this->argument_size());
  return this->_data[i+1];
}

 
template<class X>
bool
operator<(const DenseDifferential<X>& x1, const DenseDifferential<X>& x2)
{
  return x1.value() < x2.value();
}

 
template<class X>
DenseDifferential<X>&
acc(DenseDifferential<X>& r, const DenseDifferential<X>& x1, const DenseDifferential<X>& x2)
{
  ARIADNE_ASSERT(r.argument_size()==x1.argument_size());
  ARIADNE_ASSERT(r.argument_size()==x2.argument_size());
  for(MultiIndex i1(x1.argument_size()); i1.degree() <= std::min(r.degree(),x1.degree()); ++i1) {
    for(MultiIndex i2(x2.argument_size()); i2.degree() <= std::min(x2.degree(),uint(r.degree()-i1.degree())); ++i2) {
      MultiIndex i0=i1+i2;
      r[i0]+=x1[i1]*x2[i2];
    }
  }
  return r;
}

 
template<class X>
DenseDifferential<X>&
acc(DenseDifferential<X>& r, const X& c, const DenseDifferential<X>& x)
{
  ARIADNE_ASSERT(r.argument_size()==x.argument_size());
  uint n=std::max(r.data().size(),x.data().size());
  for(uint i=0; i!=n; ++i) {
    r.data()[i]+=c*x.data()[i];
  }
  return r;
}


 
template<class X>
DenseDifferential<X>&
DenseDifferential<X>::assign(const DenseDifferential<X>& x)
{
  ARIADNE_ASSERT(this->argument_size()==x.argument_size());
  ARIADNE_ASSERT(this->degree()>=x.degree());
  for(uint i=0; i!=x.data().size(); ++i) {
    this->_data[i]=x.data()[i]; 
  }
  return *this;
}


 
template<class X>
DenseDifferential<X>&
operator+=(DenseDifferential<X>& r, const DenseDifferential<X>& x)
{
  ARIADNE_ASSERT(r.argument_size()==x.argument_size());
  ARIADNE_ASSERT(r.degree()==x.degree());
  for(uint i=0; i!=r.data().size(); ++i) {
    r.data()[i]+=x.data()[i];
  }
  //reinterpret_cast<LinearAlgebra::Vector&>(r._data)
  //  += reinterpret_cast<LinearAlgebra::Vectorconst&>(x._data);
  return r;
}


template<class X>
void 
compute_composition(DenseDifferential<X>& z, 
                    const Series<X>& y, 
                    const DenseDifferential<X>& x)
{
  uint as=x.argument_size();
  uint d=z.degree();

  DenseDifferential<X> w=x;
  w.value()=0;
  DenseDifferential<X> t(as,d);
  t.value()=y[d];
  for(uint n=1; n<=d; ++n) {
    DenseDifferential<X> u(as,d);
    acc(u,t,w);
    t=u; t+=y[d-n];
  }
  z=t;
  return;
}


template<class X>
DenseDifferential<X> 
compose(const Series<X>& y, const DenseDifferential<X>& x)
{
  DenseDifferential<X> z(x.argument_size(),std::min(x.degree(),y.degree()));
  compute_composition(z,y,x);
  return z;
}


template<class X>
DenseDifferential<X> 
reduce(const DenseDifferential<X>& x, const uint& d)
{
  assert(x.degree()>=d);
  DenseDifferential<X> r(x.argument_size(),d);
  for(MultiIndex i(x.argument_size()); i.degree() <= x.degree(); ++i) {
    r[i]=x[i];
  }
  return r;
}


template<class X>
DenseDifferential<X> scalar_constant(uint as,uint d,const X& x) 
{
  return DenseDifferential<X>::constant(as,d,x);
}

template<class X>
DenseDifferential<X> scalar_variable(uint as,uint d,const X& x, uint i) 
{
  return DenseDifferential<X>::variable(as,d,x,i);
}





 
template<class X>
std::ostream& 
operator<<(std::ostream& os, const DenseDifferential<X>& x) {
  //  return os << "Taylor<X>( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
  //os << "Taylor<X>(";
  os << "D";
  uint degree=0;
  for(MultiIndex i(x.argument_size()); i.degree()<=x.degree(); ++i) {
    if(i.degree()==0) {
      os << '[';
    } else if(i.degree()==degree) {
      os << ',';
    } else {
      degree=i.degree();
      os << ';';
    }
    os << x[i];
  }
  os << ']';
  //os << ")";
  return os;

//  return os << "Taylor<X>( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
}









template<class X>
void 
compute_composition(DenseDifferentialVector<X>& z, 
                    const DenseDifferentialVector<X>& y, 
                    const DenseDifferentialVector<X>& x)
{
  ARIADNE_ASSERT(y.argument_size()==x.result_size());
  //std::cerr<<"y="<<y<<std::endl;
  //std::cerr<<"x="<<x<<std::endl;

  uint d=std::min(x.degree(),y.degree());
  uint rs=y.result_size();
  uint ms=x.result_size();
  uint as=x.argument_size();
  
  Vector<X> xv=x.value();
  DenseDifferentialVector<X>& xnc=const_cast< DenseDifferentialVector<X>& >(x);
  for(uint i=0; i!=x.result_size(); ++i) { xnc[i]._data[0]=0; }

  DenseDifferentialVector<X> r(rs,as,d);
  DenseDifferential<X> t(as,d);


  // Use inefficient brute-force approach with lots of storage...
  array< array< DenseDifferential<X> > > val(ms, array< DenseDifferential<X> >(d+1));
  for(uint j=0; j!=ms; ++j) {
    val[j][0]=DenseDifferential<X>(as,d,1.0);
    for(uint k=1; k<=d; ++k) {
      val[j][k]=val[j][k-1]*x[j];
    }
  }

  for(MultiIndex j(ms); j.degree()<=d; ++j) {
    t=DenseDifferential<X>(as,d,1.0);
    for(uint k=0; k!=ms; ++k) {
      t=t*val[k][j[k]];
    }
    for(uint i=0; i!=rs; ++i) {
      r[i]+=y[i][j]*t;
    }
  }


  z=r;
 
  for(uint i=0; i!=x.result_size(); ++i) { xnc[i]._data[0]=xv[i]; }
  //std::cerr<<"z="<<z<<std::endl;
}


 



















 
template<class X>
DenseDifferentialVector<X> 
compose(const DenseDifferentialVector<X>& y, const DenseDifferentialVector<X>& x)
{
  assert(x.result_size()==y.argument_size());
  DenseDifferentialVector<X> z(y.result_size(),x.argument_size(),std::min(x.degree(),y.degree()));
  compute_composition(z,y,x);
  return z;
}
 

template<class X>
X pow(const Vector<X>& x, const MultiIndex& a)
{
  X r=1;
  for(uint i=0; i!=x.size(); ++i) {
    r*=pow(x[i],a[i]);
  }
  return r;
}

template<class X>
DenseDifferentialVector<X> 
operator-(const DenseDifferentialVector<X>& x, const DenseDifferentialVector<X>& y)
{
  return Vector< DenseDifferential<X> >(static_cast<const Vector< DenseDifferential<X> >&>(x)-static_cast<const Vector<DenseDifferential<X> >&>(y));
}

template<class X>
DenseDifferentialVector<X> 
operator+(const DenseDifferentialVector<X>& x, const Vector<X>& c)
{
  DenseDifferentialVector<X> r(x); r+=c; return r;
  //return DenseDifferentialVector<X>(x)+=c;
}

template<class X>
DenseDifferentialVector<X> 
operator-(const DenseDifferentialVector<X>& x, const Vector<X>& c)
{
  DenseDifferentialVector<X> r(x); r-=c; return r;
}

template<class X>
DenseDifferentialVector<X> 
operator*(const Matrix<X>& A, const DenseDifferentialVector<X>& x)
{
  DenseDifferentialVector<X> r(A.row_size(), x.argument_size(), x.degree());
  for(uint i=0; i!=A.row_size(); ++i) {
    for(uint j=0; j!=A.column_size(); ++j) {
      if(A[i][j]!=0) {
        r[i]+=A[i][j]*x[j];
      }
    }
  }
  return r;
}


template<class X>
DenseDifferentialVector<X> 
vector_variable(uint n, uint d, const Vector<X>& x) 
{
  assert(x.size()==n);
  DenseDifferentialVector<X> r(n,n,d);
  for(uint i=0; i!=n; ++i) {
    r[i].data()[0]=x[i];
    r[i].data()[1+i]=1.0;
  }
  return r;
}


template<class X>
Vector<X> 
evaluate(const DenseDifferentialVector<X>& x, const Vector<X>& y)
{
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  Vector<X> r(x.result_size());
  for(MultiIndex j(x.argument_size()); j.degree()<=x.degree(); ++j) {
    X s=pow(y,j);
    for(uint i=0; i!=x.result_size(); ++i) {
      r[i]+=s*x[i][j];
    }
  }
  return r;
}
            
template<class X, class Y>
Vector<Y> 
evaluate(const DenseDifferentialVector<X>& x, const Vector<Y>& y)
{
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  //std::cerr<<"x="<<x<<"\ny="<<y<<std::endl;
  assert(x.argument_size()==y.size());

  Y zero=y[0]; zero*=0; 
  Y one=zero; one+=1.0;
  Vector<Y> r(x.result_size(),zero);
  for(MultiIndex j(x.argument_size()); j.degree()<=x.degree(); ++j) {
    //Y s=pow(y,j);
    Y s=one;
    for(uint k=0; k!=j.size(); ++k) {
      for(uint l=0; l!=j[k]; ++l) {
        s*=y[k]; 
      }
    }
    for(uint i=0; i!=x.result_size(); ++i) {
      Y sxij=s; 
      sxij*=x[i][j];
      r[i]+=sxij;
    }
  }
  //std::cerr<<" r="<<r<<std::endl;
  return r;
}
 

template<class X>
DenseDifferentialVector<X> 
project(const DenseDifferentialVector<X>& x, ublas::range rng)
{
  DenseDifferentialVector<X> r(rng.size(),x.argument_size(),x.degree());
  for(uint i=0; i!=rng.size(); ++i) {
    r[i]=x[i+rng.start()];
  }
  return r;
}


template<class X>
DenseDifferentialVector<X> 
translate(const DenseDifferentialVector<X>& x, const Vector<X>& c)
{
  assert(x.argument_size()==c.size());
  DenseDifferentialVector<X> r(x);
  for(uint k=0; k!=x.argument_size(); ++k) {
    if(c[k]!=0) { 
      DenseDifferentialVector<X> y(r);
      r*=0.0;
      MultiIndex e(y.argument_size());
      e.set(k,1);
      for(MultiIndex j(y.argument_size()); j.degree()<=y.degree(); ++j) {
        MultiIndex tj=j;
        for(uint l=0; l<=y.degree()-j.degree(); ++l) {
          assert(tj.degree()<=y.degree());
          for(uint i=0; i!=x.result_size(); ++i) {
            r[i][j]+=y[i][tj]*pow(c[k],l)*bin(tj[k],l);
          }
          tj+=e;
        }
      }
    }
  }
  return r;
}
 

template<class X>
DenseDifferentialVector<X> 
restrict(const DenseDifferentialVector<X>& x, const array<uint>& p)
{
  uint d=x.degree();
  uint rs=x.result_size();
  DenseDifferentialVector<X> r(x.result_size(),p.size(),x.degree());
  MultiIndex rj(r.argument_size());
  MultiIndex xj(x.argument_size());
  for( ; rj.degree()<=d; ++rj) {
    for(uint k=0; k!=p.size(); ++k) {
      xj.set(p[k],rj[k]);
    }
    for(uint i=0; i!=rs; ++i) {
      r[i][rj]=x[i][xj];
    }
  }
  return r;
}

template<class X>
DenseDifferentialVector<X> 
expand(const DenseDifferentialVector<X>& x, uint as, const array<uint>& p)
{
  uint d=x.degree();
  uint rs=x.result_size();
  DenseDifferentialVector<X> r(rs,as,d);
  MultiIndex rj(r.argument_size());
  MultiIndex xj(x.argument_size());
  for( ; xj.degree()<=d; ++xj) {
    for(uint k=0; k!=p.size(); ++k) {
      rj.set(p[k],xj[k]);
    }
    for(uint i=0; i!=rs; ++i) {
      r[i][rj]=x[i][xj];
    }
  }
  return r;
}


template<class X>
DenseDifferentialVector<X> 
join(const DenseDifferentialVector<X>& x1, const DenseDifferentialVector<X>& x2)
{
  assert(x1.argument_size()==x2.argument_size());
  assert(x1.degree()==x2.degree());
  DenseDifferentialVector<X> r(x1.result_size()+x2.argument_size(),x1.argument_size(),std::min(x1.degree(),x2.degree()));
  for(uint i=0; i!=x1.result_size(); ++i) {
    r[i]=x1[i];
  }
  for(uint i=0; i!=x2.result_size(); ++i) {
    r[i+x1.result_size()]=x2[i];
  }
  return r;
}


template<class X>
DenseDifferentialVector<X> 
join(const DenseDifferentialVector<X>& x1, const DenseDifferential<X>& x2)
{
  assert(x1.argument_size()==x2.argument_size());
  assert(x1.degree()==x2.degree());
  DenseDifferentialVector<X> r(x1.result_size()+1u,x1.argument_size(),std::min(x1.degree(),x2.degree()));
  for(uint i=0; i!=x1.result_size(); ++i) {
    r[i]=x1[i];
  }
  r[x1.result_size()]=x2;
  return r;
}


template<class X>
DenseDifferentialVector<X> 
derivative(const DenseDifferentialVector<X>& x, uint k)
{
  DenseDifferential<X> r(x.result_size(),x.argument_size(),x.degree()-1);
  MultiIndex e(x.argument_size());
  e.set(k,1);
  for(MultiIndex j(r.argument_size()); j.degree() <= r.degree(); ++j) {
    for(uint i=0; i!=r.result_size(); ++i) {
      r[i][j]=(j[k]+1)*x[i][j+e];
    }
  }
  return r;
}

template<class X>
DenseDifferentialVector<X>
antiderivative(const DenseDifferentialVector<X>& x, uint k)
{
  DenseDifferentialVector<X> r(x.result_size(),x.argument_size(),x.degree()+1);
  MultiIndex e(x.argument_size());
  e.set(k,1);
  for(MultiIndex j(r.argument_size()); j.degree() <= r.degree(); ++j) {
    for(uint i=0; i!=r.result_size(); ++i) {
      r[i][j+e]=x[i][j]/(j[k]+1);
    }
  }
  return r;
}

template<class X> 
DenseDifferentialVector<X> 
inverse(const DenseDifferentialVector<X>& x)
{
  assert(x.result_size()==x.argument_size());
  uint n=x.result_size();
  uint d=x.degree();
  Vector<X> z(n,0);
  Matrix<X> J=inverse(x.jacobian());
  Matrix<X> I=identity_matrix<X>(n);

  DenseDifferentialVector<X> y(n,n,d);
  DenseDifferentialVector<X> id(z,I,d);
  
  y.set_jacobian(J);
  for(uint i=2; i<=d; ++i) {
    y=y-J*(compose(x,y)-id);
  }
  return y;
}




// Computes the vector such that x(w,y(w))=c; y(0)=0
template<class X> 
DenseDifferentialVector<X> 
implicit(const DenseDifferentialVector<X>& x)
{
  assert(x.result_size()<=x.argument_size());
  //std::cerr<<"  x="<<x<<std::endl;

  uint rs=x.result_size();
  uint as=x.argument_size();
  uint d=x.degree();

  Matrix<X> A1=x.jacobian();
  Matrix<X> A2=project(A1,range(0,rs),range(as-rs,as));
  
  Matrix<X> J(as,rs);
  project(J,range(as-rs,as),range(0,rs))=inverse(A2);
  //std::cerr<<"  J="<<J<<std::endl;

  DenseDifferentialVector<X> y(as,as-rs,d);
  for(uint i=0; i!=as-rs; ++i) {
    y[i]=DenseDifferential<X>::variable(as-rs,d,0.0,i);
  }
  for(uint i=0; i!=d; ++i) {
    //std::cerr<<"  y="<<y<<std::endl;
    DenseDifferentialVector<X> z=compose(x,y);
    //std::cerr<<"  z="<<z<<std::endl;
    y=y-J*z;
  }
  //std::cerr<<"  y="<<y<<std::endl;

  DenseDifferentialVector<X> r(rs,as-rs,d);
  for(uint i=0; i!=rs; ++i) {
    r[i]=y[as-rs+i];
  }
  //std::cerr<<"  r="<<r<<std::endl;
  //DenseDifferentialVector<X> r=Vector< DenseDifferential<X> >(project(y,range(as-rs,as)));
  return r;
}





template<class X> 
Vector< Series< DenseDifferential<X> > >
flow(const DenseDifferentialVector<X>& vf, const Vector<X>& x, uint ox)
{
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  ARIADNE_ASSERT(vf.result_size()==vf.argument_size());
  ARIADNE_ASSERT(vf.argument_size()==x.dimension());
  uint n=x.dimension();
  uint ot=vf.degree();
  Vector< Series< DenseDifferential<X> > > y(n);
  Vector< Series< DenseDifferential<X> > > yp(n);
  for(uint i=0; i!=n; ++i) {
    y[i]=Series< DenseDifferential<X> >(0);
    y[i][0]=DenseDifferential<X>(n,ox,x[i],i);
  }
  for(uint j=0; j<ot; ++j) {
    yp=evaluate(vf,y);
    for(uint i=0; i!=n; ++i) {  
      y[i]=antiderivative(yp[i],y[i][0]);
    }
  } 
  return y;
}

} // namespace Ariadne




#endif /* ARIADNE_DIFFERENTIAL_H */

