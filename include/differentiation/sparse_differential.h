/***************************************************************************
 *            sparse_differential.h
 *
 *  Copyright 2008  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
/*! \file sparse_differential.h
 *  \brief Derivatives of scalar functions of many variables.
 */
 
#ifndef ARIADNE_SPARSE_DIFFERENTIAL_H
#define ARIADNE_SPARSE_DIFFERENTIAL_H

#include <cassert>
#include <map>

#include "base/array.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "differentiation/multi_index.h"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;

template<class X> class SparseSeries;
template<class X> class SparseDifferential;
template<class X> class SparseDifferentialVector;

template<class X> SparseDifferential<X> operator-(const SparseDifferential<X>& x);
template<class X> SparseDifferential<X> operator+(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
template<class X> SparseDifferential<X> operator-(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
template<class X> SparseDifferential<X> operator*(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
template<class X> SparseDifferential<X> operator/(const SparseDifferential<X>& x, const SparseDifferential<X>& y);

template<class X, class Y> Y evaluate(const SparseDifferential<X>& y, const Vector<Y>& z);
template<class X> SparseDifferential<X> compose(const SparseSeries<X>& x, const SparseDifferential<X>& y);
template<class X> SparseDifferential<X> derivative(const SparseDifferential<X>& x, uint i);
template<class X> SparseDifferential<X> antiderivative(const SparseDifferential<X>& x, uint i);



template<class X> SparseDifferentialVector<X> join(const SparseDifferentialVector<X>& x1, const SparseDifferentialVector<X>& x2);
template<class X> SparseDifferentialVector<X> join(const SparseDifferentialVector<X>& x1, const SparseDifferential<X>& x2);
template<class X> SparseDifferentialVector<X> project(const SparseDifferentialVector<X>& y, Range rng);

template<class X, class Y> Vector<Y> evaluate(const SparseDifferentialVector<X>& x, const Vector<Y>& y);
template<class X> SparseDifferentialVector<X> translate(const SparseDifferentialVector<X>& y, const Vector<X>& c);
template<class X> SparseDifferentialVector<X> compose(const SparseDifferentialVector<X>& x, const SparseDifferentialVector<X>& y);
template<class X> SparseDifferentialVector<X> inverse(const SparseDifferentialVector<X>& x);
template<class X> SparseDifferentialVector<X> implicit(const SparseDifferentialVector<X>& x);
template<class X> SparseDifferentialVector<X> derivative(const SparseDifferentialVector<X>& x, uint i);
template<class X> SparseDifferentialVector<X> antiderivative(const SparseDifferentialVector<X>& x, uint j);
template<class X> SparseDifferentialVector<X> flow(const SparseDifferentialVector<X>& vf, ushort ox);



template<class X>
class SparseSeries
{
 public:
  typedef typename std::map<uint,X>::const_iterator const_iterator;
  SparseSeries() : _data() { _data[0u]; }
  SparseSeries(uint d) : _data() { _data[d]; }
  const_iterator begin() const { return this->_data.begin(); }
  const_iterator end() const { return this->_data.end(); }
  ushort degree() const { return (--this->_data.end())->first; }
  X& operator[](const uint& a) { return this->_data[a]; }
  X operator[](const uint& a) const { const_iterator iter=this->_data.find(a); 
    if(iter==this->_data.end()) { return X(); } else { return iter->second; } }
  static SparseSeries<X> rec(ushort d, const X& c);
  static SparseSeries<X> pow(ushort d, const X& c, int k);
 private:
  std::map<uint,X> _data;
};

template<class X> SparseSeries<X> SparseSeries<X>::rec(ushort d, const X& c) {
  SparseSeries<X> y; X mr = (-1)/c; 
  for(uint i=0; i<=d; ++i) {
    y[i]=-::pow(mr,i+1u); }
  return y;
}

template<class X> SparseSeries<X> SparseSeries<X>::pow(ushort d, const X& c, int k) {
  uint n=k; SparseSeries<X> y;
  for(uint i=0; i<=std::min(uint(d),n); ++i) {
    uint j=n-i; y[i]=X(bin(n,j))*::pow(c,j); }
  return y;
}

template<class X> SparseSeries<X> antiderivative(const SparseSeries<X>& x) {
  SparseSeries<X> r(x.degree()+1);
  for(uint i=0; i!=x.degree(); ++i) { r[i+1]=x[i]/(i+1); }
  return r;
}

template<class X> std::ostream& operator<<(std::ostream& os, const SparseSeries<X>& x) {
  for(uint i=0; i<=x.degree(); ++i) { os << (i==0 ? "S[" : ",") << x[i] << std::flush; } return os << "]"; }






template<class X>
class SparseDifferential
{
 public:
  typedef typename std::map<MultiIndex,X>::iterator iterator;
  typedef typename std::map<MultiIndex,X>::const_iterator const_iterator;

  SparseDifferential() : _as(1), _deg(0), _data() { }
  SparseDifferential(uint as, uint deg) : _as(as), _deg(deg), _data() { _data[MultiIndex(as)]=0; }
  SparseDifferential(uint as, uint deg, const X& c) : _as(as), _deg(deg), _data() { _data[MultiIndex(as)]=c; }
  SparseDifferential(uint as, uint deg, const X& c, uint j) : _as(as), _deg(deg), _data() { _data[MultiIndex(as)]=c; _data[MultiIndex(as,j)]=1; }
  
  SparseDifferential<X>& operator=(const X& c) { this->_data.clear(); this->_data[MultiIndex(this->argument_size())]=c; return *this; }

  SparseDifferential<X>& operator+=(const SparseDifferential<X>& x);
  SparseDifferential<X>& operator-=(const SparseDifferential<X>& x);
  SparseDifferential<X>& operator+=(const X& c);
  SparseDifferential<X>& operator*=(const X& c);

  void set_degree(ushort d) { this->_deg = d; }

  X& operator[](const uint& j) { return this->_data[MultiIndex(this->_as,j)]; }
  X& operator[](const MultiIndex& a) { return this->_data[a]; }
  X& value() { return this->operator[](MultiIndex(this->_as)); }
  X& gradient(uint j) { return this->operator[](MultiIndex(this->_as,j)); }

  const_iterator begin() const { return this->_data.begin(); }
  const_iterator end() const { return this->_data.end(); }
  ushort argument_size() const { return this->_as; }
  ushort degree() const { return this->_deg; }
  X operator[](const MultiIndex& a) const { const_iterator iter=this->_data.find(a); 
    if(iter==this->_data.end()) { return X(0); } else { return iter->second; } }
  X value() const { return this->operator[](MultiIndex(this->_as)); }
  X gradient(uint j) const { return this->operator[](MultiIndex(this->_as,j)); }

  static SparseDifferential<X> constant(uint as, ushort d, const X& x) {
    return SparseDifferential<X>(as,d,x); }
  static SparseDifferential<X> variable(uint as, ushort d, const X& x, uint i) {
    return SparseDifferential<X>(as,d,x,i); }

  friend SparseDifferential<X> operator-<>(const SparseDifferential<X>& x);
  friend SparseDifferential<X> operator+<>(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
  friend SparseDifferential<X> operator-<>(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
  friend SparseDifferential<X> operator*<>(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
  friend SparseDifferential<X> operator/<>(const SparseDifferential<X>& x, const SparseDifferential<X>& y);
  friend SparseDifferential<X> compose<>(const SparseSeries<X>& x, const SparseDifferential<X>& y);
  friend SparseDifferential<X> derivative<>(const SparseDifferential<X>& x, uint i);
  friend SparseDifferential<X> antiderivative<>(const SparseDifferential<X>& x, uint i);
 public:
  void cleanup() { 
    for(typename std::map<MultiIndex,X>::iterator iter=this->_data.begin(); iter!=this->_data.end(); ++iter) { 
      if(iter->second==0) { this->_data.erase(iter); } } 
    this->_data[MultiIndex(this->_as)]; }
 private:
  uint _as;
  uint _deg;
  std::map<MultiIndex,X> _data;
};

template<class X>
SparseDifferential<X>& SparseDifferential<X>::operator+=(const SparseDifferential<X>& x)
{
  for(const_iterator iter=x._data.begin(); iter!=x._data.end(); ++iter) {
    this->_data[iter->first]+=iter->second;
  }
  return *this;
}

template<class X>
SparseDifferential<X>& SparseDifferential<X>::operator-=(const SparseDifferential<X>& x)
{
  for(const_iterator iter=x._data.begin(); iter!=x._data.end(); ++iter) {
    this->_data[iter->first]-=iter->second;
  }
  return *this;
}

template<class X>
SparseDifferential<X>& SparseDifferential<X>::operator+=(const X& c)
{
  this->_data[MultiIndex(this->_as)]+=c; return *this;
}

template<class X>
SparseDifferential<X>& SparseDifferential<X>::operator*=(const X& c)
{
  for(iterator iter=this->_data.begin(); iter!=this->_data.end(); ++iter) {
    iter->second*=c;
  }
  return *this;
}



template<class X>
SparseDifferential<X> operator-(const SparseDifferential<X>& x)
{
  SparseDifferential<X> r(x.argument_size(),x.degree()); r-=x; return r; 
}


template<class X>
SparseDifferential<X> operator+(const X& c, const SparseDifferential<X>& x)
{
  SparseDifferential<X> r(x); r+=c; return r; 
}

template<class X>
SparseDifferential<X> operator+(const SparseDifferential<X>& x, const X& c)
{
  SparseDifferential<X> r(x); r+=c; return r; 
}

template<class X>
SparseDifferential<X> operator-(const SparseDifferential<X>& x, const X& c)
{
  SparseDifferential<X> r(x); r+=(-c); return r; 
}

template<class X>
SparseDifferential<X> operator-(const X& c, const SparseDifferential<X>& x)
{
  SparseDifferential<X> r(-x); r+=(c); return r; 
}

template<class X>
SparseDifferential<X> operator*(const X& c, const SparseDifferential<X>& x)
{
  SparseDifferential<X> r(x); r*=c; return r; 
}

template<class X>
SparseDifferential<X> operator+(const SparseDifferential<X>& x, const SparseDifferential<X>& y)
{
  SparseDifferential<X> r(x); r+=y; r.cleanup(); return r; 
}

template<class X>
SparseDifferential<X> operator-(const SparseDifferential<X>& x, const SparseDifferential<X>& y)
{
  SparseDifferential<X> r(x); r-=y; r.cleanup(); return r; 
}


template<class X>
SparseDifferential<X> operator*(const SparseDifferential<X>& x, const SparseDifferential<X>& y)
{
  typedef typename SparseDifferential<X>::const_iterator const_iterator;
  assert(x._as==y._as);
  SparseDifferential<X> r(x._as,std::min(x._deg,y._deg));
  for(const_iterator xiter=x._data.begin(); xiter!=x._data.end(); ++xiter) {
    if(xiter->first.degree()>=r.degree()) { break; }
    for(const_iterator yiter=y._data.begin(); yiter!=y._data.end(); ++yiter) {
      if(xiter->first.degree()+yiter->first.degree()>=r.degree()) { break; }
      r._data[xiter->first+yiter->first]+=(xiter->second*yiter->second);
    }
  }
  r.cleanup();
  return r;
}

template<class X>
SparseDifferential<X> rec(const SparseDifferential<X>& x)
{
  return compose(SparseSeries<X>::rec(x.degree(),x.value()),x);
}

template<class X>
SparseDifferential<X> compose(const SparseSeries<X>& x, const SparseDifferential<X>& y)
{
  uint as=y.argument_size();
  uint d=std::min(x.degree(),y.degree());

  SparseDifferential<X> w=y;
  w[MultiIndex(as)]=0;
  SparseDifferential<X> r(as,d);
  r[MultiIndex(as)]=x[d];
  for(uint n=1; n<=d; ++n) {
    r=r*w; r+=x[d-n];
  }
  return r;
}


template<class X>
SparseDifferential<X> derivative(const SparseDifferential<X>& x, uint i)
{
  if(x.degree()==0) { return SparseDifferential<X>(x.argument_size(),0u); }
  SparseDifferential<X> r(x.argument_size(), x.degree()-1); 
  MultiIndex da(x.argument_size(),i); MultiIndex ai(x.argument_size(),i);
  for(typename SparseDifferential<X>::const_iterator iter=x.begin(); iter!=x.end(); ++iter) 
  {
    const MultiIndex& a=iter->first;
    if(a[i]!=0) { 
      da=a-ai;
      r[da]=x[a]*a[i];
    }
  }
  return r;
}

template<class X>
SparseDifferential<X> antiderivative(const SparseDifferential<X>& x, uint i)
{
  SparseDifferential<X> r(x.argument_size(), x.degree()+1); 
  MultiIndex da(x.argument_size()); MultiIndex ai(x.argument_size(),i);
  for(typename SparseDifferential<X>::const_iterator iter=x.begin(); iter!=x.end(); ++iter) 
  {
    const MultiIndex& a=iter->first;
    da=a+ai;
    r[da]=x[a]/da[i];
  }
  return r;
}


template<class X>
std::ostream& operator<<(std::ostream& os, const SparseDifferential<X>& x)
{
  for(typename SparseDifferential<X>::const_iterator iter=x.begin(); iter!=x.end(); ++iter) {
    if(iter==x.begin()) { os << "D{"; } else { os << ","; }
    os << "[" << iter->first << "]:" << iter->second ;
  }
  return os << "}";
}





template<class X>
class SparseDifferentialVector
  : public Vector< SparseDifferential<X> >
{
 public:
  SparseDifferentialVector() 
    : Vector< SparseDifferential<X> >(0,SparseDifferential<X>()) { }
  SparseDifferentialVector(uint rs, uint as, ushort d) 
    : Vector< SparseDifferential<X> >(rs,SparseDifferential<X>(as,d)) { }

  uint result_size() const { return this->Vector< SparseDifferential<X> >::size(); }
  uint argument_size() const { return (*this)[0].argument_size(); }
  uint degree() const { return (*this)[0].degree(); }

  Vector<X> value() const { 
    Vector<X> r(this->result_size()); for(uint i=0; i!=r.size(); ++i) { r[i]=(*this)[i].value(); } return r; }
  Matrix<X> jacobian() const { Matrix<X> r(this->result_size(),this->argument_size()); 
    for(uint i=0; i!=r.row_size(); ++i) { for(uint j=0; j!=r.column_size(); ++j) { r[i][j]=(*this)[i].gradient(j); } } return r; }

  static SparseDifferentialVector<X> constant(uint rs, uint as, ushort d, const Vector<X>& c) {
    ARIADNE_ASSERT(c.size()==rs);
    SparseDifferentialVector<X> result(rs,as,d);
    for(uint i=0; i!=rs; ++i) { result[i]=c[i]; }
    return result;
  }
};


template<class X>
SparseDifferentialVector<X> 
operator*(const Matrix<X>& A, const SparseDifferentialVector<X>& x)
{  
  assert(A.column_size()==x.result_size());
  SparseDifferentialVector<X> r(A.row_size(),x.argument_size(),x.degree());
  for(uint i=0; i!=A.row_size();++i) {
    for(uint j=0; j!=A.column_size();++j) {
      r[i]+=A[i][j]*x[j];
    }
  }
  return r;
}

template<class X>
SparseDifferentialVector<X> 
operator+(const SparseDifferentialVector<X>& x, const Vector<X>& c)
{  
  assert(x.result_size()==c.size());
  SparseDifferentialVector<X> r(x);
  for(uint i=0; i!=c.size();++i) {
    r[i]+=c[i];
  }
  return r;
}

  

template<class X>
SparseDifferentialVector<X> 
evaluate(const SparseDifferentialVector<X>& y, 
         const SparseDifferentialVector<X>& x)
{  
  using namespace std;
  assert(y.argument_size()==x.result_size());
  SparseDifferentialVector<X> z(y.result_size(),x.argument_size(),std::min(y.degree(),x.degree()));
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;

  ushort d=std::min(x.degree(),y.degree());
  uint rs=y.result_size();
  uint ms=x.result_size();
  uint as=x.argument_size();
  
  SparseDifferentialVector<X> r(rs,as,d);
  SparseDifferential<X> t(as,d);

  // Use inefficient brute-force approach with lots of storage...
  array< array< SparseDifferential<X> > > val(ms, array< SparseDifferential<X> >(d+1));
  for(uint j=0; j!=ms; ++j) {
    val[j][0]=SparseDifferential<X>(as,d,1.0);
    for(uint k=1; k<=d; ++k) {
      val[j][k]=val[j][k-1]*x[j];
    }
  }
  for(MultiIndex j(ms); j.degree()<=d; ++j) {
    t=SparseDifferential<X>(as,d,1.0);
    for(uint k=0; k!=ms; ++k) {
      t=t*val[k][j[k]];
    }
    for(uint i=0; i!=rs; ++i) {
      r[i]+=y[i][j]*t;
    }
  }
  
  for(uint i=0; i!=r.result_size(); ++i) { r[i].cleanup(); }
  return r;
}


template<class X>
SparseDifferentialVector<X> 
compose(const SparseDifferentialVector<X>& y, 
        const SparseDifferentialVector<X>& x)
{  
  SparseDifferentialVector<X> xnc(x); 
  for(uint i=0; i!=xnc.result_size(); ++i) { xnc[i].value()=0; }
  return evaluate(y,xnc);
}

template<class X> 
SparseDifferentialVector<X> 
implicit(const SparseDifferentialVector<X>& x, const Vector<X>& c)
{
  assert(x.result_size()<=x.argument_size());
  assert(c.size()==x.result_size());
  //std::cerr << "x=" << x << std::endl;
  
  uint rs=x.result_size();
  uint as=x.argument_size()-x.result_size();
  uint d=x.degree();

  Matrix<X> A1(as,rs);
  for(uint i=0; i!=as; ++i) {
    for(uint j=0; j!=rs; ++j) {
      A1(i,j)=x[i].gradient(j);
    }
  }
  
  Matrix<X> A2(rs,rs);
  for(uint i=0; i!=rs; ++i) {
    for(uint j=0; j!=rs; ++j) {
      A2(i,j)=x[i].gradient(as+j);
    }
  }
  
  Matrix<X> J(as+rs,rs);
  //  J(slice(as,rs),slice(0,rs))=inverse(A2);

  SparseDifferentialVector<X> y(as+rs,as,d);
  for(uint i=0; i!=as; ++i) {
    y[i]=SparseDifferential<X>(as,d,1.0,i);
  }
  for(uint i=0; i!=rs; ++i) {
    // y[as+i]=TaylorVariable<X>::constant(as,d,0.0);
  }

  for(uint i=0; i!=d; ++i) {
    SparseDifferentialVector<X> z=compose(x,y);
    y-=J*z;
  }

  SparseDifferentialVector<X> r(rs,as,d);
  for(uint i=0; i!=rs; ++i) {
    r[i]=y[as+i];
  }
  return r;
}



template<class X> 
SparseDifferentialVector<X> 
inverse(const SparseDifferentialVector<X>& x, const Vector<X>& c)
{
  using namespace std;
  assert(x.result_size()==x.argument_size());
  assert(x.result_size()==c.size());
  //std::cerr << "x=" << x << std::endl;
  uint n=x.result_size();
  uint d=x.degree();
  Vector<X> z(n,0);
  Matrix<X> J=inverse(x.jacobian());

  SparseDifferentialVector<X> y(n,n,d);
  SparseDifferentialVector<X> id(n,n,d);
  for(uint i=0; i!=n; ++i) { id[i][MultiIndex(n,i)]=1.0; }
  
  for(uint i=0; i!=n; ++i) { 
    y[i][MultiIndex(n)]=c[i]; 
    for(uint j=0; j!=n; ++j) { 
      y[i][MultiIndex(n,j)]=J[i][j];
    }
  }

  for(uint i=2; i<=d; ++i) {
    SparseDifferentialVector<X> z=compose(x,y);
    z-=id;
    y-=J*z;
  }
  return y;
}



template<class X> 
SparseDifferentialVector<X>
flow1(const SparseDifferentialVector<X>& f, const Vector<X>& x, uint to, uint so)
{
  //Change f to give a timed vector field with function derivative 1
  assert(f.result_size()==f.argument_size());
  assert(f.result_size()==x.size());
  uint n=x.size();
  uint d=f.degree();

  Matrix<X> Q(n+1,n); for(uint i=0; i!=n; ++i) { Q[i][i]=1; }
  
  Vector<X> tx(n+1);
  for(uint i=0; i!=n; ++i) { tx[i]=x[i]; } 
  tx[n]=0.0;

  SparseDifferentialVector<X> tf(n+1,n+1,d);
  MultiIndex ta(n+1);
  for(uint i=0; i!=n; ++i) {
    for(typename SparseDifferential<X>::const_iterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
      const MultiIndex& a=iter->first;
      for(uint k=0; k!=n; ++k) { ta.set(k,a[k]); }
      tf[i][ta]=iter->second;
    }
  }
  tf[n].value()=1;
  
  std::cerr << "f=" << f << std::endl;
  std::cerr << "tf=" << tf << std::endl << std::endl;

  SparseDifferentialVector<X> y(n+1,n+1,d);
  for(uint i=0; i!=n+1; ++i) { y[i].gradient(i)=1; }
  std::cerr << "y=" << y << std::endl;
  SparseDifferentialVector<X> yp(n+1,n+1,d);
  for(uint j=0; j<d; ++j) {
    yp=compose(tf,y);
    std::cerr << "yp=" << yp << std::endl;
    for(uint i=0; i!=n+1; ++i) {  
      y[i]=antiderivative(yp[i],n);
      y[i].value()=0;
      y[i].gradient(i)=1;
    }
    std::cerr << "y=" << y << std::endl << std::endl;
  } 
  for(uint i=0; i!=n; ++i) { y[i].value()=x[i]; }
  return y;
}



template<class X> 
SparseDifferentialVector<X>
flow(const SparseDifferentialVector<X>& f, const Vector<X>& x, uint to, uint so)
{
  return flow1(f,x,to,so);
}



//! Compute the flow map to the crossing set g under the vector field \a vf
template<class X> 
SparseDifferentialVector<X>  
hitting(const SparseDifferentialVector<X>& vf, const SparseDifferential<X>& g)
{
}




} //namespace Ariadne

#endif /* ARIADNE_SPARSE_DIFFERENTIAL_H */
