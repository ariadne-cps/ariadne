/***************************************************************************
 *            sparse_differential_vector.h
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

#ifndef ARIADNE_SPARSE_DIFFERENTIAL_VECTOR_H
#define ARIADNE_SPARSE_DIFFERENTIAL_VECTOR_H

#include <map>

#include "macros/throw.h"
#include "base/array.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "differentiation/multi_index.h"
#include "differentiation/power_series.h"
#include "differentiation/function_series.h"
#include "differentiation/sparse_differential.h"

#include "linear_algebra/matrix.code.h"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class PowerSeries;

template<class X> struct AffineTransformation;
template<class X> class SparseDifferential;

template<class X> class SparseDifferentialVector;

template<class X, class Y> Vector<Y> evaluate(const SparseDifferentialVector<X>& x, const Vector<Y>& y);
template<class X> SparseDifferentialVector<X> embed(const SparseDifferentialVector<X>& x, uint size, uint start);
template<class X> SparseDifferentialVector<X> translate(const SparseDifferentialVector<X>& x, const Vector<X>& c);
template<class X> SparseDifferentialVector<X> compose(const SparseDifferentialVector<X>& x, const SparseDifferentialVector<X>& y);
template<class X> SparseDifferentialVector<X> inverse(const SparseDifferentialVector<X>& x);
template<class X> SparseDifferentialVector<X> implicit(const SparseDifferentialVector<X>& x);
template<class X> SparseDifferentialVector<X> derivative(const SparseDifferentialVector<X>& x, uint j);
template<class X> SparseDifferentialVector<X> antiderivative(const SparseDifferentialVector<X>& x, uint j);
template<class X> SparseDifferentialVector<X> flow(const SparseDifferentialVector<X>& vf);
template<class X> SparseDifferentialVector<X> hitting(const SparseDifferentialVector<X>& vf, const SparseDifferentialVector<X>& g);



/*! \brief A class representing the derivatives of a vector quantity depending on multiple arguments. */
template<class X>
class SparseDifferentialVector
  : public Vector< SparseDifferential<X> >
{
  //BOOST_CONCEPT_ASSERT((DifferentialVectorConcept<SparseDifferentialVector<X> >));
 public:
  SparseDifferentialVector() 
    : Vector< SparseDifferential<X> >(0,SparseDifferential<X>()) { }
  SparseDifferentialVector(uint rs) 
    : Vector< SparseDifferential<X> >(rs,SparseDifferential<X>()) { }
  SparseDifferentialVector(uint rs, uint as, ushort d) 
    : Vector< SparseDifferential<X> >(rs,SparseDifferential<X>(as,d)) { }
  SparseDifferentialVector(const Vector< SparseDifferential<X> >& vsd) 
    : Vector< SparseDifferential<X> >(vsd) { }
  template<class XX> SparseDifferentialVector(uint rs, uint as, ushort d, const XX* ptr) 
    : Vector< SparseDifferential<X> >(rs,SparseDifferential<X>(as,d)) 
  { 
    for(uint i=0; i!=rs; ++i) { for(MultiIndex j(as); j.degree()<=d; ++j) {
        if(*ptr!=0) { (*this)[i][j]=*ptr; } ++ptr; } } 
  }
  SparseDifferentialVector(uint rs, uint as, ushort d, 
                           const Vector<X>& v, const Matrix<X>& A)
    :  Vector< SparseDifferential<X> >(rs,SparseDifferential<X>(as,d)) {
    ARIADNE_ASSERT(rs==v.size());
    ARIADNE_ASSERT(rs==A.row_size());
    ARIADNE_ASSERT(as==A.column_size());
    for(uint i=0; i!=this->result_size(); ++i) { 
      (*this)[i]=v[i]; 
      for(uint j=0; j!=this->argument_size(); ++j) { 
        const X& x=A[i][j];
        if(x!=0) { (*this)[i][j]=x; } } } 
  }

  uint result_size() const { return this->Vector< SparseDifferential<X> >::size(); }
  uint argument_size() const { return (*this)[0].argument_size(); }
  ushort degree() const { return (*this)[0].degree(); }

  Vector<X> value() const { 
    Vector<X> r(this->result_size()); for(uint i=0; i!=r.size(); ++i) { r[i]=(*this)[i].value(); } return r; }
  Matrix<X> jacobian() const { Matrix<X> r(this->result_size(),this->argument_size()); 
    for(uint i=0; i!=r.row_size(); ++i) { for(uint j=0; j!=r.column_size(); ++j) { r[i][j]=(*this)[i].gradient(j); } } return r; }

  void set_value(const Vector<X>& c) {
    ARIADNE_ASSERT(this->result_size()==c.size());
    for(uint i=0; i!=c.size(); ++i) { (*this)[i].set_value(c[i]); } }

  static SparseDifferentialVector<X> zero(uint rs, uint as, ushort d) {
    return SparseDifferentialVector<X>(rs,as,d);
  }

  static SparseDifferentialVector<X> constant(uint rs, uint as, ushort d, const Vector<X>& c) {
    ARIADNE_ASSERT(c.size()==rs);
    SparseDifferentialVector<X> result(rs,as,d);
    for(uint i=0; i!=rs; ++i) { result[i]=c[i]; }
    return result;
  }

  static SparseDifferentialVector<X> variable(uint rs, uint as, ushort d, const Vector<X>& x) {
    ARIADNE_ASSERT(x.size()==rs);
    SparseDifferentialVector<X> result(rs,as,d);
    for(uint i=0; i!=rs; ++i) { result[i]=x[i]; result[i][i]=X(1.0); }
    return result;
  }

  static SparseDifferentialVector<X> affine(uint rs, uint as, ushort d, const Vector<X>& b, const Matrix<X>& A) {
    ARIADNE_ASSERT(b.size()==rs);
    ARIADNE_ASSERT(A.row_size()==rs);
    ARIADNE_ASSERT(A.column_size()==as);
    ARIADNE_ASSERT(d>=1);
    SparseDifferentialVector<X> result(rs,as,d);
    for(uint i=0; i!=rs; ++i) { 
      result[i]=b[i]; 
      for(uint j=0; j!=as; ++j) { 
        result[i][j]=A[i][j]; 
      } 
    }
    return result;
  }

};


template<class X>
SparseDifferentialVector<X>&
operator+=(SparseDifferentialVector<X>& x, const Vector<X>& c)
{  
  assert(x.result_size()==c.size());
  for(uint i=0; i!=c.size();++i) {
    x[i]+=c[i];
  }
  return x;
}

  
template<class X>
SparseDifferentialVector<X>&
operator-=(SparseDifferentialVector<X>& x, const Vector<X>& c)
{  
  assert(x.result_size()==c.size());
  for(uint i=0; i!=c.size();++i) {
    x[i]-=c[i];
  }
  return x;
}

  
template<class X>
SparseDifferentialVector<X> 
operator+(const SparseDifferentialVector<X>& x, const Vector<X>& c)
{  
  SparseDifferentialVector<X> r(x);
  return r+=c;
}


template<class X>
SparseDifferentialVector<X> 
operator-(const SparseDifferentialVector<X>& x, const Vector<X>& c)
{  
  SparseDifferentialVector<X> r(x);
  return r-=c;
}


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
operator+(const SparseDifferentialVector<X>& x, const SparseDifferentialVector<X>& y)
{  
  SparseDifferentialVector<X> r(x);
  return r+=y;
}


template<class X>
SparseDifferentialVector<X> 
operator-(const SparseDifferentialVector<X>& x, const SparseDifferentialVector<X>& y)
{  
  SparseDifferentialVector<X> r(x);
  return r-=y;
}


template<class X>
SparseDifferentialVector<X>
join(const SparseDifferentialVector<X>& f, const SparseDifferentialVector<X>& g)
{
  ARIADNE_ASSERT(f.argument_size()==g.argument_size());
  SparseDifferentialVector<X> h(f.result_size()+g.result_size(),f.argument_size(),std::max(f.degree(),g.degree()));
  for(uint i=0; i!=f.result_size(); ++i) {
    h[i]=f[i];
  }
  for(uint i=0; i!=g.result_size(); ++i) {
    h[i+f.result_size()]=g[i];
  }
  return h;
}


template<class X>
SparseDifferentialVector<X>
join(const SparseDifferentialVector<X>& f, const SparseDifferential<X>& g)
{
  ARIADNE_ASSERT(f.argument_size()==g.argument_size());
  SparseDifferentialVector<X> h(f.result_size()+1u,f.argument_size(),std::max(f.degree(),g.degree()));
  for(uint i=0; i!=f.result_size(); ++i) {
    h[i]=f[i];
  }
  h[f.result_size()]=g;
  return h;
}


template<class X, class Y>
Y
evaluate(const SparseDifferential<X>& x, 
         const Vector<Y>& y)
{  
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  using namespace std;
  assert(x.argument_size()==y.size());
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;

  ushort d=x.degree();
  uint s=y.size();

  Y zero=y[0]; zero*=0;
  Y one=zero; one+=1;

  Y r=zero;

  //std::cerr << "zero="<<zero<<std::endl;
  //std::cerr << "one="<<one<<std::endl;
  // Use inefficient brute-force approach with lots of storage...
  array< array< Y > > val(s, array< Y >(d+1,zero));
  for(uint j=0; j!=s; ++j) {
    val[j][0]=one;
    for(ushort k=1; k<=d; ++k) {
      val[j][k]=val[j][k-1]*y[j];
    }
  }
  //std::cerr << "val="<<val<<std::endl;
  for(MultiIndex j(s); j.degree()<=d; ++j) {
    Y t=one;
    for(ushort k=0; k!=s; ++k) {
      t=t*val[k][j[k]];
    }
    r+=x[j]*t;
  }
  
  return r;
}


template<class X, class Y>
Vector<Y>
evaluate(const SparseDifferentialVector<X>& x, 
         const Vector<Y>& y)
{  
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  using namespace std;
  assert(x.argument_size()==y.size());
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;

  ushort d=x.degree();
  uint rs=x.result_size();
  uint as=x.argument_size();

  Y zero=y[0]; zero*=0;
  Y one=zero; one+=1;

  Vector<Y> r(rs,zero);

  // Use inefficient brute-force approach with lots of storage...
  array< array< Y > > val(as, array< Y >(d+1,zero));
  for(uint j=0; j!=as; ++j) {
    val[j][0]=one;
    for(ushort k=1; k<=d; ++k) {
      val[j][k]=val[j][k-1]*y[j];
    }
  }
  for(MultiIndex j(as); j.degree()<=d; ++j) {
    Y t=one;
    for(uint k=0; k!=as; ++k) {
      t=t*val[k][j[k]];
    }
    for(uint i=0; i!=rs; ++i) {
      //const X& xij=x[i][j];
      Y& ri=r[i];
      X xij=x[i][j];
      Y txij=xij*t;
      
      ri+=txij;
    }
  }
  
  return r;
}



template<class X>
SparseDifferentialVector<X> 
evaluate(const SparseDifferentialVector<X>& x, 
         const SparseDifferentialVector<X>& y)
{  
  assert(x.argument_size()==y.result_size());
  SparseDifferentialVector<X> r;
  
  static_cast<Vector< SparseDifferential<X> >&>(r) = 
    evaluate(x,static_cast<const Vector< SparseDifferential<X> >&>(y));
  for(uint i=0; i!=r.result_size(); ++i) { r[i].cleanup(); }
  return r;
}


template<class X>
SparseDifferentialVector<X> 
embed(const SparseDifferentialVector<X>& x, 
      uint size, uint start)
{  
  assert(start+x.argument_size()<=size);
  SparseDifferentialVector<X> r(x.result_size(),size,x.degree());
  for(uint i=0; i!=x.result_size(); ++i) { r[i]=embed(x[i],size,start); }
  return r;
}


/*! \brief Compose the series of \a x and the series of \a y, assuming that \a x is centred at the value of \a y. The value of \a y is therefore unused by this method. */
template<class X>
SparseDifferential<X> 
compose(const SparseDifferential<X>& x, 
        const SparseDifferentialVector<X>& y)
{  
  Vector<X> yv=y.value();
  SparseDifferentialVector<X>& ync=const_cast<SparseDifferentialVector<X>&>(y); 
  for(uint i=0; i!=ync.result_size(); ++i) { ync[i].value()=0; }
  SparseDifferential<X> r=evaluate(x,ync);
  ync+=yv;
  return r;
}


/*! \brief Compose the series of \a x and the series of \a y, assuming that \a x is centred at the value of \a y. The value of \a y is therefore unused by this method. */
template<class X>
SparseDifferentialVector<X> 
compose(const SparseDifferentialVector<X>& x, 
        const SparseDifferentialVector<X>& y)
{  
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  Vector<X> yv=y.value();
  //std::cerr<<"yv="<<yv<<std::endl;
  SparseDifferentialVector<X>& ync=const_cast<SparseDifferentialVector<X>&>(y); 
  for(uint i=0; i!=ync.result_size(); ++i) { ync[i].value()=0; }
  //std::cerr<<"ync="<<ync<<std::endl;
  SparseDifferentialVector<X> r=evaluate(x,ync);
  //std::cerr<<"r="<<r<<std::endl;
  ync+=yv;
  //std::cerr<<"ync="<<ync<<std::endl;
  return r;
}

template<class X> 
SparseDifferentialVector<X> 
implicit(const SparseDifferentialVector<X>& x)
{
  assert(x.result_size()<=x.argument_size());
  //std::cerr << "x=" << x << std::endl;
  
  uint rs=x.result_size();
  uint xas=x.argument_size();
  uint zas=x.argument_size()-x.result_size();
  ushort d=x.degree();

  Matrix<X> A1(rs,zas);
  for(uint i=0; i!=rs; ++i) {
    for(uint j=0; j!=zas; ++j) {
      A1(i,j)=x[i].gradient(j);
    }
  }
  
  Matrix<X> A2(rs,rs);
  for(uint i=0; i!=rs; ++i) {
    for(uint j=0; j!=rs; ++j) {
      A2(i,j)=x[i].gradient(zas+j);
    }
  }

  //std::cerr << "inverse(A2)="<<inverse(A2)<<std::endl;

  Matrix<X> A2inv=inverse(A2);
  Matrix<X> J(xas,rs);
  project(J,range(zas,zas+rs),range(0,rs)) = A2inv;
  //std::cerr << "J="<<J<<std::endl;

  SparseDifferentialVector<X> y(xas,zas,d);
  for(uint i=0; i!=zas; ++i) {
    y[i]=SparseDifferential<X>::variable(zas,d,1.0,i);
  }
  for(uint i=0; i!=rs; ++i) {
    // y[as+i]=Differential<X>::constant(as,d,0.0);
  }
  //std::cerr << "y="<<y<<std::endl;

  for(ushort i=0; i!=d; ++i) {
    SparseDifferentialVector<X> z=compose(x,y);
    //std::cerr << "z="<<z<<std::endl;
    y-=J*z;
    //std::cerr << "y="<<y<<std::endl;
  }

  SparseDifferentialVector<X> r(rs,zas,d);
  for(uint i=0; i!=rs; ++i) {
    r[i]=y[zas+i];
  }
  return r;
}



template<class X> 
SparseDifferentialVector<X> 
inverse(const SparseDifferentialVector<X>& x)
{
  return inverse(x,Vector<X>(x.argument_size()));
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
  ushort d=x.degree();
  Vector<X> z(n,0);
  Matrix<X> J=inverse(x.jacobian());

  SparseDifferentialVector<X> y(n,n,d);
  SparseDifferentialVector<X> id(n,n,d);
  for(uint i=0; i!=n; ++i) { id[i][i]=1.0; }
  
  for(uint i=0; i!=n; ++i) { 
    y[i].value()=c[i]; 
    for(uint j=0; j!=n; ++j) { 
      y[i][j]=J[i][j];
    }
  }

  for(ushort i=2; i<=d; ++i) {
    SparseDifferentialVector<X> z=compose(x,y);
    z-=id;
    y-=J*z;
  }
  return y;
}


template<class X> 
SparseDifferentialVector<X> 
derivative(const SparseDifferentialVector<X>& x, uint j)
{
  ushort d=std::min(uint(x.degree()),1u)-1u;
  SparseDifferentialVector<X> r(x.result_size(),x.argument_size(), d);
  for(uint i=0; i!=x.result_size(); ++i) {
    r[i]=derivative(x[i],j);
  }
  return r;
}


template<class X> 
SparseDifferentialVector<X> 
antiderivative(const SparseDifferentialVector<X>& x, uint j)
{
  SparseDifferentialVector<X> r(x.result_size(),x.argument_size(), x.degree()+1u);
  for(uint i=0; i!=x.result_size(); ++i) {
    r[i]=antiderivative(x[i],j);
  }
  return r;
}


template<class X> 
SparseDifferentialVector<X>
flow1(const SparseDifferentialVector<X>& f, const Vector<X>& x, ushort to, ushort so)
{
  // f is an untimed vector field
  assert(f.result_size()==f.argument_size());
  assert(f.result_size()==x.size());
  uint n=x.size();
  ushort d=f.degree();
  Matrix<X> Q(n+1,n); for(uint i=0; i!=n; ++i) { Q[i][i]=1; }
  
  Vector<X> tx(n+1);
  for(uint i=0; i!=n; ++i) { tx[i]=x[i]; } 
  tx[n]=0.0;

  /*
    SparseDifferentialVector<X> tf(n,n+1,d);
  MultiIndex ta(n+1);
  for(uint i=0; i!=n; ++i) {
    for(typename SparseDifferential<X>::const_iterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
      const MultiIndex& a=iter->first;
      for(uint k=0; k!=n; ++k) { ta.set(k,a[k]); }
      tf[i][ta]=iter->second;
    }
  }
  */
  //std::cerr << "f=" << f << std::endl;
  //std::cerr << "tf=" << tf << std::endl << std::endl;

  SparseDifferentialVector<X> y(n,n+1,1u);
  for(uint i=0; i!=n; ++i) { y[i].gradient(i)=1; }
  //std::cerr << "y=" << y << std::endl;
  SparseDifferentialVector<X> yp(n,n+1,1u);
  for(ushort j=0; j<d; ++j) {
    yp=compose(f,y);
    //std::cerr << "yp=" << yp << std::endl;
    for(uint i=0; i!=n; ++i) {  
      y[i]=antiderivative(yp[i],n);
      y[i].value()=0;
      y[i].gradient(i)=1;
    }
    //std::cerr << "y=" << y << std::endl << std::endl;
  } 
  for(uint i=0; i!=n; ++i) { y[i].value()=x[i]; }
  return y;
}



template<class X> 
SparseDifferentialVector<X>
flow(const SparseDifferentialVector<X>& f, const Vector<X>& x, ushort to, ushort so)
{
  return flow1(f,x,to,so);
}



//! Compute the flow map to the crossing set g under the vector field \a vf
template<class X> 
SparseDifferentialVector<X>  
hitting(const SparseDifferentialVector<X>& vf, const SparseDifferential<X>& g)
{
}


//! Translate the polynomial given by \a x to one with centre \a v.
template<class X> 
SparseDifferentialVector<X>  
translate(const SparseDifferentialVector<X>& x, const Vector<X>& v)
{
  uint as=v.size();
  ushort d=x.degree();
  SparseDifferentialVector<X> t=SparseDifferentialVector<X>::variable(as,as,d,v);
  return evaluate(x,t);
}


//! Scale the polynomial given by \a x by the values in the array \a s. 
template<class X> 
SparseDifferentialVector<X>  
scale(const SparseDifferentialVector<X>& x, const Vector<X>& s)
{
  uint as=s.size();
  ushort d=x.degree();
  SparseDifferentialVector<X> t(as,as,d);
  for(uint i=0; i!=as; ++i) { t[i][i]=s[i]; } 
  return evaluate(x,t);
}


template<class X>
std::ostream& operator<<(std::ostream& os, const SparseDifferentialVector<X>& x)
{
  const Vector< SparseDifferential<X> >& xv=x;
  return os << xv;
}


} //namespace Ariadne

#endif /* ARIADNE_SPARSE_DIFFERENTIAL_VECTOR_H */
