/***************************************************************************
 *            taylor_derivative.template.h
 *
 *  Copyright 2007  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

namespace Ariadne {

template<class X> template<class V> inline
Function::TaylorDerivative<X> 
Function::TaylorDerivative<X>::constant(size_type r, size_type a, smoothness_type d, const V& c) 
{
  TaylorDerivative<X> result(r,a,d);
  ARIADNE_ASSERT(c.size()==r);
  for(size_type i=0; i!=r; ++i) {
    result._variables[i].value()=c[i];
  }
  return result;
}

template<class X> template<class V> inline
Function::TaylorDerivative<X> 
Function::TaylorDerivative<X>::variable(size_type r, size_type a, smoothness_type d, const V& x) 
{
  ARIADNE_ASSERT(a==r);
  ARIADNE_ASSERT(x.size()==r);
  TaylorDerivative<X> result(r,a,d);
  //size_type inc=compute_polynomial_data_size(1u,a,d);
  for(size_type i=0; i!=r; ++i) {
    result._variables[i].value()=x[i];
    result._variables[i].data()[i+1u]=1;
  }
  return result;
}

template<class X> template<class V> inline
Function::TaylorDerivative<X> 
Function::TaylorDerivative<X>::variable(const V& x, smoothness_type d) 
{
  TaylorDerivative<X> result(x.size(),x.size(),d);
  //size_type inc=compute_polynomial_data_size(1u,a,d);
  for(size_type i=0; i!=x.size(); ++i) {
    result._variables[i].value()=x[i];
    result._variables[i].data()[i+1u]=1;
  }
  return result;
}



template<class P, class M> 
void
Function::evaluate_polynomial(M& r, const P& p, const M& x)
{
  ARIADNE_ASSERT(p.argument_size()==x.size());
  ARIADNE_ASSERT(p.result_size()==r.size());

  typedef typename M::value_type R;
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;
  size_type d=p.degree();
  size_type rs=r.size();
  size_type ms=x.size();
  assert(d>=1);

  R zero = x[0]; zero*=0;
  R one = zero; one+=1;

  // Use inefficient brute-force approach with lots of storage...

  // Set up an array of powers of elements
  array< array< R > > val(ms, array< R >(d+1));
  for(uint j=0; j!=ms; ++j) {
    val[j][0]=one;
    if(d>=1) {
      val[j][1]=x[j];
      if(d>=2) {
        val[j][1]=pow(x[j],2);
        for(uint k=3; k<=d; ++k) {
          val[j][k]=val[j][k-2]*val[j][2];
        }
      }
    }
  }

  // Set all elements of the result to zero
  for(uint i=0; i!=rs; ++i) { r[i]=zero; }

  // Add powers of the argument to the result
  for(MultiIndex j(ms); j.degree()<=d; ++j) {
    R t=one;
    for(uint k=0; k!=ms; ++k) {
      t=t*val[k][j[k]];
    }
    for(uint i=0; i!=rs; ++i) {
      R tf=t; tf*=x[i][j];
      r[i]+=tf;
    }
  }
  return r;
}


template<class X, class R> 
array<R>
Function::evaluate(const TaylorDerivative<X>& y, const array<R>& x)
{
  using namespace std;
  ARIADNE_ASSERT(y.argument_size()==x.size());
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;
  size_type d=y.degree();
  size_type rs=y.result_size();
  size_type ms=x.size();
  assert(d>=1);

  R zero = x[0]; zero*=0;
  R one = zero; one+=1;

  // Use inefficient brute-force approach with lots of storage...
  array< array< R > > val(ms, array< R >(d+1));
  for(uint j=0; j!=ms; ++j) {
    val[j][0]=one;
    val[j][1]=x[j];
    for(uint k=2; k<=d; ++k) {
      val[j][k]=val[j][k-1]*x[j];
    }
  }

  array<R> r(rs,zero);
  for(MultiIndex j(ms); j.degree()<=d; ++j) {
    R t=one;
    for(uint k=0; k!=ms; ++k) {
      t=t*val[k][j[k]];
    }
    for(uint i=0; i!=rs; ++i) {
      R tf=t; tf*=y[i][j];
      r[i]+=tf;
    }
  }
  return r;
}


} //namespace Ariadne
