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
 

namespace Ariadne {

template<class X0, class X1, class X2> inline
void
Function::compute_product(TaylorDerivative<X0>& x0, const TaylorDerivative<X1>& x1, const TaylorDerivative<X2>& x2)
{
  assert(x0.argument_size()==x1.argument_size());
  assert(x0.argument_size()==x2.argument_size());
  for(MultiIndex i1(x1.argument_size()); i1.degree() <= x1.degree(); ++i1) {
    for(MultiIndex i2(x2.argument_size()); i2.degree() <= std::min(x2.degree(),x0.degree()-i1.degree()); ++i2) {
      MultiIndex i0=i1+i2;
      //std::cout << "i0=" << i0 << ", i1=" << i1 << ", i2=" << i2 << std::endl;
      // FIXME: Use Integer
      //Numeric::Integer c=i0.factorial()/(i1.factorial()*i2.factorial());
      unsigned int c=choose(i0,i1);
      x0[i0]+=X0(c)*x1[i1]*x2[i2];
    }
  }
}

template<class X> 
Function::TaylorDerivative<X> 
Function::reduce(const TaylorDerivative<X>& x, const size_type& d)
{
  assert(x.degree()>=d);
  Function::TaylorDerivative<X> r(x.argument_size(),d);
  for(MultiIndex i(x.argument_size()); i.degree() <= x.degree(); ++i) {
    r[i]=x[i];
  }
}


template<class X> 
Function::TaylorDerivative<X> 
Function::derivative(TaylorDerivative<X>& x, const size_type& k)
{
  if(x.degree()==0) {
    return TaylorDerivative<X>(x.argument_size(),0);
  } 
  TaylorDerivative<X> r(x.argument_size(),x.degree()-1);
  MultiIndex e(x.argument_size());
  e.set(k,1);
  for(MultiIndex j(r.argument_size()); j.degree() <= r.degree(); ++j) {
    r[j]=x[j+e];
  }
}






template<class X0, class X1, class X2> 
void 
Function::compute_composition(TaylorDerivative<X0>& z, const ScalarDerivative<X1>& y, const TaylorDerivative<X2>& x)
{
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "z=" << z << std::endl;
  assert(z.degree()==x.degree());
  assert(z.degree()==y.degree());
  size_type d=z.degree();
  TaylorDerivative<X2> w=x;
  w.value()=0;
  //std::cerr << "w=" << w << std::endl;
  TaylorDerivative<X0> t(x.argument_size(),0,y[d]);
  //std::cerr << "t[0]=" << t << std::endl;
  for(uint n=1; n<=d; ++n) {
    TaylorDerivative<X0> u(x.argument_size(),n);
    compute_product(u,t,w);
    u.value()=y[d-n];
    t=u;
    //std::cerr << "t[" << n << "]=" << t << std::endl;
  }
  z=t;
}

template<class X> inline
std::ostream& 
Function::operator<<(std::ostream& os, const TaylorDerivative<X>& x) {
  //  return os << "TaylorDerivative( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
  size_type degree=0;
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
  return os;

//  return os << "TaylorDerivative( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
}

} //namespace Ariadne
