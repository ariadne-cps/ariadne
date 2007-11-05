/***************************************************************************
 *            multivariable_derivative.template.h
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
 
#include "taylor_derivative.h"

namespace Ariadne {


template<class X> 
Function::MultivariableDerivative<X> 
Function::reduce(const MultivariableDerivative<X>& x, const uint& d)
{
  assert(x.degree()>=d);
  Function::MultivariableDerivative<X> r(x.argument_size(),d);
  for(uint i=0; i!=x.result_size(); ++i) {
    for(MultiIndex j(x.argument_size()); j.degree() <= x.degree(); ++j) {
      r.set(i,j,x.get(i,j));
    }
  }
}







template<class X0, class X1, class X2> 
void 
Function::compute_composition(MultivariableDerivative<X0>& z, const MultivariableDerivative<X1>& y, const MultivariableDerivative<X2>& x)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  std::cerr << "y=" << y << std::endl;
  std::cerr << "z=" << z << std::endl;
  assert(z.degree()==x.degree());
  assert(z.degree()==y.degree());
  uint d=z.degree();
  MultivariableDerivative<X2> w=x;
  MultiIndex j(x.argument_size());
  for(uint i=0; i!=x.result_size(); ++i) {
    w.set(i,j,0);
  }
  std::cerr << "w=" << w << std::endl;
  MultivariableDerivative<X0> t(y.result_size(),x.argument_size(),0);
  std::cerr << "t[0]=" << t << std::endl;
  for(uint n=1; n<=d; ++n) {
    MultivariableDerivative<X0> u(y.result_size(),x.argument_size(),n);
    t=u;
    std::cerr << "t[" << n << "]=" << t << std::endl;
  }
  z=t;
}

template<class X> inline
std::ostream& 
Function::operator<<(std::ostream& os, const MultivariableDerivative<X>& x) {
  //  return os << "MultivariableDerivative( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
  uint degree=0;
  for(uint i=0; i!=x.result_size(); ++i) {
    os << (i==0 ? '[' : ';');
    for(MultiIndex j(x.argument_size()); j.degree()<=x.degree(); ++j) {
      if(j.degree()==0) {
        os << '[';
      } else if(j.degree()==degree) {
        os << ',';
      } else {
        degree=j.degree();
        os << ';';
      }
      os << x(i,j);
    }
    os << ']';
  }
  os << ']';
  return os;

//  return os << "MultivariableDerivative( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
}

} //namespace Ariadne
