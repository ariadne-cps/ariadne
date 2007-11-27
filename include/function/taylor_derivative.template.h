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

template<class X>
LinearAlgebra::Vector<X>
Function::TaylorDerivative<X>::value() const
{
  return LinearAlgebra::Vector<X>(this->result_size(),this->data().begin(),this->_increment());
}

template<class X>
LinearAlgebra::Matrix<X>
Function::TaylorDerivative<X>::jacobian() const
{
  return LinearAlgebra::Matrix<X>(this->result_size(),this->argument_size(),this->data().begin()+1u,this->_increment(),1u);
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
Function::compute_composition(TaylorVariable<X0>& z, const TaylorVariable<X1>& y, const TaylorDerivative<X2>& x)
{
  // FIXME: Rewrite this function
}

template<class X0, class X1, class X2> 
void 
Function::compute_composition(TaylorDerivative<X0>& z, const TaylorDerivative<X1>& y, const TaylorDerivative<X2>& x)
{
  // FIXME: Rewrite this function
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "z=" << z << std::endl;
  assert(z.degree()==x.degree());
  assert(z.degree()==y.degree());
  size_type d=z.degree();
  TaylorDerivative<X2> w=x;
  //w.value()=0;
  //std::cerr << "w=" << w << std::endl;
  TaylorDerivative<X0> t(y.result_size(),x.argument_size(),0);
  //std::cerr << "t[0]=" << t << std::endl;
  for(uint n=1; n<=d; ++n) {
    TaylorDerivative<X0> u(y.result_size(),x.argument_size(),n);
    //compute_product(u,t,w);
    //u.value()=y[d-n];
    t=u;
    //std::cerr << "t[" << n << "]=" << t << std::endl;
  }
  z=t;
}

template<class X> inline
std::ostream& 
Function::operator<<(std::ostream& os, const TaylorDerivative<X>& x) {
  //  return os << "TaylorDerivative( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
  os << '[';
  for(size_type i=0; i!=x.result_size(); ++i) {
    size_type degree=0;
    for(MultiIndex j(x.argument_size()); j.degree()<=x.degree(); ++j) {
      if(j.degree()==0) {
        os << '[';
      } else if(j.degree()==degree) {
        os << ',';
      } else {
        degree=j.degree();
        os << ';';
      }
      os << x.get(i,j);
    }
    os << ']';
  }
  os << ']';
  return os;

//  return os << "TaylorDerivative( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
}

} //namespace Ariadne
