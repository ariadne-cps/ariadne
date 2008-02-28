/***************************************************************************
 *            build_function.h
 *
 *  Copyright 2007  Pieter Collins
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
 
/*! \file build_function.h
 *  \brief Macro to build an %Ariadne function from a raw function template.
 */

#ifndef ARIADNE_BUILD_FUNCTION_H
#define ARIADNE_BUILD_FUNCTION_H

#include "linear_algebra/vector.h"
#include "linear_algebra/covector.h"
#include "linear_algebra/matrix.h"

#include "function/taylor_variable.h"
#include "function/taylor_derivative.h"

#define ARIADNE_BUILD_FUNCTION(Nm,f,rs,as,np,sm)   \
  template<class R> \
  class Nm \
    : public Function::FunctionInterface<R> \
  { \
   private: \
    typedef typename Numeric::traits<R>::arithmetic_type X; \
   public: \
    template<class P> explicit Nm(const P& p) : _p(p) { } \
    Nm<R>* clone() const { return new Nm<R>(*this); } \
    virtual LinearAlgebra::Vector<X> evaluate(const LinearAlgebra::Vector<X>& x) const { \
      LinearAlgebra::Vector<X> res(rs); \
      f(res,x,this->_p); \
      return res; \
    } \
    virtual LinearAlgebra::Matrix<X> jacobian(const LinearAlgebra::Vector<X>& x) const { \
      LinearAlgebra::Matrix<X> r(rs,as); \
      const LinearAlgebra::Vector<X>& p=this->_p; \
      Function::TaylorDerivative<X> dr(as,as,1u); \
      Function::TaylorDerivative<X> dx=Function::TaylorDerivative<X>::variable(x,1u); \
      f(dr,dx,p); \
      return dr.jacobian(); \
    } \
    virtual Function::TaylorDerivative<X> derivative(const LinearAlgebra::Vector<X>& x, const smoothness_type& s) const { \
      const LinearAlgebra::Vector<X>& p=this->_p; \
      Function::TaylorDerivative<X> dx(as,as,s); \
      Function::TaylorDerivative<X> dr(rs,as,s); \
      for(uint i=0; i!=as; ++i) { dx[i]=x[i]; } \
      f(dr,dx,p);      \
      return dr; \
    } \
    size_type result_size() const { return rs; } \
    size_type argument_size() const { return as; } \
    size_type number_of_parameters() const { return np; } \
    smoothness_type smoothness() const { return sm; } \
    virtual std::string name() const { return "Nm"; } \
    virtual std::ostream& write(std::ostream& os) const  { return os << "Nm"; } \
   private: \
    LinearAlgebra::Vector<X> _p; \
  }; \


#endif /* ARIADNE_BUILD_FUNCTION_H */
