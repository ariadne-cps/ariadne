/***************************************************************************
 *            build_vector_field.h
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
 
/*! \file build_vector_field.h
 *  \brief Macro to build a vector_field from a raw function template.
 */

#ifndef ARIADNE_BUILD_VECTOR_FIELD_H
#define ARIADNE_BUILD_VECTOR_FIELD_H

#include "numeric/traits.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/covector.h"
#include "linear_algebra/matrix.h"

#include "geometry/point.h"

#include "function/affine_derivative.h"
#include "function/taylor_variable.h"
#include "function/taylor_derivative.h"

#include "system/vector_field_interface.h"

#define ARIADNE_BUILD_VECTOR_FIELD(Nm,f,d,np,sm)      \
  template<class R> \
  class Nm \
    : public System::VectorFieldInterface<R> \
  { \
   private: \
    typedef typename Numeric::traits<R>::arithmetic_type X; \
   public: \
    template<class P> explicit Nm(const P& p) : _p(p) { } \
    Nm<R>* clone() const { return new Nm<R>(*this); }                     \
    const Geometry::Point<R>& parameters() const { return this->_p; } \
    const R& parameter(dimension_type j) const { return this->_p[j]; } \
    virtual LinearAlgebra::Vector<X> image(const Geometry::Point<X>& x) const { \
      LinearAlgebra::Vector<X> res(d); \
      f(res,x,this->_p); \
      return res; \
    } \
    virtual LinearAlgebra::Matrix<X> jacobian(const Geometry::Point<X>& x) const { \
      LinearAlgebra::Matrix<X> r(d,d); \
      const Geometry::Point<R>& p=this->_p; \
      Function::AffineDerivative<X> dr(d,d);  \
      Function::AffineDerivative<X> dx(d,d);                              \
      for(uint i=0; i!=d; ++i) { dx[i]=Function::AffineVariable<X>::variable(d,i,x[i]); } \
      f(dr,dx,p); \
      for(uint i=0; i!=d; ++i) { \
        for(uint j=0; j!=d; ++j) { \
          r(i,j)=dr[i].derivative(j); \
        } \
      }  \
      return r; \
    } \
    virtual Function::TaylorDerivative<X> derivative(const Geometry::Point<X>& x, const smoothness_type s) const { \
      const Geometry::Point<R>& p=this->_p; \
      Function::TaylorDerivative<X> dr(d,d,s);   \
      Function::TaylorDerivative<X> dx(d,d,s);   \
      for(uint i=0; i!=d; ++i) { dx[i]=Function::TaylorVariable<X>::variable(d,s,i,x[i]); } \
      f(dr,dx,p); \
      return dr; \
    } \
    dimension_type dimension() const { return d; } \
    size_type number_of_parameters() const { return np; } \
    smoothness_type smoothness() const { return sm; } \
    virtual std::string name() const { return "Nm"; } \
   private: \
    Geometry::Point<R> _p; \
  }; \


#endif /* ARIADNE_BUILD_VECTOR_FIELD_H */
