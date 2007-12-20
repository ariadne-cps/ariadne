/***************************************************************************
 *            build_map.h
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
 
/*! \file build_map.h
 *  \brief Macro to build a map from a raw function template.
 */

#ifndef ARIADNE_BUILD_MAP_H
#define ARIADNE_BUILD_MAP_H

#include "linear_algebra/vector.h"
#include "linear_algebra/covector.h"
#include "linear_algebra/matrix.h"

#include "geometry/point.h"

#include "function/affine_variable.h"
#include "function/affine_derivative.h"
#include "function/taylor_variable.h"
#include "function/taylor_derivative.h"

#define ARIADNE_BUILD_MAP(Nm,f,rd,ad,np,sm)   \
  template<class R> \
  class Nm \
    : public System::MapInterface<R> \
  { \
   private: \
    typedef typename Numeric::traits<R>::arithmetic_type X; \
   public: \
    template<class P> explicit Nm(const P& p) : _p(p) { } \
    Nm<R>* clone() const { return new Nm<R>(*this); } \
    const Geometry::Point<R>& parameters() const { return this->_p; } \
    const R& parameter(dimension_type i) const { return this->_p[i]; } \
    virtual Geometry::Point<X> image(const Geometry::Point<X>& x) const { \
      Geometry::Point<X> res(rd); \
      f(res,x,this->_p); \
      return res; \
    } \
    virtual LinearAlgebra::Matrix<X> jacobian(const Geometry::Point<X>& x) const { \
      LinearAlgebra::Matrix<X> r(rd,ad); \
      const Geometry::Point<R>& p=this->_p; \
      Function::AffineDerivative<X> dr(ad,ad); \
      Function::AffineDerivative<X> dx(rd,ad);                           \
      for(uint i=0; i!=ad; ++i) { dx[i]=Function::AffineVariable<X>::variable(ad,i,x[i]); } \
      f(dr,dx,p); \
      for(uint i=0; i!=rd; ++i) { \
        for(uint j=0; j!=ad; ++j) { \
          r(i,j)=dr[i].derivative(j); \
        } \
      }  \
      return r; \
    } \
    virtual Function::TaylorDerivative<X> derivative(const Geometry::Point<X>& x, const smoothness_type s) const { \
      const Geometry::Point<R>& p=this->_p; \
      Function::TaylorDerivative<X> dx(ad,ad,s); \
      Function::TaylorDerivative<X> dr(rd,ad,s); \
      for(uint i=0; i!=ad; ++i) { dx[i]=x[i]; } \
      f(dr,dx,p);      \
      return dr; \
    } \
    dimension_type result_dimension() const { return rd; } \
    dimension_type argument_dimension() const { return ad; } \
    size_type number_of_parameters() const { return np; } \
    smoothness_type smoothness() const { return sm; } \
    virtual std::string name() const { return "Nm"; } \
   private: \
    Geometry::Point<R> _p; \
  }; \


#endif /* ARIADNE_BUILD_MAP_H */
