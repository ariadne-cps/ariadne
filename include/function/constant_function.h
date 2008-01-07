/***************************************************************************
 *            function/constant_function.h
 *
 *  Copyright  2007  Pieter Collins
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
 
/*! \file function/constant_function.h
 *  \brief Functions of constant form \f$x\rightarrow c\f$.
 */

#ifndef ARIADNE_CONSTANT_FUNCTION_H
#define ARIADNE_CONSTANT_FUNCTION_H

#include "base/types.h"

#include "numeric/declarations.h"
#include "numeric/traits.h"

#include "linear_algebra/declarations.h"
#include "linear_algebra/vector.h"

#include "geometry/declarations.h"

#include "function/function_interface.h"


namespace Ariadne {
  namespace Function {

    /*! \brief An constant function \f$f(x)=c\f$ on Euclidean space. 
     *  \ingroup FunctionTypes
     */
    template<class R>
    class ConstantFunction
      : public FunctionInterface<R> 
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef typename Numeric::traits<R>::interval_type I;
     private:
      size_type _as;
      LinearAlgebra::Vector<F> _c;
     public:
      /*! \brief The type of denotable real number used to describe the system. */
      typedef R real_type;
      /*! \brief The type of denotable state the system acts on. */
      typedef Geometry::Point<R> state_type;
      
      /*! \brief Default constructor constructs a function on a zero-dimensional space. */
      explicit ConstantFunction() {}
      /*! \brief Construct from the result size \a rs, argument size \a as, and a C array \a ptr. */
      template<class RR> explicit ConstantFunction(size_type rs, size_type as, const RR* ptr) 
        : _as(as), _c(rs,ptr) { }
      
      /*! \brief Copy constructor. */
      ConstantFunction(const ConstantFunction<R>& f)
        : _as(f._as), _c(f._c) { }
      /*! \brief Assignment operator. */
      ConstantFunction<R>& operator=(const ConstantFunction<R>& f) {
        this->_as=f._as; this->_c=f._c; return *this; }
      /*! \brief Returns a pointer to a dynamically-allocated copy of the function. */
      ConstantFunction<R>* clone() const { 
        return new ConstantFunction<R>(*this); }

      
      /*! \brief  An approximation to the image of an approximate point. */
      LinearAlgebra::Vector<F> evaluate(const LinearAlgebra::Vector<F>& x) const;
      /*! \brief The Jacobian derivative matrix at a point. */
      LinearAlgebra::Matrix<F> jacobian(const LinearAlgebra::Vector<F>& x) const;
      /*! \brief All the derivative values up to degree \a s. */
      TaylorDerivative<F> derivative(const LinearAlgebra::Vector<F>& x, const smoothness_type& s) const;

           
      /*! \brief  The linear transformation of the function. */
      const LinearAlgebra::Vector<F>& c() const { return _c; }
      
      /*! \brief The size of the result. */
      virtual size_type result_size() const {
        return _c.size(); }

      /*! \brief  The size of the argument. */
      virtual size_type argument_size() const {
        return _as; }
      
      /*! \brief The smoothness of the function. */
      virtual smoothness_type smoothness() const { 
        return std::numeric_limits<smoothness_type>::max(); }
      

      /*! \brief  The name of the system. */
      std::string name() const { return "ConstantFunction"; }
      
      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;
    };


  }
}


#endif /* ARIADNE_CONSTANT_FUNCTION_H */
