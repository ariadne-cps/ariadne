/***************************************************************************
 *            vector_field.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file vector_field.h
 *  \brief Vector_type field interface.
 */
 
#ifndef ARIADNE_VECTOR_FIELD_H
#define ARIADNE_VECTOR_FIELD_H

#include <limits>

#include "base/types.h"
#include "numeric/declarations.h"
#include "numeric/traits.h"
#include "linear_algebra/declarations.h"
#include "function/declarations.h"
#include "geometry/declarations.h"

#include <boost/shared_ptr.hpp>
#include "function/function_interface.h"
#include "geometry/euclidean_space.h"


namespace Ariadne {

  namespace Function { 
    template<class X> class TaylorSeriesAffineVariable; 
    template<class X> class TaylorSeriesTaylorVariable; 
  }

  namespace System {

    /*!\ingroup System
     * \ingroup ContinuousTime
     * \brief Abstract base class for (differentiable) vector fields.
     * 
     * The system is specified by the method operator()(const Geometry::Point<F>& pt) const,
     * This method should compute an interval vector \f$v=\overline{f}(A)\f$ with the
     * following properties:
     *   -# \f$f(p)\subset\overline{f}(A)\f$,
     *   -# If \f$A_1\subset A_0\f$, then \f$\overline{f}(A_1)\subset 
     *       \overline{f}(A_0)\f$, and
     *   -# If \f$\bigcap_{n\in\mathbb{N}}A_n=\{x\}\f$, then 
     *       \f$\bigcap_{n\in\mathbb{N}}\overline{f}(A_n)=\{f(x)\}\f$.
     *
     * More succinctly, we say that \f$\overline{f}(A_n)\f$ converges monotonically 
     * as \f$A_n\f$ tends to a point.
     *
     * Additional accuracy can be obtained be using derivatives.
     * The method jacobian(const Geometry::Point<F>& pt) const computes the derivative matrix at/over the point \a pt.
     */
    template<class R>
    class VectorField {
     protected:
      typedef typename Numeric::traits<R>::arithmetic_type F; 
      typedef typename Numeric::traits<R>::interval_type I; 
     public:
      /*! \brief The type used to represent time. */
      typedef Numeric::Rational time_type;
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type used to describe the state space. */
      typedef Geometry::EuclideanSpace state_space_type;
      /*! \brief The type of denotable state the system acts on. */
      typedef Geometry::Point<R> state_type;
      
      /*! \brief Destructor. */
      ~VectorField();

      /*! \brief Construct from a function interface and parameters. */
      VectorField(const Function::FunctionInterface<R>& f);
      /*! \brief Make a copy (clone) of the vector field. */
      VectorField<R>* clone() const { return new VectorField<R>(*this); }
     
      /*! \brief The function defining the vector field. */
      Function::FunctionInterface<R>& function() const { return *this->_function_ptr; }

      /*! \brief An approximation to the vector field at a point. */
      LinearAlgebra::Vector<F> operator() (const Geometry::Point<F>& x) const { return this->evaluate(x); }

      /*! \brief An approximation to the vector field at a point. */
      LinearAlgebra::Vector<F> evaluate(const Geometry::Point<F>& x) const;

      /*! \brief An approximation to the Jacobian derivative at a point. */
      LinearAlgebra::Matrix<F> jacobian(const Geometry::Point<F>& x) const;
    
       /*! \brief An approximation to the vector field at a point. */
      Function::TaylorDerivative<F> derivative(const Geometry::Point<F>& x, const smoothness_type& s) const;

      // Used in integration method
      void compute(Function::TaylorSeriesAffineVariable<F>*, const Function::TaylorSeriesAffineVariable<F>* x) const { };
      // Used in integration method
      void compute(Function::TaylorSeriesTaylorVariable<F>*, const Function::TaylorSeriesTaylorVariable<F>* x) const { };

      /*! \brief The state space of the vector field. */
      Geometry::EuclideanSpace state_space() const { return Geometry::EuclideanSpace(this->_function_ptr->result_size()); }
      /*! \brief The dimension of the space the vector field lives in. */
      dimension_type dimension() const { return this->_function_ptr->result_size(); }
      /*! \brief The degree of differentiability of the vector field. */
      smoothness_type smoothness() const { return this->_function_ptr->smoothness(); }

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
     private:
      boost::shared_ptr< Function::FunctionInterface<R> > _function_ptr;
    };
   
    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const VectorField<R>& vf) {
      return vf.write(os);
    }
    
  }
}

#endif /* ARIADNE_VECTOR_FIELD_H */
