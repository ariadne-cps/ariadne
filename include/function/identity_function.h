/***************************************************************************
 *            identity_function.h
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
 
/*! \file function/identity_function.h
 *  \brief The identity function
 */

#ifndef ARIADNE_IDENTITY_FUNCTION_H
#define ARIADNE_IDENTITY_FUNCTION_H

#include "base/types.h"

#include "numeric/declarations.h"
#include "numeric/traits.h"

#include "linear_algebra/declarations.h"
#include "linear_algebra/vector.h"

#include "function/function_interface.h"


namespace Ariadne {

    /*! \brief The identity function on Euclidean space. 
     *  \ingroup FunctionTypes
     */
    template<class R>
    class IdentityFunction
      : public FunctionInterface<R> 
    {
      typedef typename traits<R>::arithmetic_type F;
      typedef typename traits<R>::interval_type I;
     public:
      /*! \brief The type of denotable real number used to describe the system. */
      typedef R real_type;
      
      /*! \brief Constructor constructs a function on a \a d-dimensional space. */
      explicit IdentityFunction(size_type d) : _dimension(d) { }
      
      /*! \brief Returns a pointer to a dynamically-allocated copy of the function. */
      IdentityFunction<R>* clone() const { return new IdentityFunction<R>(*this); }

      
      /*! \brief  An approximation to the image of an approximate point. */
      Vector<F> evaluate(const Vector<F>& x) const;
      /*! \brief The Jacobian derivative matrix at a point. */
      Matrix<F> jacobian(const Vector<F>& x) const;
      /*! \brief The Jacobian derivative matrix at a point. */
      TaylorDerivative<F> derivative(const Vector<F>& x, const smoothness_type& s) const;

           
      /*! \brief The size of the result. */
      virtual size_type result_size() const { return this->_dimension; }

      /*! \brief  The size of the argument. */
      virtual size_type argument_size() const { return this->_dimension; }
      
      /*! \brief The smoothness of the result. */
      virtual smoothness_type smoothness() const { return std::numeric_limits<smoothness_type>::max(); }
      

      /*! \brief  The name of the system. */
      std::string name() const { return "IdentityFunction"; }
      
      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;
     private:
      size_type _dimension;
    };


  
} // namespace Ariadne


#endif /* ARIADNE_IDENTITY_FUNCTION_H */
