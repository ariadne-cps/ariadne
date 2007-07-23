/***************************************************************************
 *            tensor_expression.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
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
 
/*! \file tensor_expression.h
 *  \brief Expression templates for tensor classes.
 */

#ifndef ARIADNE_TENSOR_EXPRESSION_H
#define ARIADNE_TENSOR_EXPRESSION_H

#include <iosfwd>

#include "../base/array.h"
#include "../numeric/numerical_traits.h"
#include "../linear_algebra/exceptions.h"
#include "../linear_algebra/multi_index.h"


namespace Ariadne {
  namespace LinearAlgebra {
    
    /*! \ingroup LinearAlgebra
     *  \brief %Base class for tensor expressions. 
     */
    template<class E> 
    class TensorExpression {
     public:
      typedef IndexArray index_array_type;
      const E& operator() () const { return static_cast<const E&>(*this); }
      E& operator() () { return static_cast<E&>(*this); }
    };
    

    /*! \ingroup LinearAlgebra
     *  \brief %Base class for symmetric tensor expressions. 
     */
    template<class E> 
    class SymmetricTensorExpression 
      : public TensorExpression<E>
    {
      typedef MultiIndex multi_index_type;
      const E& operator() () const { return static_cast<const E&>(*this); }
      E& operator() () { return static_cast<E&>(*this); }
    };
    

    /*! \ingroup LinearAlgebra
     *  \brief %Base class for derivative tensor expressions. 
     */
    template<class E> 
    class DerivativeTensorExpression 
      : public TensorExpression<E>
    {
      typedef MultiIndex multi_index_type;
      const E& operator() () const { return static_cast<const E&>(*this); }
      E& operator() () { return static_cast<E&>(*this); }
    };
    
  
  }
}

#endif /* ARIADNE_TENSOR_EXPRESSION_H */
