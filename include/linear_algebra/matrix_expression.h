/***************************************************************************
 *            matrix_expression.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file matrix_expressions.h
 *  \brief Matrix expression templates.
 */

#ifndef ARIADNE_MATRIX_EXPRESSION_H
#define ARIADNE_MATRIX_EXPRESSION_H

#include <iosfwd>

#include "../base/types.h"
#include "../base/array.h"
#include "../numeric/traits.h"
#include "../numeric/integer.h"
#include "../numeric/interval.h"

#include "../linear_algebra/exceptions.h"
#include "../linear_algebra/vector_expression.h"

namespace Ariadne {
  namespace LinearAlgebra {

    /*!\brief %Base class for all matrix expressions. */
    template<class E>
    class MatrixExpression 
    {
     public:
      /*!\brief Convert \a *this to a reference to E. */
      E& operator() () { return static_cast<E&>(*this); }
      /*!\brief Convert \a *this to a constant reference to E. */
      const E& operator() () const { return static_cast<const E&>(*this); }
    };
    


    /* Proxy for a matrix row in operator[] */
    template<class Mx>
    class MatrixRow
    {
     public:
      MatrixRow(Mx& A, const size_type& i) : _mx(A), _i(i) { }
      typename Mx::value_type& operator[](const size_type& j) { return _mx(_i,j); }
     private:
      Mx& _mx; const size_type _i;
    };
    

  }
}


#endif /* ARIADNE_MATRIX_EXPRESSION_H */
