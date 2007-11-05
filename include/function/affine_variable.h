/***************************************************************************
 *            affine_variable.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
 
/*! \file system/affine_variable.h
 *  \brief An affine approximation to a function.
 */
 
#ifndef ARIADNE_AFFINE_VARIABLE_H
#define ARIADNE_AFFINE_VARIABLE_H

#include <vector>
#include <string>
#include <exception>
#include <stdexcept>

#include "../base/types.h"
#include "../base/array.h"
#include "../numeric/numerical_traits.h"
#include "../linear_algebra/declarations.h"
#include "../output/logging.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/covector.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {
  namespace Geometry {
    template<class X> class Point;
    template<class X> class Rectangle;
  }

  namespace Function {
      
    template<class R0, class R1> class AffineVariable;
    template<class R0, class R1> void neg(AffineVariable<R0,R1>&, const AffineVariable<R0,R1>&);
    template<class R0, class R1> void rec(AffineVariable<R0,R1>&, const AffineVariable<R0,R1>&);
    template<class R0, class R1> void add(AffineVariable<R0,R1>&, const AffineVariable<R0,R1>&, const AffineVariable<R0,R1>&);
    template<class R0, class R1> void sub(AffineVariable<R0,R1>&, const AffineVariable<R0,R1>&, const AffineVariable<R0,R1>&);
    template<class R0, class R1> void mul(AffineVariable<R0,R1>&, const AffineVariable<R0,R1>&, const AffineVariable<R0,R1>&);
    template<class R0, class R1> void div(AffineVariable<R0,R1>&, const AffineVariable<R0,R1>&, const AffineVariable<R0,R1>&);
    template<class R0, class R1> void compose(AffineVariable<R0,R1>&, const AffineVariable<R0,R1>&, const AffineVariable<R0,R1>&);
    template<class R0, class R1> void reduce(AffineVariable<R0,R1>&, const AffineVariable<R0,R1>&, smoothness_type s);
    template<class R0, class R1> std::ostream& operator<<(std::ostream&, const AffineVariable<R0,R1>&);



    /*!\ingroup FunctionVariable
     * \brief Concrete class for functions.
     */
    template<class X, class CV>
    class AffineVariable
    {
     public:
      /*! \brief Destructor. */
      ~AffineVariable();
     
      /*! \brief Default constructor constructs an affine variable with value zero and no arguments. */
      AffineVariable();
     
      /*! \brief Construct an affine variable with \a d arguments. */
      AffineVariable(dimension_type d);
     
      /*! \brief Construct an affine variable with \a ad arguments based on the array \a a1. */
      template<class XT, class CVXT> AffineVariable(const dimension_type& d, const XT& x, const CVXT* dxp);
     
      /*! \brief Construct an affine variable with value \a x and derivatives \a dx. */
      template<class XT, class CVT> AffineVariable(const XT& x, const CVT& dx);
     
      /*! \brief Construct an affine variable with \a ad arguments based on the array \a a1. */
      AffineVariable(const X& x, const CV& dx);
     
      /*! \brief Copy constructor. */
      AffineVariable(const AffineVariable<X,CV>& av);
     
      /*! \brief Assignement operator. */
      AffineVariable<X,CV>& operator=(const AffineVariable<X,CV>& av);

      /*! \brief Resize to a variable in \a ad independent quantities. */
      void resize(const dimension_type& ad);

      /*! \brief Set an element in degree 0 to x. */
      template<class XT> void set(const XT& x);

      /*! \brief Set an element in degree 1 to x. */
      template<class XT> void set(dimension_type j, const XT& x);

      /*! \brief The smoothness of the model. */
      smoothness_type smoothness() const;

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;

      template<class X0, class X1> friend void neg(AffineVariable<X0,X1>&, const AffineVariable<X0,X1>&);
      template<class X0, class X1> friend void rec(AffineVariable<X0,X1>&, const AffineVariable<X0,X1>&);
      template<class X0, class X1> friend void add(AffineVariable<X0,X1>&, const AffineVariable<X0,X1>&, const AffineVariable<X0,X1>&);
      template<class X0, class X1> friend void sub(AffineVariable<X0,X1>&, const AffineVariable<X0,X1>&, const AffineVariable<X0,X1>&);
      template<class X0, class X1> friend void mul(AffineVariable<X0,X1>&, const AffineVariable<X0,X1>&, const AffineVariable<X0,X1>&);
      template<class X0, class X1> friend void div(AffineVariable<X0,X1>&, const AffineVariable<X0,X1>&, const AffineVariable<X0,X1>&);
      template<class X0, class X1> friend void compose(AffineVariable<X0,X1>&, const AffineVariable<X0,X1>&, const AffineVariable<X0,X1>&);

     private:
      static void instantiate();

     private:
      X _x;
      CV _dx;
    };
   
  }
}

#include "affine_variable.inline.h"

#endif /* ARIADNE_AFFINE_VARIABLE_H */
