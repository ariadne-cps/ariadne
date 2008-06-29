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

#include "base/types.h"
#include "base/array.h"
#include "numeric/traits.h"
#include "linear_algebra/declarations.h"
#include "output/logging.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/covector.h"
#include "linear_algebra/matrix.h"

#include "differentiation/taylor_variable.h"

namespace Ariadne {

      


    /*!\ingroup Differentiation
     * \brief Concrete class for functions.
     *
     * \deprecated This class has unknown bugs, and is deprecated in favour of TaylorVariable.
     */
    template<class X>
    class AffineVariable
      : public TaylorVariable<X>
    {
     public:
      /*! \brief Default constructor constructs an affine variable with value zero and no arguments. */
      explicit AffineVariable() : TaylorVariable<X>(0u,1u) { }
     
      /*! \brief Construct an affine variable with \a d arguments. */
      explicit AffineVariable(const size_type& as) : TaylorVariable<X>(as,1u) { } ;
     
      /*! \brief Construct an affine variable with \a d arguments based on the array \a ary. */
      template<class XX> explicit AffineVariable(const size_type& as, const XX* ary) : TaylorVariable<X>(as,1u,ary) { }
           
      /*! \brief Upcast from a Taylor variable. */
      AffineVariable(const TaylorVariable<X>& tv)
        : TaylorVariable<X>(tv) { ARIADNE_ASSERT(tv.degree()==1u); }
    
      // Dispatch assignment operator
      template<class XX> AffineVariable<X> operator=(const XX& x) { 
        this->TaylorVariable<X>::operator=(x); return *this; }

           
      /*! \brief Construct a constant variable with respect to \a as variables and value \a c. */
      static AffineVariable<X> constant(uint as, const X& c) { 
        AffineVariable<X> r(as); r.value()=c; return r; }
      /*! \brief Construct the variable of degree \a d at value \a x with respect to the \a i<sup>th</sup> variable of \a as. */
      static AffineVariable<X> variable(uint as, const X& x, uint i) { 
        AffineVariable<X> r(as); r.data()[i+1u]=x; return r; }

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
     private:
      static void instantiate();
    };
   
    template<class X>
    AffineVariable<X>& operator+=(AffineVariable<X>& r, const AffineVariable<X>& x) {
      static_cast<TaylorVariable<X>&>(r)+=static_cast<const TaylorVariable<X>&>(x); return r; }

    template<class X>
    AffineVariable<X>& operator-=(AffineVariable<X>& r, const AffineVariable<X>& x) {
      static_cast<TaylorVariable<X>&>(r)-=static_cast<const TaylorVariable<X>&>(x); return r; }


    template<class X>
    AffineVariable<X> operator+(const AffineVariable<X>& x1, const AffineVariable<X>& x2) {
      return static_cast<const TaylorVariable<X>&>(x1)+static_cast<const TaylorVariable<X>&>(x2); }

    template<class X>
    AffineVariable<X> operator-(const AffineVariable<X>& x1, const AffineVariable<X>& x2) {
      return static_cast<const TaylorVariable<X>&>(x1)-static_cast<const TaylorVariable<X>&>(x2); }

    template<class X>
    AffineVariable<X> operator*(const AffineVariable<X>& x1, const AffineVariable<X>& x2) {
      return static_cast<const TaylorVariable<X>&>(x1)*static_cast<const TaylorVariable<X>&>(x2); }

    template<class X>
    AffineVariable<X> operator/(const AffineVariable<X>& x1, const AffineVariable<X>& x2) {
      return static_cast<const TaylorVariable<X>&>(x1)/static_cast<const TaylorVariable<X>&>(x2); }

  
    // FIXME: We need to provide explicit operations here to prevent interception by 
    // templated Taylor variable operators.
    template<class X>
    AffineVariable<X> operator+(const TaylorVariable<X>& x1, const AffineVariable<X>& x2) {
      return x1+static_cast<const TaylorVariable<X>&>(x2); }

    template<class X>
    AffineVariable<X> operator-(const TaylorVariable<X>& x1, const AffineVariable<X>& x2) {
      return x1-static_cast<const TaylorVariable<X>&>(x2); }

    template<class X>
    AffineVariable<X> operator*(const TaylorVariable<X>& x1, const AffineVariable<X>& x2) {
      return x1*static_cast<const TaylorVariable<X>&>(x2); }



  
} // namespace Ariadne

#include "affine_variable.inline.h"

#endif /* ARIADNE_AFFINE_VARIABLE_H */
