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

namespace Ariadne {
  namespace Geometry {
    template<class X> class Point;
    template<class X> class Rectangle;
  }

  namespace Function {
      
    template<class X> class AffineVariable;
    template<class X> void neg(AffineVariable<X>&, const AffineVariable<X>&);
    template<class X> void rec(AffineVariable<X>&, const AffineVariable<X>&);
    template<class X> void add(AffineVariable<X>&, const AffineVariable<X>&, const AffineVariable<X>&);
    template<class X> void sub(AffineVariable<X>&, const AffineVariable<X>&, const AffineVariable<X>&);
    template<class X> void mul(AffineVariable<X>&, const AffineVariable<X>&, const AffineVariable<X>&);
    template<class X> void div(AffineVariable<X>&, const AffineVariable<X>&, const AffineVariable<X>&);
    template<class X> void compose(AffineVariable<X>&, const AffineVariable<X>&, const AffineVariable<X>&);
    template<class X> void reduce(AffineVariable<X>&, const AffineVariable<X>&, smoothness_type s);
    template<class X> std::ostream& operator<<(std::ostream&, const AffineVariable<X>&);



    /*!\ingroup FunctionVariable
     * \brief Concrete class for functions.
     */
    template<class X>
    class AffineVariable
    {
     public:
      /*! The type used to represent numbers. */
      typedef X value_type;

      /*! \brief Destructor. */
      ~AffineVariable();
     
      /*! \brief Default constructor constructs an affine variable with value zero and no arguments. */
      AffineVariable();
     
      /*! \brief Construct an affine variable with \a d arguments. */
      AffineVariable(const dimension_type& d);
     
      /*! \brief Construct an affine variable with \a d arguments based on the array \a ary. */
      template<class XX> AffineVariable(const dimension_type& d, const XX* ary);
      template<class XX> AffineVariable(const uint& d, const XX* ary);
     
      /*! \brief Construct an affine variable with \a d arguments based on the array \a a1. */
      template<class XT, class CVXT> AffineVariable(const dimension_type& d, const XT& x, const CVXT* dxp);
     
      /*! \brief Construct an affine variable with value \a x and derivatives \a dx. */
      AffineVariable(const X& x, const LinearAlgebra::Covector<X>& dx);
     
      /*! \brief Copy constructor. */
      AffineVariable(const AffineVariable<X>& av);
     
      /*! \brief Assignment operator. */
      AffineVariable<X>& operator=(const AffineVariable<X>& av);

      /*! \brief Assign a constant. */
      template<class XX> AffineVariable<X>& operator=(const XX& c);

      /*! \brief Construct a constant variable with respect to \a as variables and value \a c. */
      static AffineVariable<X> constant(uint as, const X& c); 
      /*! \brief Construct the variable of degree \a d at value \a value with respect to the \a i<sup>th</sup> variable of \a as. */
      static AffineVariable<X> variable(uint as, const X& value, uint i);

      //@{
      //! \name Data access
      /*! \brief The number of variables of the argument. */
      size_type argument_size() const; 
      /*! \brief The degree of the model (number of derivatives computed). Equal to one. */
      smoothness_type degree() const;
      /*! \brief The value of the variable. */
      const X& value() const { return this->_x; }
      /*! \brief The differential of the variable with respect to the indepenent variables. */
      const LinearAlgebra::Covector<X>& derivative() const { return this->_dx; }
      /*! \brief The differential of the variable with respect to the \a j th indepenent variable. (Deprecated)  */
      const X& derivative(size_type j) const { return this->_dx[j]; }
      //@}

      //@{
      //! \name Modifying operations.
      /*! \brief Resize to a variable in \a ad independent quantities. */
      void resize(const dimension_type& ad);

      /*! \brief Set an element in degree 0 to x. */
      template<class XT> void set(const XT& x);

      /*! \brief Set an element in degree 1 to x. */
      template<class XT> void set(dimension_type j, const XT& x);
      //@}
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;

      //@{
      //! \name Inplace arithmetic operations
      AffineVariable<X>& operator+=(const AffineVariable<X>&);
      //! \name Inplace arithmetic operations
      AffineVariable<X>& operator*=(const X&);

 #ifdef DOXYGEN
      //@{
      //! \name Comparison operators.
      /*! \brief Equality operator. */
      friend template<class X11, class X22> tribool operator==(const AffineVariable<X11>& x1, const AffineVariable<X22>& x2); 
      /*! \brief Inequality operator. */
      friend template<class X11, class X22> tribool operator!=(const AffineVariable<X11>& x1, const AffineVariable<X22>& x2); 
      //@}
     
      //@{
      //! \name Arithmetic operations
      
      /*! \brief The differential with the minimal value from \a x1 or \a x2. */
      friend AffineVariable<X> min(const AffineVariable<X>& x1, const AffineVariable<X>& x2);
      /*! \brief The differential with the maximal value from \a x1 or \a x2. */
      friend AffineVariable<X> max(const AffineVariable<X>& x1, const AffineVariable<X>& x2);
      /*! \brief Absolute value function. Returns \f$(-x,-\dot{x})\f$ if \f$x<0\f$ and \f$(x,\dot{x})\f$ if \f$x>0\f$. 
       *  If interval arithmetic is used, returns \f$(0,[-1,1]\cdot\dot{x})\f$ at zero. */
      friend AffineVariable<X> abs(const AffineVariable<X>& x);
      
      /*! \brief The positive value of \a x. Returns a copy \f$(x,\; \dot{x})\f$. */
      friend AffineVariable<X> operator+(const AffineVariable<X>& x);
      /*! \brief Negation of a differential. Returns \f$(-x,\; -\dot{x})\f$. */
      friend AffineVariable<X> operator-(const AffineVariable<X>& x);
      /*! \brief  Addition of differentials. Returns \f$(x_1+x_2,\; \dot{x}_1+\dot{x_2})\f$. Addition of a constant is also available. */
      friend AffineVariable<X> operator+(const AffineVariable<X>& x1, const AffineVariable<X>& x2);
      /*! \brief  Subtraction of differentials. Returns \f$(x_1-x_2,\; \dot{x}_1-\dot{x_2})\f$. Subtraction of/from a constant is also available. */
      friend AffineVariable<X> operator-(const AffineVariable<X>& x1, const AffineVariable<X>& x2);
      /*! \brief Multiplication of differentials. Returns \f$(x_1x_2,\; x_2\dot{x}_1+x_1\dot{x_2})\f$. Multiplication by a constant is also available. */
      friend AffineVariable<X> operator*(const AffineVariable<X>& x1, const AffineVariable<X>& x2);
      /*! \brief Division of differentials. Returns \f$(x_1/x_2,\; \dot{x}_1/x2-x_1\dot{x_2}/x_2^2)\f$. Division of/by a constant is also available. */
      friend AffineVariable<X> operator/(const AffineVariable<X>& x1, const AffineVariable<X>& x2);
      /*! \brief Integer power of a differential. Returns \f$(x^n,\; nx^{n-1}\dot{x})\f$. */
      friend template<class N> AffineVariable<X> pow(const AffineVariable<X>& x, const N& n);
      //@}
      
      //@{
      //! \name Algebraic and transcendental operations
      /*! \brief Square root. Returns \f$(\sqrt{x},\; \dot{x}/2\sqrt{x})\f$. */
      friend AffineVariable<X> sqrt(const AffineVariable<X>& x); 
      /*! \brief Exponential. Returns \f$(e^x,\; \dot{x}e^x)\f$. */
      friend AffineVariable<X> exp(const AffineVariable<X>& x); 
      /*! \brief Natural logarithm. Returns \f$(\log x,\; \dot{x}/x)\f$. */
      friend AffineVariable<X> log(const AffineVariable<X>& x); 
      /*! \brief Sine function. Returns \f$(\sin x,\; \dot{x}\cos x)\f$. */
      friend AffineVariable<X> sin(const AffineVariable<X>& x); 
      /*! \brief Cosine function. Returns \f$(\cos x,\; -\dot{x}\sin x)\f$. */
      friend AffineVariable<X> cos(const AffineVariable<X>& x); 
      /*! \brief Tangent function. Returns \f$(\tan x,\; \dot{x}\sec^2 x)\f$. */
      friend AffineVariable<X> tan(const AffineVariable<X>& x); 
      /*! \brief Inverse sine function. Returns \f$(\sin^{-1} x,\; \dot{x}/\sqrt{1-x^2})\f$. */
      friend AffineVariable<X> asin(const AffineVariable<X>& x); 
      /*! \brief Inverse cosine function. Returns \f$(\cos^{-1} x,\; \dot{x}/\sqrt{1-x^2})\f$. */
      friend AffineVariable<X> acos(const AffineVariable<X>& x); 
      /*! \brief Inverse tangent function. Returns \f$(\tan^{-1} x,\; \dot{x}/(1+x^2)\f$. */
      friend AffineVariable<X> atan(const AffineVariable<X>& x); 
      //@}
      
      //@{
      //! \name Input/output operators.
      /*! \brief Stream insertion operator. */
      friend std::ostream& operator<<(std::ostream& os, const AffineVariable<X>& x); 
     //@}
#endif

     private:
      friend void neg<>(AffineVariable<X>&, const AffineVariable<X>&);
      friend void rec<>(AffineVariable<X>&, const AffineVariable<X>&);
      friend void add<>(AffineVariable<X>&, const AffineVariable<X>&, const AffineVariable<X>&);
      friend void sub<>(AffineVariable<X>&, const AffineVariable<X>&, const AffineVariable<X>&);
      friend void mul<>(AffineVariable<X>&, const AffineVariable<X>&, const AffineVariable<X>&);
      friend void div<>(AffineVariable<X>&, const AffineVariable<X>&, const AffineVariable<X>&);
      friend void compose<>(AffineVariable<X>&, const AffineVariable<X>&, const AffineVariable<X>&);
      friend std::ostream& operator<< <>(std::ostream&, const AffineVariable<X>&);
     private:
      static void instantiate();
     private:
      X _x;
      LinearAlgebra::Covector<X> _dx;
    };
   
    template<class X11, class X22> bool operator==(const AffineVariable<X11>& x1, const AffineVariable<X22>& x2); 
    template<class X11, class X22> bool operator!=(const AffineVariable<X11>& x1, const AffineVariable<X22>& x2); 

    template<class X> std::ostream& operator<<(std::ostream& os, const AffineVariable<X>& x);

    template<class X> AffineVariable<X> operator+(const AffineVariable<X>& x);
    template<class X> AffineVariable<X> operator-(const AffineVariable<X>& x);
    template<class X> AffineVariable<X> operator+(const AffineVariable<X>& x1, const AffineVariable<X>& x2); 
    template<class X> AffineVariable<X> operator-(const AffineVariable<X>& x1, const AffineVariable<X>& x2); 
    template<class X> AffineVariable<X> operator*(const AffineVariable<X>& x1, const AffineVariable<X>& x2); 
    template<class X> AffineVariable<X> operator/(const AffineVariable<X>& x1, const AffineVariable<X>& x2);

    template<class C, class X> AffineVariable<X> operator+(const C& c, const AffineVariable<X>& x); 
    template<class C, class X> AffineVariable<X> operator+(const AffineVariable<X>& x, const C& c); 
    template<class C, class X> AffineVariable<X> operator-(const C& c, const AffineVariable<X>& x); 
    template<class C, class X> AffineVariable<X> operator-(const AffineVariable<X>& x, const C& c); 
    template<class C, class X> AffineVariable<X> operator*(const C& c, const AffineVariable<X>& x); 
    template<class C, class X> AffineVariable<X> operator*(const AffineVariable<X>& x, const C& c); 
    template<class C, class X> AffineVariable<X> operator/(const C& c, const AffineVariable<X>& x); 
    template<class C, class X> AffineVariable<X> operator/(const AffineVariable<X>& x, const C& c); 

    template<class X> AffineVariable<X> operator/(const AffineVariable<X>& x, const X& c); 

    template<class X> AffineVariable<X> min(const AffineVariable<X>& x1,const AffineVariable<X>& x2); 
    template<class X> AffineVariable<X> max(const AffineVariable<X>& x1,const AffineVariable<X>& x2); 
    template<class X> AffineVariable<X> abs(const AffineVariable<X>& x); 

    template<class X> AffineVariable<X> pos(const AffineVariable<X>& x); 
    template<class X> AffineVariable<X> neg(const AffineVariable<X>& x); 
    template<class X> AffineVariable<X> add(const AffineVariable<X>& x1, const AffineVariable<X>& x2); 
    template<class X> AffineVariable<X> sub(const AffineVariable<X>& x1, const AffineVariable<X>& x2); 
    template<class X> AffineVariable<X> mul(const AffineVariable<X>& x1, const AffineVariable<X>& x2); 
    template<class X> AffineVariable<X> div(const AffineVariable<X>& x1, const AffineVariable<X>& x2); 
    template<class X> AffineVariable<X> pow(const AffineVariable<X>& x, int n); 
    template<class X> AffineVariable<X> sqrt(const AffineVariable<X>& x);
    template<class X> AffineVariable<X> exp(const AffineVariable<X>& x); 
    template<class X> AffineVariable<X> log(const AffineVariable<X>& x); 
    template<class X> AffineVariable<X> sin(const AffineVariable<X>& x); 
    template<class X> AffineVariable<X> cos(const AffineVariable<X>& x); 
    template<class X> AffineVariable<X> tan(const AffineVariable<X>& x); 
    template<class X> AffineVariable<X> asin(const AffineVariable<X>& x); 
    template<class X> AffineVariable<X> acos(const AffineVariable<X>& x); 
    template<class X> AffineVariable<X> atan(const AffineVariable<X>& x); 

  }
}

#include "affine_variable.inline.h"

#endif /* ARIADNE_AFFINE_VARIABLE_H */
