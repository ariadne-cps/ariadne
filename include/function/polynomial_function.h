/***************************************************************************
 *            function/polynomial_function.h
 *
 *  Copyright  2007 Pieter Collins
 *  pieter.collins@cwi.nl
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
 
/*! \file function/polynomial_function.h
 *  \brief Polynomial functions.
 */

#ifndef ARIADNE_POLYNOMIAL_FUNCTION_H
#define ARIADNE_POLYNOMIAL_FUNCTION_H

#include <iosfwd>
#include <string>
#include <sstream>

#include "base/types.h"
#include "base/array.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "differentiation/sorted_index.h"
#include "differentiation/multi_index.h"


namespace Ariadne {
  
   class latexstream;


  
    template<class R> class PolynomialFunction;

    template<class R0, class R1, class R2> void add(PolynomialFunction<R0>&, const PolynomialFunction<R1>&, const PolynomialFunction<R2>&);
    template<class R0, class R1, class R2> void sub(PolynomialFunction<R0>&, const PolynomialFunction<R1>&, const PolynomialFunction<R2>&);
    template<class R0, class R1, class R2> void mul(PolynomialFunction<R0>&, const PolynomialFunction<R1>&, const PolynomialFunction<R2>&);
    template<class R0, class R1> void pow(PolynomialFunction<R0>&, const PolynomialFunction<R1>&, const unsigned int&);

    template<class R0, class R1> void scale(PolynomialFunction<R0>&, const R1&);

    template<class R0, class R1> void derivative(PolynomialFunction<R0>&, const PolynomialFunction<R1>&, const size_type&);
    template<class R0, class R1, class R2> void compose(PolynomialFunction<R0>&, const PolynomialFunction<R1>&, const PolynomialFunction<R2>&);
    template<class R0, class R1, class R2> void combine(PolynomialFunction<R0>&, const PolynomialFunction<R1>&, const PolynomialFunction<R2>&);
    template<class R0, class R1, class R2> void join(PolynomialFunction<R0>&, const PolynomialFunction<R1>&, const PolynomialFunction<R2>&);


    template<class R1, class R2> PolynomialFunction<typename traits<R1,R2>::arithmetic_type> operator+(const PolynomialFunction<R1>&, const PolynomialFunction<R2>&);
    template<class R1, class R2> PolynomialFunction<typename traits<R1,R2>::arithmetic_type> operator-(const PolynomialFunction<R1>&, const PolynomialFunction<R2>&);
    template<class R1, class R2> PolynomialFunction<typename traits<R1,R2>::arithmetic_type> operator*(const PolynomialFunction<R1>&, const PolynomialFunction<R2>&);

    template<class R1, class R2> PolynomialFunction<typename traits<R1,R2>::arithmetic_type> operator*(const R1&, const PolynomialFunction<R2>&);
    template<class R1, class R2> PolynomialFunction<typename traits<R1,R2>::arithmetic_type> operator*(const PolynomialFunction<R1>&, const R2&);
    template<class R1, class R2> PolynomialFunction<typename traits<R1,R2>::arithmetic_type> operator/(const PolynomialFunction<R1>&, const R2&);

    template<class R1, class R2> PolynomialFunction<typename traits<R1,R2>::arithmetic_type> compose(const PolynomialFunction<R1>&, const PolynomialFunction<R2>&);
    template<class R> PolynomialFunction<typename traits<R>::arithmetic_type> pow(const PolynomialFunction<R>& p, const unsigned int& n);
    template<class R> PolynomialFunction<typename traits<R>::arithmetic_type> derivative(const PolynomialFunction<R>&, const size_type& k);

    template<class R> std::ostream& operator<<(std::ostream&, const PolynomialFunction<R>&);
    template<class R> std::istream& operator>>(std::istream&, PolynomialFunction<R>&);
  

    /*! \brief A polynomial with multivalued output, using a den.
     *  \ingroup FunctionTypes
     */
    template<class R>
    class PolynomialFunction {
      typedef typename traits<R>::arithmetic_type F;
      typedef typename traits<R>::interval_type I;
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
     public:
      /*! \brief Default constructor constructs a polynomial of degree zero with no arguments and no result variables. */
      PolynomialFunction();
      /*! \brief Construct from a string literal. */
      PolynomialFunction(const std::string& s);
      /*! \brief The zero polynomial in \a n variables with size \a m image and degree d. */
      PolynomialFunction(const size_type& rs, const size_type& as, const size_type& d);
      /*! \brief Construct a polynomial in \a n variables with size \a m image and degree \a d from the data array \a a. */
      template<class RR> PolynomialFunction(const size_type& rs, const size_type& as, const size_type& d, const RR* a);
      /*! \brief Copy constructor. */
      PolynomialFunction(const PolynomialFunction<R>& p);
      /*! \brief Copy assignment. */
      PolynomialFunction<R>& operator=(const PolynomialFunction<R>& p);
      /*! \brief Conversion constructor. */
      template<class RR> explicit PolynomialFunction(const PolynomialFunction<RR>& p);
      /*! \brief Conversion assignment. */
      template<class RR> PolynomialFunction<R>& operator=(const PolynomialFunction<RR>& p);
        
      /*! \brief Equality operator. */
      bool operator==(const PolynomialFunction<R>& p) const;
      /*! \brief Inequality operator. */
      bool operator!=(const PolynomialFunction<R>& p) const;

      /*! \brief The size of the argument. */
      size_type argument_size() const;
      /*! \brief The size of the result. */
      size_type result_size() const;
      /*! \brief The degree of the polynomial. */
      size_type degree() const;
      /*! \brief The smoothness of the function. */
      size_type smoothness() const;
      /*! \brief The data used to define the polynomial. */
      const array<R>& data() const;
      
      /*! \brief Set the \a j th value of the \a i th component to \a x. */
      void set(const size_type& i, const MultiIndex& j, const R& x); 
      /*! \brief A reference to \a j th value of the i th component. */
      R& at(const size_type& i, const MultiIndex& j);
      /*! \brief The \a j th value of the i th component. */
      const R& get(const size_type& i, const MultiIndex& j) const;
      /*! \brief Resize to a polynomial in \a n variables with size \a m image and degree d. */
      void resize(const size_type& rs, const size_type& as, const size_type& d);

      /*! \brief The \a i th component polynomial. */
      PolynomialFunction<R> component(const size_type& i) const;

      /*! \brief Evaluate the polynomial at the point \a x. */
      Vector<F> evaluate(const Vector<F>& x) const;
      
      /*! \brief Compute the derivate of the map at a point. */
      Matrix<F> jacobian(const Vector<F>& s) const;

      /*! \brief The zero polynomial with result size \a rs and argument size \a as. */
      static PolynomialFunction<R> zero(const size_type& rs, const size_type& as);
      /*! \brief The unit polynomial with result size 1 and argument size \a as. */
      static PolynomialFunction<R> one(const size_type& as);
      /*! \brief The constant polynomial with result size 1 and argument size \a as. */
      static PolynomialFunction<R> constant(const size_type& as, const R& c);

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);

#ifdef DOXYGEN
      /*! \brief Addition. */
      friend template<class R1, class R2> PolynomialFunction<typename traits<R1,R2>::arithmetic_type> operator+(const PolynomialFunction<R1>&, const PolynomialFunction<R2>&);
      /*! \brief Subtraction. */
      friend template<class R1, class R2> PolynomialFunction<typename traits<R1,R2>::arithmetic_type> operator-(const PolynomialFunction<R1>&, const PolynomialFunction<R2>&);
      /*! \brief Multiplication. At least one argument must be scalar-valued. */
      friend template<class R1, class R2> PolynomialFunction<typename traits<R1,R2>::arithmetic_type> operator*(const PolynomialFunction<R1>&, const PolynomialFunction<R2>&);

      /*! \brief Multiplication by a scalar. */
      friend template<class R1, class R2> PolynomialFunction<typename traits<R1,R2>::arithmetic_type> operator*(const R1&, const PolynomialFunction<R2>&);
      /*! \brief Multiplication by a scalar. */
      friend template<class R1, class R2> PolynomialFunction<typename traits<R1,R2>::arithmetic_type> operator*(const PolynomialFunction<R1>&, const R2&);
      /*! \brief Division by a scalar. */
      friend template<class R1, class R2> PolynomialFunction<typename traits<R1,R2>::arithmetic_type> operator/(const PolynomialFunction<R1>&, const R2&);

      /*! \brief Composition \f$p\circ q(x)=p(q(x))\f$. */
      friend template<class R1, class R2> PolynomialFunction<typename traits<R1,R2>::arithmetic_type> compose(const PolynomialFunction<R1>&, const PolynomialFunction<R2>&);
      /*! \brief Power of a scalar polynomial. */
      friend template<class R> PolynomialFunction<typename traits<R>::arithmetic_type> pow(const PolynomialFunction<R>& p, const unsigned int& n);
      /*! \brief Derivative with respect to variable \a k. */
      friend template<class R> PolynomialFunction<typename traits<R>::arithmetic_type> derivative(const PolynomialFunction<R>&, const size_type& k);
#endif

     private:
      static void instantiate();
      void _compute_jacobian() const;
      void _set_argument_size(const size_type& n);
      size_type _compute_maximum_component_size() const;
     private:
      template<class R0,class R1,class R2> friend void add(PolynomialFunction<R0>&,const PolynomialFunction<R1>&,const PolynomialFunction<R2>&);
      template<class R0,class R1,class R2> friend void sub(PolynomialFunction<R0>&,const PolynomialFunction<R1>&,const PolynomialFunction<R2>&);
      template<class R0,class R1,class R2> friend void mul(PolynomialFunction<R0>&,const PolynomialFunction<R1>&,const PolynomialFunction<R2>&);
      template<class R0,class R1,class R2> friend void div(PolynomialFunction<R0>&,const PolynomialFunction<R1>&,const PolynomialFunction<R2>&);
      template<class R0,class R1,class R2> friend void compose(PolynomialFunction<R0>&,const PolynomialFunction<R1>&,const PolynomialFunction<R2>&);
      template<class R0,class R1> friend void scale(PolynomialFunction<R0>&,const R1&);
     private:
      /* Components of the map. */
      size_type _result_size; 
      size_type _argument_size;
      size_type _degree; 
      array<R> _data;
    };
    
  
} // namespace Ariadne

#include "polynomial_function.inline.h"

#endif /* ARIADNE_POLYNOMIAL_FUNCTION_H */
