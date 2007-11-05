/***************************************************************************
 *            polynomial.h
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
 
/*! \file polynomial.h
 *  \brief Polynomials.
 */

#ifndef ARIADNE_POLYNOMIAL_H
#define ARIADNE_POLYNOMIAL_H

#include <iosfwd>
#include <string>
#include <sstream>

#include "../base/types.h"
#include "../base/array.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "sorted_index.h"
#include "multi_index.h"


namespace Ariadne {
  
  namespace Geometry { template<class R> class Rectangle; }
  namespace Output { class latexstream; }


  namespace Function {
  
    template<class R> class Polynomial;

    template<class R0, class R1, class R2> void add(Polynomial<R0>&, const Polynomial<R1>&, const Polynomial<R2>&);
    template<class R0, class R1, class R2> void sub(Polynomial<R0>&, const Polynomial<R1>&, const Polynomial<R2>&);
    template<class R0, class R1, class R2> void mul(Polynomial<R0>&, const Polynomial<R1>&, const Polynomial<R2>&);
    template<class R0, class R1> void pow(Polynomial<R0>&, const Polynomial<R1>&, const unsigned int&);

    template<class R0, class R1> void scale(Polynomial<R0>&, const R1&);

    template<class R0, class R1> void derivative(Polynomial<R0>&, const Polynomial<R1>&, const size_type&);
    template<class R0, class R1, class R2> void compose(Polynomial<R0>&, const Polynomial<R1>&, const Polynomial<R2>&);
    template<class R0, class R1, class R2> void combine(Polynomial<R0>&, const Polynomial<R1>&, const Polynomial<R2>&);
    template<class R0, class R1, class R2> void join(Polynomial<R0>&, const Polynomial<R1>&, const Polynomial<R2>&);

    template<class R0, class R1, class R2> void truncate(Polynomial<R0>&, const Polynomial<R1>&, const Geometry::Rectangle<R1>&, const size_type&, const size_type&);


    template<class R1, class R2> Polynomial<typename Numeric::traits<R1,R2>::arithmetic_type> operator+(const Polynomial<R1>&, const Polynomial<R2>&);
    template<class R1, class R2> Polynomial<typename Numeric::traits<R1,R2>::arithmetic_type> operator-(const Polynomial<R1>&, const Polynomial<R2>&);
    template<class R1, class R2> Polynomial<typename Numeric::traits<R1,R2>::arithmetic_type> operator*(const Polynomial<R1>&, const Polynomial<R2>&);

    template<class R1, class R2> Polynomial<typename Numeric::traits<R1,R2>::arithmetic_type> operator*(const R1&, const Polynomial<R2>&);
    template<class R1, class R2> Polynomial<typename Numeric::traits<R1,R2>::arithmetic_type> operator*(const Polynomial<R1>&, const R2&);
    template<class R1, class R2> Polynomial<typename Numeric::traits<R1,R2>::arithmetic_type> operator/(const Polynomial<R1>&, const R2&);

    template<class R1, class R2> Polynomial<typename Numeric::traits<R1,R2>::arithmetic_type> compose(const Polynomial<R1>&, const Polynomial<R2>&);
    template<class R> Polynomial<typename Numeric::traits<R>::arithmetic_type> pow(const Polynomial<R>& p, const unsigned int& n);
    template<class R> Polynomial<typename Numeric::traits<R>::arithmetic_type> derivative(const Polynomial<R>&, const size_type& k);
    template<class R> Polynomial<typename Numeric::traits<R>::interval_type> truncate(const Polynomial<R>&, const Geometry::Rectangle<typename Numeric::traits<R>::number_type>&, const size_type&, const size_type&);

    template<class R> std::ostream& operator<<(std::ostream&, const Polynomial<R>&);
    template<class R> std::istream& operator>>(std::istream&, Polynomial<R>&);
  

    /*! \brief A polynomial with multivalued output, using a den.
     *  \ingroup FunctionTypes
     */
    template<class R>
    class Polynomial {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef typename Numeric::traits<R>::interval_type I;
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
     public:
      /*! \brief Default constructor constructs a polynomial of degree zero with no arguments and no result variables. */
      Polynomial();
      /*! \brief Construct from a string literal. */
      Polynomial(const std::string& s);
      /*! \brief The zero polynomial in \a n variables with size \a m image and degree d. */
      Polynomial(const size_type& rs, const size_type& as, const size_type& d);
      /*! \brief Construct a polynomial in \a n variables with size \a m image and degree \a d from the data array \a a. */
      template<class RR> Polynomial(const size_type& rs, const size_type& as, const size_type& d, const RR* a);
      /*! \brief Copy constructor. */
      Polynomial(const Polynomial<R>& p);
      /*! \brief Copy assignment. */
      Polynomial<R>& operator=(const Polynomial<R>& p);
      /*! \brief Conversion constructor. */
      template<class RR> explicit Polynomial(const Polynomial<RR>& p);
      /*! \brief Conversion assignment. */
      template<class RR> Polynomial<R>& operator=(const Polynomial<RR>& p);
        
      /*! \brief Equality operator. */
      bool operator==(const Polynomial<R>& p) const;
      /*! \brief Inequality operator. */
      bool operator!=(const Polynomial<R>& p) const;

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
      Polynomial<R> component(const size_type& i) const;

      /*! \brief Truncate the polynomial to a polynomial of degree \a d and smoothness \a s within the domain \a domain. */
      Polynomial<I> truncate(const size_type& degree, const size_type& smoothness, const Geometry::Rectangle<R>& domain) const;

      /*! \brief Evaluate the polynomial at the point \a x. */
      LinearAlgebra::Vector<F> evaluate(const LinearAlgebra::Vector<F>& x) const;
      
      /*! \brief Compute the derivate of the map at a point. */
      LinearAlgebra::Matrix<F> jacobian(const LinearAlgebra::Vector<F>& s) const;

      /*! \brief Compute the derivate of the function with respect to the \a k<sup>th</sup> variable. */
      Polynomial<F> derivative(const size_type& k) const;

      /*! \brief The zero polynomial with result size \a rs and argument size \a as. */
      static Polynomial<R> zero(const size_type& rs, const size_type& as);
      /*! \brief The unit polynomial with result size 1 and argument size \a as. */
      static Polynomial<R> one(const size_type& as);
      /*! \brief The constant polynomial with result size 1 and argument size \a as. */
      static Polynomial<R> constant(const size_type& as, const R& c);

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);

#ifdef DOXYGEN
      /*! \brief Addition. */
      friend template<class R1, class R2> Polynomial<typename Numeric::traits<R1,R2>::arithmetic_type> operator+(const Polynomial<R1>&, const Polynomial<R2>&);
      /*! \brief Subtraction. */
      friend template<class R1, class R2> Polynomial<typename Numeric::traits<R1,R2>::arithmetic_type> operator-(const Polynomial<R1>&, const Polynomial<R2>&);
      /*! \brief Multiplication. At least one argument must be scalar-valued. */
      friend template<class R1, class R2> Polynomial<typename Numeric::traits<R1,R2>::arithmetic_type> operator*(const Polynomial<R1>&, const Polynomial<R2>&);

      /*! \brief Multiplication by a scalar. */
      friend template<class R1, class R2> Polynomial<typename Numeric::traits<R1,R2>::arithmetic_type> operator*(const R1&, const Polynomial<R2>&);
      /*! \brief Multiplication by a scalar. */
      friend template<class R1, class R2> Polynomial<typename Numeric::traits<R1,R2>::arithmetic_type> operator*(const Polynomial<R1>&, const R2&);
      /*! \brief Division by a scalar. */
      friend template<class R1, class R2> Polynomial<typename Numeric::traits<R1,R2>::arithmetic_type> operator/(const Polynomial<R1>&, const R2&);

      /*! \brief Composition \f$p\circ q(x)=p(q(x))\f$. */
      friend template<class R1, class R2> Polynomial<typename Numeric::traits<R1,R2>::arithmetic_type> compose(const Polynomial<R1>&, const Polynomial<R2>&);
      /*! \brief Power of a scalar polynomial. */
      friend template<class R> Polynomial<typename Numeric::traits<R>::arithmetic_type> pow(const Polynomial<R>& p, const unsigned int& n);
      /*! \brief Derivative with respect to variable \a k. */
      friend template<class R> Polynomial<typename Numeric::traits<R>::arithmetic_type> derivative(const Polynomial<R>&, const size_type& k);
      /*! \brief Truncate within \a r to a polynomial of degree at most \a d, putting the error into terms of degree \a s. */
      friend template<class R> Polynomial<typename Numeric::traits<R>::arithmetic_type> truncate(const Polynomial<R>& p, const Rectangle<R>& bb, const size_type& d, const size_type& s);
#endif
     private:
      static void instantiate();
      void _compute_jacobian() const;
      void _set_argument_size(const size_type& n);
      size_type _compute_maximum_component_size() const;
     private:
      template<class R0,class R1,class R2> friend void add(Polynomial<R0>&,const Polynomial<R1>&,const Polynomial<R2>&);
      template<class R0,class R1,class R2> friend void sub(Polynomial<R0>&,const Polynomial<R1>&,const Polynomial<R2>&);
      template<class R0,class R1,class R2> friend void mul(Polynomial<R0>&,const Polynomial<R1>&,const Polynomial<R2>&);
      template<class R0,class R1,class R2> friend void div(Polynomial<R0>&,const Polynomial<R1>&,const Polynomial<R2>&);
      template<class R0,class R1,class R2> friend void compose(Polynomial<R0>&,const Polynomial<R1>&,const Polynomial<R2>&);
      template<class R0,class R1> friend void scale(Polynomial<R0>&,const R1&);
     private:
      /* Components of the map. */
      size_type _result_size; 
      size_type _argument_size;
      size_type _degree; 
      array<R> _data;
    };
    
  }
}

#include "polynomial.inline.h"

#endif /* ARIADNE_POLYNOMIAL_H */
