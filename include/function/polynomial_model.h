/***************************************************************************
 *            polynomial_model.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file polynomial_model.h
 *  \brief PolynomialFunction models for maps.
 */

#ifndef ARIADNE_POLYNOMIAL_MODEL_H
#define ARIADNE_POLYNOMIAL_MODEL_H

#include <iosfwd>
#include <string>
#include <sstream>

#include "base/array.h"
#include "base/types.h"
#include "linear_algebra/declarations.h"

namespace Ariadne {
  namespace Function {
    template<class R> class Monomial;
    template<class R> class PolynomialFunction;
    template<class R> class PolynomialModel;

    /*! \brief Graded lexicographical ordering. */
    template<class R> bool operator<(const Monomial<R>&, const Monomial<R>&);
  }

  namespace Numeric {
    template<class R> struct traits< Function::PolynomialFunction<R> > {
      typedef Function::PolynomialFunction<typename traits<R>::arithmetic_type> arithmetic_type;
    };
  }
  
  namespace Function {
    
    /*! \brief A monomial in several variables. */
    template<class R>
    class Monomial {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef typename Numeric::traits<R>::arithmetic_type result_type;
     public:
      /*! \brief The type of denotable real number used for the coefficients. */
      typedef R real_type;
     public:
      /*! \brief Construct from a string literal. 
       *
       *  The literal must be of the form <tt>c x_0^a0 x_1^a1 ... x_n^an</tt>,
       *  where \c c is a denotable real number and \c a0,...\c an are positive integers.
       *  Terms may be omitted.
       */
      Monomial(const std::string& s);
      /*! \brief Default constructor. */
      Monomial() : _coefficient(0), _multi_index() { }
      /*! \brief The zero monomial in \a n variables. */
      Monomial(const size_type& n) : _coefficient(0), _multi_index(n,0) { }
      /*! \brief Construct a monomial from a real number and an array of indices. */
      Monomial(const real_type& a, const array<size_type>& i) : _coefficient(a), _multi_index(i) { }
      /*! \brief Convert from a monomial with a different real type. */
      template<class Rl>
      Monomial(const Monomial<Rl>& t) : _coefficient(t.coefficient()), _multi_index(t.multi_index()) { }
     
      /*! \brief The size of the argument.
       *
       *  Includes variables with index 0.
       */
      size_type argument_size() const { return _multi_index.size(); }
      /*! \brief The coefficient of the monomial. */
      const real_type& coefficient() const { return _coefficient; }
      /*! \brief The index of the \a j th variable. */
      const size_type& index(const size_type& j) const { return _multi_index[j]; }
      /*! \brief The multi index of the monomial. */
      const array<size_type>& multi_index() const { return _multi_index; }
      /*! \brief The degree of the monomial, equal to the sum of the indices of the multi-index. */
      size_type degree() const;
      
      /*! \brief Compute the image of a point under the monomial. */
      F evaluate(const LinearAlgebra::Vector<F>& pt) const;
 
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an intput stream. */
      std::istream& read(std::istream& is);
    private:
      friend class PolynomialFunction<R>;
      friend class PolynomialModel<R>;
     private:
      void _resize(const size_type& n);
      real_type _coefficient;
      array<size_type> _multi_index;
    };
    
    /*! \brief Unary minus. */
    template<class R>
    inline Monomial<R> operator-(const Monomial<R>& m) 
    {
      return Monomial<R>(-m.coefficient(), m.multi_index());
    }
    
    /*! \brief A polynomial in several variables. */
    template<class R>
    class PolynomialFunction {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef typename Numeric::traits<R>::arithmetic_type result_type;
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
     public:
      /*! \brief Construct from a string literal. */
      PolynomialFunction(const std::string &);
      /*! \brief Default constructor. */
      PolynomialFunction() : _argument_size(0), _terms() { }
      /*! \brief The zero polynomial in \a n variables. */
      PolynomialFunction(const size_type& m) : _argument_size(m), _terms() { }

      /*! \brief The size of the argument. */
      size_type argument_size() const { return _argument_size; }
      /*! \brief The \a j th term. */
      const Monomial<R>& term(size_type j) const { return _terms[j]; }
      /*! \brief The total number of (nonzero) terms. */
      size_type number_of_terms() const { return _terms.size(); }
      
      /*! \brief Evaluate the polynomial at \a s. */
      F evaluate(const LinearAlgebra::Vector<F>& s) const;

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an intput stream. */
      std::istream& read(std::istream& is);
     private:
      void _sort();
      void _set_argument_size(const size_type& n);
      size_type _compute_maximum_term_size() const;
     private:
      friend class PolynomialModel<R>;
      friend class PolynomialModel<typename Numeric::traits<R>::number_type>;
     private:
      /* Polynomials's terms. */
      size_type _argument_size;
      std::vector< Monomial<R> > _terms;
    };

    /*! \brief A polynomial map with multivalued output.
     *  \ingroup FunctionModel
     */
    template<class R>
    class PolynomialModel
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
     public:
      /*! \brief Construct from a string literal. */
      PolynomialModel(const std::string& s);
      /*! \brief The zero polynomial map in \a n variables with \a m sizeal image. */
      PolynomialModel(const size_type& m, const size_type& n)
        : _argument_size(n), _components(m,PolynomialFunction<R>(n)) { }
      /*! \brief Construct from an array of polynomials. */
      PolynomialModel(const array< PolynomialFunction<R> >& c) : _components(c) { 
        this->_set_argument_size(this->_compute_maximum_component_size()); }
      /*! \brief Copy constructor. */
      PolynomialModel(const PolynomialModel<R>& pm)
        : _argument_size(pm._argument_size), _components(pm._components) { }
        
      /*! \brief The \a i th component polynomial. */
      const PolynomialFunction<R>& component(size_type i) const { return _components[i]; }
      /*! \brief The size of the argument. */
      size_type argument_size() const { return _argument_size; }
      /*! \brief The size of the result. */
      size_type result_size() const { return _components.size(); }
      /*! \brief The size of the result. */
      smoothness_type smoothness() const { return std::numeric_limits<smoothness_type>::max(); }
      
      /*! \brief Compute the image of a point under the polynomial map. */
      LinearAlgebra::Vector<F> evaluate(const LinearAlgebra::Vector<F>& s) const;
      
      /*! \brief Compute the derivate of the map at a point. */
      LinearAlgebra::Matrix<F> jacobian(const LinearAlgebra::Vector<F>& s) const;

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an intput stream. */
      std::istream& read(std::istream& is);
     private:
      void _compute_jacobian() const;
      void _set_argument_size(const size_type& n);
      size_type _compute_maximum_component_size() const;
     private:
      /* Components of the map. */
      size_type _argument_size;
      array< PolynomialFunction<R> > _components;
    };
    
  
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const Monomial<R>& m) {
      return m.write(os);
    }

    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const PolynomialFunction<R>& p) {
      return p.write(os);
    }

    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const PolynomialModel<R>& pm) {
      return pm.write(os);
    }

    template<class R> inline
    std::istream& operator>>(std::istream& is, Monomial<R>& m) {
      return m.read(is);
    }

    template<class R> inline
    std::istream& operator>>(std::istream& is, PolynomialFunction<R>& p) {
      return p.read(is);
    }

    template<class R> inline
    std::istream& operator>>(std::istream& is, PolynomialModel<R>& pm) {
      return pm.read(is);
    }



  }
}

#endif /* ARIADNE_POLYNOMIAL_MODEL_H */
