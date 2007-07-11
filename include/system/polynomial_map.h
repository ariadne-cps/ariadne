/***************************************************************************
 *            polynomial_map.h
 *
 *  17 January 2006
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file polynomial_map.h
 *  \brief Polynomial maps.
 */

#ifndef ARIADNE_POLYNOMIAL_MAP_H
#define ARIADNE_POLYNOMIAL_MAP_H

#include <iosfwd>
#include <string>
#include <sstream>

#include "../base/array.h"
#include "../linear_algebra/matrix.h"
#include "../system/map.h"

namespace Ariadne {
  namespace System {
    template<class R> class Monomial;
    template<class R> class Polynomial;
    template<class R> class PolynomialMap;
    template<class R> class PolynomialMatrix;

    template<class R> std::ostream& operator<<(std::ostream&, const Monomial<R>&);
    template<class R> std::ostream& operator<<(std::ostream&, const Polynomial<R>&);
    template<class R> std::ostream& operator<<(std::ostream&, const PolynomialMap<R>&);
    template<class R> std::ostream& operator<<(std::ostream&, const PolynomialMatrix<R>&);
   
    template<class R> std::istream& operator>>(std::istream&, Monomial<R>&);
    template<class R> std::istream& operator>>(std::istream&, Polynomial<R>&);
    template<class R> std::istream& operator>>(std::istream&, PolynomialMap<R>&);
   
    /*! \brief Graded lexicographical ordering. */
    template<class R> bool operator<(const Monomial<R>&, const Monomial<R>&);
  }

  namespace Numeric {
    template<class R> struct traits< System::Polynomial<R> > {
      typedef System::Polynomial<typename traits<R>::arithmetic_type> arithmetic_type;
    };
  }
  
  namespace System {
    
    /*! \brief A monomial in several variables. */
    template<class R>
    class Monomial {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef typename Numeric::traits<R>::arithmetic_type result_type;
     public:
      /*! \brief The type of denotable real number used for the coefficients. */
      typedef R real_type;
      /*! \brief The type of denotable state accepted as argument. */
      typedef Geometry::Point<R> state_type;
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
     
      /*! \brief The dimension of the argument.
       *
       *  Includes variables with index 0.
       */
      dimension_type argument_dimension() const { return _multi_index.size(); }
      /*! \brief The coefficient of the monomial. */
      const real_type& coefficient() const { return _coefficient; }
      /*! \brief The index of the \a j th variable. */
      const size_type& index(const dimension_type& j) const { return _multi_index[j]; }
      /*! \brief The multi index of the monomial. */
      const array<size_type>& multi_index() const { return _multi_index; }
      /*! \brief The degree of the monomial, equal to the sum of the indices of the multi-index. */
      size_type degree() const;
      
      /*! \brief Compute the image of a point under the monomial. */
      F apply(const Geometry::Point<F>& pt) const;
     private:
      friend std::istream& operator>> <> (std::istream&, Monomial<R>&);
      friend class Polynomial<R>;
      friend class PolynomialMap<R>;
     private:
      void _resize(const dimension_type& n);
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
    class Polynomial {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef typename Numeric::traits<R>::arithmetic_type result_type;
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type of denotable state contained by the simplex. */
      typedef Geometry::Point<R> state_type;
     public:
      /*! \brief Construct from a string literal. */
      Polynomial(const std::string &);
      /*! \brief Default constructor. */
      Polynomial() : _argument_dimension(0), _terms() { }
      /*! \brief The zero polynomial in \a n variables. */
      Polynomial(const size_type& m) : _argument_dimension(m), _terms() { }

      /*! \brief The dimension of the argument. */
      dimension_type argument_dimension() const { return _argument_dimension; }
      /*! \brief The \a j th term. */
      const Monomial<R>& term(size_type j) const { return _terms[j]; }
      /*! \brief The total number of (nonzero) terms. */
      size_type number_of_terms() const { return _terms.size(); }
      
      /*! \brief Compute the image of a point under the polynomial. */
      F image(const Geometry::Point<F>& s) const;
     private:
      void _sort();
      void _set_argument_dimension(const dimension_type& n);
      dimension_type _compute_maximum_term_dimension() const;
     private:
      friend class PolynomialMap<R>;
      friend class PolynomialMap<typename Numeric::traits<R>::number_type>;
      friend std::istream& operator>> <> (std::istream&, Polynomial<R>&);
     private:
      /* Polynomials's terms. */
      dimension_type _argument_dimension;
      std::vector< Monomial<R> > _terms;
    };

    /*! \brief A polynomial map with multivalued output.
     *  \ingroup DiscreteTime
     */
    template<class R>
    class PolynomialMap : public MapInterface<R> {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef Geometry::Point<F> result_type;
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type of denotable state contained by the simplex. */
      typedef Geometry::Point<R> state_type;
     public:
      /*! \brief Construct from a string literal. */
      PolynomialMap(const std::string& s);
      /*! \brief The zero polynomial map in \a n variables with \a m dimensional image. */
      PolynomialMap(const dimension_type& m, const dimension_type& n)
        : _argument_dimension(n), _components(m,Polynomial<R>(n)) { }
      /*! \brief Construct from an array of polynomials. */
      PolynomialMap(const array< Polynomial<R> >& c) : _components(c) { 
        this->_set_argument_dimension(this->_compute_maximum_component_dimension()); }
      /*! \brief Copy constructor. */
      PolynomialMap(const PolynomialMap<R>& pm)
        : _argument_dimension(pm._argument_dimension), _components(pm._components) { }
        
      /*! \brief Returns a pointer to a dynamically-allocated copy of the map. */
      PolynomialMap<R>* clone() const { return new PolynomialMap<R>(*this); }
      
      /*! \brief The \a i th component polynomial. */
      const Polynomial<R>& component(size_type i) const { return _components[i]; }
      /*! \brief The dimension of the argument. */
      dimension_type argument_dimension() const { return _argument_dimension; }
      /*! \brief The dimension of the result. */
      dimension_type result_dimension() const { return _components.size(); }
      /*! \brief The dimension of the result. */
      size_type smoothness() const { return (size_type) -1;; }
      
      /*! \brief Compute the image of a point under the polynomial map. */
      virtual Geometry::Point<F> image(const Geometry::Point<F>& s) const;
      
      /*! \brief Compute the derivate of the map at a point. */
      virtual LinearAlgebra::Matrix<F> jacobian(const Geometry::Point<F>& s) const;
     
      /*! \brief Compute a closed form for the derivative of the map. */
      const PolynomialMatrix<F>& jacobian() const;
      
      std::string name() const { return "PolynomialMap"; }
     private:
      void _compute_jacobian() const;
      void _set_argument_dimension(const dimension_type& n);
      dimension_type _compute_maximum_component_dimension() const;
     private:
      friend std::istream& operator>> <> (std::istream&, PolynomialMap<R>&);
      friend std::ostream& operator<< <> (std::ostream&, const PolynomialMap<R>&);
     private:
      /* Components of the map. */
      dimension_type _argument_dimension;
      array< Polynomial<R> > _components;
      mutable PolynomialMatrix<F> _jacobian;
    };
    
    /*! \brief A matrix with polynomial entries. */
    template<class R>
    class PolynomialMatrix {
      typedef typename Numeric::traits<R>::arithmetic_type F;
     public:
      /*! \brief Default constructor creates a 0 by 0 matrix. */
      PolynomialMatrix() : _matrix() { }
      /*! \brief Construct an \a m by \a n matrix, all of whose entries are zero. */
      PolynomialMatrix(const size_type& r, const size_type& c) : _matrix(r,c) { }
     
      /*! \brief The number of rows of the matrix. */
      size_type number_of_rows() const { return _matrix.number_of_rows(); }
      /*! \brief The number of columns of the matrix. */
      size_type number_of_columns() const { return _matrix.number_of_columns(); }

      /*! \brief A reference to the (\a i,\a j)th entry. */
      Polynomial<R>& operator() (const size_type& i, const size_type& j) { return this->_matrix(i,j); }
      /*! \brief A constant reference to the (\a i,\a j)th entry. */
      const Polynomial<R>& operator() (const size_type& i, const size_type& j) const { return this->_matrix(i,j); }
     private:     
      void set(const size_type& i, const size_type& j, const Polynomial<R>& p) { this->_matrix(i,j)=p; }
      const Polynomial<R>& get(const size_type& i, const size_type& j) const { return _matrix(i,j); }
     private:
      friend class PolynomialMap<R>;
     private:
      LinearAlgebra::Matrix< Polynomial<R> > _matrix;
    };
    

  }
}

#endif /* ARIADNE_POLYNOMIAL_MAP_H */
