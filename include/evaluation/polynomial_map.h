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

#ifndef _ARIADNE_POLYNOMIAL_MAP_H
#define _ARIADNE_POLYNOMIAL_MAP_H

#include <iosfwd>
#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "linear_algebra/linear_algebra.h"

#include "base/utility.h"
#include "base/array.h"
#include "base/interval.h"

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"

namespace Ariadne {
  namespace Evaluation {
    template<typename R> class PolynomialMap;

    template<typename R> std::ostream& operator<<(std::ostream&, const PolynomialMap<R>&);
    
    /*! \brief A monomial in several variables. */
    template<typename R>
    class Monomial {
     public:
      /*! \brief The unsigned integer type used to denote the array positions. */
      typedef size_t size_type;
      /*! \brief The type of denotable real number used for the coefficients. */
      typedef R Real;
      /*! \brief The type of denotable state accepted as argument. */
      typedef Geometry::Point<R> State;
     public:
      Monomial(size_type n) : _coefficient(0), _multi_index(n,0) { }
      Monomial(Real a, array<size_type>& i) : _coefficient(a), _multi_index(i) { }
     
      size_type argument_dimension() const { return _multi_index.size(); }
      const Real& coefficient() const { return _coefficient; }
      const array<size_type>& multi_index() const { return _multi_index; }
      const size_type& index(size_type n) const { return _multi_index[n]; }
    
      Real operator() (const State& s) const {
        assert(s.size() == this->dimension());
        Real result=_coefficient;
        for(size_type k=0; k!=this->dimension(); ++k) {
          result *= s[k]^_multi_index[k];
        }
        return result;
      }
      
     private:
      Real _coefficient;
      array<size_type> _multi_index;
    };
     
    /*! \brief A polynomial in several variables. */
    template<typename R>
    class Polynomial {
     public:
      /*! \brief The unsigned integer type used to denote the array positions. */
      typedef size_t size_type;
      /*! \brief The type of denotable real number used for the corners. */
      typedef R Real;
      /*! \brief The type of denotable state contained by the simplex. */
      typedef Geometry::Point<R> State;
     private:
      /* Simplex's vertices. */
      std::vector< Monomial<Real> > _terms;
     public:
      Polynomial() : _terms(1,Monomial<R>(0)) { }
      Polynomial(size_type m) : _terms(1,Monomial<R>(m)) { }
      Polynomial(const std::vector< Monomial<Real> >& t) : _terms(t) { }
      Polynomial(const std::vector<Real>& , std::vector< array<size_type> >&);

      size_type argument_dimension() const { return _terms[0].argument_dimension(); }
      const Monomial<R>& term(size_type j) const { return _terms[j]; }
      size_type number_of_terms() const { return _terms.size(); }
      
      Real operator() (const State& s) const {
        assert(s.size() == this->argument_dimension());
        Real result=0;
        for(size_type j=0; j!=term.size(); ++j) {
          result += _terms[j](s);
        }
        return result;
      }
    };

    /*! \brief A polynomial map with multivalued output. */
    template <typename R>
    class PolynomialMap {
     public:
      /*! \brief The unsigned integer type used to denote the array positions. */
      typedef size_t size_type;
      /*! \brief The type of denotable real number used for the corners. */
      typedef R Real;
      /*! \brief The type of denotable state contained by the simplex. */
      typedef Geometry::Point<R> State;
     private:
      /* Components of the map. */
      array< Polynomial<R> > _components;
     public:
      PolynomialMap(size_type m, size_type n) : _components(n) { _components[0] = Polynomial<R>(m); }
      PolynomialMap(const array< Polynomial<R> >& c) : _components(c) { }
      
      const Polynomial<R>& component(size_type i) const { return _components[i]; }
      const array< Polynomial<R> >& components() const { return _components; }
      size_type argument_dimension() const { return _components[0].argument_dimension(); }
      size_type result_dimension() const { return _components.size(); }
      
      State operator() (const State& s) const {
        assert(s.dimension() == this->argument_dimension());
        State result(this->result_dimension());
        for(size_type i=0; i!=this->result_dimension(); ++i) {
          result[i] = _components[i](s);
        }
        return result;
      }
    };
    
    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Monomial<R>& m) 
    {
      return os;
    }

    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Polynomial<R>& p) 
    {
      for(size_type j=0; j!=p.number_of_terms(); ++j) {
        const Monomial<R>& m = p.term(j);
        os << " + " << m.coefficient();
        for(size_type k=0; k!=m.argument_dimension(); ++k) {
          os << " * x" << k << "^" << m.index(k);
        }
      }
      return os;
    }

    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const PolynomialMap<R>& p) 
    {
      return os << p.components(); 
    }
    
  }
}

#endif /* _ARIADNE_POLYNOMIAL_MAP_H */
