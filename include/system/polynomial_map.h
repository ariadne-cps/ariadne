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

#include "base/array.h"
#include "linear_algebra/matrix.h"
#include "function/polynomial_function.h"
#include "system/map.h"

namespace Ariadne {
  namespace System {
    

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
        : _polynomial(m,n,0) { }
      /*! \brief Construct from an array of polynomials. */
      PolynomialMap(const Function::PolynomialFunction<R>& p) : _polynomial(p) { }
      /*! \brief Copy constructor. */
      PolynomialMap(const PolynomialMap<R>& pm)
        : _polynomial(pm._polynomial) { }
        
      /*! \brief Returns a pointer to a dynamically-allocated copy of the map. */
      virtual PolynomialMap<R>* clone() const { return new PolynomialMap<R>(*this); }
      
      /*! \brief The \a i th component polynomial. */
      Function::PolynomialFunction<R> component(size_type i) const { return _polynomial.component(i); }
      /*! \brief The dimension of the argument. */
      virtual dimension_type argument_dimension() const { return _polynomial.argument_size(); }
      /*! \brief The dimension of the result. */
      virtual dimension_type result_dimension() const { return _polynomial.result_size(); }
      /*! \brief The dimension of the result. */
      virtual smoothness_type smoothness() const { return std::numeric_limits<smoothness_type>::max(); }
      
      /*! \brief Compute the image of a point under the polynomial map. */
      virtual Geometry::Point<F> image(const Geometry::Point<F>& s) const;
      
      /*! \brief Compute the derivate of the map at a point. */
      virtual LinearAlgebra::Matrix<F> jacobian(const Geometry::Point<F>& s) const;
     
       /*! \brief The name of the map class. */
      virtual std::string name() const { return "PolynomialMap"; }

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an intput stream. */
      virtual std::istream& read(std::istream& is);
     private:
      /* Components of the map. */
      Function::PolynomialFunction<R> _polynomial;
    };
    

  }
}

#endif /* ARIADNE_POLYNOMIAL_MAP_H */
