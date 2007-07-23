/***************************************************************************
 *            simplex.h
 *
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
 
/*! \file simplex.h
 *  \brief Simplices.
 */

#ifndef ARIADNE_SIMPLEX_H
#define ARIADNE_SIMPLEX_H

#include <iosfwd>

#include "../base/array.h"
#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/polytope.h"

namespace Ariadne {
  namespace Geometry {

    class basic_set_tag;

    /*!\ingroup BasicSet
     * \brief A simplex of arbitrary dimension.
     *
     * A simplex is represented by its corner vertices. (We could instead define a
     * simplex by the inequalities describing its faces, but this is less 
     * convenient for simplicial complexes.) Since the representation is the
     * same as that of a polyhedron, the %Simplex class is a subclass of Polytope.
     *
     * The contains(const Point<R>&) operation is much more efficient for simplices
     * than general polytopes, since it can be performed by a simple matrix inversion
     * without the need to solve a linear programming problem.
     */
    template<class R>
    class Simplex : public Polytope<R>
    {
     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the simplex. */
      typedef Point<R> state_type;
     
     public:
      //@{
      //! \name Constructors
      /*! \brief Default constructor constructs simplex of dimension \a 0. */
      Simplex();
    
      /*! \brief Construct from matrix giving the vertices in column form. */
      explicit Simplex(const LinearAlgebra::Matrix<R>& A);
      
      /*! \brief Construct from list of vertices. */
      explicit Simplex(const PointList<R>& v);
      
      /*! \brief Copy constructor. */
      Simplex(const Simplex<R>& s);

      /*! \brief Copy assignment operator. */
      Simplex<R>& operator=(const Simplex<R>& s);
      //@}
    
      //! \name Geometric predicates
      //! \brief Specialized containment predicate
      tribool contains(const Point<R>& pt) const;
      //@}
      
      //@{ 
      //! \name Input/output operations
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      //@}
     public:
      LinearAlgebra::Vector<typename Numeric::traits<R>::arithmetic_type> coordinates(const Point<R>& s) const;
    };
  
    template<class R> std::ostream& operator<<(std::ostream& os, const Simplex<R>& s);
    template<class R> std::istream& operator>>(std::istream& is, Simplex<R>& s);
    
  }
}

#include "simplex.inline.h"

#endif /* ARIADNE_SIMPLEX_H */
