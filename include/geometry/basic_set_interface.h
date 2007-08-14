/***************************************************************************
 *            basic_set_interface.h
 *
 *  Copyright 2007  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file basic_set_interface.h
 *  \brief Abstract interface for basic set classes.
 */


namespace Ariadne {
  namespace Geometry {


    /*!\ingroup BasicSet
     * \brief Interface for basic sets to be used in polymorphic algorithms.
     * 
     * Basic sets implementing the interface are described by their relation to rectangles 
     * in the state space.
     */
    template<class R>
    class BasicSetInterface
    {
     public:
      //!\name Constructors and assignment operators
    
      /*!\brief Cloning operation. */
      virtual BasicSetInterface<R>* clone() const = 0;
    
      //!\name Geometric operations and predicates
      /*!\brief The dimension of the Euclidean space the set lies in. */
      virtual dimension_type dimension() const = 0;

      /*!\brief A rectangle containing the given set */
      virtual Rectangle<R> bounding_box() const = 0;

      /*!\brief Tests if a point is an element of the set. */
      virtual tribool contains(const Point<R>& pt) const = 0;
      
      /*!\brief Tests if a rectangle is disjoint from the set. */
      virtual tribool disjoint(const Rectangle<R>& r) const = 0;
      
      /*!\brief Tests if a point is an element of the set. */
      virtual tribool intersects(const Rectangle<R>& r) const = 0;
      
      /*!\brief Tests if a point is an element of the set. */
      virtual tribool subset(const Rectangle<R>& r) const = 0;
      
      /*!\brief Tests if a point is an element of the set. */
      virtual tribool superset(const Rectangle<R>& r) const = 0;
      
      //!\name Input/output operations
      /*!\brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const = 0;
    };
  }
}
