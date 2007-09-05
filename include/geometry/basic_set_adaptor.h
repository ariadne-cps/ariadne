/***************************************************************************
 *            basic_set_adaptor.h
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
 
/*! \file basic_set_adaptor.h
 *  \brief Adaptor for concrete sets to the basic set interface.
 */

#include "basic_set_interface.h"

namespace Ariadne {
  namespace Geometry {

    template<class R> inline
    dimension_type dimension(const Rectangle<R>& r) {
      return r.dimension();
    }

    template<class R> inline
    tribool contains(const Rectangle<R>& r, const Point<R>& pt) {
      return r.contains(pt);
    }

    template<class R> inline
    Rectangle<R> bounding_box(const Rectangle<R>& r) {
      return r;
    }

    template<class R0, class R1> inline
    dimension_type dimension(const Zonotope<R0,R1>& z) {
      return z.dimension();
    }

    template<class R0, class R1> inline
    Rectangle<typename Zonotope<R0,R1>::real_type> bounding_box(const Zonotope<R0,R1>& z) {
      return z.bounding_box();
    }


    
    /*! \ingroup BasicSet
     *  \brief Adaptor for concrete sets to the basic set interface.
     *
     *  Basic set classes in %Ariadne are implemented as concrete types, 
     *  and operations on these types are mostly provided by non-member functions.
     *  Sometimes, it is useful to give an interface without specifying the type 
     *  of basic set used. This adaptor takes a set conforming to the BasicSetConcept
     *  and makes it satisfy the BasicSetInterface. 
     */
    template<class BS>
    class BasicSetAdaptor
      : public BasicSetInterface<typename BS::real_type>,
        public BS // use inheritance so we can use whenever the base type is needed.
    {
      // Private convenience typedef
      typedef typename BS::real_type R;
     public:
      /*! \brief The type of real number used to describle the set. */
      typedef typename BS::real_type real_type;
      /*! \brief The type of point that the set contains. */
      typedef typename BS::state_type state_type;
     public:
      //@{
      //! \name Constructors and assignment operators
    
      /*! \brief Construct from a basic set. */
      BasicSetAdaptor(const BS& bs) : BS(bs) { }
    
      /*! \brief Copy constructor. */
      BasicSetAdaptor(const BasicSetAdaptor<BS>& bs) : BS(bs) { }
    
      /*! \brief Cloning operation. */
      virtual BasicSetAdaptor<BS>* clone() const { return new BasicSetAdaptor<BS>(*this); }
    
      //@}
      
      
      //@{
      //! \name Required geometric operations
      /*! \brief The dimension of the Euclidean space the set lies in. */
      virtual dimension_type dimension() const { return Geometry::dimension(*this); }
      
      /*! \brief A rectangle containing the given set */
      virtual Rectangle<R> bounding_box() const { return Geometry::bounding_box(*this); };

      /*! \brief Tests if a point is an element of the set. */
      virtual tribool contains(const state_type& pt) const { return Geometry::contains(*this,pt); }
      
      /*! \brief Tests if a rectangle is disjoint from the set. */
      virtual tribool disjoint(const Rectangle<R>& r) const { return Geometry::disjoint(*this,r); }
      
      /*! \brief Tests if a point is an element of the set. */
      virtual tribool intersects(const Rectangle<R>& r) const { return !Geometry::disjoint(*this,r); }
      
      /*! \brief Tests if a point is an element of the set. */
      virtual tribool subset(const Rectangle<R>& r) const { return Geometry::subset(*this,r); }
      
      /*! \brief Tests if a point is an element of the set. */
      virtual tribool superset(const Rectangle<R>& r) const { return Geometry::subset(r,*this); }
      
      //@}
      //@{ 
      //! \name Input/output operations
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const { return Geometry::operator<<(os,*this); }
      //@}
    };
  }
}
