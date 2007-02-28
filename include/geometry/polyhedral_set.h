/***************************************************************************
 *            polyhedral_set.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file polyhedral_set.h
 *  \brief Wrapper for polyhedron class conforming to the set interface.
 */

#ifndef ARIADNE_POLYHEDRAL_SET_H
#define ARIADNE_POLYHEDRAL_SET_H

#include <iosfwd>

#include "../declarations.h"

#include "../base/tribool.h"

#include "geometry/set.h"
#include "geometry/polyhedron.h"

namespace Ariadne {
  namespace Geometry {
    
 
    
    //! \ingroup ExactSet
    /*! \brief An adaptor for the Polyhedron class conforming to the Set interface. */
    template<class R>
    class PolyhedralSet : public Set<R>, public Polyhedron<R>
    {
     public:
      PolyhedralSet(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b)
        : Set<R>(), Polyhedron<R>(A,b) { }
      PolyhedralSet(const Rectangle<R>& r)
        : Set<R>(), Polyhedron<R>(r) { }
      PolyhedralSet(const Polyhedron<R>& ph)
        : Set<R>(), Polyhedron<R>(ph) { }
      
      virtual ~PolyhedralSet<R>() { }
      virtual PolyhedralSet<R>* clone() const { return new PolyhedralSet<R>(*this); }
      virtual dimension_type dimension() const { return Polyhedron<R>::dimension(); }
      virtual tribool contains(const Point<R>& pt) const { return Polyhedron<R>::contains(pt); }
      virtual Rectangle<R> bounding_box() const { return Polyhedron<R>::bounding_box(); }      
      virtual tribool disjoint(const Rectangle<R>& r) const { 
        return Geometry::disjoint(r,static_cast<const Polyhedron<R>&>(*this)); }
      virtual tribool superset(const Rectangle<R>& r) const { 
        return Geometry::subset(r,static_cast<const Polyhedron<R>&>(*this)); }
      virtual tribool subset(const Rectangle<R>& r) const { 
        return Geometry::subset(static_cast<const Polyhedron<R>&>(*this),r); }
      virtual std::ostream& write(std::ostream& os) const {
        return Polyhedron<R>::write(os);
      }
    };
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const PolyhedralSet<R>& pset) {
      return pset.write(os);
    }
    
  }
}

#endif /* ARIADNE_POLYHEDRAL_SET_H */
