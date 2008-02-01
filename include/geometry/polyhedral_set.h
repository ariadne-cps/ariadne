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

#include "base/tribool.h"

#include "geometry/set_interface.h"
#include "geometry/polyhedron.h"
#include "geometry/constraint_set.h"
#include "function/affine_function.h"

namespace Ariadne {
  namespace Geometry {
    
 
    
    //! \ingroup ExactSet
    /*! \brief An adaptor for the Polyhedron class conforming to the SetInterface interface. */
    template<class R>
    class PolyhedralSet
      : public ConstraintSet<R>
    {
     public:
      /*! \brief */
      PolyhedralSet(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b)
        : ConstraintSet<R>(make_constraint_set(Polyhedron<R>(A,b))), _polyhedron(A,b) { }
      /*! \brief */
      PolyhedralSet(const std::string& str)
        : ConstraintSet<R>(make_constraint_set(Polyhedron<R>(str))), _polyhedron(str) { }
      /*! \brief */
      PolyhedralSet(const Box<R>& r)
        : ConstraintSet<R>(make_constraint_set(Polyhedron<R>(r))), _polyhedron(r) { }
      /*! \brief */
      PolyhedralSet(const Polytope<R>& p)
        : ConstraintSet<R>(make_constraint_set(Polyhedron<R>(p))), _polyhedron(p) { }
      /*! \brief */
      PolyhedralSet(const Polyhedron<R>& ph)
        : ConstraintSet<R>(make_constraint_set(ph)), _polyhedron(ph) { }
      /*! \brief */
      PolyhedralSet(const PolyhedralSet<R>& ph)
        : ConstraintSet<R>(ph), _polyhedron(static_cast<const Polyhedron<R>&>(ph)) { }
      /*! \brief */
      operator const Polyhedron<R>& () const { return this->_polyhedron; }

      /*! \brief */
      virtual ~PolyhedralSet<R>() { }
      /*! \brief */
      virtual PolyhedralSet<R>* clone() const { return new PolyhedralSet<R>(this->_polyhedron); }
      /*! \brief */
      virtual dimension_type dimension() const { return this->_polyhedron.dimension(); }
      /*! \brief */
      virtual tribool contains(const Point<R>& pt) const { return this->_polyhedron.contains(pt); }
      /*! \brief */
      virtual tribool superset(const Box<R>& r) const { 
        return Geometry::subset(r,this->_polyhedron); }
      /*! \brief */
      virtual tribool intersects(const Box<R>& r) const { 
        return !Geometry::disjoint(r,this->_polyhedron); }
      /*! \brief */
      virtual tribool disjoint(const Box<R>& r) const { 
        return Geometry::disjoint(r,this->_polyhedron); }
      /*! \brief */
      virtual tribool subset(const Box<R>& r) const { 
        return Geometry::subset(this->_polyhedron,r); }
      /*! \brief */
      virtual tribool bounded() const { return this->_polyhedron.bounded(); }      
      /*! \brief */
      virtual Box<R> bounding_box() const { return this->_polyhedron.bounding_box(); }      
      /*! \brief */
      virtual std::ostream& write(std::ostream& os) const {
        return this->_polyhedron.write(os);
      }
     private:
      static ConstraintSet<R> make_constraint_set(const Polyhedron<R>& ply) {
        return ConstraintSet<R>(Function::AffineFunction<R>(-ply.A(),ply.b()),Box<R>::positive_orthant(ply.number_of_constraints())); }
     private:
      Polyhedron<R> _polyhedron;
    };
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const PolyhedralSet<R>& pset) {
      return pset.write(os);
    }
    
    
  }
}

#endif /* ARIADNE_POLYHEDRAL_SET_H */
