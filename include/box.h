/***************************************************************************
 *            box.h
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
 *
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

/*! \file box.h
 *  \brief Coordinate-aligned boxes in Euclidean space.
 */

#ifndef ARIADNE_BOX_H
#define ARIADNE_BOX_H

#include "numeric.h"
#include "vector.h"

#include "set_interface.h"
#include "graphics_interface.h"
#include "point.h"

namespace Ariadne {

class Point;
class Box;

Box widen(const Box& bx);

//! \ingroup BasicSetSubModule GeometryModule
//! \brief A box in Euclidean space.
class Box
    : public SetInterface,
      public DrawableInterface,
      public Vector<Interval>
{
  public:
    typedef Float real_type;

    //! Construct a singleton point in zero dimensions.
    Box() : Vector<Interval>() { }
    //! Construct from an integer giving the dimension and a list of floating-point values
    //! giving alternately lower and upper bounds.
    explicit Box(uint d, double x0l, double x0u, ...);

    // Templated constructor; useful for automatic conversion from vector expressions
    template<class T> Box(const T& t) : Vector<Interval>(t) { }
    template<class T1, class T2> Box(const T1& t1, const T2& t2) : Vector<Interval>(t1,t2) { }

    //! Construct from a string literal of the form "[a1,b1]x[a2,b2]x...x[ad,bd]".
    explicit Box(const std::string& str);

    //! The unit box \f$[-1,1]^n\f$ in \a n dimensions.
    static Box unit_box(uint n) {
        return Box(n,Interval(-1,1));
    }

    //! The upper quadrant box \f$[0,infty]^n\f$ in \a n dimensions.
    static Box upper_quadrant(uint n) {
        return Box(n,Interval(0,inf<Float>()));
    }

    //! An explicit case to an interval vector. Useful to prevent ambiguous function overloads.
    const Vector<Interval>& vector() const { return *this; }

    //! The set of vertices of the box.
    std::vector<Point> vertices() const;

    //! An approximation to the centre of the box.
    Point centre() const {
        return Point(midpoint(static_cast<const Vector<Interval>&>(*this)));
    }

    //! The radius of the box in the supremum norm.
    Float radius() const {
        Float dmax=0;
        for(uint i=0; i!=this->size(); ++i) {
            dmax = max( dmax, (*this)[i].width() );
        }
        return up(dmax/2);
    }

    //! An approximation to the Lesbegue measure (area, volume) of the box.
    Float measure() const {
        Float meas=1;
        for(uint i=0; i!=this->size(); ++i) {
            meas *= (*this)[i].width();
        }
        return meas;
    }

    //! \brief Test if the box is empty.
    bool empty() const {
        return Ariadne::empty(*this);
    }

    //! \brief Test if the box is bounded.
    bool bounded() const {
        for(uint i=0; i!=this->Vector<Interval>::size(); ++i) {
            if(!Ariadne::bounded((*this)[i])) {
                return false;
            }
        }
        return true;
    }

    //! \brief Make a dynamically-allocated copy.
    virtual Box* clone() const {
        return new Box(*this);
    }

    //! \brief Test if the box is a superset of another box.
    bool contains(const Point& pt) const {
        return Ariadne::contains(this->vector(),pt.vector());
    }

    //! \brief Test if the box is a subset of another box.
    bool subset(const Box& bx) const {
        return Ariadne::subset(this->vector(),bx.vector());
    }

    //! \brief Test if the box is a superset of another box.
    bool superset(const Box& bx) const {
        return Ariadne::subset(bx.vector(),this->vector());
    }

    //! \brief Test if the box intersects another box. Returns \a true even
    //! if the intersection occurs only on the boundary of the boxes.
    //! Use Box::overlaps(const Box& bx) to test robust intersection
    //! of the interiors.
    bool intersects(const Box& bx) const {
        return Ariadne::intersect(this->vector(),bx.vector());
    }

    //! \brief Test if the box intersects another box. Returns \a false even
    //! if the intersection only occurs only on the boundary of the boxes.
    bool disjoint(const Box& bx) const {
        return Ariadne::disjoint(this->vector(),bx.vector());
    }

    //! \brief The dimension of the space the box lies in.
    virtual uint dimension() const {
        return this->size();
    }

    //! \brief Tests if the box is disjoint from another box.
    //! Only returns true if the closures are disjoint.
    virtual tribool separated(const Box& other) const {
        return Ariadne::disjoint(this->vector(), other.vector());
    }

    //! \brief Tests if the box overlaps another box.
    //! Only returns true if the interior of one box intersects the other.
    virtual tribool overlaps(const Box& other) const {
        return Ariadne::overlap(this->vector(), other.vector());
    }

    //! \brief Tests if the box covers another box.
    //! Only returns true if the interior of the box is a superset
    //! of the closure of the other.
    virtual tribool covers(const Box& other) const {
        return Ariadne::inside(other.vector(), this->vector());
    }

    //! \brief Tests if the box covers another box.
    //! Only returns true if the closure of the box is a subset
    //! of the interior of the other.
    virtual tribool inside(const Box& other) const {
        return Ariadne::inside(this->vector(), other.vector());
    }

    //! \brief Returns an enclosing bounding box for the set.
    //! The result is guaranteed to have nonempty interior, and floating-point
    //! boundary coefficients so the centre and radius are exactly computable.
    virtual Box bounding_box() const {
        return Ariadne::widen(*this);
    }

    //! \brief Widens the box by the minimal floating-point increment.
    //! The result is guaranteed to contain the box in its interior.
    Void widen() {
        Box& result=*this;
        for(uint i=0; i!=result.dimension(); ++i) {
            result[i]=Ariadne::widen((*this)[i]);
        }
    }

    //! \brief Narrows the box by the minimal floating-point increment.
    //! The result is guaranteed to contain the box in its interior.
    Void narrow() {
        Box& result=*this;
        for(uint i=0; i!=result.dimension(); ++i) {
            result[i]=Ariadne::narrow((*this)[i]);
        }
    }

    //! \brief Split into two along the largest side.
    std::pair<Box,Box> split() const { return Ariadne::split(*this); }
    //! \brief Split into two along side with index \a i.
    std::pair<Box,Box> split(uint i) const { return Ariadne::split(*this,i); };

    //! \brief Draw on a canvas.
    virtual void draw(CanvasInterface& c, const Projection2d& p) const;

    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const {
        return os << *static_cast<const Vector<Interval>*>(this);
    }
};

//! \relates Box \brief The cartesian product of two boxes.
Box product(const Box& bx1, const Box& bx2);
//! \relates Box \brief The smallest box containing the two boxes.
Box hull(const Box& bx1, const Box& bx2);
//! \relates Box \brief The intersection of the two boxes.
Box intersection(const Box& bx1, const Box& bx2);
//! \relates Box \brief A box which is wider than the input, and has single-precision values.
Box widen(const Box& bx);
//! \relates Box \brief A box which is narrower than the input, and has single-precision values.
Box narrow(const Box& bx);

Box make_box(const std::string& str);

} // namespace Ariadne

#endif /* ARIADNE_BOX_H */
