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

#include "utility/container.h"

#include "utility/declarations.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "geometry/interval.h"

#include "geometry/set_interface.h"
#include "output/graphics_interface.h"
#include "geometry/point.h"

namespace Ariadne {

typedef OutputStream OutputStream;
typedef Nat Nat;

class IntervalSet;
class BoxSet;

template<class X> class Point;
typedef Point<ExactNumber> ExactPoint;

template<class IVL> class Box;
typedef Box<ExactInterval> ExactBox;
typedef Box<UpperInterval> UpperBox;
typedef Box<ApproximateInterval> ApproximateBox;

typedef ExactBox BoxDomainType;
typedef UpperBox BoundingBoxType;

Bool element(const Vector<ExactNumber>& v1, const Vector<ExactInterval>& v2);
Bool element(const Vector<ValidatedNumber>& v1, const Vector<ExactInterval>& v2);
Bool element(const Vector<ApproximateNumber>& v1, const Vector<ExactInterval>& v2);
Bool contains(const Vector<ExactInterval>& v1, const Vector<ExactNumber>& v2);
Bool contains(const Vector<ExactInterval>& v1, const Vector<ValidatedNumber>& v2);
Bool contains(const Vector<ExactInterval>& v1, const Vector<ApproximateNumber>& v2);

Bool intersect(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2);
Bool subset(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2);

Bool disjoint(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2);
Bool overlap(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2);
Bool inside(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2);
Bool covers(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2);
Bool empty(const Vector<ExactInterval>& v);

Vector<ExactInterval> split(const Vector<ExactInterval>& v, Nat k, Tribool lr);
Vector<ExactInterval> split(const Vector<ExactInterval>& v, Tribool lr);
Pair< Vector<ExactInterval>, Vector<ExactInterval> > split(const Vector<ExactInterval>& v);
Pair< Vector<ExactInterval>, Vector<ExactInterval> > split(const Vector<ExactInterval>& v, Nat k);

Vector<ExactInterval> hull(const Vector<ExactNumber>& v1, const Vector<ExactNumber>& v2);
Vector<ExactInterval> hull(const Vector<ExactInterval>& v1, const Vector<ExactNumber>& v2);
Vector<ExactInterval> hull(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2);
Vector<ExactInterval> intersection(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2);
Vector<ExactNumber> midpoint(const Vector<ExactInterval>& v);
Vector<ExactFloat> lower_bounds(const Vector<ExactInterval>& v);
Vector<ExactFloat> upper_bounds(const Vector<ExactInterval>& v);
PositiveUpperNumber radius(const Vector<ExactInterval>& z);
PositiveUpperNumber volume(const Vector<ExactInterval>& z);

Vector<UpperInterval> hull(const Vector<UpperInterval>& v1, const Vector<UpperInterval>& v2);
Vector<UpperInterval> intersection(const Vector<UpperInterval>& v1, const Vector<UpperInterval>& v2);
Tribool subset(const Vector<UpperInterval>& bx1, const Vector<ExactInterval>& bx2);
Tribool disjoint(const Vector<UpperInterval>& bx1, const Vector<UpperInterval>& bx2);
Tribool inside(const Vector<UpperInterval>& v1, const Vector<ExactInterval>& v2);




//! \ingroup GeometryModule ExactSetSubModule
//! \brief An exact interval in \f$\mathbb{R}\f$.
class IntervalSet {
    Real _lower, _upper;
  public:
    IntervalSet() : _lower(-1), _upper(+1) { }
    IntervalSet(const Real& l, const Real& u) : _lower(l), _upper(u) { }
    const Real& lower() const { return _lower; }
    const Real& upper() const { return _upper; }
    const Real midpoint() const { return (_lower+_upper)/2; }
    const Real radius() const { return (_upper-_lower)/2; }
};
inline OutputStream& operator<<(OutputStream& os, const IntervalSet& ivl) {
    return os << "{" << ivl.lower() << ":" << ivl.upper() << "}";
}
inline ExactInterval under_approximation(const IntervalSet& rivl) {
    return ExactInterval(rivl.lower().upper().raw(),rivl.upper().lower().raw());
}
inline ExactInterval over_approximation(const IntervalSet& rivl) {
    return ExactInterval(rivl.lower().lower().raw(),rivl.upper().upper().raw());
}
inline ExactInterval approximation(const IntervalSet& rivl) {
    return ExactInterval(rivl.lower().approx().raw(),rivl.upper().approx().raw());
}


//! \ingroup GeometryModule ExactSetSubModule
//! \brief An exact coordinate-aligned box in \f$\mathbb{R}^n\f$.
class BoxSet {
    Array<IntervalSet> _ary;
  public:
    BoxSet() : _ary() { }
    explicit BoxSet(const Vector<ExactInterval>& iv);
    explicit BoxSet(Nat n) : _ary(n,IntervalSet()) { }
    BoxSet(const List<IntervalSet>& t) : _ary(t.begin(),t.end()) { }
    BoxSet(Nat n, const IntervalSet& ivl) : _ary(n,ivl) { }
    Nat size() const { return _ary.size(); }
    Nat dimension() const { return _ary.size(); }
    IntervalSet const& operator[](Nat i) const { return _ary[i]; }
    IntervalSet& operator[](Nat i) { return _ary[i]; }
    friend OutputStream& operator<<(OutputStream& os, const BoxSet& bx) { return os << bx._ary; }
};
ExactBox under_approximation(const BoxSet& rbx);
ExactBox over_approximation(const BoxSet& rbx);
ExactBox approximation(const BoxSet& rbx);

ExactBox widen(const ExactBox& bx);

//! \ingroup BasicSetSubModule GeometryModule
//! \brief A box in Euclidean space.
template<> class Box<ExactInterval>
    : public SetInterface,
      public DrawableInterface,
      public Vector<ExactInterval>
{
  public:
    //! Construct a singleton point in zero dimensions.
    Box<ExactInterval>() : Vector<ExactInterval>() { }
    //! Construct an empty box in \a d dimensions.
    explicit Box<ExactInterval>(Nat d) : Vector<ExactInterval>(d) { }
    //! Construct from an initializer list of pairs of floating-point values
    //! giving lower and upper bounds.
    Box<ExactInterval>(InitializerList<ExactInterval> lst);

    explicit Box<ExactInterval>(Nat d, ExactInterval ivl) : Vector<ExactInterval>(d,ivl) { }
    Box<ExactInterval>(const Vector<ExactFloat>& pvec) : Vector<ExactInterval>(pvec) { }
    Box<ExactInterval>(const Vector<ExactInterval>& ivec) : Vector<ExactInterval>(ivec) { }
    Box<ExactInterval>(const List<ExactInterval>& ilst) : Vector<ExactInterval>(ilst) { }

    //! Construct from a string literal of the form "[a1,b1]x[a2,b2]x...x[ad,bd]".
    explicit Box<ExactInterval>(const StringType& str);

    //! The unit box \f$[-1,1]^n\f$ in \a n dimensions.
    static ExactBox unit_box(Nat n) {
        return ExactBox(n,ExactInterval(-1,1));
    }

    //! The upper quadrant box \f$[0,\infty]^n\f$ in \a n dimensions.
    static ExactBox upper_quadrant(Nat n) {
        return ExactBox(n,ExactInterval(0,inf));
    }

    //! An explicit case to an interval vector. Useful to prevent ambiguous function overloads.
    const Vector<ExactInterval>& vector() const { return *this; }

    //! The set of vertices of the box.
    std::vector<ExactPoint> vertices() const;

    //! An approximation to the centre of the box.
    ExactPoint centre() const {
        return ExactPoint(midpoint(static_cast<const Vector<ExactInterval>&>(*this)));
    }

    //! An over-approximation to radius of the box in the supremum norm.
    ErrorType radius() const {
        ErrorType dmax=0u;
        for(Nat i=0; i!=this->size(); ++i) {
            dmax = max( dmax, (*this)[i].width() );
        }
        return half(dmax);
    }

    //! An approximation to the Lesbegue measure (area, volume) of the box.
    ApproximateNumber measure() const {
        ApproximateNumber meas=1u;
        for(Nat i=0; i!=this->size(); ++i) {
            meas *= (*this)[i].width();
        }
        return meas;
    }

    //! \brief Test if the box is empty.
    Bool empty() const {
        return Ariadne::empty(*this);
    }

    //! \brief Test if the box is bounded.
    Bool bounded() const {
        for(Nat i=0; i!=this->Vector<ExactInterval>::size(); ++i) {
            if(!Ariadne::bounded((*this)[i])) {
                return false;
            }
        }
        return true;
    }

    //! \brief Make a dynamically-allocated copy.
    virtual ExactBox* clone() const {
        return new ExactBox(*this);
    }

    //! \brief Test if the box is a superset of another box.
    Bool contains(const ExactPoint& pt) const {
        return Ariadne::contains(this->vector(),pt.vector());
    }

    //! \brief Test if the box is a subset of another box.
    Bool subset(const ExactBox& bx) const {
        return Ariadne::subset(this->vector(),bx.vector());
    }

    //! \brief Test if the box is a superset of another box.
    Bool superset(const ExactBox& bx) const {
        return Ariadne::subset(bx.vector(),this->vector());
    }

    //! \brief Test if the box intersects another box. Returns \a true even
    //! if the intersection occurs only on the boundary of the boxes.
    //! Use ExactBox::overlaps(const ExactBox& bx) to test robust intersection
    //! of the interiors.
    Bool intersects(const ExactBox& bx) const {
        return Ariadne::intersect(this->vector(),bx.vector());
    }

    //! \brief Test if the box intersects another box. Returns \a false even
    //! if the intersection only occurs only on the boundary of the boxes.
    Bool disjoint(const ExactBox& bx) const {
        return Ariadne::disjoint(this->vector(),bx.vector());
    }

    //! \brief The dimension of the space the box lies in.
    virtual Nat dimension() const {
        return this->size();
    }

    //! \brief Tests if the box is disjoint from another box.
    //! Only returns true if the closures are disjoint.
    virtual Tribool separated(const ExactBox& other) const {
        return Ariadne::disjoint(this->vector(), other.vector());
    }

    //! \brief Tests if the box overlaps another box.
    //! Only returns true if the interior of one box intersects the other.
    virtual Tribool overlaps(const ExactBox& other) const {
        return Ariadne::overlap(this->vector(), other.vector());
    }

    //! \brief Tests if the box covers another box.
    //! Only returns true if the interior of the box is a superset
    //! of the closure of the other.
    virtual Tribool covers(const ExactBox& other) const {
        return Ariadne::inside(other.vector(), this->vector());
    }

    //! \brief Tests if the box covers another box.
    //! Only returns true if the closure of the box is a subset
    //! of the interior of the other.
    virtual Tribool inside(const ExactBox& other) const {
        return Ariadne::inside(this->vector(), other.vector());
    }

    //! \brief Returns an enclosing bounding box for the set.
    //! The result is guaranteed to have nonempty interior, and floating-point
    //! boundary coefficients so the centre and radius are exactly computable.
    virtual inline UpperBox bounding_box() const;

    //! \brief Widens the box by the minimal floating-point increment.
    //! The result is guaranteed to contain the box in its interior.
    Void widen() {
        ExactBox& result=*this;
        for(Nat i=0; i!=result.dimension(); ++i) {
            result[i]=Ariadne::widen((*this)[i]);
        }
    }

    //! \brief Narrows the box by the minimal floating-point increment.
    //! The result is guaranteed to contain the box in its interior.
    Void narrow() {
        ExactBox& result=*this;
        for(Nat i=0; i!=result.dimension(); ++i) {
            result[i]=Ariadne::narrow((*this)[i]);
        }
    }

    //! \brief Split into two along the largest side.
    Pair<ExactBox,ExactBox> split() const { return Ariadne::split(*this); }
    //! \brief Split into two along side with index \a i.
    Pair<ExactBox,ExactBox> split(Nat i) const { return Ariadne::split(*this,i); };

    //! \brief Draw on a canvas.
    virtual Void draw(CanvasInterface& c, const Projection2d& p) const;

    //! \brief Write to an output stream.
    virtual OutputStream& write(OutputStream& os) const {
        return os << *static_cast<const Vector<ExactInterval>*>(this);
    }

    operator UpperBox & () { return reinterpret_cast<UpperBox&>(static_cast<Vector<ExactInterval>&>(*this)); }
    operator UpperBox const& () const { return reinterpret_cast<UpperBox const&>(static_cast<Vector<ExactInterval>const&>(*this)); }
    operator Vector<UpperInterval> const& () const { return reinterpret_cast<Vector<UpperInterval>const&>(static_cast<Vector<ExactInterval>const&>(*this)); }
};

//! \relates ExactBox \brief The cartesian product of two boxes.
ExactBox product(const ExactBox& bx1, const ExactBox& bx2);
ExactBox product(const ExactBox& bx1, const ExactInterval& ivl2);
ExactBox product(const ExactInterval& ivl1, const ExactBox& bx2);
//! \relates ExactBox \brief The smallest box containing the two boxes.
ExactBox hull(const ExactBox& bx1, const ExactBox& bx2);
ExactBox hull(const ExactBox& bx1, const ExactPoint& pt2);
ExactBox hull(const ExactPoint& pt1, const ExactPoint& pt2);
//! \relates ExactBox \brief The intersection of the two boxes.
ExactBox intersection(const ExactBox& bx1, const ExactBox& bx2);
//! \relates ExactBox \brief A box which is wider than the input, and has single-precision values.
ExactBox widen(const ExactBox& bx);
//! \relates ExactBox \brief A box which is narrower than the input, and has single-precision values.
ExactBox narrow(const ExactBox& bx);

ExactBox make_box(const StringType& str);


//! \ingroup BasicSetSubModule GeometryModule
//! \brief A box in Euclidean space.
template<> class Box<UpperInterval>
    : public Vector<UpperInterval>
{
  public:
    //! Construct a singleton point in zero dimensions.
    Box<UpperInterval>() : Vector<UpperInterval>() { }
    //! Construct an empty box in \a d dimensions.
    explicit Box<UpperInterval>(Nat d) : Vector<UpperInterval>(d,UpperInterval()) { }
    explicit Box<UpperInterval>(Nat d, UpperInterval ivl) : Vector<UpperInterval>(d,ivl) { }
    //! Construct from an initializer list of pairs of floating-point values
    //! giving lower and upper bounds.
    Box<UpperInterval>(InitializerList<UpperInterval> lst);

    Box<UpperInterval>(Vector<ExactInterval>const& vec) : Vector<UpperInterval>(vec) { }
    Box<UpperInterval>(Vector<UpperInterval>const& vec) : Vector<UpperInterval>(vec) { }
    template<class E> Box<UpperInterval>(const VectorExpression<E>& t) : Vector<UpperInterval>(t) { }

    Box<UpperInterval>(const Vector<ValidatedFloat>& vec) : Vector<UpperInterval>(vec) { }

    //! The unit box \f$[-1,1]^n\f$ in \a n dimensions.
    static UpperBox unit_box(Nat n) {
        return UpperBox(n,UpperInterval(-1.0,1.0));
    }

    //! The upper quadrant box \f$[0,\infty]^n\f$ in \a n dimensions.
    static UpperBox upper_quadrant(Nat n) {
        return UpperBox(n,UpperInterval(0,inf));
    }

    //! An explicit case to an interval vector. Useful to prevent ambiguous function overloads.
    const Vector<UpperInterval>& vector() const { return *this; }

    //! An approximation to the centre of the box.
    ApproximatePoint centre() const {
        ApproximatePoint r(this->dimension());
        for(Nat i=0; i!=this->dimension(); ++i) { r[i]=midpoint((*this)[i]); }
        return r;
    }

    //! An over-approximation to radius of the box in the supremum norm.
    ErrorType radius() const {
        ErrorType dmax=0u;
        for(Nat i=0; i!=this->size(); ++i) { dmax = max( dmax, (*this)[i].width() ); }
        return half(dmax);
    }

    //! An approximation to the Lesbegue measure (area, volume) of the box.
    ApproximateNumber measure() const {
        ApproximateNumber meas=1u;
        for(Nat i=0; i!=this->size(); ++i) { meas *= (*this)[i].width(); }
        return meas;
    }

    //! \brief Test if the box is bounded.
    Bool bounded() const {
        for(Nat i=0; i!=dimension(); ++i) {
            if(!Ariadne::bounded((*this)[i])) { return false; } }
        return true;
    }

    //! \brief Test if the box is a subset of another box.
    Bool refines(const UpperBox& bx) const {
        for(Nat i=0; i!=dimension(); ++i) {
            if(!Ariadne::refines((*this)[i],bx[i])) { return false; } }
        return true;
    }

    friend Bool refines(const UpperBox& bx1, const UpperBox& bx2) {
        return bx1.refines(bx2);
    }

    friend UpperBox intersection(const UpperBox& bx1, const UpperBox& bx2) {
        assert(bx1.dimension()==bx2.dimension());
        UpperBox rbx(bx1.dimension());
        for(Nat i=0; i!=rbx.dimension(); ++i) { rbx[i]=intersection(bx1[i],bx2[i]); }
        return rbx;
    }

    friend UpperBox hull(const UpperBox& bx1, const UpperBox& bx2) {
        assert(bx1.dimension()==bx2.dimension());
        UpperBox rbx(bx1.dimension());
        for(Nat i=0; i!=rbx.dimension(); ++i) { rbx[i]=hull(bx1[i],bx2[i]); }
        return rbx;
    }

    //! \brief Test if the box is a subset of another box.
    Tribool inside(const ExactBox& bx) const {
        return Ariadne::inside(this->vector(),bx.vector());
    }

    //! \brief Test if the box is a subset of another box.
    Tribool subset(const ExactBox& bx) const {
        return Ariadne::subset(this->vector(),bx.vector());
    }

    //! \brief Test if the box intersects another box. Returns \a false even
    //! if the intersection only occurs only on the boundary of the boxes.
    Tribool disjoint(const UpperBox& bx) const {
        return Ariadne::disjoint(this->vector(),bx.vector());
    }

    //! \brief The dimension of the space the box lies in.
    Nat dimension() const {
        return this->size();
    }

    //! \brief Tests if the box is disjoint from another box.
    //! Only returns true if the closures are disjoint.
    Tribool empty() const{
        for(Nat i=0; i!=this->dimension(); ++i) {
            if(definitely((*this)[i].empty())) { return true; } }
        return indeterminate;
    }

    //! \brief Tests if the box is disjoint from another box.
    //! Only returns true if the closures are disjoint.
    Tribool separated(const UpperBox& other) const {
        return this->disjoint(other);
    }

    //! \brief Returns an enclosing bounding box for the set.
    //! The result is guaranteed to have nonempty interior, and floating-point
    //! boundary coefficients so the centre and radius are exactly computable.
    UpperBox bounding_box() const {
        UpperBox result=*this; result.widen(); return result;
    }

    //! \brief Widens the box by the minimal floating-point increment.
    //! The result is guaranteed to contain the box in its interior.
    Void widen() {
        UpperBox& result=*this;
        for(Nat i=0; i!=result.dimension(); ++i) {
            result[i]=Ariadne::widen((*this)[i]);
        }
    }

    //! \brief Split into two along the largest side.
    Pair<UpperBox,UpperBox> split() const {
        return ExactBox(reinterpret_cast<Vector<ExactInterval>const&>(*this)).split(); }

    //! \brief Split into two along side with index \a i.
    Pair<UpperBox,UpperBox> split(Nat i) const;

    //! \brief Write to an output stream.
    OutputStream& write(OutputStream& os) const {
        return os << *static_cast<const Vector<UpperInterval>*>(this);
    }

    Void draw(CanvasInterface& c, const Projection2d& p) const;
};

inline UpperBox product(const Vector<UpperInterval>& bx1, const Vector<UpperInterval>& bx2) {
    return UpperBox(join(bx1,bx2));
}

inline UpperBox product(const Vector<UpperInterval>& bx1, const UpperInterval& ivl2) {
    return UpperBox(join(bx1,ivl2));
}


//! \brief An over-approximation to an interval set.
template<> class Box<ApproximateInterval>
    : public Vector<ApproximateInterval>
{
  public:
    explicit Box<ApproximateInterval>(Nat d) : ApproximateBox(d,ApproximateInterval()) { }
    explicit Box<ApproximateInterval>(Nat d, ApproximateInterval const& ivl) : Vector<ApproximateInterval>(d,ivl) { }
    explicit Box<ApproximateInterval>(Vector<ApproximateInterval> const& vec) : Vector<ApproximateInterval>(vec) { }
    Box<ApproximateInterval>(ExactBox const& bx) : Vector<ApproximateInterval>(bx) { }
    Box<ApproximateInterval>(UpperBox const& bx) : Vector<ApproximateInterval>(bx) { }

    Nat dimension() const { return this->Vector<ApproximateInterval>::size(); }
};

inline ExactBox make_exact_box(UpperBox const& ubx) {
    return ExactBox(reinterpret_cast<Vector<ExactInterval>const&>(static_cast<Vector<UpperInterval>const&>(ubx)));
}


inline Bool subset(const ExactBox& bx1, const ExactBox& bx2) {
    return bx1.subset(bx2);
}

inline UpperBox ExactBox::bounding_box() const {
    return Ariadne::widen(*this);
}


inline UpperBox widen(UpperBox bx) {
    bx.widen(); return bx;
}

inline ApproximatePoint midpoint(const UpperBox& bx) {
    return midpoint(static_cast<Vector<UpperInterval>const&>(bx));
}

inline PositiveUpperFloat radius(const UpperBox& bx) {
    return bx.radius();
}

inline Void UpperBox::draw(CanvasInterface& c, const Projection2d& p) const {
    make_exact_box(*this).draw(c,p);
}


inline Bool contains(const ApproximateBox& abx, ApproximatePoint apt) {
    assert(abx.dimension()==apt.dimension());
    for(Nat i=0; i!=abx.dimension(); ++i) {
        if(!contains(abx[i],apt[i])) { return false; }
    }
    return true;
}

inline ApproximateBox widen(const UpperBox& bx, ApproximateFloat eps) {
    ApproximateBox r(bx.dimension());
    for(Nat i=0; i!=bx.size(); ++i) {
        r[i]=ApproximateInterval(bx[i].lower()-eps,bx[i].upper()+eps);
    }
    return r;
}

} // namespace Ariadne

#endif /* ARIADNE_BOX_H */
