/***************************************************************************
 *            box.h
 *
 *  Copyright 2008-13  Alberto Casagrande, Pieter Collins
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

#include "numeric/logical.h"
#include "numeric/float64.h"
#include "geometry/interval.h"
#include "geometry/point.h"
#include "geometry/set_interface.h"

#include "box.decl.h"

namespace Ariadne {

extern const ExactFloat64 infty;

struct DoubleInput {
    double _d;
    DoubleInput(double d) : _d(d) { }
    operator double () const { return _d; }
};



//! \ingroup BasicSetSubModule GeometryModule
//! \brief A box in Euclidean space.
template<class I>
class Box
    : public Vector<I>
{
    typedef typename I::UpperBoundType U;
    typedef typename I::LowerBoundType L;
    typedef typename I::MidpointType M;
    typedef typename I::RadiusType R;
    typedef typename I::WidthType W;
    typedef decltype(max(declval<L>(),declval<U>())) V;
  public:
    //! The type used for a component interval.
    typedef I IntervalType;
    //! The type used for the lower bound of a component interval.
    typedef L LowerBoundType;
    //! The type used for the upper bound of a component interval.
    typedef U UpperBoundType;
    //! The type used for the midpoint of a component interval.
    typedef M MidpointValueType;
    //! A type which can be used for the both the lower and upper bound of a component interval.
    typedef V VertexBoundType;
    //! A type which can be used for the radius or length of a side.
    typedef R RadiusType;
  public:
    //! A type which can be used for the vertices of the box.
    typedef Point<V> VertexType;
    //! A type which can be used for the midpoint of the box.
    typedef Point<M> MidpointType;
  public:
    Box(); // DEPRECATED
    explicit Box(SizeType n);
    Box(SizeType n, IntervalType ivl);
    Box(const InitializerList<IntervalType>& lst);

    template<class II> explicit Box(const Array<II>& bx);

    // Non-templated constructor to allow conversion from VectorExpression
    Box(const Vector<I>& bx) : Vector<I>(bx) { }

    //! Convert from a box of a different type
    template<class II, EnableIf<IsConvertible<II,I>> = dummy>
        Box(const Vector<II>& bx) : Vector<I>(bx) { }

    //! Convert from a box of a different type
    template<class II, EnableIf<And<IsConstructible<I,II>,Not<IsConvertible<II,I>>>> = dummy>
        explicit Box(const Vector<II>& bx) : Vector<I>(bx) { }

    //! The unit box \f$[-1,1]^n\f$ in \a n dimensions.
    static Box<IntervalType> unit_box(SizeType n);
    //! The upper quadrant box \f$[0,\infty]^n\f$ in \a n dimensions.
    static Box<IntervalType> upper_quadrant(SizeType n);

    //! The dimension of the set.
    SizeType dimension() const;
    //! Indexing the bounds in each dimension yields an exact interval.
    const IntervalType& operator[](SizeType i) const;
    IntervalType& operator[](SizeType i);

    //! The midpoint of the box.
    MidpointType midpoint() const;
    MidpointType centre() const;
    //! The lower corner of the box.
    VertexType lower_bounds() const;
    //! The upper corner of the box.
    VertexType  upper_bounds() const;
    //! The vertices of the box.
    List<VertexType> vertices() const;
    //! The bounding box.
    Box<IntervalType> bounding_box() const;

    //! \brief Test if the box is empty.
    auto is_empty() const -> decltype(declval<IntervalType>().is_empty());
    auto empty() const -> decltype(declval<IntervalType>().empty());
    //! \brief Test if the box is bounded.
    auto is_bounded() const -> decltype(declval<IntervalType>().is_bounded());
    auto bounded() const -> decltype(declval<IntervalType>().bounded());

    //! Splits the box along coordinate \a k and takes the lower, middle, or upper part as given by \a lmu.
    Box<IntervalType> split(SizeType k, SplitPart lmu) const;
    //! Determines an coordinate \a k and splits the box along coordinate \a k and takes the lower, middle, or upper part as given by \a lmu.
    Box<IntervalType> split(SplitPart lmu) const;
    //! Splits the box along coordinate \a k and takes the lower and upper parts.
    Pair< Box<IntervalType>, Box<IntervalType> > split(SizeType k) const;
    //! Determines an coordinate \a k, and splits the box along coordinate \a k and takes the lower and upper parts.
    Pair< Box<IntervalType>, Box<IntervalType> > split() const;

    //! \brief Test if the box contains the point \a pt.
    template<class X> auto contains(const Point<X>& pt) const -> decltype(Ariadne::contains(declval<I>(),declval<X>()));

    template<class II> auto intersects(const Box<II>& bx) const -> decltype(Ariadne::intersect(declval<I>(),declval<II>()));
    template<class II> auto disjoint(const Box<II>& bx) const -> decltype(Ariadne::disjoint(declval<I>(),declval<II>()));
    template<class II> auto subset(const Box<II>& bx) const -> decltype(Ariadne::subset(declval<I>(),declval<II>()));
    template<class II> auto superset(const Box<II>& bx) const -> decltype(Ariadne::superset(declval<I>(),declval<II>()));
    template<class II> auto covers(const Box<II>& bx) const -> decltype(Ariadne::covers(declval<I>(),declval<II>()));
    template<class II> auto overlaps(const Box<II>& bx) const -> decltype(Ariadne::overlap(declval<I>(),declval<II>()));
    template<class II> auto inside(const Box<II>& bx) const -> decltype(Ariadne::inside(declval<I>(),declval<II>()));
    template<class II> auto separated(const Box<II>& bx) const -> decltype(Ariadne::separated(declval<I>(),declval<II>()));

    //! The area or volume of the box.
    RadiusType measure() const;
    //! The sum of the lengths of the sides.
    RadiusType lengths() const;
    //! Half the length of the longest side.
    RadiusType radius() const;

    Void draw(CanvasInterface& c, const Projection2d& p) const;
  public:
    static Box<I> _product(const Box<I>& bx1, const Box<I>& bx2, const Box<I>& bx3);
    static Box<I> _product(const Box<I>& bx1, const Box<I>& bx2);
    static Box<I> _product(const Box<I>& bx1, const I& ivl2);

    static Box<I> _project(const Box<I>& bx, const Array<SizeType>& rng);
};

template<class I> inline OutputStream& operator<<(OutputStream& os, const Box<I>& bx) {
    return os << static_cast<const Vector<I>&>(bx);
}

template<class I> template<class II> inline Box<I>::Box(const Array<II>& bx) : Vector<I>(bx) { }
template<class I> inline Box<I>::Box(InitializerList<I> const& lst) : Vector<I>(lst) { }

template<class I> decltype(declval<Box<I>>().empty()) empty(const Array<I>& bx) { return static_cast<Box<I>const&>(bx).empty(); }
template<class I> decltype(declval<Box<I>>().bounded()) bounded(const Array<I>& bx) { return static_cast<Box<I>const&>(bx).bounded(); }

template<class I> decltype(declval<Box<I>>().radius()) radius(const Array<I>& bx) { return static_cast<Box<I>const&>(bx).radius(); }
template<class I> decltype(declval<Box<I>>().widths()) widths(const Array<I>& bx) { return static_cast<Box<I>const&>(bx).widths(); }
template<class I> decltype(declval<Box<I>>().measure()) measure(const Array<I>& bx) { return static_cast<Box<I>const&>(bx).measure(); }

template<class I> inline Box<I> split(const Vector<I>& bx, SizeType k, SplitPart lmu) { return static_cast<Box<I>const&>(bx).split(k,lmu); }
template<class I> inline Box<I> split(const Vector<I>& bx, SplitPart lmu) { return static_cast<Box<I>const&>(bx).split(lmu); }

template<class I> inline Pair<Box<I>,Box<I>> split(const Vector<I>& bx) { return static_cast<Box<I>const&>(bx).split(); }
template<class I> inline Pair<Box<I>,Box<I>> split(const Vector<I>& bx, SizeType k) { return static_cast<Box<I>const&>(bx).split(k); }

//! \relates ExactFloatBox \brief The cartesian product of two boxes.
template<class I> inline Box<I> product(const Box<I>& bx1, const Box<I>& bx2) { return Box<I>::_product(bx1,bx2); }
template<class I> inline Box<I> product(const Box<I>& bx1, const I& ivl2) { return Box<I>::_product(bx1,ivl2); }
template<class I> inline Box<I> product(const Box<I>& bx1, const Box<I>& bx2, const Box<I>& bx3) { return Box<I>::_product(bx1,bx2,bx3); }

template<class I> inline decltype(refines(declval<I>(),declval<I>())) refines(const Box<typename Box<I>::IntervalType>& bx1, const Box<I>& bx2) {
    ARIADNE_ASSERT(bx1.size()==bx2.size());
    for(SizeType i=0; i!=bx1.size(); ++i) {
        if(!refines(bx1[i],bx2[i])) { return false; }
    }
    return true;
}

UpperInterval apply(ScalarFunction<ValidatedTag>const& f, const Box<UpperInterval>& x);
Box<UpperInterval> apply(VectorFunction<ValidatedTag>const& f, const Box<UpperInterval>& x);

//! Project onto the variables \a rng.
//template<class I> inline Box<I> project(const Box<I> & bx, Range rng) { return Box<I>::_project(bx,rng); }

template<class I> auto is_empty(Box<I> const& bx) -> decltype(bx.is_empty()) { return bx.is_empty(); }

// Provide the following avoiding the Point class
// Use decltype to allow action on arbitrary classes suppoerting midpoint etc.
template<class I> auto midpoint(const Box<I>& bx) -> typename Box<I>::MidpointType {
    return bx.midpoint();
}

template<class I> auto lower_bounds(const Box<I>& bx) -> typename Box<I>::VertexType {
    return bx.lower_bounds();
}

template<class I> auto upper_bounds(const Box<I>& bx) -> typename Box<I>::VertexType {
    return bx.upper_bounds();
}


//! \relates EBox \brief The smallest box containing the two boxes.
template<class I1, class I2> Box<decltype(hull(declval<I1>(),declval<I2>()))> hull(const Box<I1>& bx1, const Box<I2>& bx2) {
    ARIADNE_ASSERT(bx1.size()==bx2.size());
    Box<decltype(hull(declval<I1>(),declval<I2>()))> r(bx1.size());
    for(SizeType i=0; i!=bx1.size(); ++i) {
        r[i]=hull(bx1[i],bx2[i]);
    }
    return r;
}

//! \relates Box \brief The smallest box containing the box and a point.
template<class I, class X> Box<I> hull(const Box<I>& bx1, const Point<X>& pt2) {
    ARIADNE_ASSERT(bx1.size()==pt2.size());
    Box<I> r(bx1.size());
    for(SizeType i=0; i!=bx1.size(); ++i) {
        r[i]=hull(bx1[i],pt2[i]);
    }
    return r;
}

//! \relates Box \brief The intersection of the two boxes.
template<class I1, class I2> Box<decltype(intersection(declval<I1>(),declval<I2>()))> intersection(const Box<I1>& bx1, const Box<I2>& bx2) {
    Box<decltype(intersection(declval<I1>(),declval<I2>()))> res(bx1.size());
    for(SizeType i=0; i!=bx1.size(); ++i) { res[i]=intersection(bx1[i],bx2[i]); }
    return res;
}

template<class I, class X> decltype(element(declval<X>(),declval<I>())) element(const Vector<X>& pt1, const Box<I>& bx2) {
    decltype(element(declval<X>(),declval<I>())) res=true;
    for(SizeType i=0; i!=pt1.size(); ++i) { res = res && element(pt1[i],bx2[i]); }
    return res;
}

template<class I, class X> decltype(contains(declval<I>(),declval<X>())) contains(const Box<I>& bx1, const Vector<X>& vec2) {
    decltype(contains(declval<I>(),declval<X>())) res=true;
    for(SizeType i=0; i!=bx1.size(); ++i) { res = res && contains(bx1[i],vec2[i]); }
    return res;
}

template<class I, class X> decltype(contains(declval<I>(),declval<X>())) contains(const Box<I>& bx1, const Point<X>& pt2) {
    decltype(contains(declval<I>(),declval<X>())) res=true;
    for(SizeType i=0; i!=bx1.size(); ++i) { res = res && contains(bx1[i],pt2[i]); }
    return res;
}

template<class I1, class I2> decltype(disjoint(declval<I1>(),declval<I2>())) disjoint(const Box<I1>& bx1, const Box<I2>& bx2) {
    decltype(disjoint(declval<I1>(),declval<I2>())) res=false;
    for(SizeType i=0; i!=bx1.size(); ++i) { res = res || disjoint(bx1[i],bx2[i]); }
    return res;
}

template<class I1, class I2> decltype(intersect(declval<I1>(),declval<I2>())) intersect(const Box<I1>& bx1, const Box<I2>& bx2) {
    decltype(intersect(declval<I1>(),declval<I2>())) res=true;
    for(SizeType i=0; i!=bx1.size(); ++i) { res = res && intersect(bx1[i],bx2[i]); }
    return res;
}

template<class I1, class I2> decltype(equal(declval<I1>(),declval<I2>())) equal(const Box<I1>& bx1, const Box<I2>& bx2) {
    decltype(equal(declval<I1>(),declval<I2>())) res=true;
    for(SizeType i=0; i!=bx1.size(); ++i) { res = res && equal(bx1[i],bx2[i]); }
    return res;
}

template<class I1, class I2> decltype(subset(declval<I1>(),declval<I2>())) subset(const Box<I1>& bx1, const Box<I2>& bx2) {
    decltype(subset(declval<I1>(),declval<I2>())) res=true;
    for(SizeType i=0; i!=bx1.size(); ++i) { res = res && subset(bx1[i],bx2[i]); }
    return res;
}

template<class I1, class I2> decltype(superset(declval<I1>(),declval<I2>())) superset(const Box<I1>& bx1, const Box<I2>& bx2) {
    decltype(superset(declval<I1>(),declval<I2>())) res=true;
    for(SizeType i=0; i!=bx1.size(); ++i) { res = res && superset(bx1[i],bx2[i]); }
    return res;
}

template<class I1, class I2> decltype(covers(declval<I1>(),declval<I2>())) covers(const Box<I1>& bx1, const Box<I2>& bx2) {
    decltype(covers(declval<I1>(),declval<I2>())) res=true;
    for(SizeType i=0; i!=bx1.size(); ++i) { res = res && covers(bx1[i],bx2[i]); }
    return res;
}

template<class I1, class I2> decltype(overlap(declval<I1>(),declval<I2>())) overlap(const Box<I1>& bx1, const Box<I2>& bx2) {
    decltype(overlap(declval<I1>(),declval<I2>())) res=true;
    for(SizeType i=0; i!=bx1.size(); ++i) { res = res && overlap(bx1[i],bx2[i]); }
    return res;
}

template<class I1, class I2> decltype(separated(declval<I1>(),declval<I2>())) separated(const Box<I1>& bx1, const Box<I2>& bx2) {
    decltype(separated(declval<I1>(),declval<I2>())) res=false;
    for(SizeType i=0; i!=bx1.size(); ++i) { res = res || separated(bx1[i],bx2[i]); }
    return res;
}

template<class I1, class I2> decltype(inside(declval<I1>(),declval<I2>())) inside(const Box<I1>& bx1, const Box<I2>& bx2) {
    decltype(inside(declval<I1>(),declval<I2>())) res=true;
    for(SizeType i=0; i!=bx1.size(); ++i) { res = res && inside(bx1[i],bx2[i]); }
    return res;
}

template<class I> template<class X> inline auto Box<I>::contains(const Point<X>& pt) const -> decltype(Ariadne::contains(declval<I>(),declval<X>())) {
    return Ariadne::contains(*this,pt); }

template<class I1> template<class I2> inline auto Box<I1>::intersects(const Box<I2>& bx) const -> decltype(Ariadne::intersect(declval<I1>(),declval<I2>())) {
    return Ariadne::intersect(*this,bx); }
template<class I1> template<class I2> inline auto Box<I1>::disjoint(const Box<I2>& bx) const -> decltype(Ariadne::disjoint(declval<I1>(),declval<I2>())) {
    return Ariadne::disjoint(*this,bx); }
template<class I1> template<class I2> inline auto Box<I1>::subset(const Box<I2>& bx) const -> decltype(Ariadne::subset(declval<I1>(),declval<I2>())) {
    return Ariadne::subset(*this,bx); }
template<class I1> template<class I2> inline auto Box<I1>::superset(const Box<I2>& bx) const -> decltype(Ariadne::superset(declval<I1>(),declval<I2>())) {
    return Ariadne::superset(*this,bx); }
template<class I1> template<class I2> inline auto Box<I1>::covers(const Box<I2>& bx) const -> decltype(Ariadne::covers(declval<I1>(),declval<I2>())) {
    return Ariadne::covers(*this,bx); }
template<class I1> template<class I2> inline auto Box<I1>::overlaps(const Box<I2>& bx) const -> decltype(Ariadne::overlap(declval<I1>(),declval<I2>())) {
    return Ariadne::overlap(*this,bx); }
template<class I1> template<class I2> inline auto Box<I1>::inside(const Box<I2>& bx) const -> decltype(Ariadne::inside(declval<I1>(),declval<I2>())) {
    return Ariadne::inside(*this,bx); }
template<class I1> template<class I2> inline auto Box<I1>::separated(const Box<I2>& bx) const -> decltype(Ariadne::separated(declval<I1>(),declval<I2>())) {
    return Ariadne::separated(*this,bx); }



// Inline functions so as not to have to instantiate box class

template<class I> Box<I>::Box()
    : Vector<I>() { }

template<class I> inline Box<I>::Box(SizeType n)
    : Box(n,IntervalType::empty_interval()) { }

template<class I> inline Box<I>::Box(SizeType n, IntervalType ivl)
    : Vector<I>(n,ivl) { }

template<class I> inline Box<I> Box<I>::unit_box(SizeType n) {
    return Box<IntervalType>(n,IntervalType(-1,1));
}

template<class I> inline SizeType Box<I>::dimension() const {
    return this->Vector<I>::size();
}
template<class I> inline const I& Box<I>::operator[](SizeType i) const {
    return this->Vector<I>::operator[](i);
}

template<class I> inline I& Box<I>::operator[](SizeType i) {
    return this->Vector<I>::operator[](i);
}

template<class I> inline auto make_singleton(Box<I> const& bx) -> Vector<decltype(make_singleton(declval<I>()))> {
    Vector<decltype(make_singleton(declval<I>()))> v(bx.dimension());
    for(SizeType i=0; i!=bx.dimension(); ++i) { v[i]=make_singleton(bx[i]); }
    return v;
}

inline ExactBox make_exact_box(ApproximateBox const& abx) {
    return ExactBox(reinterpret_cast<ExactBox const&>(abx));
}

inline ExactBox make_exact_box(Vector<BoundedFloat64> const& bv) {
    return ExactBox(reinterpret_cast<Vector<ExactInterval>const&>(bv));
}

inline UpperBox widen(const UpperBox& bx, UpperFloat64 eps) {
    UpperBox r(bx.dimension());
    for(Nat i=0; i!=bx.size(); ++i) {
        r[i]=widen(bx[i],eps);
    }
    return r;
}

inline UpperBox widen(const UpperBox& bx) {
    return widen(bx,UpperFloat64(Float64::min()));
}

inline ApproximateBox widen(const ApproximateBox& bx, Float64 e) {
    ApproximateFloat64 eps(e);
    ApproximateBox r(bx.dimension());
    for(Nat i=0; i!=bx.size(); ++i) {
        r[i]=ApproximateInterval(bx[i].lower()-eps,bx[i].upper()+eps);
    }
    return r;
}

inline LowerBox narrow(const LowerBox& bx, UpperFloat64 eps) {
    LowerBox r(bx.dimension());
    for(Nat i=0; i!=bx.size(); ++i) {
        r[i]=narrow(bx[i],eps);
    }
    return r;
}

inline LowerBox narrow(const LowerBox& bx) {
    return narrow(bx,UpperFloat64(Float64::min()));
}



Void draw(CanvasInterface& c, Projection2d const& p, ApproximateFloatBox const& bx);

template<class I> inline Void Box<I>::draw(CanvasInterface& c, const Projection2d& p) const {
    return Ariadne::draw(c,p,ApproximateFloatBox(*this)); }

ExactBox make_box(const String& str);

class ExactBoxSet
    : public virtual SetInterface,
      public virtual DrawableInterface,
      public Box<ExactInterval>
{

  public:
    using Box<ExactInterval>::Box;
    ExactBoxSet(Box<ExactInterval>const& bx) : Box<ExactInterval>(bx) { }

    virtual ExactBoxSet* clone() const { return new ExactBoxSet(*this); }
    virtual SizeType dimension() const final { return this->ExactBox::size(); }
    virtual Tribool separated(const ExactBox& other) const final { return this->ExactBox::separated(other); }
    virtual Tribool overlaps(const ExactBox& other) const final { return this->ExactBox::overlaps(other); }
    virtual Tribool covers(const ExactBox& other) const final { return this->ExactBox::covers(other); }
    virtual Tribool inside(const ExactBox& other) const final { return this->ExactBox::inside(other); }
    virtual UpperBox bounding_box() const final { return this->ExactBox::bounding_box(); }
    virtual Void draw(CanvasInterface& c, const Projection2d& p) const final { return Ariadne::draw(c,p,*this); }
    virtual OutputStream& write(OutputStream& os) const final { return os << static_cast<const ExactBox&>(*this); }
};

class ApproximateBoxSet
    : public virtual DrawableInterface,
      public Box<ApproximateInterval>
{

  public:
    using Box<ApproximateInterval>::Box;
    ApproximateBoxSet(Box<ApproximateInterval>const& bx) : Box<ApproximateInterval>(bx) { }

    virtual ApproximateBoxSet* clone() const { return new ApproximateBoxSet(*this); }
    virtual SizeType dimension() const final { return this->ApproximateBox::size(); }
    virtual Void draw(CanvasInterface& c, const Projection2d& p) const final { return Ariadne::draw(c,p,*this); }
    virtual OutputStream& write(OutputStream& os) const final { return os << static_cast<const ApproximateBox&>(*this); }
};

//inline bool separated(const ExactFloatBox& bx1, const ExactFloatBox& bx2) { return bx1.separated(bx2); }
//inline bool overlap(const ExactFloatBox& bx1, const ExactFloatBox& bx2) { return bx1.overlaps(bx2); }
//inline bool inside(const ExactFloatBox& bx1, const ExactFloatBox& bx2) { return bx1.inside(bx2); }
//inline bool covers(const ExactFloatBox& bx1, const ExactFloatBox& bx2) { return bx1.covers(bx2); }
/*
Sierpinski overlap(const Vector<LowerFloatInterval>& bx1, const Vector<LowerFloatInterval>& bx2);
Sierpinski separated(const Vector<UpperFloatInterval>& bx1, const Vector<UpperFloatInterval>& bx2);
Sierpinski inside(const Vector<UpperFloatInterval>& bx1, const Vector<LowerFloatInterval>& bx2);
Sierpinski covers(const Vector<LowerFloatInterval>& bx1, const Vector<UpperFloatInterval>& bx2);

Box<UpperFloatInterval> over_approximation(const Box<RealInterval>& bx);
Box<LowerFloatInterval> under_approximation(const Box<RealInterval>& bx);
Box<ApproximateFloatInterval> approximation(const Box<RealInterval>& bx);

Box<UpperFloatInterval> widen(Box<UpperFloatInterval> bx, UpperFloat e);

Box<UpperFloatInterval> widen(Box<UpperFloatInterval> bx);
Box<LowerFloatInterval> narrow(Box<LowerFloatInterval> bx);

Vector<ValidatedFloat> make_singleton(const Vector<UpperFloatInterval>& bx);

inline ExactFloatBox const& make_exact(const Vector<ExactFloatInterval>& bx) { return reinterpret_cast<ExactFloatBox const&>(bx); }
inline ExactFloatBox const& make_exact(const Vector<UpperFloatInterval>& bx) { return reinterpret_cast<ExactFloatBox const&>(bx); }
inline ExactFloatBox const& make_exact(const Vector<LowerFloatInterval>& bx) { return reinterpret_cast<ExactFloatBox const&>(bx); }
inline ExactFloatBox const& make_exact(const Vector<ApproximateFloatInterval>& bx) { return reinterpret_cast<ExactFloatBox const&>(bx); }
*/

} // namespace Ariadne


#endif /* ARIADNE_BOX_H */
