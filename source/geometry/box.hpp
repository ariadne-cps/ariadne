/***************************************************************************
 *            box.hpp
 *
 *  Copyright 2008-17  Alberto Casagrande, Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file box.hpp
 *  \brief Coordinate-aligned boxes in Euclidean space.
 */

#ifndef ARIADNE_BOX_HPP
#define ARIADNE_BOX_HPP

#include "../utility/container.hpp"

#include "../numeric/logical.hpp"
#include "../numeric/floatdp.hpp"
#include "../geometry/interval.hpp"
#include "../geometry/point.hpp"
#include "../geometry/set_interface.hpp"

#include "box.decl.hpp"

namespace Ariadne {

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
    typedef typename I::CentreType C;
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
    typedef M CentreValueType;
    //! The type used for the midpoint of a component interval.
    typedef C MidpointValueType;
    //! A type which can be used for the both the lower and upper bound of a component interval.
    typedef V VertexBoundType;
    //! A type which can be used for the radius or length of a side.
    typedef R RadiusType;
  public:
    //! A type which can be used for the vertices of the box.
    typedef Point<V> VertexType;
    //! A type which can be used for the midpoint of the box.
    typedef Point<M> MidpointType;
    //! A type which can be used for the exact centre of the box.
    typedef Point<C> CentreType;
  public:
    Box(); // DEPRECATED
    explicit Box(SizeType n); // DEPRECATED
    Box(SizeType n, IntervalType ivl);
    Box(const InitializerList<IntervalType>& lst);

    template<class II> explicit Box(const Array<II>& bx);

    // Non-templated constructor to allow conversion from VectorExpression
    Box(const Vector<I>& bx) : Vector<I>(bx) { }

    //! Convert from a box of a different type
    template<class II, EnableIf<IsConvertible<II,I>> = dummy>
        Box(const Vector<II>& bx) : Vector<I>(bx) { }

    //! Construct from a box of a different type
    template<class II, EnableIf<And<IsConstructible<I,II>,Not<IsConvertible<II,I>>>> = dummy>
        explicit Box(const Vector<II>& bx) : Vector<I>(bx) { }

    template<class II, class PR, EnableIf<IsConstructible<I,II,PR>> = dummy>
        Box(const Vector<II>& bx, PR pr) : Vector<I>(bx,pr) { }

    //! \brief Construct from an interval of a different type using a default precision.
    template<class II, EnableIf<IsConstructibleGivenDefaultPrecision<I,II>> =dummy, DisableIf<IsConstructible<I,II>> =dummy>
        explicit Box(Box<II> const& x) : Box(x,PrecisionType<II>()) { }

    //! The unit box \f$[-1,1]^n\f$ in \a n dimensions.
    static Box<IntervalType> unit_box(SizeType n);
    //! The upper quadrant box \f$[0,\infty]^n\f$ in \a n dimensions.
    static Box<IntervalType> upper_orthant(SizeType n);

    //! The dimension of the set.
    DimensionType dimension() const;
    //! Indexing the bounds in each dimension yields an exact interval.
    const IntervalType& operator[](SizeType i) const;
    IntervalType& operator[](SizeType i);

    const Box<IntervalType> operator[](Range is) const;

    //! The (mid)point of the box.
    MidpointType midpoint() const;
    //! The exact centre of the box.
    CentreType centre() const;
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
    //! \brief Test if the box is singleton.
    auto is_bounded() const -> decltype(declval<IntervalType>().is_bounded());

    //! Splits the box along coordinate \a k and takes the lower, middle, or upper part as given by \a lmu.
    Box<IntervalType> split(SizeType k, SplitPart lmu) const;
    //! Splits the box along the longest direction \a k and takes the lower, middle, or upper part as given by \a lmu.
    Box<IntervalType> split(SplitPart lmu) const;
    //! Splits the box along widest coordinate \a k and takes the lower and upper parts.
    Pair< Box<IntervalType>, Box<IntervalType> > split(SizeType k) const;
    //! Splits the box along widest coordinate \a k and takes the lower and upper parts.
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

template<class I> template<class II> inline Box<I>::Box(const Array<II>& ary) : Vector<I>(ary) { }
template<class I> inline Box<I>::Box(InitializerList<I> const& lst) : Vector<I>(lst) { }

template<class I> decltype(declval<Box<I>>().is_empty()) is_empty(const Vector<I>& bx) { return static_cast<Box<I>const&>(bx).is_empty(); }
template<class I> decltype(declval<Box<I>>().is_bounded()) is_bounded(const Vector<I>& bx) { return static_cast<Box<I>const&>(bx).is_bounded(); }

template<class I> decltype(declval<Box<I>>().radius()) radius(const Vector<I>& bx) { return static_cast<Box<I>const&>(bx).radius(); }
template<class I> decltype(declval<Box<I>>().widths()) widths(const Vector<I>& bx) { return static_cast<Box<I>const&>(bx).widths(); }
template<class I> decltype(declval<Box<I>>().measure()) measure(const Vector<I>& bx) { return static_cast<Box<I>const&>(bx).measure(); }

template<class I> inline Box<I> split(const Vector<I>& bx, SizeType k, SplitPart lmu) { return static_cast<Box<I>const&>(bx).split(k,lmu); }
template<class I> inline Box<I> split(const Vector<I>& bx, SplitPart lmu) { return static_cast<Box<I>const&>(bx).split(lmu); }
template<class I> inline Pair<Box<I>,Box<I>> split(const Vector<I>& bx) { return static_cast<Box<I>const&>(bx).split(); }
template<class I> inline Pair<Box<I>,Box<I>> split(const Vector<I>& bx, SizeType k) { return static_cast<Box<I>const&>(bx).split(k); }

//! \relates FloatDPExactBox \brief The cartesian product of two boxes.
template<class I> inline Box<I> product(const Box<I>& bx1, const Box<I>& bx2) { return Box<I>::_product(bx1,bx2); }
template<class I> inline Box<I> product(const Box<I>& bx1, const I& ivl2) { return Box<I>::_product(bx1,ivl2); }
template<class I> inline Box<I> product(const Box<I>& bx1, const Box<I>& bx2, const Box<I>& bx3) { return Box<I>::_product(bx1,bx2,bx3); }

template<class S1, class S2, class S3> inline decltype(auto) product(S1 const& s1, S2 const& s2, S3 const& s3) { return product(product(s1,s2),s3); }

template<class I> inline Box<I> remove(const Box<I>& bx, SizeType k) {
    Box<I> rbx(bx.dimension()-1); for(SizeType i=0; i!=k; ++i) { rbx[i]=bx[i]; } for(SizeType i=k; i!=rbx.dimension(); ++i) { rbx[i]=bx[i+1]; } return rbx; }

UpperIntervalType apply(ScalarFunction<ValidatedTag>const& f, const Box<UpperIntervalType>& x);
Box<UpperIntervalType> apply(VectorFunction<ValidatedTag>const& f, const Box<UpperIntervalType>& x);
UpperBoxType image(UpperBoxType bx, ValidatedVectorFunction const& f);

//! \relates Box \brief Project onto the variables \a rng.
template<class I> inline Box<I> project(const Box<I> & bx, Array<SizeType> const& rng) { return Box<I>::_project(bx,rng); }


template<class I> auto midpoint(const Box<I>& bx) -> typename Box<I>::MidpointType {
    return bx.midpoint();
}

template<class I> auto lower_bounds(const Box<I>& bx) -> typename Box<I>::VertexType {
    return bx.lower_bounds();
}

template<class I> auto upper_bounds(const Box<I>& bx) -> typename Box<I>::VertexType {
    return bx.upper_bounds();
}


template<class I> Bool same(const Box<I>& bx1, const Box<I>& bx2) {
    if(bx1.size()!=bx2.size()) { return false; }
    for(SizeType i=0; i!=bx1.size(); ++i) {
        if(not same(bx1[i],bx2[i])) { return false; }
    }
    return true;
}

//! \relates EBoxType \brief The smallest box containing the two boxes.
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

//! \relates Box \brief The smallest box containing the box and a point.
template<class I, class X> Box<I> hull(const Point<X>& pt1, const Box<I>& bx2) {
    ARIADNE_ASSERT(pt1.size()==bx2.size());
    Box<I> r(bx2.size());
    for(SizeType i=0; i!=bx2.size(); ++i) {
        r[i]=hull(pt1[i],bx2[i]);
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


template<class I> inline decltype(refines(declval<I>(),declval<I>())) refines(const Box<typename Box<I>::IntervalType>& bx1, const Box<I>& bx2) {
    ARIADNE_ASSERT(bx1.size()==bx2.size());
    for(SizeType i=0; i!=bx1.size(); ++i) {
        if(!refines(bx1[i],bx2[i])) { return false; }
    }
    return true;
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



// Declare related functions
Void draw(CanvasInterface& c, Projection2d const& p, ApproximateBoxType const& bx);
ExactBoxType make_box(const String& str);


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

template<class I> inline const Box<I> Box<I>::operator[](Range rng) const {
    return Box<I>(Vector<I>(project(*this,rng)));
}

template<class I> inline auto cast_singleton(Box<I> const& bx) -> Vector<decltype(cast_singleton(declval<I>()))> {
    Vector<decltype(cast_singleton(declval<I>()))> v(bx.dimension());
    for(SizeType i=0; i!=bx.dimension(); ++i) { v[i]=cast_singleton(bx[i]); }
    return v;
}

template<class I> inline Void Box<I>::draw(CanvasInterface& c, const Projection2d& p) const {
    Ariadne::draw(c,p,ApproximateBoxType(*this)); }

template<> inline Void Box<Interval<Real>>::draw(CanvasInterface& c, const Projection2d& p) const {
    Ariadne::draw(c,p,ApproximateBoxType(*this,dp)); }

inline FloatDPExactBox cast_exact_box(FloatDPApproximateBox const& abx) {
    return FloatDPExactBox(reinterpret_cast<FloatDPExactBox const&>(abx));
}

inline Box<FloatDPExactInterval> cast_exact_box(Vector<FloatDPBounds> const& bv) {
    return Box<FloatDPExactInterval>(reinterpret_cast<Vector<FloatDPExactInterval>const&>(bv));
}


inline FloatDPApproximateBox widen(const FloatDPApproximateBox& bx, FloatDPApproximation e) {
    FloatDPApproximation eps(e);
    FloatDPApproximateBox r(bx.dimension());
    for(Nat i=0; i!=bx.size(); ++i) {
        r[i]=ApproximateIntervalType(bx[i].lower()-eps,bx[i].upper()+eps);
    }
    return r;
}

inline FloatDPUpperBox widen(const FloatDPUpperBox& bx, FloatDPUpperBound eps) {
    FloatDPUpperBox r(bx.dimension());
    for(Nat i=0; i!=bx.size(); ++i) {
        r[i]=widen(bx[i],eps);
    }
    return r;
}
// TODO: Add widen for other generic values
inline FloatDPUpperBox widen(const FloatDPUpperBox& bx, ValidatedUpperNumber eps) {
    return widen(bx,eps.get(UpperTag(),dp));
}

inline FloatDPUpperBox widen(const FloatDPExactBox& bx, FloatDPValue eps) {
    return widen(reinterpret_cast<const FloatDPUpperBox&>(bx),FloatDPUpperBound(eps));
}

inline FloatDPUpperBox widen(const FloatDPUpperBox& bx) {
    FloatDPUpperBox r(bx.dimension());
    for(Nat i=0; i!=bx.size(); ++i) {
        r[i]=widen(bx[i]);
    }
    return r;
}

inline FloatDPExactBox widen_domain(const FloatDPUpperBox& bx) {
    FloatDPExactBox r(bx.dimension());
    for(Nat i=0; i!=bx.size(); ++i) {
        r[i]=widen_domain(bx[i]);
    }
    return r;
}

inline FloatDPLowerBox narrow(const FloatDPLowerBox& bx, FloatDPUpperBound eps) {
    FloatDPLowerBox r(bx.dimension());
    for(Nat i=0; i!=bx.size(); ++i) {
        r[i]=narrow(bx[i],eps);
    }
    return r;
}

inline FloatDPLowerBox narrow(const FloatDPLowerBox& bx) {
    FloatDPLowerBox r(bx.dimension());
    for(Nat i=0; i!=bx.size(); ++i) {
        r[i]=narrow(bx[i]);
    }
    return r;
}


template<class IVL> class BoxSet;

template<> class BoxSet<ExactIntervalType>
    : public virtual RegularLocatedSetInterface,
      public virtual DrawableInterface,
      public Box<ExactIntervalType>
{

  public:
    using Box<ExactIntervalType>::Box;
    BoxSet<ExactIntervalType>(Box<ExactIntervalType>const& bx) : Box<ExactIntervalType>(bx) { }

    virtual ExactBoxSet* clone() const { return new ExactBoxSet(*this); }
    virtual DimensionType dimension() const final { return this->ExactBoxType::size(); }
    virtual LowerKleenean separated(const ExactBoxType& other) const final { return this->ExactBoxType::separated(other); }
    virtual LowerKleenean overlaps(const ExactBoxType& other) const final { return this->ExactBoxType::overlaps(other); }
    virtual LowerKleenean covers(const ExactBoxType& other) const final { return this->ExactBoxType::covers(other); }
    virtual LowerKleenean inside(const ExactBoxType& other) const final { return this->ExactBoxType::inside(other); }
    virtual ValidatedLowerKleenean separated(const ExactBoxType& other, Effort eff) const final { return this->ExactBoxType::separated(other); }
    virtual ValidatedLowerKleenean overlaps(const ExactBoxType& other, Effort eff) const final { return this->ExactBoxType::overlaps(other); }
    virtual ValidatedLowerKleenean covers(const ExactBoxType& other, Effort eff) const final { return this->ExactBoxType::covers(other); }
    virtual ValidatedLowerKleenean inside(const ExactBoxType& other, Effort eff) const final { return this->ExactBoxType::inside(other); }
    virtual UpperBoxType bounding_box() const final { return this->ExactBoxType::bounding_box(); }
    virtual Void draw(CanvasInterface& c, const Projection2d& p) const final { return Ariadne::draw(c,p,*this); }
    virtual OutputStream& write(OutputStream& os) const final { return os << static_cast<const ExactBoxType&>(*this); }
};

template<> class BoxSet<ApproximateIntervalType>
    : public virtual DrawableInterface,
      public Box<ApproximateIntervalType>
{

  public:
    using Box<ApproximateIntervalType>::Box;
    BoxSet<ApproximateIntervalType>(Box<ApproximateIntervalType>const& bx) : Box<ApproximateIntervalType>(bx) { }

    virtual ApproximateBoxSet* clone() const { return new ApproximateBoxSet(*this); }
    virtual DimensionType dimension() const final { return this->ApproximateBoxType::size(); }
    virtual Void draw(CanvasInterface& c, const Projection2d& p) const final { return Ariadne::draw(c,p,*this); }
    virtual OutputStream& write(OutputStream& os) const final { return os << static_cast<const ApproximateBoxType&>(*this); }
};

template<> class BoxSet<RealInterval>
    : public virtual DrawableInterface,
      public Box<RealInterval>
{

  public:
    using Box<RealInterval>::Box;
    BoxSet<RealInterval>(Box<RealInterval>const& bx) : Box<RealInterval>(bx) { }

    virtual RealBoxSet* clone() const { return new RealBoxSet(*this); }
    virtual DimensionType dimension() const final { return this->Box<RealInterval>::size(); }
    virtual Void draw(CanvasInterface& c, const Projection2d& p) const final { return this->Box<RealInterval>::draw(c,p); }
    virtual OutputStream& write(OutputStream& os) const final { return os << static_cast<const Box<RealInterval>&>(*this); }
};


} // namespace Ariadne


#endif /* ARIADNE_BOX_HPP */
