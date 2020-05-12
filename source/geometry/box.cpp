/***************************************************************************
 *            geometry/box.cpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file geometry/box.cpp
 *  \brief
 */


#include "../numeric/module.hpp"
#include "../function/function.hpp"
#include "../function/taylor_model.hpp"

#include "box.hpp"
#include "box.tpl.hpp"

#include "../function/formula.hpp"

#include "../algebra/algebra.hpp"

namespace Ariadne {

typedef FloatDPBounds ValidatedNumericType;
typedef Interval<FloatDPUpperBound> UpperIntervalType;

UpperIntervalType apply(ValidatedScalarUnivariateFunction const& f, UpperIntervalType const& ivl) {
    return static_cast<UpperIntervalType>(f(reinterpret_cast<ValidatedNumericType const&>(ivl))); }
UpperBoxType apply(ValidatedVectorUnivariateFunction const& f, UpperIntervalType const& ivl) {
    return static_cast<UpperBoxType>(f(reinterpret_cast<ValidatedNumericType const&>(ivl))); }
UpperIntervalType apply(ValidatedScalarMultivariateFunction const& f, UpperBoxType const& bx) {
    return static_cast<UpperIntervalType>(f(reinterpret_cast<Vector<ValidatedNumericType>const&>(bx))); }
UpperBoxType apply(ValidatedVectorMultivariateFunction const& f, UpperBoxType const& bx) {
    return static_cast<UpperBoxType>(f(reinterpret_cast<Vector<ValidatedNumericType>const&>(bx))); }

UpperIntervalType image(UpperIntervalType const& ivl, ValidatedScalarUnivariateFunction const& f) {
    return static_cast<UpperIntervalType>(f(reinterpret_cast<ValidatedNumericType const&>(ivl))); }
UpperBoxType image(UpperIntervalType const& ivl, ValidatedVectorUnivariateFunction const& f) {
    return static_cast<UpperBoxType>(f(reinterpret_cast<ValidatedNumericType const&>(ivl))); }
UpperIntervalType image(UpperBoxType const& bx, ValidatedScalarMultivariateFunction const& f) {
    return static_cast<UpperIntervalType>(f(reinterpret_cast<Vector<ValidatedNumericType>const&>(bx))); }
UpperBoxType image(UpperBoxType const& bx, ValidatedVectorMultivariateFunction const& f) {
    return static_cast<UpperBoxType>(f(reinterpret_cast<Vector<ValidatedNumericType>const&>(bx))); }

template class Box<Interval<Real>>;
template class Box<Interval<FloatDPValue>>;
template class Box<Interval<FloatDPUpperBound>>;
template class Box<Interval<FloatDPApproximation>>;

Void draw(CanvasInterface& c, Projection2d const& p, ApproximateBoxType const& bx) {
    Nat ix=p.x_coordinate(); Nat iy=p.y_coordinate();
    ApproximateIntervalType x=bx[ix]; ApproximateIntervalType y=bx[iy];
    c.move_to(numeric_cast<double>(x.lower()),numeric_cast<double>(y.lower()));
    c.line_to(numeric_cast<double>(x.upper()),numeric_cast<double>(y.lower()));
    c.line_to(numeric_cast<double>(x.upper()),numeric_cast<double>(y.upper()));
    c.line_to(numeric_cast<double>(x.lower()),numeric_cast<double>(y.upper()));
    c.line_to(numeric_cast<double>(x.lower()),numeric_cast<double>(y.lower()));
    c.fill();
}

ExactBoxType make_box(const String& str)
{
    // Representation as a literal
    //   "[a1,b1]x[a2,b2]x...x[an,bn]"

    StringStream ss(str);
    std::vector<ExactIntervalType> vec;
    ExactIntervalType ivl;
    char c;

    c='x';
    while(c=='x') {
        ss >> ivl;
        vec.push_back(ivl);
        c=' ';
        while( ss && c==' ') {
            ss >> c;
        }
    }
    if(ss) {
        ss.putback(c);
    }

    ExactBoxType bx(vec.size());
    for(Nat i=0; i!=bx.dimension(); ++i) {
        bx[i]=vec[i];
    }
    return bx;
}


template<class BX> using VertexType = typename BX::VertexType;

template<class BX> void make_vertices_down(const BX& bx, SizeType i, SizeType n, VertexType<BX>& pt, std::vector<VertexType<BX>>& v);

template<class BX> void make_vertices_up(const BX& bx, SizeType i, SizeType n, VertexType<BX>& pt, std::vector<VertexType<BX>>& v) {
    ARIADNE_ASSERT(i <= n);
    if(i == n) {    // base case: we are at the last dimension of the box
        pt[i] = bx[i].lower();
        v.push_back(pt);
        pt[i] = bx[i].upper();
        v.push_back(pt);
    } else {        // recursive case: we are still scanning dimensions
        pt[i] = bx[i].lower();
        make_vertices_up(bx, i+1, n, pt, v);
        pt[i] = bx[i].upper();
        make_vertices_down(bx, i+1, n, pt, v);
    }
}

template<class BX> void make_vertices_down(const BX& bx, SizeType i, SizeType n, VertexType<BX>& pt, std::vector<VertexType<BX>>& v) {
    ARIADNE_ASSERT(i <= n);
    if(i == n) {    // base case: we are at the last dimension of the box
        pt[i] = bx[i].upper();
        v.push_back(pt);
        pt[i] = bx[i].lower();
        v.push_back(pt);
    } else {        // recursive case: we are still scanning dimensions
        pt[i] = bx[i].upper();
        make_vertices_up(bx, i+1, n, pt, v);
        pt[i] = bx[i].lower();
        make_vertices_down(bx, i+1, n, pt, v);
    }
}

template<class I> List<typename Box<I>::VertexType> Box<I>::vertices() const {
    std::vector<VertexType> v;
    SizeType n = this->dimension();
    if(n > 0) {
        VertexType pt(n);
        make_vertices_up(*this, 0, n-1, pt, v);
    }
    return v;
}

template List<typename FloatDPExactBox::VertexType> Box<FloatDPExactInterval>::vertices() const;

#warning Fix code below
Int abslog2floor(FloatDP x) { return -53; }
Int abslog2floor(FloatMP const& x) { return x.exponent(); }

#warning Should convert to ValidatedReal
/*
Interval<FloatMPBounds> compute(const Interval<Real>& rivl, Effort eff) {
    MultiplePrecision pr(eff.work());
    return Interval<FloatMPBounds>(rivl.lower().get(pr),rivl.upper().get(pr));
}

Box<Interval<FloatMPBounds>> compute(const Box<Interval<Real>>& rbx, Effort eff) {
    return Box<Interval<FloatMPBounds>>(Array<Interval<FloatMPBounds>>(rbx.dimension(),[&rbx,eff](SizeType i){return compute(rbx[i],eff);}));
}
*/

ValidatedReal operator-(ValidatedReal const&, ValidatedReal const&);
ValidatedReal max(ValidatedReal const&, ValidatedReal const&);
Positive<ValidatedReal> cast_positive(ValidatedReal const&);
ValidatedNegatedSierpinskian operator==(ValidatedReal const&, ValidatedReal const&);
ValidatedKleenean operator<(ValidatedReal const& r1, ValidatedReal const& r2) {
    return r1.get() < r2.get(); }
ValidatedKleenean operator>(ValidatedReal const& r1, ValidatedReal const& r2) {
    return r1.get() > r2.get(); }
ValidatedKleenean operator<(ValidatedReal const& r1, FloatDPValue const& x2) {
    return r1.get() < Dyadic(x2); }
ValidatedKleenean operator>(ValidatedReal const& r1, FloatDPValue const& x2) {
    return r1.get() > Dyadic(x2); }

template<> struct MidpointTrait<DyadicBounds> { typedef Dyadic Type; };

Interval<FloatDPLowerBound> compute_under_approximation(const Interval<Real>& rivl, Effort eff) {
//    Interval<ValidatedReal> vrivl(rivl.lower().compute(eff),rivl.upper().compute(eff));
    DoublePrecision pr; return Interval<FloatDPLowerBound>(rivl.lower().get(pr),rivl.upper().get(pr));
}
Box<Interval<FloatDPLowerBound>> compute_under_approximation(const Box<Interval<Real>>& rbx, Effort eff) {
    return Box<Interval<FloatDPLowerBound>>(Array<Interval<FloatDPLowerBound>>(rbx.dimension(),[&rbx,eff](SizeType i){return compute_under_approximation(rbx[i],eff);}));
}
/*
Interval<FloatDPUpperBound> compute_over_approximation(const Interval<Real>& rivl, Effort eff) {
    DoublePrecision pr; return Interval<FloatDPUpperBound>(rivl.lower().get(pr),rivl.upper().get(pr));
}
Box<Interval<FloatDPUpperBound>> compute_over_approximation(const Box<Interval<Real>>& rbx, Effort eff) {
        return Box<Interval<FloatDPUpperBound>>(Array<Interval<FloatDPUpperBound>>(rbx.dimension(),[&rbx,eff](SizeType i){return compute_over_approximation(rbx[i],eff);}));
}
*/
Interval<ValidatedReal> compute(const Interval<Real>& rivl, Effort eff) {
    return Interval<ValidatedReal>(rivl.lower().compute(eff),rivl.upper().compute(eff));
}
Box<Interval<ValidatedReal>> compute(const Box<Interval<Real>>& rbx, Effort eff) {
    return Box<Interval<ValidatedReal>>(Array<Interval<ValidatedReal>>(rbx.dimension(),[&rbx,eff](SizeType i){return compute(rbx[i],eff);}));
}

Interval<DyadicBounds> compute_get(const Interval<Real>& rivl, Effort eff) {
    return Interval<DyadicBounds>(rivl.lower().compute(eff).get(),rivl.upper().compute(eff).get());
}
Box<Interval<DyadicBounds>> compute_get(const Box<Interval<Real>>& rbx, Effort eff) {
    return Box<Interval<DyadicBounds>>(Array<Interval<DyadicBounds>>(rbx.dimension(),[&rbx,eff](SizeType i){return compute_get(rbx[i],eff);}));
}

ValidatedLowerKleenean BoxSet<Interval<Real>>::separated(const ExactBoxType& other) const {
    Effort eff(static_cast<Nat>(std::max(0,-abslog2floor(other.radius().upper().raw()))));
    return compute_get(*this,eff).separated(other);
}
ValidatedLowerKleenean BoxSet<Interval<Real>>::overlaps(const ExactBoxType& other) const {
    Effort eff(static_cast<Nat>(std::max(0,-abslog2floor(other.radius().upper().raw()))));
    return compute_get(*this,eff).overlaps(other);
}
ValidatedLowerKleenean BoxSet<Interval<Real>>::covers(const ExactBoxType& other) const {
    Effort eff(static_cast<Nat>(std::max(0,-abslog2floor(other.radius().upper().raw()))));
    return compute(*this,eff).covers(other);
}
ValidatedLowerKleenean BoxSet<Interval<Real>>::inside(const ExactBoxType& other) const {
    Effort eff(static_cast<Nat>(std::max(0,-abslog2floor(other.radius().upper().raw()))));
    return compute(*this,eff).inside(other);
}

/*
FloatDPLowerBox under_approximation(const RealBox& rbx) {
    FloatDPLowerBox bx(rbx.size());
    for(SizeType i=0; i!=bx.size(); ++i) {
        bx[i]=under_approximation(rbx[i]);
    }
    return bx;
}

FloatDPUpperBox over_approximation(const Box& rbx) {
    FloatDPUpperBox bx(rbx.size(),FloatDPUpperInterval(-inf,+inf));
    for(SizeType i=0; i!=bx.size(); ++i) {
        bx[i]=over_approximation(rbx[i]);
    }
    return bx;
}

FloatDPApproximationBox approximation(const RealBox& rbx) {
    FloatDPApproximationBox bx(rbx.size(),FloatDPApproximationInterval(-inf,+inf));
    for(SizeType i=0; i!=bx.size(); ++i) {
        bx[i]=approximation(rbx[i]);
    }
    return bx;
}
*/


/*
Vector<FloatBounds> cast_singleton(const Vector<FloatDPUpperInterval>& bx) {
    Vector<FloatBounds> r(bx.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=cast_singleton(bx[i]); }
    return r;
}

template class Box<RealInterval>;
template class Box<FloatDPExactInterval>;
template class Box<FloatDPUpperInterval>;
template class Box<FloatDPApproximationInterval>;
*/

} //namespace Ariadne
