/***************************************************************************
 *            box.cc
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file box.cc
 *  \brief
 */



#include "numeric/module.h"
#include "function/function.h"
#include "function/taylor_model.h"

#include "box.h"
#include "box.tpl.h"

#include "function/formula.h"

#include "algebra/algebra.h"

namespace Ariadne {

UpperIntervalType apply(ScalarFunction<ValidatedTag>const& f, const Box<UpperIntervalType>& x) {
    return static_cast<UpperIntervalType>(f(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
Box<UpperIntervalType> apply(VectorFunction<ValidatedTag>const& f, const Box<UpperIntervalType>& x) {
    return static_cast<Box<UpperIntervalType>>(f(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }

template class Box<Interval<Real>>;
template class Box<Interval<Float64Value>>;
template class Box<Interval<Float64UpperBound>>;
template class Box<Interval<Float64Approximation>>;

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

template List<typename Float64ExactBox::VertexType> Box<Float64ExactInterval>::vertices() const;

/*
Float64LowerBox under_approximation(const RealBox& rbx) {
    Float64LowerBox bx(rbx.size());
    for(SizeType i=0; i!=bx.size(); ++i) {
        bx[i]=under_approximation(rbx[i]);
    }
    return bx;
}

Float64UpperBox over_approximation(const RealBox& rbx) {
    Float64UpperBox bx(rbx.size(),Float64UpperInterval(-inf,+inf));
    for(SizeType i=0; i!=bx.size(); ++i) {
        bx[i]=over_approximation(rbx[i]);
    }
    return bx;
}

Float64ApproximationBox approximation(const RealBox& rbx) {
    Float64ApproximationBox bx(rbx.size(),Float64ApproximationInterval(-inf,+inf));
    for(SizeType i=0; i!=bx.size(); ++i) {
        bx[i]=approximation(rbx[i]);
    }
    return bx;
}
*/


/*
Vector<FloatBounds> cast_singleton(const Vector<Float64UpperInterval>& bx) {
    Vector<FloatBounds> r(bx.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=cast_singleton(bx[i]); }
    return r;
}

template class Box<RealInterval>;
template class Box<Float64ExactInterval>;
template class Box<Float64UpperInterval>;
template class Box<Float64ApproximationInterval>;
*/

} //namespace Ariadne
