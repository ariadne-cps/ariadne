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

#include "box.h"
#include "box.tcc"



namespace Ariadne {

UpperIntervalType apply(ScalarFunction<ValidatedTag>const& f, const Box<UpperIntervalType>& x) {
    return static_cast<UpperIntervalType>(f(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
Box<UpperIntervalType> apply(VectorFunction<ValidatedTag>const& f, const Box<UpperIntervalType>& x) {
    return static_cast<Box<UpperIntervalType>>(f(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }

template class Box<Interval<Real>>;
template class Box<Interval<ExactFloat64>>;
template class Box<Interval<UpperFloat64>>;
template class Box<Interval<ApproximateFloat64>>;

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
/*
void make_vertices_down(const ExactFloat64Box& bx, SizeType i, SizeType n, ExactFloatPoint& pt, std::vector<ExactFloatPoint>& v);

void make_vertices_up(const ExactFloat64Box& bx, SizeType i, SizeType n, ExactFloatPoint& pt, std::vector<ExactFloatPoint>& v) {
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

void make_vertices_down(const ExactFloat64Box& bx, SizeType i, SizeType n, ExactFloatPoint& pt, std::vector<ExactFloatPoint>& v) {
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
*/


/*
LowerFloat64Box under_approximation(const RealBox& rbx) {
    LowerFloat64Box bx(rbx.size());
    for(SizeType i=0; i!=bx.size(); ++i) {
        bx[i]=under_approximation(rbx[i]);
    }
    return bx;
}

UpperFloat64Box over_approximation(const RealBox& rbx) {
    UpperFloat64Box bx(rbx.size(),UpperFloat64Interval(-inf,+inf));
    for(SizeType i=0; i!=bx.size(); ++i) {
        bx[i]=over_approximation(rbx[i]);
    }
    return bx;
}

ApproximateFloat64Box approximation(const RealBox& rbx) {
    ApproximateFloat64Box bx(rbx.size(),ApproximateFloat64Interval(-inf,+inf));
    for(SizeType i=0; i!=bx.size(); ++i) {
        bx[i]=approximation(rbx[i]);
    }
    return bx;
}
*/



/*
List<Point> Box<I>::vertices() const {
    std::vector<Point> v;
    SizeType n = this->dimension();
    if(n > 0) {
        Point pt(n);
        make_vertices_up(*this, 0, n-1, pt, v);
    }
    return v;
}
*/


/*
Vector<BoundedFloat> cast_singleton(const Vector<UpperFloat64Interval>& bx) {
    Vector<BoundedFloat> r(bx.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=cast_singleton(bx[i]); }
    return r;
}

template class Box<RealInterval>;
template class Box<ExactFloat64Interval>;
template class Box<UpperFloat64Interval>;
template class Box<ApproximateFloat64Interval>;
*/

} //namespace Ariadne
