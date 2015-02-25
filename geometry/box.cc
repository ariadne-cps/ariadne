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

UpperInterval apply(ScalarFunction<ValidatedTag>const& f, const Box<UpperInterval>& x) {
    return static_cast<UpperInterval>(f(reinterpret_cast<Vector<ValidatedNumber>const&>(x))); }
Box<UpperInterval> apply(VectorFunction<ValidatedTag>const& f, const Box<UpperInterval>& x) {
    return static_cast<Box<UpperInterval>>(f(reinterpret_cast<Vector<ValidatedNumber>const&>(x))); }

template class Box<Interval<Real>>;
template class Box<Interval<ExactFloat64>>;
template class Box<Interval<UpperFloat64>>;
template class Box<Interval<ApproximateFloat64>>;

Void draw(CanvasInterface& c, Projection2d const& p, ApproximateFloatBox const& bx) {
    Nat ix=p.x_coordinate(); Nat iy=p.y_coordinate();
    ApproximateFloatInterval x=bx[ix]; ApproximateFloatInterval y=bx[iy];
    c.move_to(numeric_cast<double>(x.lower()),numeric_cast<double>(y.lower()));
    c.line_to(numeric_cast<double>(x.upper()),numeric_cast<double>(y.lower()));
    c.line_to(numeric_cast<double>(x.upper()),numeric_cast<double>(y.upper()));
    c.line_to(numeric_cast<double>(x.lower()),numeric_cast<double>(y.upper()));
    c.line_to(numeric_cast<double>(x.lower()),numeric_cast<double>(y.lower()));
    c.fill();
}

ExactBox make_box(const String& str)
{
    // Representation as a literal
    //   "[a1,b1]x[a2,b2]x...x[an,bn]"

    StringStream ss(str);
    std::vector<ExactInterval> vec;
    ExactInterval ivl;
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

    ExactBox bx(vec.size());
    for(Nat i=0; i!=bx.dimension(); ++i) {
        bx[i]=vec[i];
    }
    return bx;
}
/*
void make_vertices_down(const ExactFloatBox& bx, SizeType i, SizeType n, ExactFloatPoint& pt, std::vector<ExactFloatPoint>& v);

void make_vertices_up(const ExactFloatBox& bx, SizeType i, SizeType n, ExactFloatPoint& pt, std::vector<ExactFloatPoint>& v) {
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

void make_vertices_down(const ExactFloatBox& bx, SizeType i, SizeType n, ExactFloatPoint& pt, std::vector<ExactFloatPoint>& v) {
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
LowerFloatBox under_approximation(const RealBox& rbx) {
    LowerFloatBox bx(rbx.size());
    for(SizeType i=0; i!=bx.size(); ++i) {
        bx[i]=under_approximation(rbx[i]);
    }
    return bx;
}

UpperFloatBox over_approximation(const RealBox& rbx) {
    UpperFloatBox bx(rbx.size(),UpperFloatInterval(-inf,+inf));
    for(SizeType i=0; i!=bx.size(); ++i) {
        bx[i]=over_approximation(rbx[i]);
    }
    return bx;
}

ApproximateFloatBox approximation(const RealBox& rbx) {
    ApproximateFloatBox bx(rbx.size(),ApproximateFloatInterval(-inf,+inf));
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
Vector<ValidatedFloat> make_singleton(const Vector<UpperFloatInterval>& bx) {
    Vector<ValidatedFloat> r(bx.size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=make_singleton(bx[i]); }
    return r;
}

template class Box<RealInterval>;
template class Box<ExactFloatInterval>;
template class Box<UpperFloatInterval>;
template class Box<ApproximateFloatInterval>;
*/

} //namespace Ariadne
