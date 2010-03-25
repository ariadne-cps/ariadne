/***************************************************************************
 *            box.cc
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

#include <sstream>
#include <string>
#include <vector>

#include <cstdarg>

#include "box.h"
#include "stlio.h"
#include "point.h"

typedef unsigned int uint;

namespace Ariadne {

//
// Helper functions needed to extract the set of vertices from a box
//
void make_vertices_down(const Box& bx, uint i, uint n, Point& pt, std::vector<Point>& v);

void make_vertices_up(const Box& bx, uint i, uint n, Point& pt, std::vector<Point>& v) {
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

void make_vertices_down(const Box& bx, uint i, uint n, Point& pt, std::vector<Point>& v) {
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


Box::Box(uint d, const Float& x0l, const Float& x0u, ...)
    : Vector<Interval>(d)
{
    assert(d>=1);
    va_list args;
    va_start(args,x0u);
    (*this)[0]=Interval(x0l,x0u);
    for(uint i=1; i!=d; ++i) {
        Float xil=va_arg(args,Float);
        Float xiu=va_arg(args,Float);
        (*this)[i]=Interval(xil,xiu);
    }
    va_end(args);
}

Box::Box(const std::string& str)
{
    *this=make_box(str);
}

std::vector<Point> Box::vertices() const {
    std::vector<Point> v;
    uint n = this->dimension();
    if(n > 0) {
        Point pt(n);
        make_vertices_up(*this, 0, n-1, pt, v);
    }
    return v;
}

Box hull(const Box& bx1, const Box& bx2) {
    return Box(Ariadne::hull(static_cast<const Vector<Interval>&>(bx1),static_cast<const Vector<Interval>&>(bx2)));
}

Box intersection(const Box& bx1, const Box& bx2) {
    return Box(Ariadne::intersection(static_cast<const Vector<Interval>&>(bx1),static_cast<const Vector<Interval>&>(bx2)));
}

Box widen(const Box& bx) {
    Box result(bx.dimension());
    for(uint i=0; i!=result.dimension(); ++i) {
        if(bx[i].lower()==bx[i].upper()) {
            result[i]=trunc(Ariadne::widen(bx[i]));
        } else {
            result[i]=trunc(bx[i]);
        }
    }
    return result;
}

void Box::draw(CanvasInterface& c) const
{
    uint ix=c.x_coordinate(); uint iy=c.y_coordinate();
    Interval x=(*this)[ix]; Interval y=(*this)[iy];
    c.move_to(x.lower(),y.lower());
    c.line_to(x.upper(),y.lower());
    c.line_to(x.upper(),y.upper());
    c.line_to(x.lower(),y.upper());
    c.line_to(x.lower(),y.lower());
    c.fill();
}

Box make_box(const std::string& str)
{
    // Representation as a literal
    //   "[a1,b1]x[a2,b2]x...x[an,bn]"

    std::stringstream ss(str);
    std::vector<Interval> vec;
    Interval ivl;
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
    return Box(vec.size(),&vec[0]);
}

} //namespace Ariadne
