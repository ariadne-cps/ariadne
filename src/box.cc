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

#include "numeric.h"
#include "config.h"

#include <sstream>
#include <string>
#include <vector>

#include <cstdarg>

#include "box.h"
#include "stlio.h"
#include "point.h"
#include "box.h"

typedef unsigned int uint;

namespace Ariadne {

bool contains(const Vector<ExactInterval>& v1, const Vector<ExactFloat>& v2);

//
// Helper functions needed to extract the set of vertices from a box
//
void make_vertices_down(const ExactBox& bx, uint i, uint n, ExactPoint& pt, std::vector<ExactPoint>& v);

void make_vertices_up(const ExactBox& bx, uint i, uint n, ExactPoint& pt, std::vector<ExactPoint>& v) {
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

void make_vertices_down(const ExactBox& bx, uint i, uint n, ExactPoint& pt, std::vector<ExactPoint>& v) {
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


ExactBox::ExactBox(std::initializer_list<ExactInterval> lst)
    : Vector<ExactInterval>(lst)
{
}

ExactBox::ExactBox(const std::string& str)
{
    *this=make_box(str);
}

std::vector<ExactPoint> ExactBox::vertices() const {
    std::vector<ExactPoint> v;
    uint n = this->dimension();
    if(n > 0) {
        ExactPoint pt(n);
        make_vertices_up(*this, 0, n-1, pt, v);
    }
    return v;
}

ExactBox product(const ExactBox& bx1, const ExactBox& bx2) {
    size_t n1=bx1.dimension();
    size_t n2=bx2.dimension();
    ExactBox r(n1+n2);
    for(size_t i=0; i!=n1; ++i) {
        r[i]=bx1[i];
    }
    for(size_t i=0; i!=n2; ++i) {
        r[n1+i]=bx2[i];
    }
    return r;
}

ExactBox hull(const ExactBox& bx1, const ExactBox& bx2) {
    ARIADNE_ASSERT(bx1.dimension()==bx2.dimension());
    ExactBox r(bx1.dimension());
    for(size_t i=0; i!=r.dimension(); ++i) {
        r[i]=hull(bx1[i],bx2[i]);
    }
    return r;
}

ExactBox hull(const ExactBox& bx1, const ExactPoint& pt2) {
    ARIADNE_ASSERT(bx1.dimension()==pt2.dimension());
    ExactBox r(bx1.dimension());
    for(size_t i=0; i!=r.dimension(); ++i) {
        r[i]=hull(bx1[i],pt2[i].raw());
    }
    return r;
}

ExactBox hull(const ExactPoint& pt1, const ExactPoint& pt2) {
    ARIADNE_ASSERT(pt1.dimension()==pt2.dimension());
    ExactBox r(pt1.dimension());
    for(size_t i=0; i!=r.dimension(); ++i) {
        r[i]=hull(pt1[i].raw(),pt2[i].raw());
    }
    return r;
}

ExactBox intersection(const ExactBox& bx1, const ExactBox& bx2) {
    ARIADNE_ASSERT(bx1.dimension()==bx2.dimension());
    ExactBox r(bx1.dimension());
    for(size_t i=0; i!=r.dimension(); ++i) {
        r[i]=intersection(bx1[i],bx2[i]);
    }
    return r;
}

ExactBox widen(const ExactBox& bx) {
    ExactBox res=bx; res.widen(); return res;
    ExactBox result(bx.dimension());
    for(uint i=0; i!=result.dimension(); ++i) {
        if(bx[i].lower()==bx[i].upper()) {
            result[i]=trunc(Ariadne::widen(bx[i]));
        } else {
            result[i]=trunc(bx[i]);
        }
    }
    return result;
}

ExactBox narrow(const ExactBox& bx) {
    ExactBox res=bx; res.narrow(); return res;
    ExactBox result(bx.dimension());
    for(uint i=0; i!=result.dimension(); ++i) {
        if(bx[i].lower()==bx[i].upper()) {
            result[i]=trunc(Ariadne::narrow(bx[i]));
        } else {
            result[i]=trunc(bx[i]);
        }
    }
    return result;
}

void ExactBox::draw(CanvasInterface& c, const Projection2d& p) const
{
    uint ix=p.x_coordinate(); uint iy=p.y_coordinate();
    ExactInterval x=(*this)[ix]; ExactInterval y=(*this)[iy];
    c.move_to(approx_cast<double>(x.lower()),approx_cast<double>(y.lower()));
    c.line_to(approx_cast<double>(x.upper()),approx_cast<double>(y.lower()));
    c.line_to(approx_cast<double>(x.upper()),approx_cast<double>(y.upper()));
    c.line_to(approx_cast<double>(x.lower()),approx_cast<double>(y.upper()));
    c.line_to(approx_cast<double>(x.lower()),approx_cast<double>(y.lower()));
    c.fill();
}

ExactBox make_box(const std::string& str)
{
    // Representation as a literal
    //   "[a1,b1]x[a2,b2]x...x[an,bn]"

    std::stringstream ss(str);
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
    for(uint i=0; i!=bx.dimension(); ++i) {
        bx[i]=vec[i];
    }
    return bx;
}

} //namespace Ariadne
