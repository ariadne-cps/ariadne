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

#include "numeric/numeric.h"
#include "config.h"

#include <sstream>
#include <string>
#include <vector>

#include <cstdarg>

#include "geometry/box.h"
#include "utility/stlio.h"
#include "numeric/logical.h"
#include "geometry/point.h"

typedef unsigned int uint;

namespace Ariadne {

Vector<ValidatedFloat>const& make_singleton(const Vector<ExactInterval>& ivlv) {
    return reinterpret_cast<Vector<ValidatedFloat>const&>(ivlv);
}

Vector<ValidatedFloat>const& make_singleton(const Vector<UpperInterval>& ivlv) {
    return reinterpret_cast<Vector<ValidatedFloat>const&>(ivlv);
}


bool element(const Vector<ExactFloat>& v1, const Vector<ExactInterval>& v2)
{
    return contains(v2,v1);
}

bool element(const Vector<ValidatedFloat>& v1, const Vector<ExactInterval>& v2)
{
    return contains(v2,v1);
}


bool contains(const Vector<ExactInterval>& v1, const Vector<ExactFloat>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!contains(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool contains(const Vector<ExactInterval>& v1, const Vector<ValidatedFloat>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!contains(v1[i],v2[i])) { return false; }
    }
    return true;
}



bool subset(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!subset(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool intersect(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!intersect(v1[i],v2[i])) { return false; }
    }
    return true;
}


bool disjoint(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(disjoint(v1[i],v2[i])) { return true; }
    }
    return false;
}

bool overlap(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!overlap(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool covers(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!covers(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool inside(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!inside(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool empty(const Vector<ExactInterval>& v)
{
    for(size_t i=0; i!=v.size(); ++i) {
        if(empty(v[i])) { return true; }
    }
    return false;
}


uint irmax(const Vector<ExactInterval>& v) {
    uint imw(0);
    Float mw=v[0].width().raw();
    for(uint i=1; i!=v.size(); ++i) {
        if(v[i].width().raw()>mw) { imw=i; mw=v[i].width().raw(); }
    }
    return imw;
}


Vector<ExactInterval> split(const Vector<ExactInterval>& v, uint k, Tribool lr) {
    ARIADNE_ASSERT(k<v.size());
    Vector<ExactInterval> r(v);
    Float c=v[k].centre().raw();
    if(definitely(lr==true)) {
        r[k].set_upper(c);
    } else if(definitely(lr==false)) {
        r[k].set_lower(c);
    } else {
        Float cl=(3*v[k].lower().raw()+v[k].upper().raw())/4;
        Float cu=(v[k].lower().raw()+3*v[k].upper().raw())/4;
        r[k].set_lower(cl);
        r[k].set_upper(cu);
    }
    return r;
}

Pair< Vector<ExactInterval>, Vector<ExactInterval> > split(const Vector<ExactInterval>& v, uint k) {
    ARIADNE_ASSERT(k<v.size());
    Pair< Vector<ExactInterval>, Vector<ExactInterval> > r(v,v);
    Float c=v[k].centre().raw();
    r.first[k].set_upper(c);
    r.second[k].set_lower(c);
    return r;
}

Vector<ExactInterval> split(const Vector<ExactInterval>& v, Tribool lr) {
    return split(v,irmax(v),lr);
}

Pair< Vector<ExactInterval>, Vector<ExactInterval> > split(const Vector<ExactInterval>& v) {
    return split(v,irmax(v));
}



Vector<ExactFloat> midpoint(const Vector<ExactInterval>& v)
{
    Vector<ExactFloat> r(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
        r[i]=v[i].centre();
    }
    return r;
}

Vector<ExactFloat> lower_bounds(const Vector<ExactInterval>& v)
{
    Vector<ExactFloat> r(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
        r[i]=v[i].lower();
    }
    return r;
}

Vector<ExactFloat> upper_bounds(const Vector<ExactInterval>& v)
{
    Vector<ExactFloat> r(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
        r[i]=v[i].upper();
    }
    return r;
}

Vector<ExactInterval> hull(const Vector<ExactFloat>& v1, const Vector<ExactFloat>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    Vector<ExactInterval> r(v1.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        r[i]=hull(v1[i],v2[i]);
    }
    return r;
}

Vector<ExactInterval> hull(const Vector<ExactInterval>& v1, const Vector<ExactFloat>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    Vector<ExactInterval> r(v1.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        r[i]=hull(v1[i],v2[i]);
    }
    return r;
}

Vector<ExactInterval> hull(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    Vector<ExactInterval> r(v1.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        r[i]=hull(v1[i],v2[i]);
    }
    return r;
}

Vector<ExactInterval> intersection(const Vector<ExactInterval>& v1, const Vector<ExactInterval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    Vector<ExactInterval> r(v1.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        r[i]=intersection(v1[i],v2[i]);
    }
    return r;
}

PositiveUpperNumber radius(const Vector<ExactInterval>& v)
{
    PositiveUpperNumber r=0;
    for(size_t i=0; i!=v.size(); ++i) {
        r=max(r,v[i].radius());
    }
    return r;
}

PositiveUpperNumber volume(const Vector<ExactInterval>& v)
{
    PositiveUpperNumber r=1u;
    for(size_t i=0; i!=v.size(); ++i) {
        r*=diam(v[i]);
    }
    return r;
}



Vector<UpperInterval> hull(const Vector<UpperInterval>& bx1, const Vector<UpperInterval>& bx2) {
    ARIADNE_ASSERT(bx1.size()==bx2.size());
    Vector<UpperInterval> r(bx1.size());
    for(size_t i=0; i!=bx1.size(); ++i) {
        r[i]=hull(bx1[i],bx2[i]);
    }
    return r;
}

Vector<UpperInterval> intersection(const Vector<UpperInterval>& bx1, const Vector<UpperInterval>& bx2) {
    ARIADNE_ASSERT(bx1.size()==bx2.size());
    Vector<UpperInterval> r(bx1.size());
    for(size_t i=0; i!=bx1.size(); ++i) {
        r[i]=intersection(bx1[i],bx2[i]);
    }
    return r;
}

Tribool disjoint(const Vector<UpperInterval>& bx1, const Vector<UpperInterval>& bx2) {
    ARIADNE_ASSERT(bx1.size()==bx2.size());
    for(uint i=0; i!=bx1.size(); ++i) {
        if(definitely(disjoint(bx1[i],bx2[i]))) {
            return true;
        }
    }
    return Tribool(indeterminate);
}

Tribool subset(const Vector<UpperInterval>& bx1, const Vector<ExactInterval>& bx2) {
    ARIADNE_ASSERT(bx1.size()==bx2.size());
    for(uint i=0; i!=bx1.size(); ++i) {
        if(definitely(subset(bx1[i],bx2[i]))) {
        } else {
            return Tribool(indeterminate);
        }
    }
    return true;
}

Tribool inside(const Vector<UpperInterval>& bx1, const Vector<ExactInterval>& bx2) {
    ARIADNE_ASSERT(bx1.size()==bx2.size());
    for(uint i=0; i!=bx1.size(); ++i) {
        if(definitely(inside(bx1[i],bx2[i]))) {
        } else {
            return Tribool(indeterminate);
        }
    }
    return true;
}

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


Box<ExactInterval>::Box(std::initializer_list<ExactInterval> lst)
    : Vector<ExactInterval>(lst)
{
}

Box<ExactInterval>::Box(const std::string& str)
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
        r[i]=hull(bx1[i],pt2[i]);
    }
    return r;
}

ExactBox hull(const ExactPoint& pt1, const ExactPoint& pt2) {
    ARIADNE_ASSERT(pt1.dimension()==pt2.dimension());
    ExactBox r(pt1.dimension());
    for(size_t i=0; i!=r.dimension(); ++i) {
        r[i]=hull(pt1[i],pt2[i]);
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
