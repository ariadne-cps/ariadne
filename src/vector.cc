/***************************************************************************
 *            vector.cc
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

#include "macros.h"
#include "vector.h"

namespace Ariadne {

Vector<ExactFloatType>const& make_exact(const Vector<ApproximateFloatType>& av) {
    return reinterpret_cast<Vector<ExactFloatType>const&>(av);
}

Vector<ValidatedFloatType> make_bounds(const Vector<ErrorFloatType>& ev) {
    Vector<ValidatedFloatType> r(ev.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=make_bounds(ev[i]);
    }
    return r;
}

Vector<ValidatedFloatType>const& make_singleton(const Vector<Interval>& ivlv) {
    return reinterpret_cast<Vector<ValidatedFloatType>const&>(ivlv);
}

Matrix<ValidatedFloatType>const& make_singleton(const Matrix<Interval>& ivlA) {
    return reinterpret_cast<Matrix<ValidatedFloatType>const&>(ivlA);
}

bool contains(const Vector<Interval>& v1, const Vector<ExactFloat>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!contains(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool contains(const Vector<Interval>& v1, const Vector<ValidatedFloat>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!contains(v1[i],v2[i])) { return false; }
    }
    return true;
}

UpperFloatType sup_error(const Vector<ValidatedFloatType>& x) {
    UpperFloatType e(0);
    for(uint i=0; i!=x.size(); ++i) {
        e=max(e,x[i].error());
    }
}

Vector<ExactFloatType> midpoint(const Vector<ValidatedFloatType>& x) {
    Vector<ExactFloatType> r(x.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=midpoint(x[i]);
    }
    return r;
}

bool models(const Vector<ValidatedFloatType>& x1, const Vector<ExactFloatType>& x2) {
    assert(x1.size()==x2.size());
    for(uint i=0; i!=x1.size(); ++i) {
        if(!models(x1[i],x2[i])) { return false; }
    }
    return true;
}

bool consistent(const Vector<ValidatedFloatType>& x1, const Vector<ValidatedFloatType>& x2) {
    assert(x1.size()==x2.size());
    for(uint i=0; i!=x1.size(); ++i) {
        if(!consistent(x1[i],x2[i])) { return false; }
    }
    return true;
}

bool inconsistent(const Vector<ValidatedFloatType>& x1, const Vector<ValidatedFloatType>& x2) {
    return !consistent(x1,x2);
}

bool refines(const Vector<ValidatedFloatType>& x1, const Vector<ValidatedFloatType>& x2) {
    assert(x1.size()==x2.size());
    for(uint i=0; i!=x1.size(); ++i) {
        if(!refines(x1[i],x2[i])) { return false; }
    }
    return true;
}

Vector<ValidatedFloatType> refinement(const Vector<ValidatedFloatType>& x1, const Vector<ValidatedFloatType>& x2) {
    assert(x1.size()==x2.size());
    Vector<ValidatedFloatType> r(x1.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=refinement(x1[i],x2[i]);
    }
    return r;
}




bool contains(const Vector<Interval>& v1, const Vector<Float>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!contains(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool subset(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!subset(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool intersect(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!intersect(v1[i],v2[i])) { return false; }
    }
    return true;
}


bool disjoint(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(disjoint(v1[i],v2[i])) { return true; }
    }
    return false;
}

bool overlap(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!overlap(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool covers(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!covers(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool inside(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        if(!inside(v1[i],v2[i])) { return false; }
    }
    return true;
}

bool empty(const Vector<Interval>& v)
{
    for(size_t i=0; i!=v.size(); ++i) {
        if(empty(v[i])) { return true; }
    }
    return false;
}


uint irmax(const Vector<Interval>& v) {
    uint imw(0);
    Float mw=v[0].width().raw();
    for(uint i=1; i!=v.size(); ++i) {
        if(v[i].width().raw()>mw) { imw=i; mw=v[i].width().raw(); }
    }
    return imw;
}


Vector<Interval> split(const Vector<Interval>& v, uint k, tribool lr) {
    ARIADNE_ASSERT(k<v.size());
    Vector<Interval> r(v);
    Float c=v[k].midpoint().raw();
    if(lr) {
        r[k].set_upper(c);
    } else if(!lr) {
        r[k].set_lower(c);
    } else {
        Float cl=(3*v[k].lower().raw()+v[k].upper().raw())/4;
        Float cu=(v[k].lower().raw()+3*v[k].upper().raw())/4;
        r[k].set_lower(cl);
        r[k].set_upper(cu);
    }
    return r;
}

std::pair< Vector<Interval>, Vector<Interval> > split(const Vector<Interval>& v, uint k) {
    ARIADNE_ASSERT(k<v.size());
    std::pair< Vector<Interval>, Vector<Interval> > r(v,v);
    Float c=v[k].midpoint().raw();
    r.first[k].set_upper(c);
    r.second[k].set_lower(c);
    return r;
}

Vector<Interval> split(const Vector<Interval>& v, tribool lr) {
    return split(v,irmax(v),lr);
}

std::pair< Vector<Interval>, Vector<Interval> > split(const Vector<Interval>& v) {
    return split(v,irmax(v));
}



Vector<ExactFloatType> midpoint(const Vector<Interval>& v)
{
    Vector<ExactFloatType> r(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
        r[i]=v[i].midpoint();
    }
    return r;
}

Vector<ExactFloatType> lower_bounds(const Vector<Interval>& v)
{
    Vector<ExactFloatType> r(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
        r[i]=v[i].lower();
    }
    return r;
}

Vector<ExactFloatType> upper_bounds(const Vector<Interval>& v)
{
    Vector<ExactFloatType> r(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
        r[i]=v[i].upper();
    }
    return r;
}

Vector<Interval> hull(const Vector<Float>& v1, const Vector<Float>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    Vector<Interval> r(v1.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        r[i]=hull(v1[i],v2[i]);
    }
    return r;
}

Vector<Interval> hull(const Vector<Interval>& v1, const Vector<Float>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    Vector<Interval> r(v1.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        r[i]=hull(v1[i],v2[i]);
    }
    return r;
}

Vector<Interval> hull(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    Vector<Interval> r(v1.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        r[i]=hull(v1[i],v2[i]);
    }
    return r;
}

Vector<Interval> intersection(const Vector<Interval>& v1, const Vector<Interval>& v2)
{
    ARIADNE_ASSERT(v1.size()==v2.size());
    Vector<Interval> r(v1.size());
    for(size_t i=0; i!=v1.size(); ++i) {
        r[i]=intersection(v1[i],v2[i]);
    }
    return r;
}

Float radius(const Vector<Interval>& v)
{
    Float r=0;
    for(size_t i=0; i!=v.size(); ++i) {
        r=Ariadne::max(r,v[i].radius().raw());
    }
    return r;
}

Float volume(const Vector<Interval>& v)
{
    Float r=1.0;
    for(size_t i=0; i!=v.size(); ++i) {
        r*=diam(v[i]);
    }
    return r;
}

} // namespace Ariadne
