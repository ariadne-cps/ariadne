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



Vector<Float> midpoint(const Vector<Interval>& v)
{
    Vector<Float> r(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
        r[i]=v[i].midpoint().raw();
    }
    return r;
}

Vector<Float> lower_bounds(const Vector<Interval>& v)
{
    Vector<Float> r(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
        r[i]=v[i].lower().raw();
    }
    return r;
}

Vector<Float> upper_bounds(const Vector<Interval>& v)
{
    Vector<Float> r(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
        r[i]=v[i].upper().raw();
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
