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
 
#include "macros.h"
#include "numeric.h"
#include "vector.h"

template class boost::numeric::ublas::vector<Ariadne::Float>;
template class boost::numeric::ublas::vector<Ariadne::Interval>;

namespace Ariadne {


template<> Vector<Float>::Vector(size_t n, const double& t0, const double& t1, ...)
    : ublas::vector<Float>(n)
{
    assert(n>=2); va_list args; va_start(args,t1);
    (*this)[0]=t0; (*this)[1]=t1; 
    for(size_t i=2; i!=n; ++i) { (*this)[i]=va_arg(args,double); } 
    va_end(args);
}

template<> Vector<Interval>::Vector(size_t n, const double& t0, const double& t1, ...)
    : ublas::vector<Interval>(n)
{
    assert(n>=1); va_list args; va_start(args,t1);
    (*this)[0]=Interval(t0,t1);
    for(size_t i=1; i!=n; ++i) {
        double l=va_arg(args,double);
        double u=va_arg(args,double);
        (*this)[i]=Interval(l,u);
    }
    va_end(args);
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

Vector<Float> midpoint(const Vector<Interval>& v) 
{
    Vector<Float> r(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
        r[i]=v[i].midpoint();
    }
    return r;
}

Vector<Float> lower(const Vector<Interval>& v) 
{
    Vector<Float> r(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
        r[i]=v[i].lower();
    }
    return r;
}

Vector<Float> upper(const Vector<Interval>& v) 
{
    Vector<Float> r(v.size());
    for(size_t i=0; i!=v.size(); ++i) {
        r[i]=v[i].upper();
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
        r=Ariadne::max(r,v[i].radius());
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

Vector<Interval> intersection(const Vector<Interval>& v1, const Vector<Interval>& v2);
Vector<Float> midpoint(const Vector<Interval>& v);

} // namespace Ariadne
