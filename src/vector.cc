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

bool contains(const Vector<Interval>& v1, const Vector<Float>& v2)
{
  ARIADNE_ASSERT(v1.size()==v2.size());
  for(size_t i=0; i!=v1.size(); ++i) {
    if(!contains(v1[i],v2[i])) { return false; }
  }
  return true;
}

bool subset(const Vector<Float>& v1, const Vector<Interval>& v2)
{
  ARIADNE_ASSERT(v1.size()==v2.size());
  for(size_t i=0; i!=v1.size(); ++i) {
    if(!subset(v1[i],v2[i])) { return false; }
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

bool disjoint(const Vector<Interval>& v1, const Vector<Interval>& v2) 
{
  ARIADNE_ASSERT(v1.size()==v2.size());
  for(size_t i=0; i!=v1.size(); ++i) {
    if(v1[i].u<v2[i].l || v1[i].l>v2[i].u) {
      return true;
    }
  }
  return false;
}

Vector<Interval> intersection(const Vector<Interval>& v1, const Vector<Interval>& v2);
Vector<Float> midpoint(const Vector<Interval>& v);

} // namespace Ariadne
