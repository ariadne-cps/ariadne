/***************************************************************************
 *            point.h
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
 
/*! \file point.h
 *  \brief Points in Euclidean space.
 */

#ifndef ARIADNE_POINT_H
#define ARIADNE_POINT_H

#include "numeric.h"
#include "vector.h"

namespace Ariadne {

class Point : public Vector<Float>
{
 public:
  typedef Float real_type;
  Point() : Vector<Float>() { }
  template<class T> Point(const T& t) : Vector<Float>(t) { }
  template<class T1, class T2> Point(const T1& t1, const T2& t2) : Vector<Float>(t1,t2) { }
  uint dimension() const { return this->size(); }
  Vector<Float> centre() const { return *this; }
};

} // namespace Ariadne

#endif // ARIADNE_POINT_H
