/***************************************************************************
 *            python/export_orbit.cc
 *
 *  Copyright  2007  Pieter Collins
 * 
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "python/float.h"

#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/orbit.h"

#include "output/epsstream.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class T, class R>
struct OrbitWrapper
{
  boost::shared_ptr< const OrbitInterface<T> > _ptr;
};

template<class T, class R>
std::ostream& operator<<(std::ostream& os, const OrbitWrapper<T,R>& orbit) {
  return os << *orbit._ptr;
}


template<class R>
void export_orbit() 
{
  typedef Zonotope<R> ZBS;
  class_< Orbit<Integer,ZBS>, boost::noncopyable > orbit_class("Orbit",init<>());
  orbit_class.def(self_ns::str(self));
}

template void export_orbit<FloatPy>();
