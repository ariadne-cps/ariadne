/***************************************************************************
 *            python/export_set.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include "geometry/set_interface.h"
#include "geometry/point.h"
#include "geometry/rectangle.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
class SetWrapper
  : public SetInterface<R>, public wrapper< SetInterface<R> >
{
 public: 
  SetInterface<R>* clone() const { return this->get_override("clone")(); }
  dimension_type dimension() const { return this->get_override("dimension")(); }
  tribool contains(const Point<R>& pt) const { return this->get_override("contains")(); }
  tribool superset(const Box<R>& r) const { return this->get_override("superset")(); }
  tribool intersects(const Box<R>& r) const { return this->get_override("intersects")(); }
  tribool disjoint(const Box<R>& r) const { return this->get_override("disjoint")(); }
  tribool subset(const Box<R>& r) const { return this->get_override("subset")(); }
  tribool bounded() const { return this->get_override("bounded")(); }
  Box<R> bounding_box() const { return this->get_override("bounding_box")(); }
  std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

template<class R>
void export_set() 
{
  class_<SetWrapper<R>, boost::noncopyable>("Set")
    .def("dimension", pure_virtual(&SetInterface<R>::dimension))
    .def("contains", pure_virtual(&SetInterface<R>::contains))
    .def("superset", pure_virtual(&SetInterface<R>::superset))
    .def("intersects", pure_virtual(&SetInterface<R>::intersects))
    .def("disjoint", pure_virtual(&SetInterface<R>::disjoint))
    .def("subset", pure_virtual(&SetInterface<R>::subset))
    .def("bounded", pure_virtual(&SetInterface<R>::bounded))
    .def("bounding_box", pure_virtual(&SetInterface<R>::bounding_box))
    .def(self_ns::str(self))
  ;
}



template void export_set<FloatPy>();
