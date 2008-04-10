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
#include "geometry/euclidean_space.h"
#include "geometry/point.h"
#include "geometry/box.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class BS>
class SetWrapper
  : public SetInterface<BS>, public wrapper< SetInterface<BS> >
{
  typedef typename BS::state_type Pt;
  typedef typename BS::space_type space_type;
 public: 
  SetInterface<BS>* clone() const { return this->get_override("clone")(); }
  space_type space() const { return this->get_override("space")(); }
  tribool contains(const Pt& pt) const { return this->get_override("contains")(); }
  tribool superset(const  BS& r) const { return this->get_override("superset")(); }
  tribool intersects(const  BS& r) const { return this->get_override("intersects")(); }
  tribool disjoint(const  BS& r) const { return this->get_override("disjoint")(); }
  tribool subset(const  BS& r) const { return this->get_override("subset")(); }
  tribool bounded() const { return this->get_override("bounded")(); }
  BS bounding_box() const { return this->get_override("bounding_box")(); }
  std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

template<class R>
void export_set() 
{
  typedef Box<R> BS;
  class_<SetWrapper<BS>, boost::noncopyable>("Set")
    .def("space", pure_virtual(&SetInterface<BS>::space))
    .def("contains", pure_virtual(&SetInterface<BS>::contains))
    .def("superset", pure_virtual(&SetInterface<BS>::superset))
    .def("intersects", pure_virtual(&SetInterface<BS>::intersects))
    .def("disjoint", pure_virtual(&SetInterface<BS>::disjoint))
    .def("subset", pure_virtual(&SetInterface<BS>::subset))
    .def("bounding_box", pure_virtual(&SetInterface<BS>::bounding_box))
    .def(self_ns::str(self))
  ;
}



template void export_set<FloatPy>();
