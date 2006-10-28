/***************************************************************************
 *            python/export_map.cc
 *
 *  13 February 2006
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

#include "real_typedef.h"

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "system/map.h"

using namespace Ariadne;
using namespace Ariadne::System;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
class MapWrapper : public Map<R>, public wrapper< Map<R> >
{
 public:
  dimension_type argument_dimension() const { return this->get_override("argument_dimension")(); }
  dimension_type result_dimension() const { return this->get_override("result_dimension")(); }
  size_type smoothness() const { return this->get_override("smoothness")(); }
  std::string name() const { return this->get_override("name")(); }
};

template<class R>
void export_map() 
{
  class_<MapWrapper<R>, boost::noncopyable>("Map")
    .def("argument_dimension", pure_virtual(&MapWrapper<R>::argument_dimension))
    .def("result_dimension", pure_virtual(&MapWrapper<R>::result_dimension))
    .def("smoothness", pure_virtual(&MapWrapper<R>::smoothness))
  ;
}

template void export_map<Real>();
