/***************************************************************************
 *            python/export_map.cc
 *
 *  13 February 2006
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can rediself_ns::stribute it and/or modify
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

#include "evaluation/map.h"

#include "python/typedefs.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

struct RMap : RMapBase, wrapper<RMapBase>
{
  dimension_type argument_dimension() const { return this->get_override("argument_dimension")(); }
  dimension_type result_dimension() const { return this->get_override("result_dimension")(); }
  std::string name() const { return this->get_override("name")(); }
};

void export_map() {
  class_<RMap, boost::noncopyable>("Map")
    .def("argument_dimension", pure_virtual(&RMapBase::argument_dimension))
    .def("result_dimension", pure_virtual(&RMapBase::result_dimension))
  ;
}
