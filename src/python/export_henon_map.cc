/***************************************************************************
 *            python/export_henon_map.cc
 *
 *  8 February 2006
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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

#include "system/henon_map.h"

#include "python/typedefs.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;
  
typedef RPoint (RHenonMap::* PointMap) (const RPoint&) const;
typedef RRectangle (RHenonMap::* RectangleMap) (const RRectangle&) const;
typedef RParallelotope (RHenonMap::* ParallelotopeMap) (const RParallelotope&) const;
typedef RMatrix (RHenonMap::* PointDerivative) (const RPoint&) const;
typedef RIntervalMatrix (RHenonMap::* RectangleDerivative) (const RRectangle&) const;

void export_henon_map() {
  class_<RHenonMap, bases<RMapBase> >("HenonMap",init<Real,Real>())
    .def("argument_dimension", &RHenonMap::argument_dimension)
    .def("result_dimension", &RHenonMap::result_dimension)
    .def("__call__", PointMap(&RHenonMap::operator()))
    .def("__call__", RectangleMap(&RHenonMap::operator()))
    .def("__call__", ParallelotopeMap(&RHenonMap::operator()))
    .def("derivative", PointDerivative(&RHenonMap::derivative))
    .def("derivative", RectangleDerivative(&RHenonMap::derivative))
    .def(self_ns::str(self))
  ;
}
