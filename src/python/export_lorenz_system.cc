/***************************************************************************
 *            python/export_lorenz_system.cc
 *
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

#include "evaluation/lorenz_system.h"

#include "python/typedefs.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

typedef RVector (RLorenzSystem::* PointMap) (const RPoint&) const;
typedef RIntervalVector (RLorenzSystem::* RectangleMap) (const RRectangle&) const;
typedef RMatrix (RLorenzSystem::* PointDerivative) (const RPoint&) const;
typedef RIntervalMatrix (RLorenzSystem::* RectangleDerivative) (const RRectangle&) const;

void export_lorenz_system() {
  class_<RLorenzSystem, bases<RVectorFieldBase> >("LorenzSystem",init<Real,Real,Real>())
    .def("dimension", &RLorenzSystem::dimension)
    .def("__call__", PointMap(&RLorenzSystem::apply))
    .def("__call__", RectangleMap(&RLorenzSystem::apply))
    .def("apply", PointMap(&RLorenzSystem::apply))
    .def("apply", RectangleMap(&RLorenzSystem::apply))
    .def("derivative", PointDerivative(&RLorenzSystem::derivative))
    .def("derivative", RectangleDerivative(&RLorenzSystem::derivative))
    .def(self_ns::str(self))    // __self_ns::str__
  ;
}
