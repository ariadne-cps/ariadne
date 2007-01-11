/***************************************************************************
 *            python/export_polyhedron.cc
 *
 *  21 October 2005
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

#include "real_typedef.h"

#include "geometry/rectangle.h"
#include "geometry/polyhedron.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_polyhedron() 
{
  typedef LinearAlgebra::Vector<R> RVector;
  typedef LinearAlgebra::Matrix<R> RMatrix;
  typedef Rectangle<R> RRectangle;
  typedef Polyhedron<R> RPolyhedron;

  def("disjoint", (tribool(*)(const RPolyhedron&, const RPolyhedron&))(&disjoint));
  def("subset", (tribool(*)(const RPolyhedron&, const RPolyhedron&))(&subset));

  class_<RPolyhedron>("Polyhedron",init<int>())
    .def(init<RMatrix,RVector>())
    .def(init<RRectangle>())
    .def(init<RPolyhedron>())
    .def("dimension", &RPolyhedron::dimension)
    .def(self_ns::str(self))
  ;
}

template void export_polyhedron<Real>();
