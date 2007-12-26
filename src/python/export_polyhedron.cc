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

#include "python/float.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/box.h"
#include "geometry/polytope.h"
#include "geometry/polyhedron.h"
#include "geometry/set_interface.h"

#include "python/utilities.h"
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;


template<class R>
void export_polyhedron() 
{
  typedef Vector<R> RVector;
  typedef Matrix<R> RMatrix;
  typedef SetInterface<R> RSetInterface;
  typedef Polyhedron<R> RPolyhedron;
  typedef Box<R> RBox;
  typedef Polytope<R> RPolytope;
  
  def("disjoint", (tribool(*)(const RPolyhedron&, const RBox&))(&disjoint));
  def("disjoint", (tribool(*)(const RBox&, const RPolyhedron&))(&disjoint));
  def("disjoint", (tribool(*)(const RPolyhedron&, const RPolyhedron&))(&disjoint));
  def("subset", (tribool(*)(const RBox&, const RPolyhedron&))(&subset));
  def("subset", (tribool(*)(const RPolytope&, const RPolyhedron&))(&subset));
  def("subset", (tribool(*)(const RPolyhedron&, const RPolyhedron&))(&subset));
  def("closed_intersection", (RPolyhedron(*)(const RPolyhedron&, const RPolyhedron&))(&closed_intersection));

  //class_< RPolyhedron,bases<RSet> >("Polyhedron",init<int>())
  class_< RPolyhedron >("Polyhedron",init<int>())
    .def(init<RMatrix,RVector>())
    .def(init<RPolyhedron>())
    .def(init<RBox>())
    .def("dimension", &RPolyhedron::dimension)
    .def(self_ns::str(self))
  ;
  
}

template void export_polyhedron<FloatPy>();
