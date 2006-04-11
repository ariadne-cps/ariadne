/***************************************************************************
 *            python/export_rectangle.cc
 *
 *  21 October 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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



#include "geometry/simplex.h"

#include "geometry/rectangle.h"
#include "geometry/polyhedron.h"

#include "python/typedefs.h"
using namespace Ariadne;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;

void export_simplex() {
  typedef bool (*SmplxSmplxBinPred) (const RSimplex&, const RSimplex&);
  typedef bool (*SmplxRectBinPred) (const RSimplex&, const RRectangle&);

  def("interiors_intersect", SmplxRectBinPred(&interiors_intersect));
  def("disjoint", SmplxSmplxBinPred(&disjoint));
  def("interiors_intersect", SmplxSmplxBinPred(&interiors_intersect));
  def("inner_subset", SmplxSmplxBinPred(&inner_subset));
  def("subset", SmplxSmplxBinPred(&subset));

  class_<RSimplex>("Simplex",init<int>())
    .def(init< std::vector<RPoint> >())
    .def(init< Ariadne::array<RPoint> >())
    .def(init<RSimplex>())
    .def(init<std::string>())
    .def("empty", &RSimplex::empty)
    .def("empty_interior", &RSimplex::empty_interior)
    .def("dimension", &RSimplex::dimension)
    .def("contains", &RSimplex::contains)
    .def("interior_contains", &RSimplex::interior_contains)
    .def(self_ns::str(self))    // __self_ns::str__
  ;
}
