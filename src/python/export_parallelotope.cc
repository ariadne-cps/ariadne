/***************************************************************************
 *            python/export_parallelotope.cc
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



#include "geometry/parallelotope.h"

#include "geometry/rectangle.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"


#include "python/typedefs.h"
#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;

void export_parallelotope() {
  typedef bool (*PltpPltpBinPred) (const RParallelotope&, const RParallelotope&);
  typedef bool (*PltpRectBinPred) (const RParallelotope&, const RRectangle&);
  typedef bool (*RectPltpBinPred) (const RRectangle&, const RParallelotope&);
  
  def("interiors_intersect", PltpPltpBinPred(&interiors_intersect));
  def("interiors_intersect", PltpRectBinPred(&interiors_intersect));
  def("interiors_intersect", RectPltpBinPred(&interiors_intersect));
  def("disjoint", PltpPltpBinPred(&disjoint));
  def("disjoint", PltpRectBinPred(&disjoint));
  def("disjoint", RectPltpBinPred(&disjoint));
  def("inner_subset", PltpPltpBinPred(&inner_subset));
  def("inner_subset", PltpRectBinPred(&inner_subset));
  def("inner_subset", RectPltpBinPred(&inner_subset));

  def("subset", PltpPltpBinPred(&subset));
  def("subset", PltpRectBinPred(&subset));
  def("subset", RectPltpBinPred(&subset));

  class_<RParallelotope>("Parallelotope",init<int>())
    .def(init<RPoint,RMatrix>())
    .def(init<RParallelotope>())
    .def(init<RRectangle>())
    .def(init<std::string>())
    .def("bounding_box", &RParallelotope::bounding_box)
    .def("subdivide", &RParallelotope::subdivide)
    .def("empty", &RParallelotope::empty)
    .def("empty_interior", &RParallelotope::empty_interior)
    .def("dimension", &RParallelotope::dimension)
    .def("contains", &RParallelotope::contains)
    .def("interior_contains", &RParallelotope::interior_contains)
    .def(self_ns::str(self))    // __self_ns::str__
  ;
}
