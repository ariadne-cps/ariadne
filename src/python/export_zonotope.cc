/***************************************************************************
 *            python/export_parallelotope.cc
 *
 *  22 March 2006
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



#include "geometry/zonotope.h"

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"


#include "python/typedefs.h"
#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;

template <template <typename> class BS>
inline
Zonotope<Real> 
touching_intersection(const Zonotope<Real> &a, 
		      const BS<Real> &b) {
	
  if (interiors_intersect(a,b))
    return a;

  return Zonotope<Real>(a.dimension());
}

template Zonotope<Real> touching_intersection(const Zonotope<Real> &,  
      	                                      const Rectangle<Real> &);
template Zonotope<Real> touching_intersection(const Zonotope<Real> &,  
               	                              const Parallelotope<Real> &);
template Zonotope<Real> touching_intersection(const Zonotope<Real> &,  
                                              const Zonotope<Real> &);


void export_zonotope() {
  typedef RZonotope (*ZntpZntpBinFunc) (const RZonotope&, const RZonotope&);
  typedef bool (*ZntpZntpBinPred) (const RZonotope&, const RZonotope&);
  typedef bool (*ZntpRectBinPred) (const RZonotope&, const RRectangle&);
  typedef bool (*RectZntpBinPred) (const RRectangle&, const RZonotope&);
  typedef bool (*PltpZntpBinPred) (const RParallelotope&, const RZonotope&);
  typedef bool (*ZntpPltpBinPred) (const RZonotope&, const RParallelotope&);

  typedef RZonotope (*ZntpRectBinFun) (const RZonotope&, const RRectangle&);
  typedef RZonotope (*ZntpPltpBinFun) (const RZonotope&, const RParallelotope&);
  typedef RZonotope (*ZntpZntpBinFun) (const RZonotope&, const RZonotope&);
  
  typedef bool (RZonotope::*RectPred)(const RRectangle &) const;
  typedef bool (RZonotope::*PointPred) (const RPoint&) const;

 
  def("interiors_intersect", ZntpZntpBinPred(&interiors_intersect));
  def("interiors_intersect", ZntpRectBinPred(&interiors_intersect));
  def("interiors_intersect", RectZntpBinPred(&interiors_intersect));
  def("interiors_intersect", PltpZntpBinPred(&interiors_intersect));
  def("interiors_intersect", ZntpPltpBinPred(&interiors_intersect));
  def("disjoint", ZntpZntpBinPred(&disjoint));
  def("disjoint", ZntpRectBinPred(&disjoint));
  def("disjoint", RectZntpBinPred(&disjoint));
  def("disjoint", PltpZntpBinPred(&disjoint));
  def("disjoint", ZntpPltpBinPred(&disjoint));
  def("inner_subset", ZntpZntpBinPred(&inner_subset));
  def("inner_subset", ZntpRectBinPred(&inner_subset));
  def("inner_subset", RectZntpBinPred(&inner_subset));
  def("inner_subset", PltpZntpBinPred(&inner_subset));
  def("inner_subset", ZntpPltpBinPred(&inner_subset));
  def("touching_intersection", ZntpRectBinFun(&touching_intersection));
  def("touching_intersection", ZntpPltpBinFun(&touching_intersection));
  def("touching_intersection", ZntpZntpBinFun(&touching_intersection));
  def("subset", ZntpZntpBinPred(&subset));
  def("subset", ZntpRectBinPred(&subset));
  def("subset", RectZntpBinPred(&subset));
  def("subset", PltpZntpBinPred(&subset));
  def("subset", ZntpPltpBinPred(&subset));
  def("minkowski_sum", ZntpZntpBinFunc(&minkowski_sum));

  class_<RZonotope>("Zonotope",init<int>())
    .def(init<RPoint,RMatrix>())
    .def(init<RZonotope>())
    .def(init<RParallelotope>())
    .def(init<RRectangle>())
    .def(init<std::string>())
    .def("bounding_box", &RZonotope::bounding_box)
    .def("subdivide", &RZonotope::subdivide)
    .def("empty", &RZonotope::empty)
    .def("empty_interior", &RZonotope::empty_interior)
    .def("dimension", &RZonotope::dimension)
    .def("contains", RectPred(&RZonotope::contains) )
    .def("contains", PointPred(&RZonotope::contains))
    //.def("==", &RZonotope::operator==)
    //.def("!=", &RZonotope::operator!=)
    .def("interior_contains", &RZonotope::interior_contains)
    .def(self_ns::str(self))    // __self_ns::str__
  ;
}
