/***************************************************************************
 *            python/export_list_set.cc
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



#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"

#include "python/typedefs.h"
#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;

inline RParallelotope plsg(const RParallelotopeListSet& s, int n) { return ::get_item(s,n); }
inline RZonotope zlsg(const RZonotopeListSet& s, int n) { return ::get_item(s,n); }
  
void export_list_set() {
  typedef bool (*RectLSBinPred) (const RRectangleListSet&, const RRectangleListSet&);
  typedef RRectangleListSet (*RectLSBinFun) (const RRectangleListSet&, const RRectangleListSet&);
  typedef void (RRectangleListSet::* RectLSadjRectLS) (const RRectangleListSet&);
  typedef void (RParallelotopeListSet::* PltpLSadjRect) (const RRectangle&);
  typedef void (RParallelotopeListSet::* PltpLSadjRectLS) (const RRectangleListSet&);
  typedef void (RParallelotopeListSet::* PltpLSadjPltp) (const RParallelotope&);
  typedef void (RParallelotopeListSet::* PltpLSadjPltpLS) (const RParallelotopeListSet&);


def("regular_intersection", RectLSBinFun(&regular_intersection));
  def("interiors_intersect", RectLSBinPred(&interiors_intersect));
  def("disjoint", RectLSBinPred(&disjoint));
  def("inner_subset", RectLSBinPred(&inner_subset));
  def("subset", RectLSBinPred(&subset));

  class_<RRectangleListSet>("RectangleListSet",init<int>())
    .def(init<RRectangle>())
    .def(init<RRectangleListSet>())
    .def(init<RGridRectangleListSet>())
    .def(init<RGridCellListSet>())
    .def(init<RGridMaskSet>())
    .def(init<RPartitionTreeSet>())
    .def("dimension", &RRectangleListSet::dimension)
    .def("push_back", &RRectangleListSet::push_back)
    .def("size", &RRectangleListSet::size)
    .def("__len__", &RRectangleListSet::size)
    .def("__getitem__", &RRectangleListSet::get, return_value_policy<copy_const_reference>())
    .def("__setitem__", &RRectangleListSet::set)
    .def("__iter__", iterator<RRectangleListSet>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;
  
  class_<RParallelotopeListSet>("ParallelotopeListSet",init<int>())
    .def(init<RParallelotope>())
    .def(init<RParallelotopeListSet>())
    .def("dimension", &RParallelotopeListSet::dimension)
    .def("push_back", &RParallelotopeListSet::push_back)
    .def("adjoin", PltpLSadjPltp(&RParallelotopeListSet::adjoin))
    .def("adjoin", PltpLSadjPltpLS(&RParallelotopeListSet::adjoin))
    .def("size", &RParallelotopeListSet::size)
    .def("__len__", &RParallelotopeListSet::size)
//    .def("__getitem__", &RParallelotopeListSet::get, return_value_policy<copy_const_reference>())
//    .def("__getitem__", &get_item<RParallelotopeListSet>)
    .def("__getitem__", &plsg)
    .def("__setitem__", &RParallelotopeListSet::set)
    .def("__iter__", iterator<RParallelotopeListSet>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;
  
  class_<RZonotopeListSet>("ZonotopeListSet",init<int>())
    .def(init<RZonotope>())
    .def(init<RZonotopeListSet>())
    .def("dimension", &RZonotopeListSet::dimension)
    .def("push_back", &RZonotopeListSet::push_back)
//    .def("adjoin", ZntpLSadjPltp(&RZonotopeListSet::adjoin))
//    .def("adjoin", ZntpLSadjPltpLS(&RZonotopeListSet::adjoin))
    .def("size", &RZonotopeListSet::size)
    .def("__len__", &RZonotopeListSet::size)
//    .def("__getitem__", &RParallelotopeListSet::get, return_value_policy<copy_const_reference>())
//    .def("__getitem__", &get_item<RParallelotopeListSet>)
    .def("__getitem__", &zlsg)
    .def("__setitem__", &RZonotopeListSet::set)
    .def("__iter__", iterator<RZonotopeListSet>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;

}
