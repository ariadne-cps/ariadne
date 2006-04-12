/***************************************************************************
 *            python/export_grid_set.cc
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
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */



#include "linear_algebra/linear_algebra.h"

#include "geometry/parallelotope.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"


#include "python/typedefs.h"
#include "python/python_utilities.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

struct RGrid : RGridBase, wrapper<RGridBase>
{
  RGrid* clone() const { return this->get_override("clone")(); }
  dimension_type dimension() const { return this->get_override("dimension")(); }
  real_type subdivision_coordinate(dimension_type d, index_type n) const { return this->get_override("subdivision_coordinate")(); }
  index_type subdivision_interval(dimension_type d, const real_type& x) const { return this->get_override("subdivision_interval")(); }
  bool bounds_enclose(const RRectangle& r) const { return this->get_override("bounds_enclose")(); }
  std::ostream& write(std::ostream& os) const { return this->get_override("write")(); }
};

inline RGridRectangle over_approximation_rectangle(const RRectangle& r, const RGridBase& g) {
  return Ariadne::Geometry::over_approximation(r,g);
}
inline RGridCellListSet over_approximation_parallelotope(const RParallelotope& p, const RGridBase& g) {
  return Ariadne::Geometry::over_approximation(p,g);
}
inline RGridMaskSet over_approximation_rectangle_list_set(const RRectangleListSet& rls, const RFiniteGrid& g) {
  return Ariadne::Geometry::over_approximation(rls,g);
}
inline RGridMaskSet over_approximation_parallelotope_list_set(const RParallelotopeListSet& pls, const RFiniteGrid& g) {
  return Ariadne::Geometry::over_approximation(pls,g);
}

inline void grid_cell_list_set_adjoin_grid_cell(RGridMaskSet& gms, const RGridCell& gc) {
  return gms.adjoin(gc);
}
inline RGridMaskSet grid_mask_set_join(const RGridMaskSet& gms1, const RGridMaskSet& gms2) {
  return Ariadne::Geometry::join(gms1,gms2);
}
inline RGridMaskSet grid_mask_set_difference(const RGridMaskSet& gms1, const RGridMaskSet& gms2) {
  return Ariadne::Geometry::difference(gms1,gms2);
}
inline RGridMaskSet grid_mask_set_regular_intersection(const RGridMaskSet& gms1, const RGridMaskSet& gms2) {
  return Ariadne::Geometry::regular_intersection(gms1,gms2);
}
inline void grid_mask_set_adjoin_grid_cell(RGridMaskSet& gms, const RGridCell& gc) {
  gms.adjoin(gc);
}
inline void grid_mask_set_adjoin_grid_rectangle(RGridMaskSet& gms, const RGridRectangle& gr) {
  gms.adjoin(gr);
}
inline void grid_mask_set_adjoin_grid_cell_list_set(RGridMaskSet& gms, const RGridCellListSet& gcls) {
  gms.adjoin(gcls);
}
inline void grid_mask_set_adjoin_grid_rectangle_list_set(RGridMaskSet& gms, const RGridRectangleListSet& grls) {
  gms.adjoin(grls);
}
inline void grid_mask_set_adjoin_grid_mask_set(RGridMaskSet& gms, const RGridMaskSet& agms) {
  gms.adjoin(agms);
}


void export_grid_set() {

  class_<RGrid, boost::noncopyable>("Grid")
    .def("dimension", pure_virtual(&RGridBase::dimension))
    .def("subdivision_coordinate", pure_virtual(&RGridBase::subdivision_coordinate))
    .def("subdivision_index", pure_virtual(&RGridBase::subdivision_index))
    .def("subdivision_lower_index", &RGridBase::subdivision_lower_index)
    .def("subdivision_upper_index", &RGridBase::subdivision_upper_index)
    .def("bounds_enclose", pure_virtual(&RGridBase::bounds_enclose))
    .def(self_ns::str(self))    // __str__
    ;

  class_< RFiniteGrid, bases<RGridBase> >("FiniteGrid",init<RRectangle,SizeArray>())
    .def(init<RRectangle,uint>())
    .def(init<RRectangleListSet>())
    .def("dimension", &RFiniteGrid::dimension)
    .def("subdivision_coordinate", &RFiniteGrid::subdivision_coordinate)
    .def("subdivision_index", &RFiniteGrid::subdivision_index)
    .def("subdivision_lower_index", &RFiniteGrid::subdivision_lower_index)
    .def("subdivision_upper_index", &RFiniteGrid::subdivision_upper_index)
    .def("bounds_enclose", &RFiniteGrid::bounds_enclose)
    .def(self_ns::str(self))    // __str__
    ;

  class_< RInfiniteGrid, bases<RGridBase> >("InfiniteGrid",init<const array<Real>&>())
    .def(init<uint,Real>())
    .def("dimension", &RInfiniteGrid::dimension)
    .def("subdivision_coordinate", &RInfiniteGrid::subdivision_coordinate)
    .def("subdivision_index", &RInfiniteGrid::subdivision_index)
    .def("subdivision_lower_index", &RInfiniteGrid::subdivision_lower_index)
    .def("subdivision_upper_index", &RInfiniteGrid::subdivision_upper_index)
    .def("bounds_enclose", &RInfiniteGrid::bounds_enclose)
    .def(self_ns::str(self))    // __str__
    ;

  class_<RGridCell>("GridCell",init<const RGrid&,LatticeCell>())
    .def("dimension", &RGridCell::dimension)
    .def(self_ns::str(self))    // __str__
    ;
  
  class_<RGridRectangle>("GridRectangle",init<const RGridBase&,LatticeRectangle>())
    .def(init<RGridCell>()) 
    .def("dimension", &RGridRectangle::dimension)
    .def("position", &RGridRectangle::position,return_value_policy<copy_const_reference>())
    .def(self_ns::str(self))    // __str__
    ;
  
  class_<RGridCellListSet>("GridCellListSet",init<const RGrid&>())
    .def(init<RGridMaskSet>())
    .def(init<RGridCellListSet>())
    .def(init<RGridRectangleListSet>())
    .def(init<RRectangleListSet>())
    .def("lattice_set", &RGridCellListSet::lattice_set,return_value_policy<copy_const_reference>())
    .def("dimension", &RGridCellListSet::dimension)
    .def("adjoin", &RGridCellListSet::adjoin)
    .def("size", &RGridCellListSet::size)
    .def("__len__", &RGridCellListSet::size)
    .def("__getitem__", &get_item<RGridCellListSet>)
    .def("__iter__", iterator<RGridCellListSet>())
    .def(self_ns::str(self))    // __str__
    ;
  
  class_<RGridRectangleListSet>("GridRectangleListSet",init<const RGrid&>())
    .def(init<RGridRectangleListSet>())
    .def(init<RRectangleListSet>())
    .def(init<RPartitionTreeSet>())
    .def("lattice_set", &RGridRectangleListSet::lattice_set,return_value_policy<copy_const_reference>())
    .def("dimension", &RGridRectangleListSet::dimension)
    .def("adjoin", &RGridRectangleListSet::adjoin)
    .def("size", &RGridRectangleListSet::size)
    .def("__len__", &RGridRectangleListSet::size)
    .def("__getitem__", &get_item<RGridRectangleListSet>)
    .def("__iter__", iterator<RGridRectangleListSet>())
    .def(self_ns::str(self))    // __str__
    ;
    
  class_<RGridMaskSet>("GridMaskSet",init<const RFiniteGrid&>())
    .def(init<const RGridBase&, LatticeRectangle>())
    .def(init<const RInfiniteGrid&, LatticeRectangle>())
    .def(init<RGridRectangleListSet>())
    .def(init<RGridCellListSet>())
    .def(init<RGridMaskSet>())
    .def(init<RRectangleListSet>())
    .def("bounding_box", &RGridMaskSet::bounding_box)
    .def("empty", &RGridMaskSet::empty)
    .def("dimension", &RGridMaskSet::dimension)
    .def("clear", &RGridMaskSet::clear)
    .def("bounds", &RGridMaskSet::bounds,return_value_policy<copy_const_reference>())
    .def("lattice_set", &RGridMaskSet::lattice_set,return_value_policy<copy_const_reference>())
    .def("adjoin", &grid_mask_set_adjoin_grid_cell)
    .def("adjoin", &grid_mask_set_adjoin_grid_rectangle)
    .def("adjoin", &grid_mask_set_adjoin_grid_cell_list_set)
    .def("adjoin", &grid_mask_set_adjoin_grid_rectangle_list_set)
    .def("adjoin", &grid_mask_set_adjoin_grid_mask_set)
    .def("neighbourhood", &RGridMaskSet::neighbourhood)
    .def("adjoining", &RGridMaskSet::adjoining)
    .def("size", &RGridMaskSet::size)
    .def("__len__", &RGridMaskSet::size)
    .def("__getitem__", &get_item<RGridMaskSet>)
    .def("__iter__", iterator<RGridMaskSet>())
    .def(self_ns::str(self))    // __str__
    ;
    
  def("join",&grid_mask_set_join);
  def("difference",&grid_mask_set_difference);
  def("regular_intersection",&grid_mask_set_regular_intersection);
    
  def("over_approximation",&over_approximation_rectangle);
  def("over_approximation",&over_approximation_parallelotope);
  def("over_approximation",&over_approximation_rectangle_list_set);
  def("over_approximation",&over_approximation_parallelotope_list_set);
}
