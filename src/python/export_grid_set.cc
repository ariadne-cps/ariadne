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

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"


#include "python/typedefs.h"
#include "python/python_utilities.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

struct RGridWrap : RGrid, wrapper<RGrid>
{
  RGrid* clone() const { return this->get_override("clone")(); }
  dimension_type dimension() const { return this->get_override("dimension")(); }
  grid_type type() const { return this->get_override("type")(); }
  real_type subdivision_coordinate(dimension_type d, index_type n) const { return this->get_override("subdivision_coordinate")(); }
  index_type subdivision_interval(dimension_type d, const real_type& x) const { return this->get_override("subdivision_interval")(); }
  bool bounds_enclose(const RRectangle& r) const { return this->get_override("bounds_enclose")(); }
  std::ostream& write(std::ostream& os) const { return this->get_override("write")(); }
  std::istream& read(std::istream& is) { return this->get_override("read")(); }
};

typedef RGridBlock (*ApprxRctGridFunc) (const RRectangle&, const RGrid&);
typedef RGridCellListSet (*ApprxZltpGridFunc) (const RZonotope&, const RGrid&);
typedef RGridCellListSet (*ApprxPltpGridFunc) (const RPolytope&, const RGrid&);
//typedef RGridCellListSet (*ApprxPlhdGridFunc) (const RPolyhedron&, const RGrid&);
typedef RGridBlock (*ApprxBlkFGridFunc) (const RRectangle&, const RFiniteGrid&);
typedef RGridMaskSet (*ApprxLSRctFGridFunc) (const RRectangleListSet&, const RFiniteGrid&);
typedef RGridMaskSet (*ApprxLSPrltpFGridFunc) (const RParallelotopeListSet&, const RFiniteGrid&);
typedef RGridMaskSet (*ApprxLSZltpFGridFunc) (const RZonotopeListSet&, const RFiniteGrid&);
typedef RGridMaskSet (*ApprxGridltpFGridFunc) (const RGridMaskSet&, const RFiniteGrid&);

typedef void (RGridMaskSet::*GMSAdjCellFunc) (const RGridCell&);
typedef void (RGridMaskSet::*GMSAdjBlkFunc) (const RGridBlock&);
typedef void (RGridMaskSet::*GMSAdjCellLSFunc) (const RGridCellListSet&);
typedef void (RGridMaskSet::*GMSAdjBlkLSFunc) (const RGridBlockListSet&);
typedef void (RGridMaskSet::*GMSAdjGMSFunc) (const RGridMaskSet&);

typedef bool (*GMSBinPred) (const RGridMaskSet&, const RGridMaskSet&);
typedef bool (*GMSRctPred) (const RGridMaskSet&, const RRectangle&);

typedef RGridMaskSet (*GMSBinFunc) (const RGridMaskSet&, const RGridMaskSet&);

typedef RGridMaskSet (*GMSRctFunc)(const RGridMaskSet&, const RRectangle&);
typedef RGridMaskSet (*GMSPrltpFunc)(const RGridMaskSet&, const RParallelotope&);
typedef RGridMaskSet (*GMSZntpFunc)(const RGridMaskSet&, const RZonotope&);
typedef RGridMaskSet (*GMSPltpFunc)(const RGridMaskSet&, const RPolytope&);
//typedef RGridMaskSet (*GMSPlhdFunc)(const RGridMaskSet&, const RPolyhedron&);
typedef RGridMaskSet (*GMSLSRctFunc)(const RGridMaskSet&, 
		                      const RRectangleListSet&);
typedef RGridMaskSet (*GMSLSPltpFunc)(const RGridMaskSet&, 
		                      const RParallelotopeListSet&);
typedef RGridMaskSet (*GMSLSZntpFunc)(const RGridMaskSet&, 
		                      const RZonotopeListSet&);
typedef RGridMaskSet (*GMSLSPlhdFunc)(const RGridMaskSet&, 
		                      const RPolytopeListSet&);

void export_grid_set() {

  class_<RGridWrap, boost::noncopyable>("Grid")
    .def("dimension", pure_virtual(&RGrid::dimension))
    .def("subdivision_coordinate", pure_virtual(&RGrid::subdivision_coordinate))
    .def("subdivision_index", pure_virtual(&RGrid::subdivision_index))
    .def("subdivision_lower_index", &RGrid::subdivision_lower_index)
    .def("subdivision_upper_index", &RGrid::subdivision_upper_index)
    .def("bounds_enclose", pure_virtual(&RGrid::bounds_enclose))
    .def(self_ns::str(self))    // __str__
    ;

  class_< RIrregularGrid, bases<RGrid> >("IrregularGrid",init<RRectangle,SizeArray>())
    .def(init<RRectangle,uint>())
    .def(init<RRectangleListSet>())
    .def("dimension", &RIrregularGrid::dimension)
    .def("subdivision_coordinate", &RIrregularGrid::subdivision_coordinate)
    .def("subdivision_index", &RIrregularGrid::subdivision_index)
    .def("subdivision_lower_index", &RIrregularGrid::subdivision_lower_index)
    .def("subdivision_upper_index", &RIrregularGrid::subdivision_upper_index)
    .def("bounds_enclose", &RIrregularGrid::bounds_enclose)
    .def(self_ns::str(self))    // __str__
    ;

  class_< RRegularGrid, bases<RGrid> >("RegularGrid",init<const array<Real>&>())
    .def(init<uint,Real>())
    .def(init<uint,double>())
    .def("dimension", &RRegularGrid::dimension)
    .def("subdivision_coordinate", &RRegularGrid::subdivision_coordinate)
    .def("subdivision_index", &RRegularGrid::subdivision_index)
    .def("subdivision_lower_index", &RRegularGrid::subdivision_lower_index)
    .def("subdivision_upper_index", &RRegularGrid::subdivision_upper_index)
    .def("bounds_enclose", &RRegularGrid::bounds_enclose)
    .def(self_ns::str(self))    // __str__
    ;

  class_< RFiniteGrid>("FiniteGrid",init<RRectangle,size_type>())
    .def(init<const RGrid&,LatticeBlock>())
    .def(init<const RGrid&,RRectangle>())
    .def("dimension", &RFiniteGrid::dimension)
    .def("subdivision_coordinate", &RFiniteGrid::subdivision_coordinate)
    .def("subdivision_lower_index", &RFiniteGrid::subdivision_lower_index)
    .def("subdivision_upper_index", &RFiniteGrid::subdivision_upper_index)
    .def("grid", &RFiniteGrid::grid,return_internal_reference<>())
    .def("bounds", &RFiniteGrid::bounds,return_value_policy<copy_const_reference>())
    .def(self_ns::str(self))    // __str__
    ;

  class_<RGridCell>("GridCell",init<const RGrid&,LatticeCell>())
    .def("dimension", &RGridCell::dimension)
    .def("lattice_set", &RGridCell::lattice_set,return_value_policy<copy_const_reference>())
    .def(self_ns::str(self))    // __str__
    ;
  
  class_<RGridBlock>("GridBlock",init<const RGrid&,LatticeBlock>())
    .def(init<RGridCell>()) 
    .def("dimension", &RGridBlock::dimension)
    .def("lattice_set", &RGridBlock::lattice_set,return_value_policy<copy_const_reference>())
    .def(self_ns::str(self))    // __str__
    ;
  
  class_<RGridCellListSet>("GridCellListSet",init<const RGrid&>())
    .def(init<RGridMaskSet>())
    .def(init<RGridCellListSet>())
    .def(init<RGridBlockListSet>())
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
  
  class_<RGridBlockListSet>("GridBlockListSet",init<const RGrid&>())
    .def(init<RGridBlockListSet>())
    .def(init<RRectangleListSet>())
    .def(init<RPartitionTreeSet>())
    .def("lattice_set", &RGridBlockListSet::lattice_set,return_value_policy<copy_const_reference>())
    .def("dimension", &RGridBlockListSet::dimension)
    .def("adjoin", &RGridBlockListSet::adjoin)
    .def("size", &RGridBlockListSet::size)
    .def("__len__", &RGridBlockListSet::size)
    .def("__getitem__", &get_item<RGridBlockListSet>)
    .def("__iter__", iterator<RGridBlockListSet>())
    .def(self_ns::str(self))    // __str__
    ;
    
  class_<RGridMaskSet>("GridMaskSet",init<const RFiniteGrid&>())
    .def(init<RGridBlockListSet>())
    .def(init<RGridCellListSet>())
    .def(init<RGridMaskSet>())
    .def(init<RRectangleListSet>())
    .def("bounding_box", &RGridMaskSet::bounding_box)
    .def("empty", &RGridMaskSet::empty)
    .def("dimension", &RGridMaskSet::dimension)
    .def("clear", &RGridMaskSet::clear)
    .def("block", &RGridMaskSet::block,return_value_policy<copy_const_reference>())
    .def("lattice_set", &RGridMaskSet::lattice_set,return_value_policy<copy_const_reference>())
    .def("adjoin", GMSAdjCellFunc(&RGridMaskSet::adjoin))
    .def("adjoin", GMSAdjBlkFunc(&RGridMaskSet::adjoin))
    .def("adjoin", GMSAdjCellLSFunc(&RGridMaskSet::adjoin))
    .def("adjoin", GMSAdjBlkLSFunc(&RGridMaskSet::adjoin))
    .def("adjoin", GMSAdjGMSFunc(&RGridMaskSet::adjoin))
    .def("neighbourhood", &RGridMaskSet::neighbourhood)
    .def("adjoining", &RGridMaskSet::adjoining)
    .def("size", &RGridMaskSet::size)
    .def("__len__", &RGridMaskSet::size)
    .def("__getitem__", &get_item<RGridMaskSet>)
    .def("__iter__", iterator<RGridMaskSet>())
    .def(self_ns::str(self))    // __str__
    ;

  def("join",GMSBinFunc(&Geometry::join));
  def("difference",GMSBinFunc(&Geometry::difference));
  def("regular_intersection",GMSBinFunc(&Geometry::regular_intersection));
  def("interiors_intersect",GMSBinPred(&Geometry::interiors_intersect));
  def("interiors_intersect",GMSRctPred(&Geometry::interiors_intersect));

  def("over_approximation",ApprxRctGridFunc(&Geometry::over_approximation));
  def("over_approximation",ApprxZltpGridFunc(&Geometry::over_approximation));
  def("over_approximation",ApprxPltpGridFunc(&Geometry::over_approximation));
  //def("over_approximation",ApprxPlhdGridFunc(&Geometry::over_approximation));
  def("over_approximation",ApprxLSRctFGridFunc(&Geometry::over_approximation));
  def("over_approximation",ApprxLSPrltpFGridFunc(&Geometry::over_approximation));
  def("over_approximation",ApprxLSZltpFGridFunc(&Geometry::over_approximation));
  def("over_approximation",ApprxGridltpFGridFunc(&Geometry::over_approximation));

  def("under_approximation",ApprxRctGridFunc(&Geometry::under_approximation));
  def("under_approximation",ApprxZltpGridFunc(&Geometry::under_approximation));
  def("under_approximation",ApprxPltpGridFunc(&Geometry::under_approximation));
  //def("under_approximation",ApprxPlhdGridFunc(&Geometry::under_approximation));
  def("under_approximation",ApprxLSRctFGridFunc(&Geometry::under_approximation));
  //def("under_approximation",ApprxLSPrltpFGridFunc(&Geometry::under_approximation));
  //def("under_approximation",ApprxLSZltpFGridFunc(&Geometry::under_approximation));
  def("under_approximation",ApprxGridltpFGridFunc(&Geometry::under_approximation));


}
