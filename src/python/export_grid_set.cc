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

#include "python/python_float.h"

#include "linear_algebra/vector.h"

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"
#include "geometry/set_interface.h"


#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Combinatoric;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;



template<class R>
void export_grid() 
{
  typedef Vector<R> RVector;
  typedef Point<R> RPoint;
  typedef Rectangle<R> RRectangle;
  typedef ListSet< Rectangle<R> > RRectangleListSet;
  typedef SetInterface<R> RSetInterface;

  typedef Grid<R> RGrid;
  typedef FiniteGrid<R> RFiniteGrid;
  
  class_< RGrid > grid_class("Grid",init<uint,R>());
  grid_class.def(init<RPoint,RVector>());
  grid_class.def(init<RVector>());
  grid_class.def(init<RRectangle,LatticeBlock>());
  grid_class.def("dimension", &RGrid::dimension);
  grid_class.def("subdivision_coordinate", &RGrid::subdivision_coordinate);
  grid_class.def("subdivision_index", &RGrid::subdivision_index);
  grid_class.def("subdivision_lower_index", &RGrid::subdivision_lower_index);
  grid_class.def("subdivision_upper_index", &RGrid::subdivision_upper_index);
  grid_class.def(self_ns::str(self));


  class_< RFiniteGrid> finite_grid_class("FiniteGrid",init<const RGrid&,LatticeBlock>());
  finite_grid_class.def(init<RRectangle,size_type>());
  finite_grid_class.def(init<const RGrid&,RRectangle>());
  finite_grid_class.def("dimension", &RFiniteGrid::dimension);
  finite_grid_class.def("subdivision_coordinate", &RFiniteGrid::subdivision_coordinate);
  finite_grid_class.def("subdivision_lower_index", &RFiniteGrid::subdivision_lower_index);
  finite_grid_class.def("subdivision_upper_index", &RFiniteGrid::subdivision_upper_index);
  finite_grid_class.def("grid", &RFiniteGrid::grid,return_internal_reference<>());
  finite_grid_class.def("lattice_block", &RFiniteGrid::lattice_block,return_value_policy<copy_const_reference>());
  finite_grid_class.def("extent", &RFiniteGrid::extent);
  finite_grid_class.def(self_ns::str(self));
}

template<class R>
void export_grid_set() 
{
  typedef Interval<R> I;

  typedef Vector<R> RVector;

  typedef Grid<R> RGrid;
  typedef FiniteGrid<R> RFiniteGrid;
  
  typedef GridCell<R> RGridCell;
  typedef GridBlock<R> RGridBlock;
  typedef GridCellListSet<R> RGridCellListSet;
  typedef GridMaskSet<R> RGridMaskSet;
  
  typedef Point<R> RPoint;
  typedef Point<I> IPoint;

  typedef SetInterface<R> RSetInterface;
  typedef Rectangle<R> RRectangle;
  typedef Polyhedron<R> RPolyhedron;
  typedef Polytope<R> RPolytope;
  typedef Zonotope<R,R> RZonotope;
  typedef Zonotope<I,R> EZonotope;
  typedef Zonotope<I,I> IZonotope;
  typedef ListSet< Rectangle<R> > RRectangleListSet;
  typedef ListSet< Zonotope<R,R> > RZonotopeListSet;
  typedef ListSet< Zonotope<I,R> > EZonotopeListSet;
  typedef ListSet< Zonotope<I,I> > IZonotopeListSet;
  typedef PartitionTreeSet<R> RPartitionTreeSet;
  

  class_<RGridCell> grid_cell_class("GridCell",init<const RGrid&,LatticeCell>());
  grid_cell_class.def("dimension", &RGridCell::dimension);
  grid_cell_class.def("lattice_set", &RGridCell::lattice_set,return_value_policy<copy_const_reference>());
  grid_cell_class.def(self_ns::str(self));
  

  class_<RGridBlock> grid_block_class("GridBlock",init<const RGrid&,LatticeBlock>());
  grid_block_class.def(init<RGridCell>());
  grid_block_class.def("dimension", &RGridBlock::dimension);
  grid_block_class.def("lattice_set", &RGridBlock::lattice_set,return_value_policy<copy_const_reference>());
  grid_block_class.def(self_ns::str(self));


  
  class_<RGridCellListSet> grid_cell_list_set_class("GridCellListSet",init<const RGrid&>());
  grid_cell_list_set_class.def(init<RGridMaskSet>());
  grid_cell_list_set_class.def(init<RGridCellListSet>());
  grid_cell_list_set_class.def("lattice_set", &RGridCellListSet::lattice_set,return_value_policy<copy_const_reference>());
  grid_cell_list_set_class.def("dimension", &RGridCellListSet::dimension);
  grid_cell_list_set_class.def("adjoin", (void(RGridCellListSet::*)(const RGridCell&))(&RGridCellListSet::adjoin));
  grid_cell_list_set_class.def("adjoin", (void(RGridCellListSet::*)(const RGridBlock&))(&RGridCellListSet::adjoin));
  grid_cell_list_set_class.def("adjoin", (void(RGridCellListSet::*)(const RGridCellListSet&))(&RGridCellListSet::adjoin));
  grid_cell_list_set_class.def("adjoin_over_approximation", &RGridCellListSet::adjoin_over_approximation);
  grid_cell_list_set_class.def("adjoin_outer_approximation", (void(RGridCellListSet::*)(const RRectangle&))(&RGridCellListSet::adjoin_outer_approximation));
  grid_cell_list_set_class.def("adjoin_outer_approximation", (void(RGridCellListSet::*)(const RPolyhedron&))(&RGridCellListSet::adjoin_outer_approximation));
  grid_cell_list_set_class.def("adjoin_outer_approximation", (void(RGridCellListSet::*)(const RPolytope&))(&RGridCellListSet::adjoin_outer_approximation));
  grid_cell_list_set_class.def("adjoin_outer_approximation", (void(RGridCellListSet::*)(const RZonotope&))(&RGridCellListSet::adjoin_outer_approximation));
  grid_cell_list_set_class.def("adjoin_outer_approximation", (void(RGridCellListSet::*)(const EZonotope&))(&RGridCellListSet::adjoin_outer_approximation));
  grid_cell_list_set_class.def("adjoin_outer_approximation", (void(RGridCellListSet::*)(const IZonotope&))(&RGridCellListSet::adjoin_outer_approximation));
  grid_cell_list_set_class.def("adjoin_outer_approximation", (void(RGridCellListSet::*)(const RSetInterface&))(&RGridCellListSet::adjoin_outer_approximation));
  grid_cell_list_set_class.def("size", &RGridCellListSet::size);
  grid_cell_list_set_class.def("__len__", &RGridCellListSet::size);
  grid_cell_list_set_class.def("__getitem__", &get_item<RGridCellListSet>);
  grid_cell_list_set_class.def("__iter__", iterator<RGridCellListSet>());
  grid_cell_list_set_class.def(self_ns::str(self));

 
  class_<RGridMaskSet, bases<RSetInterface> > grid_mask_set_class("GridMaskSet",init<const RFiniteGrid&>());
  grid_mask_set_class.def(init<const RGrid&,LatticeBlock>());
  grid_mask_set_class.def(init<const RGrid&,RRectangle>());
  grid_mask_set_class.def(init<RGridMaskSet>());
  grid_mask_set_class.def("bounding_box", &RGridMaskSet::bounding_box);
  grid_mask_set_class.def("empty", &RGridMaskSet::empty);
  grid_mask_set_class.def("dimension", &RGridMaskSet::dimension);
  grid_mask_set_class.def("clear", &RGridMaskSet::clear);
  grid_mask_set_class.def("block", &RGridMaskSet::block,return_value_policy<copy_const_reference>());
  grid_mask_set_class.def("lattice_set", &RGridMaskSet::lattice_set,return_value_policy<copy_const_reference>());
  grid_mask_set_class.def("adjoin", (void(RGridMaskSet::*)(const RGridCell&))(&RGridMaskSet::adjoin));
  grid_mask_set_class.def("adjoin", (void(RGridMaskSet::*)(const RGridBlock&))(&RGridMaskSet::adjoin));
  grid_mask_set_class.def("adjoin", (void(RGridMaskSet::*)(const RGridCellListSet&))(&RGridMaskSet::adjoin));
  grid_mask_set_class.def("adjoin", (void(RGridMaskSet::*)(const RGridMaskSet&))(&RGridMaskSet::adjoin));
  grid_mask_set_class.def("restrict", (void(RGridMaskSet::*)(const RGridCellListSet&))(&RGridMaskSet::restrict));
  grid_mask_set_class.def("restrict", (void(RGridMaskSet::*)(const RGridMaskSet&))(&RGridMaskSet::restrict));
  grid_mask_set_class.def("adjoin_over_approximation", &RGridMaskSet::adjoin_over_approximation);
  grid_mask_set_class.def("adjoin_outer_approximation", (void(RGridMaskSet::*)(const RRectangle&))(&RGridMaskSet::adjoin_outer_approximation));
  grid_mask_set_class.def("adjoin_outer_approximation", (void(RGridMaskSet::*)(const RPolyhedron&))(&RGridMaskSet::adjoin_outer_approximation));
  grid_mask_set_class.def("adjoin_outer_approximation", (void(RGridMaskSet::*)(const RPolytope&))(&RGridMaskSet::adjoin_outer_approximation));
  grid_mask_set_class.def("adjoin_outer_approximation", (void(RGridMaskSet::*)(const RZonotope&))(&RGridMaskSet::adjoin_outer_approximation));
  grid_mask_set_class.def("adjoin_outer_approximation", (void(RGridMaskSet::*)(const EZonotope&))(&RGridMaskSet::adjoin_outer_approximation));
  grid_mask_set_class.def("adjoin_outer_approximation", (void(RGridMaskSet::*)(const IZonotope&))(&RGridMaskSet::adjoin_outer_approximation));
  grid_mask_set_class.def("adjoin_outer_approximation", (void(RGridMaskSet::*)(const RSetInterface&))(&RGridMaskSet::adjoin_outer_approximation));
  grid_mask_set_class.def("adjoin_over_approximation", &RGridMaskSet::adjoin_over_approximation);
  grid_mask_set_class.def("neighbourhood", &RGridMaskSet::neighbourhood);
  grid_mask_set_class.def("adjoining", &RGridMaskSet::adjoining);
  grid_mask_set_class.def("size", &RGridMaskSet::size);
  grid_mask_set_class.def("capacity", &RGridMaskSet::capacity);
  grid_mask_set_class.def("__len__", &RGridMaskSet::size);
  grid_mask_set_class.def("__getitem__", &get_item<RGridMaskSet>);
  grid_mask_set_class.def("__iter__", iterator<RGridMaskSet>());
  grid_mask_set_class.def(self_ns::str(self));


  def("join",(RGridMaskSet(*)(const RGridMaskSet&,const RGridMaskSet&))(&Geometry::join));
  def("difference",(RGridMaskSet(*)(const RGridMaskSet&,const RGridMaskSet&))(&Geometry::difference));
  def("regular_intersection",(RGridMaskSet(*)(const RGridMaskSet&,const RGridMaskSet&))(&Geometry::regular_intersection));
  def("overlap",(tribool(*)(const RGridMaskSet&,const RGridMaskSet&))(&Geometry::overlap));
  def("subset",(tribool(*)(const RGridMaskSet&,const RGridMaskSet&))(&Geometry::subset));
  def("subset",(tribool(*)(const RGridCellListSet&,const RGridMaskSet&))(&Geometry::subset));

  def("outer_approximation",(RGridBlock(*)(const IPoint&,const RGrid&))(&Geometry::outer_approximation));

  def("over_approximation",(RGridBlock(*)(const RRectangle&,const RGrid&))(&Geometry::over_approximation));
  def("under_approximation",(RGridBlock(*)(const RRectangle&,const RGrid&))(&Geometry::under_approximation));
  def("outer_approximation",(RGridBlock(*)(const RRectangle&,const RGrid&))(&Geometry::outer_approximation));
  def("inner_approximation",(RGridBlock(*)(const RRectangle&,const RGrid&))(&Geometry::inner_approximation));

  def("outer_approximation",(RGridCellListSet(*)(const RPolyhedron&,const RGrid&))(&Geometry::outer_approximation));
  def("outer_approximation",(RGridCellListSet(*)(const RPolytope&,const RGrid&))(&Geometry::outer_approximation));
  def("outer_approximation",(RGridCellListSet(*)(const RZonotope&,const RGrid&))(&Geometry::outer_approximation));
  def("outer_approximation",(RGridCellListSet(*)(const EZonotope&,const RGrid&))(&Geometry::outer_approximation));
  def("outer_approximation",(RGridCellListSet(*)(const IZonotope&,const RGrid&))(&Geometry::outer_approximation));
  def("outer_approximation",(RGridMaskSet(*)(const RRectangleListSet&,const RFiniteGrid&))(&Geometry::outer_approximation));
  def("outer_approximation",(RGridMaskSet(*)(const RZonotopeListSet&,const RFiniteGrid&))(&Geometry::outer_approximation));
  def("outer_approximation",(RGridMaskSet(*)(const EZonotopeListSet&,const RFiniteGrid&))(&Geometry::outer_approximation));
  def("outer_approximation",(RGridMaskSet(*)(const IZonotopeListSet&,const RFiniteGrid&))(&Geometry::outer_approximation));
  def("outer_approximation",(RGridMaskSet(*)(const RSetInterface&,const RFiniteGrid&))(&Geometry::outer_approximation));

  def("inner_approximation",(RGridCellListSet(*)(const RPolytope&,const RGrid&))(&Geometry::inner_approximation));
  def("inner_approximation",(RGridCellListSet(*)(const RPolyhedron&,const RGrid&))(&Geometry::inner_approximation));
  def("inner_approximation",(RGridCellListSet(*)(const RZonotope&,const RGrid&))(&Geometry::inner_approximation));
  def("inner_approximation",(RGridMaskSet(*)(const RSetInterface&,const RFiniteGrid&))(&Geometry::inner_approximation));

}

template void export_grid<Float>();
template void export_grid_set<Float>();
