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


#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Combinatoric;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;


template<class PS, class S>
void
adjoin_over_approximation(PS& ps, const S& s) 
{
  ps.adjoin_over_approximation(s);
}


template<class R>
void export_grid() 
{
  typedef Vector<R> RVector;
  typedef Point<R> RPoint;
  typedef Rectangle<R> RRectangle;
  typedef ListSet< Rectangle<R> > RRectangleListSet;

  typedef Grid<R> RGrid;
  typedef FiniteGrid<R> RFiniteGrid;
  
  class_< RGrid >("Grid",init<uint,R>())
    .def(init<RPoint,RVector>())
    .def(init<RVector>())
    .def(init<RRectangle,LatticeBlock>())
    .def("dimension", &RGrid::dimension)
    .def("subdivision_coordinate", &RGrid::subdivision_coordinate)
    .def("subdivision_index", &RGrid::subdivision_index)
    .def("subdivision_lower_index", &RGrid::subdivision_lower_index)
    .def("subdivision_upper_index", &RGrid::subdivision_upper_index)
    .def(self_ns::str(self))    // __str__
    ;


  class_< RFiniteGrid>("FiniteGrid",init<const RGrid&,LatticeBlock>())
    .def(init<RRectangle,size_type>())
    .def(init<const RGrid&,RRectangle>())
    .def("dimension", &RFiniteGrid::dimension)
    .def("subdivision_coordinate", &RFiniteGrid::subdivision_coordinate)
    .def("subdivision_lower_index", &RFiniteGrid::subdivision_lower_index)
    .def("subdivision_upper_index", &RFiniteGrid::subdivision_upper_index)
    .def("grid", &RFiniteGrid::grid,return_internal_reference<>())
    .def("lattice_block", &RFiniteGrid::lattice_block,return_value_policy<copy_const_reference>())
    .def("extent", &RFiniteGrid::extent)
    .def(self_ns::str(self))    // __str__
    ;
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
  typedef Rectangle<R> RRectangle;
  typedef Parallelotope<R> RParallelotope;
  typedef Zonotope<R> RZonotope;
  typedef Polytope<R> RPolytope;
  typedef Zonotope<I> IZonotope;
  typedef ListSet< Rectangle<R> > RRectangleListSet;
  typedef ListSet< Parallelotope<R> > RParallelotopeListSet;
  typedef ListSet< Zonotope<R> > RZonotopeListSet;
  typedef ListSet< Zonotope<I> > IZonotopeListSet;
  typedef PartitionTreeSet<R> RPartitionTreeSet;
  

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
    .def("lattice_set", &RGridCellListSet::lattice_set,return_value_policy<copy_const_reference>())
    .def("dimension", &RGridCellListSet::dimension)
    .def("adjoin", (void(RGridCellListSet::*)(const RGridCell&))(&RGridCellListSet::adjoin))
    .def("adjoin", (void(RGridCellListSet::*)(const RGridBlock&))(&RGridCellListSet::adjoin))
    .def("adjoin", (void(RGridCellListSet::*)(const RGridCellListSet&))(&RGridCellListSet::adjoin))
    .def("adjoin_over_approximation", &adjoin_over_approximation<RGridCellListSet,RRectangle>)
    .def("adjoin_over_approximation", &adjoin_over_approximation<RGridCellListSet,RZonotope>)
    .def("adjoin_over_approximation", &adjoin_over_approximation<RGridCellListSet,IZonotope>)
    .def("adjoin_over_approximation", &adjoin_over_approximation<RGridCellListSet,RPolytope>)
    .def("size", &RGridCellListSet::size)
    .def("__len__", &RGridCellListSet::size)
    .def("__getitem__", &get_item<RGridCellListSet>)
    .def("__iter__", iterator<RGridCellListSet>())
    .def(self_ns::str(self))    // __str__
    ;
  
  class_<RGridMaskSet, bases<RSetInterface> >("GridMaskSet",init<const RFiniteGrid&>())
    .def(init<const RGrid&,LatticeBlock>())
    .def(init<const RGrid&,RRectangle>())
    .def(init<RGridMaskSet>())
    .def("bounding_box", &RGridMaskSet::bounding_box)
    .def("empty", &RGridMaskSet::empty)
    .def("dimension", &RGridMaskSet::dimension)
    .def("clear", &RGridMaskSet::clear)
    .def("block", &RGridMaskSet::block,return_value_policy<copy_const_reference>())
    .def("lattice_set", &RGridMaskSet::lattice_set,return_value_policy<copy_const_reference>())
    .def("adjoin", (void(RGridMaskSet::*)(const RGridCell&))(&RGridMaskSet::adjoin))
    .def("adjoin", (void(RGridMaskSet::*)(const RGridBlock&))(&RGridMaskSet::adjoin))
    .def("adjoin", (void(RGridMaskSet::*)(const RGridCellListSet&))(&RGridMaskSet::adjoin))
    .def("adjoin", (void(RGridMaskSet::*)(const RGridMaskSet&))(&RGridMaskSet::adjoin))
    .def("restrict", (void(RGridMaskSet::*)(const RGridCellListSet&))(&RGridMaskSet::restrict))
    .def("restrict", (void(RGridMaskSet::*)(const RGridMaskSet&))(&RGridMaskSet::restrict))
    .def("adjoin_over_approximation", &adjoin_over_approximation<RGridMaskSet,RRectangle>)
    .def("adjoin_over_approximation", &adjoin_over_approximation<RGridMaskSet,RZonotope>)
    .def("adjoin_over_approximation", &adjoin_over_approximation<RGridMaskSet,IZonotope>)
    .def("adjoin_over_approximation", &adjoin_over_approximation<RGridMaskSet,RPolytope>)
    //.def("adjoin_over_approximation", &adjoin_over_approximation<RGridMaskSet,RRectangleListSet>)
    //.def("adjoin_over_approximation", &adjoin_over_approximation<RGridMaskSet,RZonotopeListSet>)
    //.def("adjoin_over_approximation", &adjoin_over_approximation<RGridMaskSet,IZonotopeListSet>)
    .def("neighbourhood", &RGridMaskSet::neighbourhood)
    .def("adjoining", &RGridMaskSet::adjoining)
    .def("size", &RGridMaskSet::size)
    .def("capacity", &RGridMaskSet::capacity)
    .def("__len__", &RGridMaskSet::size)
    .def("__getitem__", &get_item<RGridMaskSet>)
    .def("__iter__", iterator<RGridMaskSet>())
    .def(self_ns::str(self))    // __str__
    ;


  def("join",(RGridMaskSet(*)(const RGridMaskSet&,const RGridMaskSet&))(&Geometry::join));
  def("difference",(RGridMaskSet(*)(const RGridMaskSet&,const RGridMaskSet&))(&Geometry::difference));
  def("regular_intersection",(RGridMaskSet(*)(const RGridMaskSet&,const RGridMaskSet&))(&Geometry::regular_intersection));
  def("overlap",(tribool(*)(const RGridMaskSet&,const RGridMaskSet&))(&Geometry::overlap));
  def("subset",(tribool(*)(const RGridMaskSet&,const RGridMaskSet&))(&Geometry::subset));
  def("subset",(tribool(*)(const RGridCellListSet&,const RGridMaskSet&))(&Geometry::subset));

  def("over_approximation",(RGridBlock(*)(const IPoint&,const RGrid&))(&Geometry::over_approximation));
  def("over_approximation",(RGridBlock(*)(const RRectangle&,const RGrid&))(&Geometry::over_approximation));
  def("over_approximation",(RGridCellListSet(*)(const RZonotope&,const RGrid&))(&Geometry::over_approximation));
  def("over_approximation",(RGridCellListSet(*)(const RPolytope&,const RGrid&))(&Geometry::over_approximation));
  def("over_approximation",(RGridCellListSet(*)(const IZonotope&,const RGrid&))(&Geometry::over_approximation));
  //def("over_approximation",(RGridCellListSet(*)(const RPolyhedron&,const RGrid&))(&Geometry::over_approximation));
  def("over_approximation",(RGridMaskSet(*)(const RRectangleListSet&,const RFiniteGrid&))(&Geometry::over_approximation));
  def("over_approximation",(RGridMaskSet(*)(const RParallelotopeListSet&,const RFiniteGrid&))(&Geometry::over_approximation));
  def("over_approximation",(RGridMaskSet(*)(const RZonotopeListSet&,const RFiniteGrid&))(&Geometry::over_approximation));
  def("over_approximation",(RGridMaskSet(*)(const RGridMaskSet&,const RFiniteGrid&))(&Geometry::over_approximation));

  def("under_approximation",(RGridBlock(*)(const RRectangle&,const RGrid&))(&Geometry::under_approximation));
  def("under_approximation",(RGridCellListSet(*)(const RZonotope&,const RGrid&))(&Geometry::under_approximation));
  def("under_approximation",(RGridCellListSet(*)(const RPolytope&,const RGrid&))(&Geometry::under_approximation));
  //def("under_approximation",(RGridCellListSet(*)(const RPolyhedron&,const RGrid&))(&Geometry::under_approximation));
  //def("under_approximation",(RGridMaskSet(*)(const RRectangleListSet&,const RFiniteGrid&))(&Geometry::under_approximation));
  //def("under_approximation",(RGridMaskSet(*)(const RParallelotopeListSet&,const RFiniteGrid&))(&Geometry::under_approximation));
  //def("under_approximation",(RGridMaskSet(*)(const RZonotopeListSet&,const RFiniteGrid&))(&Geometry::under_approximation));
  def("under_approximation",(RGridMaskSet(*)(const RGridMaskSet&,const RFiniteGrid&))(&Geometry::under_approximation));


}

template void export_grid<Float>();
template void export_grid_set<Float>();
