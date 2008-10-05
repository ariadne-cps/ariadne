#include "numeric.h"
#include "paving.h"

using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;



void export_paving() 
{
  typedef Vector<Float> RVector;
  typedef Vector<Interval> IVector;
  
  class_< Grid > grid_class("Grid",init<uint,Float>());
  grid_class.def(init<RVector,RVector>());
  grid_class.def("dimension", &Grid::dimension);
  grid_class.def("origin", &Grid::origin);
  grid_class.def("lengths", &Grid::lengths);
  grid_class.def(self_ns::str(self));


  class_<Cell> cell_class("Cell",no_init);
  cell_class.def("dimension", &GridCell::dimension);
  cell_class.def("box", &GridCell::box);
  cell_class.def(self_ns::str(self));
  

  class_<RGridBlock> grid_block_class("GridBlock",init<const RGrid&,LatticeBlock>());
  grid_block_class.def(init<RGridCell>());
  grid_block_class.def("dimension", &RGridBlock::dimension);
  grid_block_class.def("lattice_set", &RGridBlock::lattice_set,return_value_policy<copy_const_reference>());
  grid_block_class.def(self_ns::str(self));


  
  class_<RGridCellListSet> grid_cell_list_set_class("GridCellListSet",init<const RGrid&>());
  grid_cell_list_set_class.def(init<RGridMaskSet>());
  grid_cell_list_set_class.def(init<RGridCellListSet>());
  grid_cell_list_set_class.def("lattice_set", &RGridCellListSet::lattice_set,return_value_policy<copy_const_reference>());
  grid_cell_list_set_class.def("clear", &RGridCellListSet::clear);
  grid_cell_list_set_class.def("dimension", &RGridCellListSet::dimension);
  grid_cell_list_set_class.def("adjoin", (void(RGridCellListSet::*)(const RGridCell&))(&RGridCellListSet::adjoin));
  grid_cell_list_set_class.def("adjoin", (void(RGridCellListSet::*)(const RGridBlock&))(&RGridCellListSet::adjoin));
  grid_cell_list_set_class.def("adjoin", (void(RGridCellListSet::*)(const RGridCellListSet&))(&RGridCellListSet::adjoin));
  grid_cell_list_set_class.def("adjoin_over_approximation", &RGridCellListSet::adjoin_over_approximation);
  grid_cell_list_set_class.def("adjoin_outer_approximation", (void(RGridCellListSet::*)(const RBox&))(&RGridCellListSet::adjoin_outer_approximation));
  grid_cell_list_set_class.def("adjoin_outer_approximation", (void(RGridCellListSet::*)(const RPolyhedron&))(&RGridCellListSet::adjoin_outer_approximation));
  grid_cell_list_set_class.def("adjoin_outer_approximation", (void(RGridCellListSet::*)(const RPolytope&))(&RGridCellListSet::adjoin_outer_approximation));
  grid_cell_list_set_class.def("adjoin_outer_approximation", (void(RGridCellListSet::*)(const RZonotope&))(&RGridCellListSet::adjoin_outer_approximation));
  grid_cell_list_set_class.def("adjoin_outer_approximation", (void(RGridCellListSet::*)(const RSetInterface&))(&RGridCellListSet::adjoin_outer_approximation));
  grid_cell_list_set_class.def("size", &RGridCellListSet::size);
  grid_cell_list_set_class.def("__len__", &RGridCellListSet::size);
  grid_cell_list_set_class.def("__getitem__", &__getitem__<RGridCellListSet>);
  grid_cell_list_set_class.def("__iter__", iterator<RGridCellListSet>());
  grid_cell_list_set_class.def("unique_sort", &RGridCellListSet::unique_sort);
  grid_cell_list_set_class.def("pop", &RGridCellListSet::pop);
  grid_cell_list_set_class.def(self_ns::str(self));

 
  class_<RGridMaskSet, bases<RSetInterface> > grid_mask_set_class("GridMaskSet",init<const RFiniteGrid&>());
  grid_mask_set_class.def(init<const RGrid&,LatticeBlock>());
  grid_mask_set_class.def(init<const RGrid&,RBox>());
  grid_mask_set_class.def(init<RGridMaskSet>());
  grid_mask_set_class.def("bounding_box", &RGridMaskSet::bounding_box);
  grid_mask_set_class.def("empty", &RGridMaskSet::empty);
  grid_mask_set_class.def("dimension", &RGridMaskSet::dimension);
  grid_mask_set_class.def("clear", &RGridMaskSet::clear);
  grid_mask_set_class.def("grid", &RGridMaskSet::grid,return_value_policy<copy_const_reference>());
  grid_mask_set_class.def("block", &RGridMaskSet::block,return_value_policy<copy_const_reference>());
  grid_mask_set_class.def("lattice_set", &RGridMaskSet::lattice_set,return_value_policy<copy_const_reference>());
  grid_mask_set_class.def("adjoin", (void(RGridMaskSet::*)(const RGridCell&))(&RGridMaskSet::adjoin));
  grid_mask_set_class.def("adjoin", (void(RGridMaskSet::*)(const RGridBlock&))(&RGridMaskSet::adjoin));
  grid_mask_set_class.def("adjoin", (void(RGridMaskSet::*)(const RGridCellListSet&))(&RGridMaskSet::adjoin));
  grid_mask_set_class.def("adjoin", (void(RGridMaskSet::*)(const RGridMaskSet&))(&RGridMaskSet::adjoin));
  grid_mask_set_class.def("restrict", (void(RGridMaskSet::*)(const RGridCellListSet&))(&RGridMaskSet::restrict));
  grid_mask_set_class.def("restrict", (void(RGridMaskSet::*)(const RGridMaskSet&))(&RGridMaskSet::restrict));
  grid_mask_set_class.def("remove", (void(RGridMaskSet::*)(const RGridMaskSet&))(&RGridMaskSet::remove));
  grid_mask_set_class.def("adjoin_over_approximation", &RGridMaskSet::adjoin_over_approximation);
  grid_mask_set_class.def("adjoin_outer_approximation", (void(RGridMaskSet::*)(const RBox&))(&RGridMaskSet::adjoin_outer_approximation));
  grid_mask_set_class.def("adjoin_outer_approximation", (void(RGridMaskSet::*)(const RPolyhedron&))(&RGridMaskSet::adjoin_outer_approximation));
  grid_mask_set_class.def("adjoin_outer_approximation", (void(RGridMaskSet::*)(const RPolytope&))(&RGridMaskSet::adjoin_outer_approximation));
  grid_mask_set_class.def("adjoin_outer_approximation", (void(RGridMaskSet::*)(const RZonotope&))(&RGridMaskSet::adjoin_outer_approximation));
  grid_mask_set_class.def("adjoin_outer_approximation", (void(RGridMaskSet::*)(const RSetInterface&))(&RGridMaskSet::adjoin_outer_approximation));
  grid_mask_set_class.def("adjoin_over_approximation", &RGridMaskSet::adjoin_over_approximation);
  grid_mask_set_class.def("neighbourhood", &RGridMaskSet::neighbourhood);
  grid_mask_set_class.def("adjoining", &RGridMaskSet::adjoining);
  grid_mask_set_class.def("size", &RGridMaskSet::size);
  grid_mask_set_class.def("capacity", &RGridMaskSet::capacity);
  grid_mask_set_class.def("__len__", &RGridMaskSet::size);
  grid_mask_set_class.def("__getitem__", &__getitem__<RGridMaskSet>);
  grid_mask_set_class.def("__iter__", iterator<RGridMaskSet>());
  grid_mask_set_class.def(self_ns::str(self));


  def("join",(RGridMaskSet(*)(const RGridMaskSet&,const RGridMaskSet&))(&join));
  def("difference",(RGridMaskSet(*)(const RGridMaskSet&,const RGridMaskSet&))(&difference));
  def("regular_intersection",(RGridMaskSet(*)(const RGridMaskSet&,const RGridMaskSet&))(&regular_intersection));
  def("overlap",(tribool(*)(const RGridMaskSet&,const RGridMaskSet&))(&overlap));
  def("subset",(tribool(*)(const RGridMaskSet&,const RGridMaskSet&))(&subset));
  def("subset",(tribool(*)(const RGridCellListSet&,const RGridMaskSet&))(&subset));

  def("outer_approximation",(RGridBlock(*)(const IPoint&,const RGrid&))(&outer_approximation));

  def("over_approximation",(RGridBlock(*)(const RBox&,const RGrid&))(&over_approximation));
  def("under_approximation",(RGridBlock(*)(const RBox&,const RGrid&))(&under_approximation));
  def("outer_approximation",(RGridBlock(*)(const RBox&,const RGrid&))(&outer_approximation));
  def("inner_approximation",(RGridBlock(*)(const RBox&,const RGrid&))(&inner_approximation));


  def("outer_approximation",(RGridCellListSet(*)(const RPolyhedron&,const RGrid&))(&outer_approximation));
  def("outer_approximation",(RGridCellListSet(*)(const RPolytope&,const RGrid&))(&outer_approximation));
  def("outer_approximation",(RGridCellListSet(*)(const RZonotope&,const RGrid&))(&outer_approximation));
  def("outer_approximation",(RGridCellListSet(*)(const RSetInterface&,const RGrid&))(&outer_approximation));
  def("outer_approximation",(RGridMaskSet(*)(const RZonotopeListSet&,const RFiniteGrid&))(&outer_approximation));
  def("outer_approximation",(RGridMaskSet(*)(const RSetInterface&,const RFiniteGrid&))(&outer_approximation));

  def("inner_approximation",(RGridCellListSet(*)(const RPolytope&,const RGrid&))(&inner_approximation));
  def("inner_approximation",(RGridCellListSet(*)(const RPolyhedron&,const RGrid&))(&inner_approximation));
  def("inner_approximation",(RGridCellListSet(*)(const RZonotope&,const RGrid&))(&inner_approximation));
  def("inner_approximation",(RGridMaskSet(*)(const RSetInterface&,const RFiniteGrid&))(&inner_approximation));

}

template void export_grid<FloatPy>();
template void export_grid_set<FloatPy>();
