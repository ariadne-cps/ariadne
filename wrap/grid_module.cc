#include "numeric.h"
#include "set.h"
#include "grid.h"

using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;



void export_grid() 
{
  typedef Vector<Float> RVector;
  typedef Vector<Interval> IVector;
  
  class Zonotope;

  class_< Grid > grid_class("Grid",init<uint>());
  grid_class.def(init<uint,Float>());
  grid_class.def(init<RVector,RVector>());
  grid_class.def("dimension", &Grid::dimension);
  grid_class.def("origin", &Grid::origin);
  grid_class.def("lengths", &Grid::lengths);
  grid_class.def(self_ns::str(self));


  class_<Cell> cell_class("Cell",no_init);
  cell_class.def("dimension", &Cell::dimension);
  cell_class.def("box", &Cell::box);
  cell_class.def(self_ns::str(self));
  

  class_<CellList> cell_list_class("CellList",init<const Grid&>());
  cell_list_class.def(init<CellList>());
  cell_list_class.def("clear", &CellList::clear);
  cell_list_class.def("dimension", &CellList::dimension);
  cell_list_class.def("size", &CellList::size);
  cell_list_class.def("__len__", &CellList::size);
  //cell_list_class.def("__getitem__", &__getitem__<CellList>);
  //cell_list_class.def("__iter__", iterator<CellList>());
  cell_list_class.def("pop", &CellList::pop);
  cell_list_class.def("unique_sort", &CellList::unique_sort);
  cell_list_class.def("adjoin", (void(CellList::*)(const Cell&))(&CellList::adjoin));
  cell_list_class.def("adjoin", (void(CellList::*)(const CellList&))(&CellList::adjoin));
  cell_list_class.def("adjoin_over_approximation", &CellList::adjoin_over_approximation);
  cell_list_class.def("adjoin_outer_approximation", &CellList::adjoin_outer_approximation);
  cell_list_class.def(self_ns::str(self));

 
  class_<GridTreeSet> grid_tree_set_class("GridTreeSet",init<Grid>());
  grid_tree_set_class.def(init<uint>());
  grid_tree_set_class.def(init<GridTreeSet>());
  grid_tree_set_class.def("bounding_box", &GridTreeSet::bounding_box);
  grid_tree_set_class.def("empty", &GridTreeSet::empty);
  grid_tree_set_class.def("dimension", &GridTreeSet::dimension);
  grid_tree_set_class.def("clear", &GridTreeSet::clear);
  grid_tree_set_class.def("grid", &GridTreeSet::grid,return_value_policy<copy_const_reference>());
  grid_tree_set_class.def("adjoin", (void(GridTreeSet::*)(const Cell&))(&GridTreeSet::adjoin));
  grid_tree_set_class.def("adjoin", (void(GridTreeSet::*)(const CellList&))(&GridTreeSet::adjoin));
  grid_tree_set_class.def("adjoin", (void(GridTreeSet::*)(const GridTreeSet&))(&GridTreeSet::adjoin));
  grid_tree_set_class.def("restrict", (void(GridTreeSet::*)(const CellList&))(&GridTreeSet::restrict));
  grid_tree_set_class.def("restrict", (void(GridTreeSet::*)(const GridTreeSet&))(&GridTreeSet::restrict));
  grid_tree_set_class.def("remove", (void(GridTreeSet::*)(const GridTreeSet&))(&GridTreeSet::remove));
  grid_tree_set_class.def("adjoin_over_approximation", &GridTreeSet::adjoin_over_approximation);
  grid_tree_set_class.def("adjoin_outer_approximation", &GridTreeSet::adjoin_outer_approximation);
  grid_tree_set_class.def("adjoin_inner_approximation", &GridTreeSet::adjoin_inner_approximation);
  grid_tree_set_class.def("neighbourhood", &GridTreeSet::neighbourhood);
  grid_tree_set_class.def("adjoining", &GridTreeSet::adjoining);
  grid_tree_set_class.def("size", &GridTreeSet::size);
  grid_tree_set_class.def("capacity", &GridTreeSet::capacity);
  grid_tree_set_class.def("__len__", &GridTreeSet::size);
  //grid_tree_set_class.def("__iter__", iterator<GridTreeSet>());
  grid_tree_set_class.def(self_ns::str(self));


  //def("join",(GridTreeSet(*)(const GridTreeSet&,const GridTreeSet&))(&join));
  //def("difference",(GridTreeSet(*)(const GridTreeSet&,const GridTreeSet&))(&difference));
  //def("regular_intersection",(GridTreeSet(*)(const GridTreeSet&,const GridTreeSet&))(&regular_intersection));
  //def("overlap",(tribool(*)(const GridTreeSet&,const GridTreeSet&))(&overlap));
  //def("subset",(tribool(*)(const GridTreeSet&,const GridTreeSet&))(&subset));
  //def("subset",(tribool(*)(const CellList&,const GridTreeSet&))(&subset));

}


BOOST_PYTHON_MODULE(grid)
{
 export_grid();
}
