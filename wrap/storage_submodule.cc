/***************************************************************************
 *            grid_submodule.cc
 *
 *  Copyright 2008  Pieter Collins
 * 
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
 

#include "numeric.h"
#include "function_set.h"
#include "grid_set.h"

using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;




void export_grid() 
{
    typedef Vector<Float> RVector;
    typedef Vector<Interval> IVector;
  
    class_< Grid > grid_class("Grid",no_init);
    grid_class.def(init<uint>());
    grid_class.def(init<uint,Float>());
    grid_class.def(init< Vector<Float>, Vector<Float> >());
    grid_class.def("dimension", &Grid::dimension);
    grid_class.def("origin", &Grid::origin, return_value_policy<copy_const_reference>());
    grid_class.def("lengths", &Grid::lengths, return_value_policy<copy_const_reference>());
    grid_class.def(self_ns::str(self));
}


void export_grid_cell() 
{
    class_<GridCell> cell_class("GridCell",no_init);
    cell_class.def("dimension", &GridCell::dimension);
    cell_class.def("box", &GridCell::box, return_value_policy<copy_const_reference>());
    cell_class.def(self_ns::str(self));
}


void export_grid_tree_set() {
    class_<GridTreeSet> grid_tree_set_class("GridTreeSet",init<Grid>());
    grid_tree_set_class.def(init<uint>());
    grid_tree_set_class.def(init<GridTreeSet>());
    grid_tree_set_class.def("bounding_box", &GridTreeSet::bounding_box);
    grid_tree_set_class.def("empty", &GridTreeSet::empty);
    grid_tree_set_class.def("dimension", &GridTreeSet::dimension);
    grid_tree_set_class.def("clear", &GridTreeSet::clear);
    grid_tree_set_class.def("grid", &GridTreeSet::grid,return_value_policy<copy_const_reference>());
    grid_tree_set_class.def("adjoin", (void(GridTreeSet::*)(const GridCell&))(&GridTreeSet::adjoin));
    grid_tree_set_class.def("adjoin", (void(GridTreeSet::*)(const GridTreeSubset&))(&GridTreeSet::adjoin));
    grid_tree_set_class.def("restrict", (void(GridTreeSet::*)(const GridTreeSubset&))(&GridTreeSet::restrict));
    grid_tree_set_class.def("remove", (void(GridTreeSet::*)(const GridTreeSubset&))(&GridTreeSet::remove));
    grid_tree_set_class.def("adjoin_over_approximation", (void(GridTreeSet::*)(const Box&,const uint)) &GridTreeSet::adjoin_over_approximation);
    grid_tree_set_class.def("adjoin_outer_approximation", (void(GridTreeSet::*)(const CompactSetInterface&,const uint)) &GridTreeSet::adjoin_outer_approximation);
    grid_tree_set_class.def("adjoin_inner_approximation", (void(GridTreeSet::*)(const OpenSetInterface&,const uint,const uint)) &GridTreeSet::adjoin_inner_approximation);
    grid_tree_set_class.def("size", &GridTreeSet::size);
    grid_tree_set_class.def("__len__", &GridTreeSet::size);
    //grid_tree_set_class.def("__iter__", boost::python::iterator<GridTreeSet>());
    //grid_tree_set_class.def("__iter__", boost::python::const_iterator<GridTreeSet>());
    grid_tree_set_class.def(self_ns::str(self));


    //def("join",(GridTreeSet(*)(const GridTreeSet&,const GridTreeSet&))(&join));
    //def("difference",(GridTreeSet(*)(const GridTreeSet&,const GridTreeSet&))(&difference));
    //def("regular_intersection",(GridTreeSet(*)(const GridTreeSet&,const GridTreeSet&))(&regular_intersection));
    //def("overlap",(tribool(*)(const GridTreeSet&,const GridTreeSet&))(&overlap));
    //def("subset",(tribool(*)(const GridTreeSet&,const GridTreeSet&))(&subset));
    //def("subset",(tribool(*)(const CellList&,const GridTreeSet&))(&subset));

}


void storage_submodule()
{
    export_grid();
    export_grid_cell();
    export_grid_tree_set();
}

