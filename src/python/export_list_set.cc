/***************************************************************************
 *            python/export_list_set.cc
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
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "real_typedef.h"


#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"

#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
inline 
bool interiors_intersect(const ListSet<R,Zonotope>& A, const Parallelotope<R>& B) 
{
  return interiors_intersect(A,ListSet<R,Zonotope>(Zonotope<R>(B)));
}

template<class R, template<class> class BS, template<class> class BS2>
inline
ListSet<R,BS> 
touching_intersection(const ListSet<R,BS>& ls, const BS2<R>& bs) 
{
  ListSet<R,BS> output(ls.dimension());
    
  for (size_t i=0; i< ls.size(); i++) {
    if (!(disjoint(ls[i],bs)))
      output.adjoin(ls[i]);
  }

  return output;
}



template<class R>
void export_list_set() 
{
  typedef Rectangle<R> RRectangle;
  typedef Parallelotope<R> RParallelotope;
  typedef Zonotope<R> RZonotope;
  typedef Polytope<R> RPolytope;
  
  typedef ListSet<R,Rectangle> RRectangleListSet;
  typedef ListSet<R,Parallelotope> RParallelotopeListSet;
  typedef ListSet<R,Zonotope> RZonotopeListSet;
  typedef ListSet<R,Polytope> RPolytopeListSet;
  
  typedef GridCellListSet<R> RGridCellListSet;
  typedef GridMaskSet<R> RGridMaskSet;
  typedef PartitionTreeSet<R> RPartitionTreeSet;
  
  def("open_intersection",(RRectangleListSet(*)(const RRectangleListSet&,const RRectangleListSet&))(&open_intersection));
  def("disjoint",(tribool(*)(const RRectangleListSet&,const RRectangleListSet&))(&disjoint));
  def("disjoint",(tribool(*)(const RParallelotopeListSet&,const RParallelotopeListSet&))(&disjoint));
  def("disjoint",(tribool(*)(const RZonotopeListSet&,const RZonotopeListSet&))(&disjoint));
  def("subset", (tribool(*)(const RRectangleListSet&,const RRectangleListSet&))(&subset));

  def("touching_intersection",(RRectangleListSet(*)(const RRectangleListSet&,const RRectangle&))(&touching_intersection));
  def("touching_intersection",(RRectangleListSet(*)(const RRectangleListSet&,const RParallelotope&))(&touching_intersection));
  def("touching_intersection",(RRectangleListSet(*)(const RRectangleListSet&,const RZonotope&))(&touching_intersection));
  def("touching_intersection",(RParallelotopeListSet(*)(const RParallelotopeListSet&,const RRectangle&))(&touching_intersection));
  def("touching_intersection",(RParallelotopeListSet(*)(const RParallelotopeListSet&,const RParallelotope&))(&touching_intersection));
  def("touching_intersection",(RParallelotopeListSet(*)(const RParallelotopeListSet&,const RZonotope&))(&touching_intersection));
  def("touching_intersection",(RZonotopeListSet(*)(const RZonotopeListSet&,const RRectangle&))(&touching_intersection));
  def("touching_intersection",(RZonotopeListSet(*)(const RZonotopeListSet&,const RParallelotope&))(&touching_intersection));
  def("touching_intersection",(RZonotopeListSet(*)(const RZonotopeListSet&,const RZonotope&))(&touching_intersection));


  class_<RRectangleListSet>("RectangleListSet",init<int>())
    .def(init<RRectangle>())
    .def(init<RRectangleListSet>())
    .def(init<RGridCellListSet>())
    .def(init<RGridMaskSet>())
    .def(init<RPartitionTreeSet>())
    .def("dimension", &RRectangleListSet::dimension)
    .def("push_back", &RRectangleListSet::push_back)
    .def("size", &RRectangleListSet::size)
    .def("empty", &RRectangleListSet::empty)
    .def("__len__", &RRectangleListSet::size)
    .def("__getitem__", &get_item<RRectangleListSet>)
//    .def("__setitem__", &set_item<RRectangleListSet>)
    .def("__iter__", iterator<RRectangleListSet>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;
  
  class_<RParallelotopeListSet>("ParallelotopeListSet",init<int>())
    .def(init<RParallelotope>())
    .def(init<RRectangleListSet>())
    .def(init<RParallelotopeListSet>())
    .def("dimension", &RParallelotopeListSet::dimension)
    .def("push_back", &RParallelotopeListSet::push_back)
    .def("adjoin", (void(RParallelotopeListSet::*)(const RParallelotope&))(&RParallelotopeListSet::adjoin))
    .def("adjoin", (void(RParallelotopeListSet::*)(const RParallelotopeListSet&))(&RParallelotopeListSet::adjoin))
    .def("size", &RParallelotopeListSet::size)
    .def("empty", &RParallelotopeListSet::empty)
    .def("__len__", &RParallelotopeListSet::size)
//    .def("__getitem__", &RParallelotopeListSet::get, return_value_policy<copy_const_reference>())
    .def("__getitem__", &get_item<RParallelotopeListSet>)
//    .def("__setitem__", &set_item<RParallelotopeListSet>)
    .def("__iter__", iterator<RParallelotopeListSet>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;
  
  class_<RZonotopeListSet>("ZonotopeListSet",init<int>())
    .def(init<RZonotope>())
    .def(init<RRectangleListSet>())
    .def(init<RParallelotopeListSet>())
    .def(init<RZonotopeListSet>())
    .def("dimension", &RZonotopeListSet::dimension)
    .def("push_back", &RZonotopeListSet::push_back)
    .def("adjoin", (void(RZonotopeListSet::*)(const RZonotope&))(&RZonotopeListSet::adjoin))
    .def("adjoin", (void(RZonotopeListSet::*)(const RZonotopeListSet&))(&RZonotopeListSet::adjoin))
    .def("size", &RZonotopeListSet::size)
    .def("empty", &RZonotopeListSet::empty)
    .def("__len__", &RZonotopeListSet::size)
//    .def("__getitem__", &RParallelotopeListSet::get, return_value_policy<copy_const_reference>())
    .def("__getitem__", &get_item<RZonotopeListSet>)
//    .def("__setitem__", &set_item<RZonotopeListSet>)
    .def("__iter__", iterator<RZonotopeListSet>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;

}

template void export_list_set<Real>();
