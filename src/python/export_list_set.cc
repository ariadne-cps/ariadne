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



#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polyhedron.h"
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

inline 
bool interiors_intersect(const RZonotopeListSet& A, const RParallelotope& B) 
{
  return interiors_intersect(A,RZonotopeListSet(RZonotope(B)));
}

template <typename R, template<typename> class BS, template<typename> class BS2>
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

template ListSet<Real,Rectangle> touching_intersection(const ListSet<Real,Rectangle>&, const Rectangle<Real> &);
template ListSet<Real,Parallelotope> touching_intersection(const ListSet<Real,Parallelotope>&, const Rectangle<Real> &);
template ListSet<Real,Zonotope> touching_intersection(const ListSet<Real,Zonotope>&, const Rectangle<Real> &);
template ListSet<Real,Rectangle> touching_intersection(const ListSet<Real,Rectangle>&, const Parallelotope<Real> &);
template ListSet<Real,Parallelotope> touching_intersection(const ListSet<Real,Parallelotope>&, const Parallelotope<Real> &);
template ListSet<Real,Zonotope> touching_intersection(const ListSet<Real,Zonotope>&, const Parallelotope<Real> &);
template ListSet<Real,Rectangle> touching_intersection(const ListSet<Real,Rectangle>&, const Zonotope<Real> &);
template ListSet<Real,Parallelotope> touching_intersection(const ListSet<Real,Parallelotope>&, const Zonotope<Real> &);
template ListSet<Real,Zonotope> touching_intersection(const ListSet<Real,Zonotope>&, const Zonotope<Real> &);


void export_list_set() {
  def("regular_intersection",(RRectangleListSet(*)(const RRectangleListSet&,const RRectangleListSet&))(&regular_intersection));
  def("interiors_intersect",(bool(*)(const RRectangleListSet&,const RRectangleListSet&))(&interiors_intersect));
  def("interiors_intersect",(bool(*)(const RParallelotopeListSet&,const RParallelotopeListSet&))(&interiors_intersect));
  def("interiors_intersect",(bool(*)(const RZonotopeListSet&,const RZonotopeListSet&))(&interiors_intersect));

  def("touching_intersection",(RRectangleListSet(*)(const RRectangleListSet&,const RRectangle&))(&touching_intersection));
  def("touching_intersection",(RRectangleListSet(*)(const RRectangleListSet&,const RParallelotope&))(&touching_intersection));
  def("touching_intersection",(RRectangleListSet(*)(const RRectangleListSet&,const RZonotope&))(&touching_intersection));
  def("touching_intersection",(RParallelotopeListSet(*)(const RParallelotopeListSet&,const RRectangle&))(&touching_intersection));
  def("touching_intersection",(RParallelotopeListSet(*)(const RParallelotopeListSet&,const RParallelotope&))(&touching_intersection));
  def("touching_intersection",(RParallelotopeListSet(*)(const RParallelotopeListSet&,const RZonotope&))(&touching_intersection));
  def("touching_intersection",(RZonotopeListSet(*)(const RZonotopeListSet&,const RRectangle&))(&touching_intersection));
  def("touching_intersection",(RZonotopeListSet(*)(const RZonotopeListSet&,const RParallelotope&))(&touching_intersection));
  def("touching_intersection",(RZonotopeListSet(*)(const RZonotopeListSet&,const RZonotope&))(&touching_intersection));

  def("disjoint", (bool(*)(const RRectangleListSet&,const RRectangleListSet&))(&disjoint));
  def("inner_subset", (bool(*)(const RRectangleListSet&,const RRectangleListSet&))(&inner_subset));
  def("subset", (bool(*)(const RRectangleListSet&,const RRectangleListSet&))(&subset));

  class_<RRectangleListSet>("RectangleListSet",init<int>())
    .def(init<RRectangle>())
    .def(init<RRectangleListSet>())
    .def(init<RGridBlockListSet>())
    .def(init<RGridCellListSet>())
    .def(init<RGridMaskSet>())
    .def(init<RPartitionTreeSet>())
    .def("dimension", &RRectangleListSet::dimension)
    .def("push_back", &RRectangleListSet::push_back)
    .def("size", &RRectangleListSet::size)
    .def("empty", &RRectangleListSet::empty)
    .def("__len__", &RRectangleListSet::size)
    .def("__getitem__", &RRectangleListSet::get, return_value_policy<copy_const_reference>())
    .def("__setitem__", &RRectangleListSet::set)
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
//    .def("__getitem__", &get_item<RParallelotopeListSet>)
    .def("__getitem__", &plsg)
    .def("__setitem__", &RParallelotopeListSet::set)
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
//    .def("__getitem__", &get_item<RParallelotopeListSet>)
    .def("__getitem__", &zlsg)
    .def("__setitem__", &RZonotopeListSet::set)
    .def("__iter__", iterator<RZonotopeListSet>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;

}
