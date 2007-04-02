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

#include "python/python_float.h"


#include "linear_algebra/matrix.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"

#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
inline 
bool interiors_intersect(const ListSet< Zonotope<R> >& A, const Parallelotope<R>& B) 
{
  return interiors_intersect(A,ListSet< Zonotope<R> >(Zonotope<R>(B)));
}

template<class BS1, class BS2>
inline
ListSet<BS1> 
touching_intersection(const ListSet<BS1>& ls, const BS2& bs) 
{
  ListSet<BS1> output(ls.dimension());
    
  for (size_t i=0; i< ls.size(); i++) {
    if (!(disjoint(ls[i],bs))) {
      output.adjoin(ls[i]);
    }
  }

  return output;
}


template<class R>
inline
ListSet< Zonotope<R> >
over_approximate_interval_zonotope_list_set(const ListSet< Zonotope<Interval<R> > >& izls)
{
  ListSet< Zonotope<R> > result(izls.dimension());
  for(typename ListSet< Zonotope<Interval<R> > >::const_iterator iz_iter=izls.begin();
      iz_iter!=izls.end(); ++iz_iter)
  {
    result.adjoin(over_approximation(*iz_iter));
  }
  return result;
}

template<class R>
inline
ListSet< Zonotope<R> >
approximate_interval_zonotope_list_set(const ListSet< Zonotope<Interval<R> > > izls)
{
  ListSet< Zonotope<R> > result(izls.dimension());
  for(typename ListSet< Zonotope<Interval<R> > >::const_iterator iz_iter=izls.begin();
      iz_iter!=izls.end(); ++iz_iter)
  {
    result.adjoin(approximation(*iz_iter));
  }
  return result;
}


template<class R>
void export_list_set() 
{
  typedef Interval<R> I;
  
  typedef Rectangle<R> RRectangle;
  typedef Parallelotope<R> RParallelotope;
  typedef Zonotope<R> RZonotope;
  typedef Polytope<R> RPolytope;
  typedef Zonotope<I> IZonotope;

  typedef ListSet< Rectangle<R> > RRectangleListSet;
  typedef ListSet< Parallelotope<R> > RParallelotopeListSet;
  typedef ListSet< Zonotope<R> > RZonotopeListSet;
  typedef ListSet< Polytope<R> > RPolytopeListSet;
  typedef ListSet< Zonotope<I> > IZonotopeListSet;
  
  typedef GridCellListSet<R> RGridCellListSet;
  typedef GridMaskSet<R> RGridMaskSet;
  typedef PartitionTreeSet<R> RPartitionTreeSet;
  
  typedef Zonotope< Interval<R> > IZonotope;
  def("open_intersection",(RRectangleListSet(*)(const RRectangleListSet&,const RRectangleListSet&))(&open_intersection));
  def("disjoint",(tribool(*)(const RRectangleListSet&,const RRectangleListSet&))(&disjoint));
  def("disjoint",(tribool(*)(const RZonotopeListSet&,const RZonotopeListSet&))(&disjoint));
  def("subset", (tribool(*)(const RRectangleListSet&,const RRectangleListSet&))(&subset));

  def("touching_intersection",(RRectangleListSet(*)(const RRectangleListSet&,const RRectangle&))(&touching_intersection));
  def("touching_intersection",(RRectangleListSet(*)(const RRectangleListSet&,const RZonotope&))(&touching_intersection));
  def("touching_intersection",(RZonotopeListSet(*)(const RZonotopeListSet&,const RRectangle&))(&touching_intersection));
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

  class_<IZonotopeListSet>("IntervalZonotopeListSet",init<int>())
    .def(init<IZonotope>())
    .def(init<RRectangleListSet>())
    .def(init<RZonotopeListSet>())
    .def(init<IZonotopeListSet>())
    .def("dimension", &IZonotopeListSet::dimension)
    .def("push_back", &IZonotopeListSet::push_back)
    .def("adjoin", (void(IZonotopeListSet::*)(const IZonotope&))(&IZonotopeListSet::adjoin))
    .def("adjoin", (void(IZonotopeListSet::*)(const IZonotopeListSet&))(&IZonotopeListSet::adjoin))
    .def("size", &IZonotopeListSet::size)
    .def("empty", &IZonotopeListSet::empty)
    .def("__len__", &IZonotopeListSet::size)
//    .def("__getitem__", &RParallelotopeListSet::get, return_value_policy<copy_const_reference>())
    .def("__getitem__", &get_item<IZonotopeListSet>)
//    .def("__setitem__", &set_item<IZonotopeListSet>)
    .def("__iter__", iterator<IZonotopeListSet>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;

  def("over_approximation",&over_approximate_interval_zonotope_list_set<R>);
  def("approximation",&approximate_interval_zonotope_list_set<R>);
}

template void export_list_set<Float>();
