/***************************************************************************
 *            python/export_list_set.cc
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *
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

#include "python/float.h"


#include "linear_algebra/matrix.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"

#include "python/utilities.h"
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

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
over_approximate_interval_zonotope_list_set(const ListSet< Zonotope<R,IntervalTag> >& izls)
{
  ListSet< Zonotope<R> > result(izls.dimension());
  Zonotope<R> z;
  for(typename ListSet< Zonotope<R,IntervalTag> >::const_iterator iz_iter=izls.begin();
      iz_iter!=izls.end(); ++iz_iter)
  {
    z=over_approximation(*iz_iter);
    result.adjoin(z);
  }
  return result;
}

template<class R>
inline
ListSet< Zonotope<R> >
approximate_interval_zonotope_list_set(const ListSet< Zonotope<R,IntervalTag> > izls)
{
  ListSet< Zonotope<R> > result(izls.dimension());
  Zonotope<R> z;
  for(typename ListSet< Zonotope<R,IntervalTag> >::const_iterator iz_iter=izls.begin();
      iz_iter!=izls.end(); ++iz_iter)
  {
    z=approximation(*iz_iter);
    result.adjoin(z);
  }
  return result;
}


template<class R>
void export_list_set() 
{
  typedef Interval<R> I;
  
  typedef Rectangle<R> RRectangle;
  typedef Polytope<R> RPolytope;
  typedef Zonotope<R,ExactTag> RZonotope;
  typedef Zonotope<R,UniformErrorTag> EZonotope;

  typedef ListSet<RRectangle> RRectangleListSet;
  typedef ListSet<RPolytope> RPolytopeListSet;
  typedef ListSet<RZonotope> RZonotopeListSet;
  typedef ListSet<EZonotope> EZonotopeListSet;

  typedef GridCellListSet<R> RGridCellListSet;
  typedef GridMaskSet<R> RGridMaskSet;
  typedef PartitionTreeSet<R> RPartitionTreeSet;
  
  def("open_intersection",(RRectangleListSet(*)(const RRectangleListSet&,const RRectangleListSet&))(&open_intersection));
  def("disjoint",(tribool(*)(const RRectangleListSet&,const RRectangleListSet&))(&disjoint));
  def("subset", (tribool(*)(const RRectangleListSet&,const RRectangleListSet&))(&subset));


  class_<RRectangleListSet>("RectangleListSet",init<int>())
    .def(init<RRectangle>())
    .def(init<RRectangleListSet>())
    .def(init<RGridCellListSet>())
    .def(init<RGridMaskSet>())
    .def(init<RPartitionTreeSet>())
    .def("dimension", &RRectangleListSet::dimension)
    .def("adjoin", (void(RRectangleListSet::*)(const RRectangle&))&RRectangleListSet::adjoin)
    .def("adjoin", (void(RRectangleListSet::*)(const RRectangleListSet&))&RRectangleListSet::adjoin)
    .def("push_back", &RRectangleListSet::push_back)
    .def("size", &RRectangleListSet::size)
    .def("empty", &RRectangleListSet::empty)
    .def("__len__", &RRectangleListSet::size)
    .def("__getitem__", &get_item<RRectangleListSet>)
//    .def("__setitem__", &set_item<RRectangleListSet>)
    .def("__iter__", iterator<RRectangleListSet>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;
  
  
  class_<RZonotopeListSet>("ZonotopeListSet",init<int>())
    .def(init<RZonotope>())
    .def(init<RRectangleListSet>())
    .def(init<RZonotopeListSet>())
    .def("dimension", &RZonotopeListSet::dimension)
    .def("push_back", &RZonotopeListSet::push_back)
    .def("adjoin", (void(RZonotopeListSet::*)(const RZonotope&))(&RZonotopeListSet::adjoin))
    .def("adjoin", (void(RZonotopeListSet::*)(const RZonotopeListSet&))(&RZonotopeListSet::adjoin))
    .def("size", &RZonotopeListSet::size)
    .def("empty", &RZonotopeListSet::empty)
    .def("__len__", &RZonotopeListSet::size)
    .def("__getitem__", &get_item<RZonotopeListSet>)
    .def("__iter__", iterator<RZonotopeListSet>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;

  class_<EZonotopeListSet>("ErrorZonotopeListSet",init<int>())
    .def(init<EZonotope>())
    .def(init<RRectangleListSet>())
    .def(init<RZonotopeListSet>())
    .def(init<EZonotopeListSet>())
    .def("dimension", &EZonotopeListSet::dimension)
    .def("push_back", &EZonotopeListSet::push_back)
    .def("adjoin", (void(EZonotopeListSet::*)(const EZonotope&))(&EZonotopeListSet::adjoin))
    .def("adjoin", (void(EZonotopeListSet::*)(const EZonotopeListSet&))(&EZonotopeListSet::adjoin))
    .def("size", &EZonotopeListSet::size)
    .def("empty", &EZonotopeListSet::empty)
    .def("__len__", &EZonotopeListSet::size)
    .def("__getitem__", &get_item<EZonotopeListSet>)
    .def("__iter__", iterator<EZonotopeListSet>())
    .def(self_ns::str(self))    // __self_ns::str__
  ;

}

template void export_list_set<FloatPy>();
