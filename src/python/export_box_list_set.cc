/***************************************************************************
 *            python/export_box_list_set.cc
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


#include "geometry/box.h"
#include "geometry/box_list_set.h"
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



template<class R>
void export_box_list_set() 
{
  typedef Box<R> RBox;
  typedef BoxListSet<R> RBoxListSet;
  
  def("open_intersection",(RBoxListSet(*)(const RBoxListSet&,const RBoxListSet&))(&open_intersection));
  def("disjoint",(tribool(*)(const RBoxListSet&,const RBoxListSet&))(&disjoint));
  def("subset", (tribool(*)(const RBoxListSet&,const RBoxListSet&))(&subset));

  class_<RBoxListSet>("BoxListSet",init<>())
    .def(init<RBoxListSet>())
    .def("dimension", &RBoxListSet::dimension)
    .def("adjoin", (void(RBoxListSet::*)(const RBox&))&RBoxListSet::adjoin)
    .def("adjoin", (void(RBoxListSet::*)(const RBoxListSet&))&RBoxListSet::adjoin)
    .def("push_back", &RBoxListSet::push_back)
    .def("size", &RBoxListSet::size)
    .def("empty", &RBoxListSet::empty)
    .def("__len__", &RBoxListSet::size)
    .def("__getitem__", &__getitem__<RBoxListSet>)
    .def("__iter__", iterator<RBoxListSet>())
    .def(self_ns::str(self))
  ;
  

}

template void export_box_list_set<FloatPy>();
