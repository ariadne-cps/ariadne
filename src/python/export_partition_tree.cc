/***************************************************************************
 *            python/export_partition_tree.cc
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

#include "base/numerical_type.h"
#include "base/binary_word.h"
#include "base/binary_tree.h"

#include "geometry/rectangle.h"
#include "geometry/partition_tree_set.h"

#include <boost/python.hpp>

#include "python/real_typedef.h"

using Ariadne::BooleanArray;
using Ariadne::BinaryWord;
using Ariadne::BinaryTree;
using Ariadne::SubdivisionSequence;

typedef Ariadne::Geometry::Point<Real> RPoint;

typedef Ariadne::Interval<Real> RInterval;
typedef Ariadne::Geometry::Rectangle<Real> RRectangle;

typedef Ariadne::Geometry::ListSet<Real,Ariadne::Geometry::Rectangle> RRectangleListSet;

typedef Ariadne::Geometry::PartitionScheme<Real> RPartitionScheme;
typedef Ariadne::Geometry::PartitionTree<Real> RPartitionTree;
typedef Ariadne::Geometry::PartitionTreeCell<Real> RPartitionTreeCell;
typedef Ariadne::Geometry::PartitionTreeSet<Real> RPartitionTreeSet;

using Ariadne::Geometry::regular_intersection;
using Ariadne::Geometry::interiors_intersect;
using Ariadne::Geometry::disjoint;
using Ariadne::Geometry::inner_subset;
using Ariadne::Geometry::subset;

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;
using boost::python::def;
using boost::python::self_ns::str;

void export_partition_tree() {
  class_<RPartitionScheme>("PartitionScheme",init<RRectangle,SubdivisionSequence>())
    .def(init<RRectangle>())
    .def("dimension", &RPartitionScheme::dimension)
    .def("bounding_box", &RPartitionScheme::bounding_box, return_value_policy<copy_const_reference>())
    .def("subdivision_coordinates", &RPartitionScheme::subdivision_coordinates, return_value_policy<copy_const_reference>())
    .def("subdivision_coordinate", &RPartitionScheme::subdivision_coordinate)
    .def(str(self))    // __str__
  ;

  class_<RPartitionTree>("PartitionTree",init<RRectangle,SubdivisionSequence,BinaryTree>())
    .def(init<RRectangle>())
    .def(init<RRectangle,BinaryTree>())
    .def(init<RPartitionScheme,BinaryTree>())
    .def("dimension", &RPartitionTree::dimension)
    .def("bounding_box", &RPartitionTree::bounding_box, return_value_policy<copy_const_reference>())
    .def("subdivision_coordinates", &RPartitionTree::subdivision_coordinates, return_value_policy<copy_const_reference>())
    .def("tree", &RPartitionTree::tree, return_value_policy<copy_const_reference>())
    .def(str(self))    // __str__
  ;

  class_<RPartitionTreeSet>("PartitionTreeSet",init<RRectangle,SubdivisionSequence,BinaryTree,BooleanArray>())
    .def(init<RPartitionScheme,BinaryTree,BooleanArray>())
    .def(init<RPartitionTree,BooleanArray>())
    .def("dimension", &RPartitionTreeSet::dimension)
    .def("bounding_box", &RPartitionTreeSet::bounding_box, return_value_policy<copy_const_reference>())
    .def("subdivision_coordinates", &RPartitionTreeSet::subdivision_coordinates, return_value_policy<copy_const_reference>())
    .def("tree", &RPartitionTreeSet::tree, return_value_policy<copy_const_reference>())
    .def("mask", &RPartitionTreeSet::mask, return_value_policy<copy_const_reference>()) 
    .def(str(self))    // __str__
  ;
  class_<RPartitionTreeCell>("PartitionTreeCell",init<RRectangle,SubdivisionSequence,BinaryWord>())
    .def(init<RPartitionScheme,BinaryWord>())
//    .def(init<RRectangle,SubdivisionSequence,BinaryWord>())
    .def("dimension", &RPartitionTreeCell::dimension)
    .def("bounding_box", &RPartitionTreeCell::bounding_box, return_value_policy<copy_const_reference>())
    .def("subdivision_coordinates", &RPartitionTreeCell::subdivision_coordinates, return_value_policy<copy_const_reference>())
    .def("word", &RPartitionTreeCell::word, return_value_policy<copy_const_reference>())
    .def(str(self))    // __str__
  ;
}
