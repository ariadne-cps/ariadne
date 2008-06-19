/***************************************************************************
 *            python/export_lattice_map.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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



#include "combinatoric/lattice_map.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

/* WARNING: for some reason, we need this helper routing to avoid segmentation faults. */
inline LatticeCellListSet apply_cell(const LatticeMultiMap& lm, const LatticeCell& lc) { 
  const LatticeCellListSet& ref=lm.apply(lc); 
  LatticeCellListSet res=ref;
  return res;
}

void export_lattice_map() {
  class_<LatticeMultiMap>("LatticeMultiMap",init<dimension_type,dimension_type>())
    .def("adjoin_to_image", (void(LatticeMultiMap::*)(const LatticeCell&,const LatticeCell&))(&LatticeMultiMap::adjoin_to_image))
    .def("adjoin_to_image", (void(LatticeMultiMap::*)(const LatticeCell&,const LatticeCellListSet&))(&LatticeMultiMap::adjoin_to_image))
    .def("__call__", &apply_cell)
    .def("__call__", (LatticeCellListSet(LatticeMultiMap::*)(const LatticeCell&)const)(&LatticeMultiMap::operator()))
    .def("__call__", (LatticeCellListSet(LatticeMultiMap::*)(const LatticeBlock&)const)(&LatticeMultiMap::operator()))
    .def("__call__", (LatticeCellListSet(LatticeMultiMap::*)(const LatticeCellListSet&)const)(&LatticeMultiMap::operator()))
//    .def("__call__", LatticeCellListSet(LatticeMultiMap::*)(const LatticeMaskSet&)const)(&LatticeMultiMap::operator()))
    .def("argument_dimension", &LatticeMultiMap::argument_dimension)
    .def("result_dimension", &LatticeMultiMap::result_dimension)
    .def(self_ns::str(self))
  ;
}
