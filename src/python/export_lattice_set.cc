/***************************************************************************
 *            python/export_lattice_set.cc
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



#include "combinatoric/lattice_set.h"

using namespace Ariadne;
using namespace Ariadne::Combinatoric;

#include <boost/python.hpp>
using namespace boost::python;

void export_lattice_set() {
  class_<LatticePoint>("LatticePoint",init<const IndexArray&>())
    .def(init<const LatticePoint&>())
    .def("__eq__",&LatticePoint::operator==)
    .def("dimension", &LatticePoint::dimension)
    .def(self_ns::str(self))
    ;

  class_<LatticeCell>("LatticeCell",init<const IndexArray&>())
    .def(init<const LatticeCell&>())
    .def("__eq__",&LatticeCell::operator==)
    .def(self < self)
    .def("dimension", &LatticeCell::dimension)
    .def(self_ns::str(self))
    ;

  class_<LatticeBlock>("LatticeBlock",init<IndexArray,IndexArray>())
    .def(init<std::string>())
    .def(init<const LatticeCell&>())
    .def(init<const LatticeBlock&>())
    .def("dimension", &LatticeBlock::dimension)
    .def(self_ns::str(self))
    ;

  class_<LatticeCellListSet>("LatticeCellListSet",init<uint>())
    .def(init<const LatticeCell&>())
    .def(init<const LatticeBlock&>())
    .def(init<const LatticeCellListSet&>())
    .def(init<const LatticeMaskSet&>())
    .def("dimension", &LatticeCellListSet::dimension)
    .def("__getitem__", &LatticeCellListSet::operator[])
    .def("adjoin", (void(LatticeCellListSet::*)(const LatticeCell&))(&LatticeCellListSet::adjoin))
    .def("adjoin", (void(LatticeCellListSet::*)(const LatticeBlock&))(&LatticeCellListSet::adjoin))
    .def("adjoin", (void(LatticeCellListSet::*)(const LatticeCellListSet&))(&LatticeCellListSet::adjoin))
    .def("adjoin", (void(LatticeCellListSet::*)(const LatticeMaskSet&))(&LatticeCellListSet::adjoin))
    .def("unique_sort",&LatticeCellListSet::unique_sort)
    .def(self_ns::str(self))
    ;
    
  class_<LatticeMaskSet>("LatticeMaskSet",init<LatticeBlock>())
    .def(init<LatticeBlock,LatticeCellListSet>())
    .def(init<LatticeMaskSet>())
    .def("dimension", &LatticeMaskSet::dimension)
    .def("adjoin", (void(LatticeMaskSet::*)(const LatticeCell&))(&LatticeMaskSet::adjoin))
    .def("adjoin", (void(LatticeMaskSet::*)(const LatticeBlock&))(&LatticeMaskSet::adjoin))
    .def("adjoin", (void(LatticeMaskSet::*)(const LatticeCellListSet&))(&LatticeMaskSet::adjoin))
    .def("adjoin", (void(LatticeMaskSet::*)(const LatticeMaskSet&))(&LatticeMaskSet::adjoin))
    .def(self_ns::str(self))
    ;

}
