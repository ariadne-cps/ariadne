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



#include "geometry/lattice_set.h"

#include "python/typedefs.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

typedef void (LatticeCellListSet::* LatCellListAdjCellFunc) (const LatticeCell&);
typedef void (LatticeCellListSet::* LatCellListAdjRectFunc) (const LatticeRectangle&);
typedef void (LatticeCellListSet::* LatCellListAdjCellListFunc) (const LatticeCellListSet&);
typedef void (LatticeCellListSet::* LatCellListAdjRectListFunc) (const LatticeRectangleListSet&);
typedef void (LatticeCellListSet::* LatCellListAdjMaskSetFunc) (const LatticeMaskSet&);

typedef void (LatticeRectangleListSet::* LatRectListAdjRectFunc) (const LatticeRectangle&);
typedef void (LatticeRectangleListSet::* LatRectListAdjRectListFunc) (const LatticeRectangleListSet&);

typedef void (LatticeMaskSet::* MaskSetadjLatCellFunc) (const LatticeCell&);
typedef void (LatticeMaskSet::* MaskSetadjLatRectFunc) (const LatticeRectangle&);
typedef void (LatticeMaskSet::* MaskSetadjLatCellListFunc) (const LatticeCellListSet&);
typedef void (LatticeMaskSet::* MaskSetadjLatRectListFunc) (const LatticeRectangleListSet&);
typedef void (LatticeMaskSet::* MaskSetadjLatMaskSetFunc) (const LatticeMaskSet&);


void export_lattice_set() {
  class_<LatticeCell>("LatticeCell",init<const IndexArray&>())
    .def(init<const LatticeCell&>())
    .def("__eq__",&LatticeCell::operator==)
    .def(self < self)
    .def("dimension", &LatticeCell::dimension)
    .def(self_ns::str(self))
    ;

  class_<LatticeRectangle>("LatticeRectangle",init<IndexArray,IndexArray>())
    .def(init<std::string>())
    .def(init<const LatticeCell&>())
    .def(init<const LatticeRectangle&>())
    .def("dimension", &LatticeRectangle::dimension)
    .def(self_ns::str(self))
    ;

  class_<LatticeCellListSet>("LatticeCellListSet",init<uint>())
    .def(init<const LatticeCell&>())
    .def(init<const LatticeRectangle&>())
    .def(init<const LatticeCellListSet&>())
    .def(init<const LatticeMaskSet&>())
    .def("dimension", &LatticeCellListSet::dimension)
    .def("__getitem__", &LatticeCellListSet::operator[])
    .def("adjoin", LatCellListAdjCellFunc(&LatticeCellListSet::adjoin))
    .def("adjoin", LatCellListAdjRectFunc(&LatticeCellListSet::adjoin))
    .def("adjoin", LatCellListAdjCellListFunc(&LatticeCellListSet::adjoin))
    .def("adjoin", LatCellListAdjRectListFunc(&LatticeCellListSet::adjoin))
    .def("adjoin", LatCellListAdjMaskSetFunc(&LatticeCellListSet::adjoin))
    .def("unique_sort",&LatticeCellListSet::unique_sort)
    .def(self_ns::str(self))
    ;
    
  class_<LatticeRectangleListSet>("LatticeRectangleListSet",init<uint>())
    .def(init<LatticeRectangleListSet>())
    .def("dimension", &LatticeRectangleListSet::dimension)
    .def("adjoin", LatRectListAdjRectFunc(&LatticeRectangleListSet::adjoin))
    .def("adjoin", LatRectListAdjRectListFunc(&LatticeRectangleListSet::adjoin))
    .def(self_ns::str(self))
    ;

  class_<LatticeMaskSet>("LatticeMaskSet",init<LatticeRectangle>())
    .def(init<LatticeRectangle,LatticeCellListSet>())
    .def(init<LatticeMaskSet>())
    .def("dimension", &LatticeMaskSet::dimension)
    .def("adjoin", MaskSetadjLatCellFunc(&LatticeMaskSet::adjoin))
    .def("adjoin", MaskSetadjLatRectFunc(&LatticeMaskSet::adjoin))
    .def("adjoin", MaskSetadjLatCellListFunc(&LatticeMaskSet::adjoin))
    .def("adjoin", MaskSetadjLatRectListFunc(&LatticeMaskSet::adjoin))
    .def("adjoin", MaskSetadjLatMaskSetFunc(&LatticeMaskSet::adjoin))
    .def(self_ns::str(self))
    ;

}
