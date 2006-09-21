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

#include "python/typedefs.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

typedef void (LatticeCellListSet::* LatCellListAdjCellFunc) (const LatticeCell&);
typedef void (LatticeCellListSet::* LatCellListAdjBlkFunc) (const LatticeBlock&);
typedef void (LatticeCellListSet::* LatCellListAdjCellListFunc) (const LatticeCellListSet&);
typedef void (LatticeCellListSet::* LatCellListAdjBlkListFunc) (const LatticeBlockListSet&);
typedef void (LatticeCellListSet::* LatCellListAdjMaskSetFunc) (const LatticeMaskSet&);

typedef void (LatticeBlockListSet::* LatBlkListAdjBlkFunc) (const LatticeBlock&);
typedef void (LatticeBlockListSet::* LatBlkListAdjBlkListFunc) (const LatticeBlockListSet&);

typedef void (LatticeMaskSet::* MaskSetadjLatCellFunc) (const LatticeCell&);
typedef void (LatticeMaskSet::* MaskSetadjLatBlkFunc) (const LatticeBlock&);
typedef void (LatticeMaskSet::* MaskSetadjLatCellListFunc) (const LatticeCellListSet&);
typedef void (LatticeMaskSet::* MaskSetadjLatBlkListFunc) (const LatticeBlockListSet&);
typedef void (LatticeMaskSet::* MaskSetadjLatMaskSetFunc) (const LatticeMaskSet&);


void export_lattice_set() {
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
    .def("adjoin", LatCellListAdjCellFunc(&LatticeCellListSet::adjoin))
    .def("adjoin", LatCellListAdjBlkFunc(&LatticeCellListSet::adjoin))
    .def("adjoin", LatCellListAdjCellListFunc(&LatticeCellListSet::adjoin))
    .def("adjoin", LatCellListAdjBlkListFunc(&LatticeCellListSet::adjoin))
    .def("adjoin", LatCellListAdjMaskSetFunc(&LatticeCellListSet::adjoin))
    .def("unique_sort",&LatticeCellListSet::unique_sort)
    .def(self_ns::str(self))
    ;
    
  class_<LatticeBlockListSet>("LatticeBlockListSet",init<uint>())
    .def(init<LatticeBlockListSet>())
    .def("dimension", &LatticeBlockListSet::dimension)
    .def("adjoin", LatBlkListAdjBlkFunc(&LatticeBlockListSet::adjoin))
    .def("adjoin", LatBlkListAdjBlkListFunc(&LatticeBlockListSet::adjoin))
    .def(self_ns::str(self))
    ;

  class_<LatticeMaskSet>("LatticeMaskSet",init<LatticeBlock>())
    .def(init<LatticeBlock,LatticeCellListSet>())
    .def(init<LatticeMaskSet>())
    .def("dimension", &LatticeMaskSet::dimension)
    .def("adjoin", MaskSetadjLatCellFunc(&LatticeMaskSet::adjoin))
    .def("adjoin", MaskSetadjLatBlkFunc(&LatticeMaskSet::adjoin))
    .def("adjoin", MaskSetadjLatCellListFunc(&LatticeMaskSet::adjoin))
    .def("adjoin", MaskSetadjLatBlkListFunc(&LatticeMaskSet::adjoin))
    .def("adjoin", MaskSetadjLatMaskSetFunc(&LatticeMaskSet::adjoin))
    .def(self_ns::str(self))
    ;

}
