/***************************************************************************
 *            python/export_text_output.cc
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins, Davide Bresolin
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl, bresolin@sci.univr.it
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

#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "geometry/polyhedron.h"
#include "geometry/list_set.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"

#include "output/textstream.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::Output;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class S> inline void write(textfstream& txt, const S& s) { txt << s; }
template<class R> inline void textfstream_open(textfstream& txt) { txt.open("Ariadne"); }
inline void textfstream_close(textfstream& txt) { txt.close(); }

void export_text_output()
{
    
  class_<textfstream, boost::noncopyable>("TextFile",init<>())
    .def("open",(void(textfstream::*)(const char* fn))&textfstream::open)
    .def("close",&textfstream_close)
    .def("write",&write< Box<FloatPy> >)
    .def("write",&write< BoxListSet<FloatPy> >)
    .def("write",&write< Zonotope<FloatPy> >)
    .def("write",&write< Polytope<FloatPy> >)
    .def("write",&write< Polyhedron<FloatPy> >)
    .def("write",&write< ListSet< Polytope<FloatPy> > >)
    .def("write",&write< ListSet< Zonotope<FloatPy> > >)
    .def("write",&write< GridCell<FloatPy> >)
    .def("write",&write< GridBlock<FloatPy> >)
    .def("write",&write< GridCellListSet<FloatPy> >)
    .def("write",&write< GridMaskSet<FloatPy> >)
    .def("write",&write< PartitionTreeSet<FloatPy> >)
  ;
  
}
