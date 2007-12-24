/***************************************************************************
 *            python/export_latex_output.cc
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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

#include "numeric/integer.h"
#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/rectangle.h"

#include "output/latexstream.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Output;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class T> inline latexfstream& write(latexfstream& txs, const T& t) { return static_cast<latexfstream&>(txs << t); }

void export_latex_output()
{

  class_<latexfstream, boost::noncopyable>("LatexFile",init<>())
    .def("open",(void(latexfstream::*)(const char*))&latexfstream::open)
    .def("open",(void(latexfstream::*)(const char*,const char*))&latexfstream::open)
    .def("close",(void(latexfstream::*)())&latexfstream::close)
    .def("write",&write< char >,return_internal_reference<1>())
    .def("write",&write< char* >,return_internal_reference<1>())
    .def("write",&write< Integer >,return_internal_reference<1>())
    .def("write",&write< Rational >,return_internal_reference<1>())
    .def("write",&write< FloatPy >,return_internal_reference<1>())
    .def("write",&write< Interval<FloatPy> >,return_internal_reference<1>())
    .def("write",&write< Vector<FloatPy> >,return_internal_reference<1>())
    .def("write",&write< Matrix<FloatPy> >,return_internal_reference<1>())
    .def("write",&write< Box<FloatPy> >,return_internal_reference<1>())
  ;
  
}
