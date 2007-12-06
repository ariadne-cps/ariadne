/***************************************************************************
 *            python/export_polynomial_map.cc
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

#include "python/float.h"

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "system/polynomial_map.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_polynomial_map() 
{
/*
  typedef FloatPy (RPolynomial::* PolyApplyPointFunc) (const RPoint&) const;
  typedef RInterval (RPolynomial::* PolyApplyRectFunc) (const RRectangle&) const;
  typedef RPoint (RPolynomialMap::* PolyMapApplyPointFunc) (const RPoint&) const;
  typedef RRectangle (RPolynomialMap::* PolyMapApplyRectFunc) (const RRectangle&) const;
  typedef RMatrix (RPolynomialMap::* PolyMapDerivPointFunc) (const RPoint&) const;
  typedef RIntervalMatrix (RPolynomialMap::* PolyMapDerivRectFunc) (const RRectangle&) const;
  typedef const RPolynomialMatrix& (RPolynomialMap::* PolyMapDerivFunc) () const;
 
  class_<RMonomial>("Monomial",init<uint>())
    .def(init<FloatPy,SizeArray>())
    .def(self < self)
    .def(self_ns::str(self))
  ;
  
  
  class_<RPolynomial>("Polynomial",init<uint,uint>())
    .def("apply", PolyApplyPointFunc(&RPolynomial::apply))
    .def("apply", PolyApplyRectFunc(&RPolynomial::apply))
    .def("argument_dimension", &RPolynomial::argument_dimension)
    .def(self_ns::str(self))
  ;
  
  class_<RPolynomialMap>("PolynomialMap",init<uint,uint,uint>())
    .def("apply", PolyMapApplyPointFunc(&RPolynomialMap::apply))
    .def("apply", PolyMapApplyRectFunc(&RPolynomialMap::apply))
    .def("derivative", PolyMapDerivPointFunc(&RPolynomialMap::derivative))
    .def("derivative", PolyMapDerivRectFunc(&RPolynomialMap::derivative))
    .def("derivative", PolyMapDerivFunc(&RPolynomialMap::derivative),return_value_policy<copy_const_reference>())
    .def("argument_dimension", &RPolynomialMap::argument_dimension)
    .def("result_dimension", &RPolynomialMap::result_dimension)
    .def(self_ns::str(self))
  ;
  
  class_<RPolynomialMatrix>("PolynomialMatrix",init<uint,uint>())
    .def(self_ns::str(self))
  ;
*/
}

template void export_polynomial_map<FloatPy>();
