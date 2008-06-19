/***************************************************************************
 *            python/export_curve.cc
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
#include "python/utilities.h"
#include "python/read_array.h"

#include "geometry/interpolated_curve.h"


using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
using namespace boost::python;


template<class R>
void export_segment() 
{
  class_< Segment<R> > segment_class(python_name<R>("Segment").c_str(),init< Point<R>, Point<R> >());
  segment_class.def("dimension", &Segment<R>::dimension);
  segment_class.def("initial_point", &Segment<R>::initial_point, return_value_policy<copy_const_reference>());
  segment_class.def("final_point", &Segment<R>::final_point, return_value_policy<copy_const_reference>());
  segment_class.def(self_ns::str(self));

}

template<class R>
void export_interpolated_curve() 
{
  class_< InterpolatedCurve<R> > interpolated_curve_class(python_name<R>("InterpolatedCurve").c_str(),init< R, Point<R> >());
  interpolated_curve_class.def(init< Point<R> >());
  interpolated_curve_class.def(init< Point<R>, Point<R> >());
  interpolated_curve_class.def(init< Segment<R> >());
  interpolated_curve_class.def("dimension", &InterpolatedCurve<R>::dimension);
  interpolated_curve_class.def("insert", (void(InterpolatedCurve<R>::*)(const double&, const Point<R>&)) &InterpolatedCurve<R>::insert);
  interpolated_curve_class.def("insert", (void(InterpolatedCurve<R>::*)(const R&, const Point<R>&)) &InterpolatedCurve<R>::insert);
  interpolated_curve_class.def(self_ns::str(self));

}

template void export_segment<Rational>();
template void export_segment<FloatPy>();

template void export_interpolated_curve<Rational>();
template void export_interpolated_curve<FloatPy>();
