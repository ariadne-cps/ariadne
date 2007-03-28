/***************************************************************************
 *            python/export_function.cc
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

#include "numeric/float64.h"
#include "numeric/mpfloat.h"

#include "numeric/function.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne::Numeric;

template<class R>
void export_function() {
  def("add_approx",(R(*)(const R&,const R&))(add_approx<R>), "approximate addition function" );
  def("sub_approx",(R(*)(const R&,const R&))(sub_approx<R>), "approximate subtraction function" );
  def("mul_approx",(R(*)(const R&,const R&))(mul_approx<R>), "approximate multiplication function" );
  def("div_approx",(R(*)(const R&,const R&))(div_approx<R>), "approximate division function" );
  def("mul_approx",(R(*)(const R&,const int&))(mul_approx<R>), "approximate multiplication function" );
  def("div_approx",(R(*)(const R&,const int&))(div_approx<R>), "approximate division function" );
  def("sqrt_approx", sqrt_approx<R>, "approximate square root function" );
  def("exp_approx", exp_approx<R>, "approximate exponential function" );
  def("log_approx", log_approx<R>, "approximate logarithm function" );
  //def("sin_approx", sin_approx<R>, "approximate cosine function" );
  //def("cos_approx", cos_approx<R>, "approximate sine function" );
  //def("tan_approx", tan_approx<R>, "approximate tangent function" );
}

template void export_function<Float64>();
template void export_function<MPFloat>();
