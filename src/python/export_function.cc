/***************************************************************************
 *            python/export_function.cc
 *
 *  21 October 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can rediself_ns::stribute it and/or modify
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

#include "numeric/function.h"

#include <boost/python.hpp>
using namespace boost::python;

#include "python/typedefs.h"
using namespace Ariadne;

void export_function() {
  def("div_approx", div_approx<Real>, "approximate division function (maximum error e)" );
  def("sqrt_approx", sqrt_approx<Real>, "approximate square root function (maximum error e)" );
  def("exp_approx", exp_approx<Real>, "approximate exponential function (maximum error e)" );
  def("cos_approx", cos_approx<Real>, "approximate sine function (maximum error e)" );
  def("sin_approx", sin_approx<Real>, "approximate cosine function (maximum error e)" );
}
