/***************************************************************************
 *            python/linear_algebra_module.cc
 *
 *  17 November 2005
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

#include <boost/python.hpp>

#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"

using namespace Ariadne::Numeric;

template<class R> void export_vector();
template<class R> void export_matrix();
template<class R> void export_tensor();
template<class R> void export_linear_program();

template<class R> void export_interval_vector();
template<class R> void export_interval_matrix();
template<class R> void export_interval_tensor();


BOOST_PYTHON_MODULE(linear_algebra)
{
  export_vector<Float64>();
  export_vector<MPFloat>();
  export_vector<Rational>();

  export_matrix<Float64>();
  export_matrix<MPFloat>();
  export_matrix<Rational>();
  
  export_tensor<Float64>();
  export_tensor<MPFloat>();
  export_tensor<Rational>();
  
  export_linear_program<Rational>();
  
  export_interval_vector<Float64>();
  export_interval_vector<MPFloat>();
  
  export_interval_matrix<Float64>();
  export_interval_matrix<MPFloat>();
  
  export_interval_tensor<Float64>();
  export_interval_tensor<MPFloat>();
}
