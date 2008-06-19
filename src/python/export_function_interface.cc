/***************************************************************************
 *            python/export_function_interface.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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

#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "function/taylor_derivative.h"
#include "function/function_interface.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;
using namespace Ariadne::Python;

template<class R>
class FunctionWrapper
  : public FunctionInterface<R>
  , public wrapper< FunctionInterface<R> >
{
  typedef typename traits<R>::arithmetic_type F;
  virtual FunctionInterface<R>* clone() const { return this->get_override("clone")(); };
  virtual std::string name() const { return this->get_override("name")(); }
  virtual size_type result_size() const { return this->get_override("result_size")(); }
  virtual size_type argument_size() const { return this->get_override("argument_size")(); }
  virtual smoothness_type smoothness() const { return this->get_override("smoothness")(); }
  virtual Vector<F> evaluate(const Vector<F>&) const { return this->get_override("evaluate")(); }
  virtual Matrix<F> jacobian(const Vector<F>&) const { return this->get_override("jacobian")(); }
  virtual TaylorDerivative<F> derivative(const Vector<F>&, const smoothness_type&) const { return this->get_override("derivative")(); }
  virtual std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};


template<class R>
void export_function_interface() 
{
  class_<FunctionWrapper<R>, boost::noncopyable> function_interface_class("FunctionInterface");
  function_interface_class.def("argument_size", pure_virtual(&FunctionInterface<R>::argument_size));
  function_interface_class.def("result_size", pure_virtual(&FunctionInterface<R>::result_size));
  function_interface_class.def("smoothness", pure_virtual(&FunctionInterface<R>::smoothness));
  function_interface_class.def("__call__",&FunctionInterface<R>::operator());
  function_interface_class.def("evaluate",pure_virtual(&FunctionInterface<R>::evaluate));
  function_interface_class.def("jacobian",pure_virtual(&FunctionInterface<R>::jacobian));
  function_interface_class.def("derivative",pure_virtual(&FunctionInterface<R>::derivative));

}

template void export_function_interface<Rational>();
template void export_function_interface<FloatPy>();
