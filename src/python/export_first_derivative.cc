/***************************************************************************
 *            python/export_first_derivative.cc
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
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 

#include "python/python_float.h"

#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "function/first_derivative.h"


#include <boost/python.hpp>
#include "python/python_utilities.h"

using namespace boost::python;
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Function;
using namespace Ariadne::Python;



template<class R>
void export_first_derivative()
{
  typedef typename Numeric::traits<R>::arithmetic_type A;
  typedef FirstDerivative< A, Vector<A> > D;

  class_<D>(python_name<R>("FirstDerivative").c_str())
    .def( init< double,Vector<R> >())
    .def( init< double,Vector<A> >())
    .def( init< R,Vector<R> >())
    .def( init< R,Vector<A> >())
    .def( init< A,Vector<R> >())
    .def( init< A,Vector<A> >())
    .def("__neg__", &Python::neg<D,D>)
    .def("__add__", &add<D,D,D>)
    .def("__add__", &add<D,D,double>)
    .def("__add__", &add<D,D,R>)
    .def("__add__", &add<D,D,A>)
    .def("__radd__", &radd<D,D,double>)
    .def("__radd__", &radd<D,D,R>)
    .def("__radd__", &radd<D,D,A>)
    .def("__sub__", &sub<D,D,D>)
    .def("__sub__", &sub<D,D,double>)
    .def("__sub__", &sub<D,D,R>)
    .def("__sub__", &sub<D,D,A>)
    .def("__rsub__", &rsub<D,D,double>)
    .def("__rsub__", &rsub<D,D,R>)
    .def("__rsub__", &rsub<D,D,A>)
    .def("__mul__", &mul<D,D,D>)
    .def("__mul__", &mul<D,D,double>)
    .def("__mul__", &mul<D,D,R>)
    .def("__mul__", &mul<D,D,A>)
    .def("__rmul__", &mul<D,D,double>)
    .def("__rmul__", &mul<D,D,R>)
    .def("__rmul__", &mul<D,D,A>)
    .def("__div__", &div<D,D,D>)
    .def("__div__", &div<D,D,double>)
    .def("__div__", &div<D,D,R>)
    .def("__div__", &div<D,D,A>)
    .def("__rdiv__", &rdiv<D,D,double>)
    .def("__rdiv__", &rdiv<D,D,R>)
    .def("__rdiv__", &rdiv<D,D,A>)
    .def("__pow__", &pow<D,D,int>)
    .def("value",(const A&(D::*)()const) &D::value, return_value_policy<copy_const_reference>())
    .def("derivative", (const Vector<A>&(D::*)()const) &D::derivative, return_value_policy<copy_const_reference>())
    .def(self_ns::str(self))
 ;
  
  def("abs",(D(*)(const D&))&Function::abs);
  def("pow",(D(*)(const D&, int))&Function::pow);

  def("sqrt", (D(*)(const D&))&Function::sqrt);
  def("exp", (D(*)(const D&))&Function::exp);
  def("log", (D(*)(const D&))&Function::log);
  def("sin", (D(*)(const D&))&Function::sin);
  def("cos", (D(*)(const D&))&Function::cos);
  def("tan", (D(*)(const D&))&Function::tan);
  def("asin", (D(*)(const D&))&Function::asin);
  def("acos", (D(*)(const D&))&Function::acos);
  def("atan", (D(*)(const D&))&Function::atan);
}

template<>
void export_first_derivative<Rational>()
{
  typedef Rational Q;
  typedef FirstDerivative< Q, Vector<Q> > D;

  class_<D>("QFirstDerivative")
    .def( init< double,Vector<Q> >())
    .def( init< Q,Vector<Q> >())
    .def("__neg__", &Python::neg<D,D>)
    .def("__add__", &add<D,D,D>)
    .def("__add__", &add<D,D,double>)
    .def("__add__", &add<D,D,Q>)
    .def("__radd__", &radd<D,D,double>)
    .def("__radd__", &radd<D,D,Q>)
    .def("__sub__", &sub<D,D,D>)
    .def("__sub__", &sub<D,D,double>)
    .def("__sub__", &sub<D,D,Q>)
    .def("__rsub__", &rsub<D,D,double>)
    .def("__rsub__", &rsub<D,D,Q>)
    .def("__mul__", &mul<D,D,D>)
    .def("__mul__", &mul<D,D,double>)
    .def("__mul__", &mul<D,D,Q>)
    .def("__rmul__", &rmul<D,D,double>)
    .def("__rmul__", &rmul<D,D,Q>)
    .def("__div__", &div<D,D,D>)
    .def("__div__", &div<D,D,double>)
    .def("__div__", &div<D,D,Q>)
    .def("__rdiv__", &rdiv<D,D,double>)
    .def("__rdiv__", &rdiv<D,D,Q>)
    .def("__pow__", &pow<D,D,int>)
    .def("value",(const Q&(D::*)()const) &D::value, return_value_policy<copy_const_reference>())
    .def("derivative", (const Vector<Q>&(D::*)()const) &D::derivative, return_value_policy<copy_const_reference>())
    .def(self_ns::str(self))
 ;
  
  def("pow",(D(*)(const D&, int))&Function::pow);
}







template void export_first_derivative<Rational>();
template void export_first_derivative<Float>();
