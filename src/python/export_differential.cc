/***************************************************************************
 *            python/export_differential.cc
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
#include "numeric/differential.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"


#include <boost/python.hpp>
#include "python/python_utilities.h"

using namespace boost::python;
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;


template<class R1, class R2> inline 
void scalar_derivative_set_item(ScalarDerivative<R1>& sd, uint i, R2 x) {
  assert(i<=sd.degree()); 
  sd[i]=x;
}

template<class R> inline 
R scalar_derivative_get_item(const ScalarDerivative<R>& sd, uint i) {
  assert(i<=sd.degree()); 
  return sd[i];
}

template<class R>
void export_differential()
{
  typedef typename Numeric::traits<R>::arithmetic_type A;
  typedef Differential< A, Vector<A> > D;

  class_<D>(python_name<R>("Differential").c_str())
    .def( init< double,Vector<R> >())
    .def( init< double,Vector<A> >())
    .def( init< R,Vector<R> >())
    .def( init< R,Vector<A> >())
    .def( init< A,Vector<R> >())
    .def( init< A,Vector<A> >())
    .def("__neg__", &Ariadne::neg<D,D>)
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
  
  def("abs",(D(*)(const D&))&Numeric::abs);
  def("pow",(D(*)(const D&, int))&Numeric::pow);

  def("sqrt", (D(*)(const D&))&Numeric::sqrt);
  def("exp", (D(*)(const D&))&Numeric::exp);
  def("log", (D(*)(const D&))&Numeric::log);
  def("sin", (D(*)(const D&))&Numeric::sin);
  def("cos", (D(*)(const D&))&Numeric::cos);
  def("tan", (D(*)(const D&))&Numeric::tan);
  def("asin", (D(*)(const D&))&Numeric::asin);
  def("acos", (D(*)(const D&))&Numeric::acos);
  def("atan", (D(*)(const D&))&Numeric::atan);
}

template<>
void export_differential<Rational>()
{
  typedef Rational Q;
  typedef Differential< Q, Vector<Q> > D;

  class_<D>("QDifferential")
    .def( init< double,Vector<Q> >())
    .def( init< Q,Vector<Q> >())
    .def("__neg__", &Ariadne::neg<D,D>)
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
  
  def("pow",(D(*)(const D&, int))&Numeric::pow);
}






template<class R>
void export_derivative()
{
  typedef typename Numeric::traits<R>::arithmetic_type A;
  typedef ScalarDerivative<A> SD;


  class_<SD>(python_name<R>("ScalarDerivative").c_str())
    .def( init< uint >())
    .def( init< uint, double >())
    .def( init< uint, double, double >())
    .def( init< uint, A >())
    .def( init< uint, A, A >())
    .def("__getitem__", &scalar_derivative_get_item<A>)
    .def("__setitem__",&scalar_derivative_set_item<A,double>)
    .def("__setitem__",&scalar_derivative_set_item<A,R>)
    .def("__setitem__",&scalar_derivative_set_item<A,A>)
    .def("__neg__", &Ariadne::neg<SD,SD>)
    .def("__add__", &add<SD,SD,SD>)
    .def("__add__", &add<SD,SD,double>)
    .def("__add__", &add<SD,SD,R>)
    .def("__add__", &add<SD,SD,A>)
    .def("__radd__", &radd<SD,SD,double>)
    .def("__radd__", &radd<SD,SD,R>)
    .def("__radd__", &radd<SD,SD,A>)
    .def("__sub__", &sub<SD,SD,SD>)
    .def("__sub__", &sub<SD,SD,double>)
    .def("__sub__", &sub<SD,SD,R>)
    .def("__sub__", &sub<SD,SD,A>)
    .def("__rsub__", &rsub<SD,SD,double>)
    .def("__rsub__", &rsub<SD,SD,R>)
    .def("__rsub__", &rsub<SD,SD,A>)
    .def("__mul__", &mul<SD,SD,SD>)
    .def("__mul__", &mul<SD,SD,double>)
    .def("__mul__", &mul<SD,SD,R>)
    .def("__mul__", &mul<SD,SD,A>)
    .def("__rmul__", &rmul<SD,SD,double>)
    .def("__rmul__", &rmul<SD,SD,R>)
    .def("__rmul__", &rmul<SD,SD,A>)
    .def("__div__", &div<SD,SD,SD>)
    .def("__div__", &div<SD,SD,double>)
    .def("__div__", &div<SD,SD,R>)
    .def("__div__", &div<SD,SD,A>)
    .def("__rdiv__", &rdiv<SD,SD,double>)
    .def("__rdiv__", &rdiv<SD,SD,R>)
    .def("__rdiv__", &rdiv<SD,SD,A>)
    .def("__pow__", &pow<SD,SD,int>)
    .def(self_ns::str(self))
 ;
  
  def("abs",(SD(*)(const SD&))&Numeric::abs);
  def("pow",(SD(*)(const SD&, int))&Numeric::pow);

  def("sqrt", (SD(*)(const SD&))&Numeric::sqrt);
  def("exp", (SD(*)(const SD&))&Numeric::exp);
  def("log", (SD(*)(const SD&))&Numeric::log);
  def("sin", (SD(*)(const SD&))&Numeric::sin);
  def("cos", (SD(*)(const SD&))&Numeric::cos);
  def("tan", (SD(*)(const SD&))&Numeric::tan);
  def("asin", (SD(*)(const SD&))&Numeric::asin);
  def("acos", (SD(*)(const SD&))&Numeric::acos);
  def("atan", (SD(*)(const SD&))&Numeric::atan);

}


template<>
void export_derivative<Rational>()
{
  typedef Rational Q;
  typedef ScalarDerivative<Q> SD;


  class_<SD>(python_name<Q>("ScalarDerivative").c_str())
    .def( init< std::string >())
    .def( init< uint >())
    .def( init< uint, double >())
    .def( init< uint, Q >())
    .def( init< uint, double, double >())
    .def( init< uint, Q, Q >())
    .def("__getitem__", &scalar_derivative_get_item<Q>)
    .def("__setitem__",&scalar_derivative_set_item<Q,double>)
    .def("__setitem__",&scalar_derivative_set_item<Q,Q>)
    .def("__neg__", &Ariadne::neg<SD,SD>)
    .def("__add__", &add<SD,SD,SD>)
    .def("__add__", &add<SD,SD,double>)
    .def("__add__", &add<SD,SD,Q>)
    .def("__radd__", &radd<SD,SD,double>)
    .def("__radd__", &radd<SD,SD,Q>)
    .def("__sub__", &sub<SD,SD,SD>)
    .def("__sub__", &sub<SD,SD,double>)
    .def("__sub__", &sub<SD,SD,Q>)
    .def("__rsub__", &rsub<SD,SD,double>)
    .def("__rsub__", &rsub<SD,SD,Q>)
    .def("__mul__", &mul<SD,SD,SD>)
    .def("__mul__", &mul<SD,SD,double>)
    .def("__mul__", &mul<SD,SD,Q>)
    .def("__rmul__", &rmul<SD,SD,double>)
    .def("__rmul__", &rmul<SD,SD,Q>)
    .def("__div__", &div<SD,SD,SD>)
    .def("__div__", &div<SD,SD,double>)
    .def("__div__", &div<SD,SD,Q>)
    .def("__rdiv__", &rdiv<SD,SD,double>)
    .def("__rdiv__", &rdiv<SD,SD,Q>)
    .def("__pow__", &pow<SD,SD,int>)
    .def(self_ns::str(self))
 ;
  
  def("abs",(SD(*)(const SD&))&Numeric::abs);
  def("pow",(SD(*)(const SD&, int))&Numeric::pow);

}




template void export_differential<Rational>();
template void export_differential<Float>();
template void export_derivative<Rational>();
template void export_derivative<Float>();
