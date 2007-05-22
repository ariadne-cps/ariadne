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

#include "numeric/differential.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"


#include <boost/python.hpp>
#include "python/python_utilities.h"

using namespace boost::python;
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;


template<class R>
void export_differential()
{
  typedef typename Numeric::traits<R>::arithmetic_type A;
  typedef Differential< A, Vector<A> > D;
  typedef Differential<A,A> SD;
  typedef SecondDifferential<A,A,A> SSD;

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


  class_<SD>(python_name<R>("ScalarDifferential").c_str())
    .def( init< double,double >())
    .def( init< R,double >())
    .def( init< A,double >())
    .def( init< R,R >())
    .def( init< A,R >())
    .def( init< A,A >())
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
    .def("value",(const A&(SD::*)()const) &SD::value, return_value_policy<copy_const_reference>())
    .def("derivative", (const A&(SD::*)()const) &SD::derivative, return_value_policy<copy_const_reference>())
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



  class_<SSD>(python_name<R>("ScalarSecondDifferential").c_str())
    .def(init<double,double,double >())
    .def(init<R,double,double >())
    .def(init<A,double,double >())
    .def(init<R,R,R>())
    .def(init<A,R,R>())  
    .def(init<A,A,A>())
    .def("__neg__", &Ariadne::neg<SSD,SSD>)
    .def("__add__", &add<SSD,SSD,SSD>)
    .def("__add__", &add<SSD,SSD,double>)
    .def("__add__", &add<SSD,SSD,R>)
    .def("__add__", &add<SSD,SSD,A>)
    .def("__radd__", &radd<SSD,SSD,double>)
    .def("__radd__", &radd<SSD,SSD,R>)
    .def("__radd__", &radd<SSD,SSD,A>)
    .def("__sub__", &sub<SSD,SSD,SSD>)
    .def("__sub__", &sub<SSD,SSD,double>)
    .def("__sub__", &sub<SSD,SSD,R>)
    .def("__sub__", &sub<SSD,SSD,A>)
    .def("__rsub__", &rsub<SSD,SSD,double>)
    .def("__rsub__", &rsub<SSD,SSD,R>)
    .def("__rsub__", &rsub<SSD,SSD,A>)
    .def("__mul__", &mul<SSD,SSD,SSD>)
    .def("__mul__", &mul<SSD,SSD,double>)
    .def("__mul__", &mul<SSD,SSD,R>)
    .def("__mul__", &mul<SSD,SSD,A>)
    .def("__rmul__", &rmul<SSD,SSD,double>)
    .def("__rmul__", &rmul<SSD,SSD,R>)
    .def("__rmul__", &rmul<SSD,SSD,A>)
    .def("__div__", &div<SSD,SSD,SSD>)
    .def("__div__", &div<SSD,SSD,double>)
    .def("__div__", &div<SSD,SSD,R>)
    .def("__div__", &div<SSD,SSD,A>)
    .def("__rdiv__", &rdiv<SSD,SSD,double>)
    .def("__rdiv__", &rdiv<SSD,SSD,R>)
    .def("__rdiv__", &rdiv<SSD,SSD,A>)
    .def("__pow__", &pow<SSD,SSD,int>)
    .def("value",(const A&(SSD::*)()const) &SSD::value, return_value_policy<copy_const_reference>())
    .def("derivative", (const A&(SSD::*)()const) &SSD::derivative, return_value_policy<copy_const_reference>())
    .def("first_derivative", (const A&(SSD::*)()const) &SSD::first_derivative, return_value_policy<copy_const_reference>())
    .def("second_derivative", (const A&(SSD::*)()const) &SSD::second_derivative, return_value_policy<copy_const_reference>())
    .def(self_ns::str(self))
    ;
  
  def("abs",(SSD(*)(const SSD&))&Numeric::abs);
  def("pow",(SSD(*)(const SSD&, int))&Numeric::pow);

  def("sqrt", (SSD(*)(const SSD&))&Numeric::sqrt);
  def("exp", (SSD(*)(const SSD&))&Numeric::exp);
  def("log", (SSD(*)(const SSD&))&Numeric::log);
  def("sin", (SSD(*)(const SSD&))&Numeric::sin);
  def("cos", (SSD(*)(const SSD&))&Numeric::cos);
  def("tan", (SSD(*)(const SSD&))&Numeric::tan);
  def("asin", (SSD(*)(const SSD&))&Numeric::asin);
  def("acos", (SSD(*)(const SSD&))&Numeric::acos);
  def("atan", (SSD(*)(const SSD&))&Numeric::atan);

}


template<>
void export_differential<Rational>()
{
  typedef Rational Q;
  typedef Differential< Q, Vector<Q> > D;
  typedef Differential<Q,Q> SD;
  typedef SecondDifferential<Q,Q,Q> SSD;

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





  class_<SD>(python_name<Q>("ScalarDifferential").c_str())
    .def(init<double,double>())
    .def(init<Q,double>())
    .def(init<Q,Q>())
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
    .def("value",(const Q&(SD::*)()const) &SD::value, return_value_policy<copy_const_reference>())
    .def("derivative", (const Q&(SD::*)()const) &SD::derivative, return_value_policy<copy_const_reference>())
    .def(self_ns::str(self))
 ;
  
  def("pow",(SD(*)(const SD&, int))&Numeric::pow);



  class_<SSD>(python_name<Q>("ScalarSecondDifferential").c_str())
    .def(init<double,double,double >())
    .def(init<Q,double,double >())
    .def(init<Q,Q,Q>())
    .def("__neg__", &Ariadne::neg<SSD,SSD>)
    .def("__add__", &add<SSD,SSD,SSD>)
    .def("__add__", &add<SSD,SSD,double>)
    .def("__add__", &add<SSD,SSD,Q>)
    .def("__radd__", &radd<SSD,SSD,double>)
    .def("__radd__", &radd<SSD,SSD,Q>)
    .def("__sub__", &sub<SSD,SSD,SSD>)
    .def("__sub__", &sub<SSD,SSD,double>)
    .def("__sub__", &sub<SSD,SSD,Q>)
    .def("__rsub__", &rsub<SSD,SSD,double>)
    .def("__rsub__", &rsub<SSD,SSD,Q>)
    .def("__mul__", &mul<SSD,SSD,SSD>)
    .def("__mul__", &mul<SSD,SSD,double>)
    .def("__mul__", &mul<SSD,SSD,Q>)
    .def("__rmul__", &rmul<SSD,SSD,double>)
    .def("__rmul__", &rmul<SSD,SSD,Q>)
    .def("__div__", &div<SSD,SSD,SSD>)
    .def("__div__", &div<SSD,SSD,double>)
    .def("__div__", &div<SSD,SSD,Q>)
    .def("__rdiv__", &div<SSD,SSD,double>)
    .def("__rdiv__", &div<SSD,SSD,Q>)
    .def("__pow__", &pow<SSD,SSD,int>)
    .def("value",(const Q&(SSD::*)()const) &SSD::value, return_value_policy<copy_const_reference>())
    .def("derivative", (const Q&(SSD::*)()const) &SSD::derivative, return_value_policy<copy_const_reference>())
    .def("first_derivative", (const Q&(SSD::*)()const) &SSD::first_derivative, return_value_policy<copy_const_reference>())
    .def("second_derivative", (const Q&(SSD::*)()const) &SSD::second_derivative, return_value_policy<copy_const_reference>())
    .def(self_ns::str(self))
 ;
  
  def("pow",(SSD(*)(const SSD&, int))&Numeric::pow);




}



template<class R>
void export_second_differential()
{
  typedef typename Numeric::traits<R>::arithmetic_type A;

}



template void export_differential<Rational>();
template void export_differential<Float>();
