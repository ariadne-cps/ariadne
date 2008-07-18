/***************************************************************************
 *            python/export_vector.cc
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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

#include "numeric/rational.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"

// Need these since we can't define v*cv using __rmul__
#include "linear_algebra/covector.h" 
#include "linear_algebra/matrix.h"

#include "python/float.h"
#include "python/name.h"
#include "python/subscripting.h"
#include "python/operators.h"
#include "python/read_array.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
using namespace boost::python;

template<class X>
std::string
__str__(const Vector<X>& v)
{
  std::stringstream ss;
  ss << v;
  return ss.str();
}

template<class X>
std::string
__repr__(const Vector<X>& v)
{
  std::stringstream ss;
  ss << "Vector(" << v << ")";
  return ss.str();
}


template<class X>
Vector<X>*
make_vector(const boost::python::object& obj) 
{
  Vector<X>* v=new Vector<X>;
  read_array(v->data(),obj);
  return v;
}




template<class R>
void export_vector()
{
  typedef typename traits<R>::approximate_arithmetic_type A;

  typedef Interval<R> I;
  typedef Vector<R> Vec;
  typedef Vector< Interval<R> > IVec;
  typedef Covector<R> Cvec;
  typedef Covector< Interval<R> > ICvec;
  typedef Matrix< Interval<R> > IMx;
  
  class_< Vector<A> > approximate_vector_class("ApproximateVector",init< Vector<R> >());
  approximate_vector_class.def("__len__", &Vector<A>::size);
  approximate_vector_class.def("__getitem__",&__getitem__< Vector<A> >);

  class_< Vector<R> > vector_class(python_name<R>("Vector").c_str(),no_init);
  vector_class.def("__init__", make_constructor(&make_vector<R>));
  vector_class.def(init<int>());
  vector_class.def(init< Vector<A> >());
  //vector_class.def(init<std::string>());
  vector_class.def(init<Vec>());
  vector_class.def("__len__", &Vec::size);
  vector_class.def("__getitem__",&__getitem__<Vec>);
  vector_class.def("__setitem__",&__setitem__<Vec,R>);
  vector_class.def("__setitem__",&__setitem__<Vec,double>);
  vector_class.def("__neg__",&__neg__<Vec,Vec>);
  vector_class.def("__add__",&__add__<IVec,Vec,Vec>);
  vector_class.def("__add__",&__add__<IVec,Vec,IVec>);
  vector_class.def("__sub__",&__sub__<IVec,Vec,Vec>);
  vector_class.def("__sub__",&__sub__<IVec,Vec,IVec>);
  vector_class.def("__rmul__",&__mul__<IVec,Vec,int,Vec,R>);
  vector_class.def("__rmul__",&__mul__<IVec,Vec,double,Vec,R>);
  vector_class.def("__rmul__",&__mul__<IVec,Vec,R>);
  vector_class.def("__rmul__",&__mul__<IVec,Vec,I>);
  vector_class.def("__mul__",&__mul__<IVec,Vec,int,Vec,R>);
  vector_class.def("__mul__",&__mul__<IVec,Vec,double,Vec,R>);
  vector_class.def("__mul__",&__mul__<IVec,Vec,R>);
  vector_class.def("__mul__",&__mul__<IVec,Vec,I>);
  vector_class.def("__mul__",&__mul__<IMx,Vec,Cvec>);
  vector_class.def("__mul__",&__mul__<IMx,Vec,ICvec>);
  vector_class.def("__div__",&__div__<IVec,Vec,int,Vec,R>);
  vector_class.def("__div__",&__div__<IVec,Vec,double,Vec,R>);
  vector_class.def("__div__",&__div__<IVec,Vec,R>);
  vector_class.def("__div__",&__div__<IVec,Vec,I>);
  vector_class.def("__str__",&__str__<R>);
  vector_class.def("__repr__",&__repr__<R>);

  def("zero_vector",&zero_vector<FloatPy>);
  def("unit_vector",&unit_vector<FloatPy>);

  def("sup_norm", (R(*)(const Vec&)) &sup_norm<R>);
  def("norm", (R(*)(const Vec&)) &sup_norm<R>);
}

template<>
void export_vector<Rational>()
{
  typedef Rational R;
  typedef Vector<R> Vec;
  typedef Covector<R> Cvec;
  typedef Matrix<R> Mx;
  
  class_< Vector<R> > vector_class(python_name<R>("Vector").c_str(),no_init);
  vector_class.def("__init__", make_constructor(&make_vector<R>) );
  vector_class.def(init<int>());
  //vector_class.def(init<std::string>());
  vector_class.def(init<Vec>());
  vector_class.def("__len__", &Vec::size);
  vector_class.def("__getitem__",&__getitem__<Vec>);
  vector_class.def("__setitem__",&__setitem__<Vec,R>);
  vector_class.def("__setitem__",&__setitem__<Vec,double>);
  vector_class.def("__neg__",&__neg__<Vec,Vec>);
  vector_class.def("__add__",&__add__<Vec,Vec,Vec>);
  vector_class.def("__sub__",&__sub__<Vec,Vec,Vec>);
  vector_class.def("__rmul__",&__mul__<Vec,Vec,int,Vec,R>);
  vector_class.def("__rmul__",&__mul__<Vec,Vec,double,Vec,R>);
  vector_class.def("__rmul__",&__mul__<Vec,Vec,R,Vec,R>);
  vector_class.def("__mul__",&__mul__<Vec,Vec,int,Vec,R>);
  vector_class.def("__mul__",&__mul__<Vec,Vec,double,Vec,R>);
  vector_class.def("__mul__",&__mul__<Vec,Vec,R,Vec,R>);
  vector_class.def("__mul__",&__mul__<Mx,Vec,Cvec>);
  vector_class.def("__div__",&__div__<Vec,Vec,int,Vec,R>);
  vector_class.def("__div__",&__div__<Vec,Vec,double,Vec,R>);
  vector_class.def("__div__",&__div__<Vec,Vec,R,Vec,R>);
  vector_class.def("__str__",&__str__<R>);
  vector_class.def("__repr__",&__repr__<R>);

  def("sup_norm", (R(*)(const Vec&)) &sup_norm<R>);
  def("norm", (R(*)(const Vec&)) &sup_norm<R>);
}


template<class R>
void export_interval_vector() {
  typedef Interval<R> I;
  typedef Vector<R> Vec;
  typedef Vector< Interval<R> > IVec;
  typedef Covector<R> Cvec;
  typedef Covector< Interval<R> > ICvec;
  typedef Matrix< Interval<R> > IMx;
  
  class_< Vector<I> > vector_class(python_name<R>("IntervalVector").c_str(),no_init);
  vector_class.def("__init__", make_constructor(&make_vector<I>) );
  vector_class.def(init<int>());;
  //vector_class.def(init<std::string>());;
  vector_class.def(init<Vec>());
  vector_class.def(init<IVec>());
  vector_class.def("__len__", &IVec::size);
  vector_class.def("__getitem__",&__getitem__<IVec>) ;
  vector_class.def("__setitem__",&__setitem__<IVec,I>);
  vector_class.def("__setitem__",&__setitem__<IVec,R>);
  vector_class.def("__setitem__",&__setitem__<IVec,double>);
  vector_class.def("__neg__",&__neg__<IVec,IVec>);
  vector_class.def("__add__",&__add__<IVec,IVec,Vec>);
  vector_class.def("__add__",&__add__<IVec,IVec,IVec>);
  vector_class.def("__sub__",&__sub__<IVec,IVec,Vec>);
  vector_class.def("__sub__",&__sub__<IVec,IVec,IVec>);
  vector_class.def("__rmul__",&__mul__<IVec,IVec,int,IVec,R>);
  vector_class.def("__rmul__",&__mul__<IVec,IVec,double,IVec,R>);
  vector_class.def("__rmul__",&__mul__<IVec,IVec,R>);
  vector_class.def("__rmul__",&__mul__<IVec,IVec,I>);
  vector_class.def("__mul__",&__mul__<IVec,IVec,int,IVec,R>);
  vector_class.def("__mul__",&__mul__<IVec,IVec,double,IVec,R>);
  vector_class.def("__mul__",&__mul__<IVec,IVec,R>);
  vector_class.def("__mul__",&__mul__<IVec,IVec,I>);
  vector_class.def("__mul__",&__mul__<IMx,IVec,Cvec>);
  vector_class.def("__mul__",&__mul__<IMx,IVec,ICvec>);
  vector_class.def("__div__",&__div__<IVec,IVec,int,IVec,R>);
  vector_class.def("__div__",&__div__<IVec,IVec,double,IVec,R>);
  vector_class.def("__div__",&__div__<IVec,IVec,R>);
  vector_class.def("__div__",&__div__<IVec,IVec,I>);
  vector_class.def("__str__",&__str__<I>);
  vector_class.def("__repr__",&__repr__<I>);
  
  def("sup_norm", (I(*)(const IVec&)) &sup_norm<I>);
  def("norm", (I(*)(const IVec&)) &sup_norm<I>);

  def("midpoint",(Vec(*)(const IVec&))&midpoint);
  def("encloses",(bool(*)(const IVec&,const Vec&))&encloses);
  def("refines",(bool(*)(const IVec&,const IVec&))&refines);
}

template void export_vector<FloatPy>();
template void export_vector<Rational>();

template void export_interval_vector<FloatPy>();
