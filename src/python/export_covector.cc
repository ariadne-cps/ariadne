/***************************************************************************
 *            python/export_covector.cc
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
#include "linear_algebra/covector.h"
#include "linear_algebra/matrix.h"

#include "python/name.h"
#include "python/subscripting.h"
#include "python/operators.h"
#include "python/float.h"
#include "python/read_array.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
using namespace boost::python;

template<class X>
std::string
__str__(const Covector<X>& cv)
{
  std::stringstream ss;
  ss << cv;
  return ss.str();
}

template<class X>
std::string
__repr__(const Covector<X>& cv)
{
  std::stringstream ss;
  ss << "Covector(" << cv << ")";
  return ss.str();
}

template<class X>
Covector<X>*
make_covector(const boost::python::object& obj) 
{
  Covector<X>* cv=new Covector<X>;
  read_array(cv->data(),obj);
  return cv;
}


template<class R> 
inline
R
covector_get_item(const Covector<R>& v, int n) {
  if(n<0) {
    n+=v.size();
  }
  assert(0<=n);
  size_t m=size_t(n);
  assert(m<v.size());
  return v(m);
}

template<class R, class A> 
inline
void
covector_set_item(Covector<R>& v, int n, const A& x) {
  if(n<0) {
    n+=v.size();
  }
  assert(0<=n);
  size_t m=size_t(n);
  assert(m<v.size());
  R& r=v(m);
  r=R(x);
}



template<class R>
void 
export_covector()
{
  typedef Interval<R> I;
  typedef Covector<R> Cvec;
  typedef Vector<R> Vec;
  typedef Matrix<R> Mx;
  typedef Vector<I> IVec;
  typedef Covector<I> ICvec;
  typedef Matrix<I> IMx;
  
  class_< Covector<R> > covector_class(python_name<R>("Covector").c_str(),no_init);
  covector_class.def("__init__", make_constructor(&make_covector<R>));
  covector_class.def(init<int>());
  covector_class.def(init<uint>());
  covector_class.def(init<Cvec>());
  covector_class.def("__len__", &Cvec::size);
  covector_class.def("__getitem__",&covector_get_item<R>);
  covector_class.def("__setitem__",&covector_set_item<R,R>);
  covector_class.def("__setitem__",&covector_set_item<R,double>);
  covector_class.def("__neg__",&__neg__<Cvec,Cvec>);
  covector_class.def("__add__",&__add__<ICvec,Cvec,Cvec>);
  covector_class.def("__add__",&__add__<ICvec,Cvec,ICvec>);
  covector_class.def("__sub__",&__sub__<ICvec,Cvec,Cvec>);
  covector_class.def("__sub__",&__sub__<ICvec,Cvec,ICvec>);
  covector_class.def("__rmul__",&__mul__<ICvec,Cvec,int,Cvec,R>);
  covector_class.def("__rmul__",&__mul__<ICvec,Cvec,double,Cvec,R>);
  covector_class.def("__rmul__",&__mul__<ICvec,Cvec,R>);
  covector_class.def("__rmul__",&__mul__<ICvec,Cvec,I>);
  //covector_class.def("__rmul__",&__mul__<IMx,Vec,Cvec>);
  //covector_class.def("__rmul__",&__mul__<IMx,IVec,Cvec>);
  covector_class.def("__mul__",&__mul__<ICvec,Cvec,int,Cvec,R>);
  covector_class.def("__mul__",&__mul__<ICvec,Cvec,double,Cvec,R>);
  covector_class.def("__mul__",&__mul__<ICvec,Cvec,R>);
  covector_class.def("__mul__",&__mul__<ICvec,Cvec,I>);
  covector_class.def("__mul__",&__mul__<I,Cvec,Vec>);
  covector_class.def("__mul__",&__mul__<I,Cvec,IVec>);
  covector_class.def("__mul__",&__mul__<ICvec,Cvec,IMx>);
  covector_class.def("__mul__",&__mul__<ICvec,Cvec,Mx>);
  covector_class.def("__div__",&__div__<ICvec,Cvec,int,Cvec,R>);
  covector_class.def("__div__",&__div__<ICvec,Cvec,double,Cvec,R>);
  covector_class.def("__div__",&__div__<ICvec,Cvec,R>);
  covector_class.def("__div__",&__div__<ICvec,Cvec,I>);
  covector_class.def("__str__",&__str__<R>);
  covector_class.def("__repr__",&__repr__<R>);

  def("zero_covector",&zero_covector<FloatPy>);
  def("unit_covector",&unit_covector<FloatPy>);

  //def("sup_norm",&sup_norm<FloatPy>);
  //def("one_norm",&one_norm<FloatPy>);
}

template<>
void 
export_covector<Rational>()
{
  typedef Rational R;
  typedef Vector<R> Vec;
  typedef Matrix<R> Mx;
  typedef Covector<R> Cvec;
  
  class_< Covector<R> > covector_class(python_name<R>("Covector").c_str(),no_init);
  covector_class.def("__init__", make_constructor(&make_covector<R>) );
  covector_class.def(init<int>());
  covector_class.def(init<Cvec>());
  covector_class.def("__len__", &Cvec::size);
  covector_class.def("__getitem__",&covector_get_item<R>);
  covector_class.def("__setitem__",&covector_set_item<R,R>);
  covector_class.def("__setitem__",&covector_set_item<R,double>);
  covector_class.def("__neg__",&__neg__<Cvec,Cvec>);
  covector_class.def("__add__",&__add__<Cvec,Cvec,Cvec>);
  covector_class.def("__sub__",&__sub__<Cvec,Cvec,Cvec>);
  covector_class.def("__rmul__",&__mul__<Cvec,Cvec,int,Cvec,R>);
  covector_class.def("__rmul__",&__mul__<Cvec,Cvec,double,Cvec,R>);
  covector_class.def("__rmul__",&__mul__<Cvec,Cvec,R,Cvec,R>);
  covector_class.def("__rmul__",&__mul__<Mx,Vec,Cvec>);
  covector_class.def("__mul__",&__mul__<Cvec,Cvec,int,Cvec,R>);
  covector_class.def("__mul__",&__mul__<Cvec,Cvec,double,Cvec,R>);
  covector_class.def("__mul__",&__mul__<Cvec,Cvec,R,Cvec,R>);
  covector_class.def("__mul__",&__mul__<R,Cvec,Vec>);
  covector_class.def("__mul__",&__mul__<Cvec,Cvec,Mx>);
  covector_class.def("__div__",&__div__<Cvec,Cvec,int,Cvec,R>);
  covector_class.def("__div__",&__div__<Cvec,Cvec,double,Cvec,R>);
  covector_class.def("__div__",&__div__<Cvec,Cvec,R,Cvec,R>);
  covector_class.def("__str__",&__str__<R>);
  covector_class.def("__repr__",&__repr__<R>);

  def("zero_covector",&zero_covector<Rational>);
  def("unit_covector",&unit_covector<Rational>);

  //def("sup_norm",&sup_norm<FloatPy>);
  //def("one_norm",&one_norm<FloatPy>);

}


template<class R>
void 
export_interval_covector() 
{
  typedef Interval<R> I;
  typedef Vector<R> Vec;
  typedef Matrix<R> Mx;
  typedef Covector<R> Cvec;
  typedef Vector<I> IVec;
  typedef Matrix<I> IMx;
  typedef Covector<I> ICvec;
  
  class_< Covector<I> > covector_class(python_name<R>("IntervalCovector").c_str(),no_init);
  covector_class.def("__init__", make_constructor(&make_covector<I>) );
  covector_class.def(init<int>());
  covector_class.def(init<Cvec>());
  covector_class.def(init<ICvec>());
  covector_class.def("__len__", &ICvec::size);
  covector_class.def("__getitem__",&covector_get_item<I>);
  covector_class.def("__setitem__",&covector_set_item<I,I>);
  covector_class.def("__setitem__",&covector_set_item<I,R>);
  covector_class.def("__setitem__",&covector_set_item<I,double>);
  covector_class.def("__neg__",&__neg__<ICvec,ICvec>);
  covector_class.def("__add__",&__add__<ICvec,ICvec,Cvec>);
  covector_class.def("__add__",&__add__<ICvec,ICvec,ICvec>);
  covector_class.def("__sub__",&__sub__<ICvec,ICvec,Cvec>);
  covector_class.def("__sub__",&__sub__<ICvec,ICvec,ICvec>);
  covector_class.def("__rmul__",&__mul__<ICvec,ICvec,int,ICvec,R>);
  covector_class.def("__rmul__",&__mul__<ICvec,ICvec,double,ICvec,R>);
  covector_class.def("__rmul__",&__mul__<ICvec,ICvec,R>);
  covector_class.def("__rmul__",&__mul__<ICvec,ICvec,I>);
  covector_class.def("__rmul__",&__mul__<IMx,Vec,ICvec>);
  covector_class.def("__rmul__",&__mul__<IMx,IVec,ICvec>);
  covector_class.def("__mul__",&__mul__<ICvec,ICvec,int,ICvec,R>);
  covector_class.def("__mul__",&__mul__<ICvec,ICvec,double,ICvec,R>);
  covector_class.def("__mul__",&__mul__<ICvec,ICvec,R>);
  covector_class.def("__mul__",&__mul__<ICvec,ICvec,I>);
  covector_class.def("__mul__",&__mul__<I,ICvec,Vec>);
  covector_class.def("__mul__",&__mul__<ICvec,ICvec,IMx>);
  covector_class.def("__mul__",&__mul__<I,ICvec,IVec>);
  covector_class.def("__mul__",&__mul__<ICvec,ICvec,Mx>);
  covector_class.def("__div__",&__div__<ICvec,ICvec,int,ICvec,R>);
  covector_class.def("__div__",&__div__<ICvec,ICvec,double,ICvec,R>);
  covector_class.def("__div__",&__div__<ICvec,ICvec,R>);
  covector_class.def("__div__",&__div__<ICvec,ICvec,I>);
  covector_class.def("__str__",&__str__<I>);
  covector_class.def("__repr__",&__repr__<I>);
 
  //def("sup_norm",&sup_norm<FloatPy>);
  //def("one_norm",&one_norm<FloatPy>);

  //def("midpoint",(Cvec(*)(const ICvec&))&midpoint);
  //def("encloses",(bool(*)(const ICvec&,const Cvec&))&encloses);
  //def("refines",(bool(*)(const ICvec&,const ICvec&))&refines);
}

template void export_covector<FloatPy>();
template void export_covector<Rational>();

template void export_interval_covector<FloatPy>();
