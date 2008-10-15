/***************************************************************************
 *            differentiation_submodule.cc
 *
 *  Copyright 2008  Pieter Collins
 * 
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
 
#include "array.h"
#include "numeric.h"
#include "dense_differential.h"
#include "sparse_differential.h"
#include "differential_vector.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;

template<class X> void read_array(array<X>&, const boost::python::object& obj) { }

template<class DIFF>
DIFF*
make_differential(const uint& as, const uint& d, const boost::python::object& obj) 
{
  typedef typename DIFF::ScalarType X;
  DIFF* result=new DIFF(as,d);
  array<X> data;
  read_array(data,obj);
  assert(data.size()==compute_polynomial_data_size(1u,as,d));
  MultiIndex i(as);
  const X* ptr=data.begin();
  while(i.degree()<=d) {
    result[i]=*ptr; ++i; ++ptr;
  }
  return result;
}


template<class DIFF>
DifferentialVector<DIFF>*
make_differential_vector(const uint& rs, const uint& as, const uint& d, const boost::python::object& obj) 
{
  typedef typename DIFF::ScalarType X;
  array<X> data;
  read_array(data,obj);
  ARIADNE_ASSERT(data.size()==compute_polynomial_data_size(rs,as,d));
  DifferentialVector<DIFF>* result=new DifferentialVector<DIFF>(rs,as,d,data.begin());
  return result;
}


template<class C, class I, class X> inline 
void set_item(C& c, const I& i, const X& x) {
  c[i]=x;
}

template<class C, class I, class J, class X> inline 
void matrix_set_item(C& c, const I& i, const J& j, const X& x) {
  c[i][j]=x;
}


template<class C, class I> inline 
typename C::ValueType 
get_item(const C& c, const I& i) {
  return c[i];
}

template<class C, class I, class J> inline 
typename C::ValueType::ValueType 
matrix_get_item(const C& c, const I& i, const J& j) {
  return c[i][j];
}


template<class DIFF>
void export_differential() 
{
  typedef typename DIFF::ScalarType X;
  typedef Vector<X> V;
  typedef Series<X> S;
  typedef DIFF D;
  typedef DifferentialVector<D> DV;


  class_<D> differential_class("DenseDifferential");
  //differential_class.def("__init__", make_constructor(&make_differential<X>) );
  differential_class.def( init< uint, uint >());
  differential_class.def("value", (const X&(D::*)()const) &D::value, return_value_policy<copy_const_reference>());
  differential_class.def("__getitem__", &get_item<D,MultiIndex>);
  differential_class.def("__setitem__",&set_item<D,MultiIndex,double>);
  differential_class.def("__setitem__",&set_item<D,MultiIndex,X>);
  differential_class.def(-self);
  differential_class.def(self+self);
  differential_class.def(self-self);
  differential_class.def(self*self);
  differential_class.def(self/self);
  differential_class.def(self+X());
  differential_class.def(self-X());
  differential_class.def(self+=X());
  differential_class.def(self-=X());
  differential_class.def(self*=X());
  differential_class.def(self/=X());
  differential_class.def(self_ns::str(self));
  
  def("constant",(D(*)(uint, ushort, const X&))&D::constant);
  def("variable",(D(*)(uint, ushort, const X&, uint))&D::variable);

  def("max",(D(*)(const D&,const D&))&max<X>);
  def("min",(D(*)(const D&,const D&))&min<X>);
  def("abs",(D(*)(const D&))&abs<X>);
  def("pos",(D(*)(const D&))&pos<X>);
  def("neg",(D(*)(const D&))&neg<X>);
  def("rec",(D(*)(const D&))&rec<X>);
  def("pow",(D(*)(const D&, int))&pow<X>);

  def("sqrt", (D(*)(const D&))&sqrt<X>);
  def("exp", (D(*)(const D&))&exp<X>);
  def("log", (D(*)(const D&))&log<X>);
  /*
  def("sin", (D(*)(const D&))&sin<X>);
  def("cos", (D(*)(const D&))&cos<X>);
  def("tan", (D(*)(const D&))&tan<X>);
  def("asin", (D(*)(const D&))&asin<X>);
  def("acos", (D(*)(const D&))&acos<X>);
  def("atan", (D(*)(const D&))&atan<X>);
  */
}

template<class DIFF> 
void
export_differential_vector()
{
  typedef typename DIFF::ScalarType X;
  typedef Vector<X> V;
  typedef Series<X> S;
  typedef DIFF D;
  typedef DifferentialVector<D> DV;

  class_<DV> differential_vector_class("DifferentialVector");
  differential_vector_class.def("__init__", make_constructor(&make_differential_vector<D>) );
  differential_vector_class.def( init< uint, uint, uint >());
  differential_vector_class.def("__getitem__", &matrix_get_item<DV,int,MultiIndex>);
  differential_vector_class.def("__setitem__",&matrix_set_item<DV,int,MultiIndex,double>);
  differential_vector_class.def("__setitem__",&set_item<DV,int,X>);
  differential_vector_class.def("__getitem__", &get_item<DV,int>);
  differential_vector_class.def("__setitem__",&set_item<DV,int,D>);
  differential_vector_class.def(-self);
  differential_vector_class.def(self+self);
  differential_vector_class.def(self-self);
  differential_vector_class.def(self+V());
  differential_vector_class.def(self-V());
  differential_vector_class.def(self*X());
  differential_vector_class.def(self+=V());
  differential_vector_class.def(self-=V());
  differential_vector_class.def(self*=X());
  differential_vector_class.def("value", &DV::value);
  differential_vector_class.def("jacobian", &DV::jacobian);
  differential_vector_class.def(self_ns::str(self));

  //def("variable",(DV(*)(const Vector<Float>&,ushort)) &DV::variable);
  //def("variable",(DV(*)(const Vector<Interval>&,ushort)) &DV::variable);

  def("evaluate",(V(*)(const DV&,const V&))&evaluate);
  //def("compose",(D(*)(const D&,const DV&))&compose);
  def("compose",(DV(*)(const DV&,const DV&))&compose);
  def("translate",(DV(*)(const DV&,const V&))&translate);
  def("inverse",(DV(*)(const DV&))&inverse);
  def("implicit",(DV(*)(const DV&))&implicit);

}

template void export_differential< DenseDifferential<Float> >();
template void export_differential< DenseDifferential<Interval> >();
template void export_differential< SparseDifferential<Float> >();
template void export_differential< SparseDifferential<Interval> >();

template void export_differential_vector< DenseDifferential<Float> >();
template void export_differential_vector< DenseDifferential<Interval> >();
template void export_differential_vector< SparseDifferential<Float> >();
template void export_differential_vector< SparseDifferential<Interval> >();

void differentiation_submodule() 
{
  export_differential< DenseDifferential<Float> >();
  export_differential< DenseDifferential<Interval> >();
  export_differential< SparseDifferential<Float> >();
  export_differential< SparseDifferential<Interval> >();

  export_differential_vector< DenseDifferential<Float> >();
  export_differential_vector< DenseDifferential<Interval> >();
  export_differential_vector< SparseDifferential<Float> >();
  export_differential_vector< SparseDifferential<Interval> >();
}



 /*
BOOST_PYTHON_MODULE(differentiation)
{
  export_differential< DenseDifferential<Float> >();
  export_differential< DenseDifferential<Interval> >();
  export_differential< SparseDifferential<Float> >();
  export_differential< SparseDifferential<Interval> >();

  export_differential_vector< DenseDifferential<Float> >();
  export_differential_vector< DenseDifferential<Interval> >();
  export_differential_vector< SparseDifferential<Float> >();
  export_differential_vector< SparseDifferential<Interval> >();
}
 */
