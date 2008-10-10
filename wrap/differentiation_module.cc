#include "array.h"
#include "numeric.h"
#include "dense_differential.h"
#include "differential_vector.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;


template<class X> void read_array(array<X>&, const boost::python::object& obj) { }

template<class X>
DenseDifferential<X>*
make_differential(const uint& as, const uint& d, const boost::python::object& obj) 
{
  DenseDifferential<X>* result=new DenseDifferential<X>(as,d);
  read_array(result->data(),obj);
  assert(result->data().size()==compute_polynomial_data_size(1u,as,d));
  return result;
}


template<class R1, class R2> inline 
void differential_set_item(DenseDifferential<R1>& td, const MultiIndex& i, R2 x) {
  assert(i.number_of_variables()==td.argument_size()); 
  assert(i.degree()<=td.degree()); 
  td[i]=x;
}

template<class R> inline 
R differential_get_item(const DenseDifferential<R>& td, const MultiIndex& i) {
  assert(i.number_of_variables()==td.argument_size()); 
  assert(i.degree()<=td.degree()); 
  return td[i];
}

template<class X, class XX> inline 
DenseDifferential<X> differential_variable(uint as, ushort d, const XX& x, uint i) {
  return DenseDifferential<X>::variable(as,d,X(x),i);
}


template<class X>
DenseDifferentialVector<X>*
make_differential_vector(const uint& rs, const uint& as, const uint& d, const boost::python::object& obj) 
{
  array<X> data;
  read_array(data,obj);
  ARIADNE_ASSERT(data.size()==compute_polynomial_data_size(rs,as,d));
  DenseDifferentialVector<X>* result=new DenseDifferentialVector<X>(rs,as,d,data.begin());
  return result;
}

template<class X1, class X2> inline 
void differential_vector_set_variable(DenseDifferentialVector<X1>& td, const uint& i, DenseDifferential<X2> x) {
  ARIADNE_ASSERT(i<td.result_size()); 
  ARIADNE_ASSERT(x.argument_size()==td.argument_size()); 
  ARIADNE_ASSERT(x.degree()<=td.degree()); 
  td[i]=x;
}

template<class X> inline 
DenseDifferential<X> differential_vector_get_variable(const DenseDifferentialVector<X>& td, const uint& i) {
  ARIADNE_ASSERT(i<td.result_size()); 
  return td[i];
}

template<class X, class XX> inline 
void differential_vector_set_item(DenseDifferentialVector<X>& td, const uint& i, const MultiIndex& j, const XX& x) {
  ARIADNE_ASSERT(i<td.result_size()); 
  ARIADNE_ASSERT(j.number_of_variables()==td.argument_size()); 
  ARIADNE_ASSERT(j.degree()<=td.degree()); 
  td[i][j]=x;
  //td.set(i,j,x);
}

template<class X> inline 
X differential_vector_get_item(const DenseDifferentialVector<X>& td, const uint& i, const MultiIndex& j) {
  ARIADNE_ASSERT(i==td.result_size()); 
  ARIADNE_ASSERT(j.number_of_variables()==td.argument_size()); 
  ARIADNE_ASSERT(j.degree()<=td.degree()); 
  return td[i][j];
  //return td.get(i,j);
}




template<class X>
void export_differential() 
{
  typedef Vector<X> V;
  typedef Series<X> S;
  typedef DenseDifferential<X> D;
  //typedef DifferentialVector< DenseDifferential<X> > DV;
  typedef DenseDifferentialVector<X> DV;


  class_<D> differential_class("DenseDifferential");
  //differential_class.def("__init__", make_constructor(&make_differential<X>) );
  differential_class.def( init< uint, uint >());
  differential_class.def("value", (const X&(D::*)()const) &D::value, return_value_policy<copy_const_reference>());
  differential_class.def("__getitem__", &differential_get_item<X>);
  differential_class.def("__setitem__",&differential_set_item<X,double>);
  differential_class.def("__setitem__",&differential_set_item<X,X>);
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
  def("variable",(D(*)(uint, ushort, const X&, uint))&differential_variable<X,X>);

  def("max",(D(*)(const D&,const D&))&max<X>);
  def("min",(D(*)(const D&,const D&))&min<X>);
  def("abs",(D(*)(const D&))&abs<X>);
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




  class_<DV> differential_vector_class("DifferentialVector");
  differential_vector_class.def("__init__", make_constructor(&make_differential_vector<X>) );
  differential_vector_class.def( init< uint, uint, uint >());
  differential_vector_class.def("__getitem__", &differential_vector_get_item<X>);
  differential_vector_class.def("__setitem__",&differential_vector_set_item<X,double>);
  differential_vector_class.def("__setitem__",&differential_vector_set_item<X,X>);
  differential_vector_class.def("__getitem__", &differential_vector_get_variable<X>);
  differential_vector_class.def("__setitem__",&differential_vector_set_variable<X,X>);
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

BOOST_PYTHON_MODULE(differentiation)
{
  export_differential<Float>();
  export_differential<Interval>();
}
