/***************************************************************************
 *            file
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
 
/*! \file orbit.h
 *  \brief Orbits of dynamic systems
 */
#include "array.h"
#include "numeric.h"
#include "vector.h"
#include "dense_differential.h"
#include "sparse_differential.h"
#include "function.h"

#include "approximate_taylor_model.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;

typedef ApproximateTaylorModel ModelType;

template<class X> 
array<X>
extract_array(const boost::python::object& obj)
{
  list elements=extract<list>(obj);
  int n=len(elements);
  array<X> result(n);
  for(int i=0; i!=n; ++i) {
    boost::python::extract<X> xv(elements[i]);
    boost::python::extract<double> dv(elements[i]);
    if(xv.check()) {
      result[i]=static_cast<X>(xv()); 
    } else if(dv.check()) {
      result[i]=static_cast<double>(dv());
    } else {
      result[i]=0;
    }
  }
  return result;
}

template<class X> void read(Vector<X>& v, const object& obj) {
  list elements=extract<list>(obj);
  ARIADNE_ASSERT(v.size()==uint(len(elements)));
  for(uint i=0; i!=v.size(); ++i) { 
    boost::python::extract<X> xv(elements[i]);
    boost::python::extract<double> dv(elements[i]);
    if(xv.check()) {
      v[i]=static_cast<X>(xv()); 
    } else if(dv.check()) {
      v[i]=static_cast<double>(dv());
    } else {
      v[i]=xv();
    }
  }
}


template<class X> 
void read(DifferentialVector< SparseDifferential<X> >& td, const object& obj) {
  list elements=extract<list>(obj);
  ARIADNE_ASSERT(td.result_size()==uint(len(elements)));
  for(uint i=0; i!=td.size(); ++i) { 
    boost::python::extract< SparseDifferential<X> > etv(elements[i]);
    boost::python::extract<double> ed(elements[i]);
    if(etv.check()) {
      td[i]=etv(); 
    } else if(ed.check()) {
      td[i]=ed();
    } else {
      td[i]=etv();
    }
  }
}

class FunctionPyWrap
  : public FunctionInterface
  , public wrapper< FunctionInterface >
{
  virtual FunctionInterface* clone() const { return this->get_override("clone")(); };
  virtual std::string name() const { return this->get_override("name")(); }
  virtual uint result_size() const { return this->get_override("result_size")(); }
  virtual uint argument_size() const { return this->get_override("argument_size")(); }
  virtual ushort smoothness() const { return this->get_override("smoothness")(); }
  virtual Vector<Interval> evaluate(const Vector<Interval>&) const { return this->get_override("evaluate")(); }
  virtual Matrix<Interval> jacobian(const Vector<Interval>&) const { return this->get_override("jacobian")(); }
  virtual DifferentialVector< SparseDifferential<Float> > expansion(const Vector<Float>&, const ushort&) const { return this->get_override("expansion")(); }
  virtual DifferentialVector< SparseDifferential<Interval> > expansion(const Vector<Interval>&, const ushort&) const { return this->get_override("expansion")(); }
  virtual std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};


class PythonFunction
  : public FunctionInterface
{
 public:
  PythonFunction(std::string& nm, uint rs, uint as, const object& pyf) : _name(nm), _result_size(rs), _argument_size(as), _pyf(pyf) { }
  PythonFunction(uint rs, uint as, const object& pyf) : _name(), _result_size(rs), _argument_size(as), _pyf(pyf) { }
  PythonFunction(const object& pyf)
    : _name(), 
      _result_size(extract<int>(pyf.attr("result_size"))), 
      _argument_size(extract<int>(pyf.attr("argument_size"))), 
      _pyf(pyf) { }

  PythonFunction* clone() const { return new PythonFunction(*this); }
  virtual uint result_size() const { return this->_result_size; }
  virtual uint argument_size() const { return this->_argument_size; }
  virtual ushort smoothness() const { return 255; }

  virtual Vector<Interval> evaluate (const Vector<Interval>& x) const { 
    Vector<Interval> r(this->_result_size); 
    read(r,this->_pyf(x)); 
    return r; }
  virtual Matrix<Interval> jacobian (const Vector<Interval>& x) const { 
    DifferentialVector< SparseDifferential<Interval> > rj(this->_result_size,this->_argument_size,1u); 
    DifferentialVector< SparseDifferential<Interval> > aj=DifferentialVector< SparseDifferential<Interval> >::variable(x.size(),x.size(),1u,x); 
    read(rj,this->_pyf(aj)); 
    return rj.jacobian(); }
  virtual DifferentialVector< SparseDifferential<Float> > expansion (const Vector<Float>& x, const ushort& d) const {  
    DifferentialVector< SparseDifferential<Float> > rd(this->_result_size,this->_argument_size,d); 
    DifferentialVector< SparseDifferential<Float> > ad=DifferentialVector< SparseDifferential<Float> >::variable(x.size(),x.size(),d,x); 
    read(rd,this->_pyf(ad)); 
    return rd; }
  virtual DifferentialVector< SparseDifferential<Interval> > expansion (const Vector<Interval>& x, const ushort& d) const {  
    DifferentialVector< SparseDifferential<Interval> > rd(this->_result_size,this->_argument_size,d); 
    DifferentialVector< SparseDifferential<Interval> > ad=DifferentialVector< SparseDifferential<Interval> >::variable(x.size(),x.size(),d,x); 
    read(rd,this->_pyf(ad)); 
    return rd; }

  virtual std::ostream& write(std::ostream& os) const { 
    os << "Function( ";
    if(this->_name.size()>0) { os << "name=" << this->_name << ", "; }
    os << "result_size="<<this->_result_size;
    os << ", argument_size="<<this->_argument_size;
    return os << " )"; }
 private:
  std::string _name;
  uint _result_size;
  uint _argument_size;
  boost::python::object _pyf;
};


typedef Vector<Float> FV;
typedef Vector<Interval> IV;
typedef DifferentialVector< SparseDifferential<Float> > FSDV;
typedef DifferentialVector< SparseDifferential<Interval> > ISDV;

void export_function_interface() 
{
  class_<FunctionPyWrap, boost::noncopyable> function_interface_class("FunctionInterface");
  function_interface_class.def("argument_size", pure_virtual(&FunctionInterface::argument_size));
  function_interface_class.def("result_size", pure_virtual(&FunctionInterface::result_size));
  function_interface_class.def("smoothness", pure_virtual(&FunctionInterface::smoothness));
  function_interface_class.def("evaluate",pure_virtual(&FunctionInterface::evaluate));
  function_interface_class.def("jacobian",pure_virtual(&FunctionInterface::jacobian));
  function_interface_class.def("expansion",pure_virtual((ISDV(FunctionInterface::*)(const IV&,const ushort&)const)&FunctionInterface::expansion));
}


void export_python_function() 
{
  class_<PythonFunction, bases< FunctionInterface > > python_function_class("AriadneFunction", init<object>());
  python_function_class.def(init<uint,uint,object>());
  python_function_class.def("argument_size", &PythonFunction::argument_size);
  python_function_class.def("result_size", &PythonFunction::result_size);
  python_function_class.def("smoothness", &PythonFunction::smoothness);
  python_function_class.def("__call__",&PythonFunction::evaluate);
  python_function_class.def("evaluate",&PythonFunction::evaluate);
  python_function_class.def("jacobian",&PythonFunction::jacobian);
  python_function_class.def("expansion",(FSDV(PythonFunction::*)(const FV&,const ushort&)const)&PythonFunction::expansion);
  python_function_class.def("expansion",(ISDV(PythonFunction::*)(const IV&,const ushort&)const)&PythonFunction::expansion);
  python_function_class.def(self_ns::str(self));
}

void export_model() 
{
  typedef ApproximateTaylorModel Model;
  typedef Float R;
  typedef Interval I;
  typedef Vector<I> Box;
  typedef Vector<R> Point;
  typedef DifferentialVector< SparseDifferential<R> > Expansion;

  class_< Model > function_model_class("ApproximateTaylorModel",init< Vector<I>, Vector<R>, const FunctionInterface&, ushort, ushort>());
  function_model_class.def(init< uint, uint, ushort, ushort >());
  //function_model_class.def(init< Vector<I> >());
  function_model_class.def(init< Vector<I>, Vector<R>, DifferentialVector< SparseDifferential<R> > >());
  function_model_class.def(init< Vector<I>, Vector<R>, DifferentialVector< SparseDifferential<R> > >());
  function_model_class.def(init< Model >());
  function_model_class.def("result_size", &Model::result_size);
  function_model_class.def("argument_size", &Model::argument_size);
  function_model_class.def("order", &Model::order);
  function_model_class.def("smoothness", &Model::smoothness);
  function_model_class.def("domain", &Model::domain);
  function_model_class.def("centre", &Model::centre);
  function_model_class.def("range", &Model::range);
  function_model_class.def("expansion", &Model::expansion, return_value_policy<copy_const_reference>());
  function_model_class.def("identity", &Model::identity);
  function_model_class.staticmethod("identity");
  //function_model_class.def("set",(void(Model::*)(uint,const MultiIndex&,const R&)) &Model::set);
  //function_model_class.def("get",(R(Model::*)(uint,const MultiIndex&)const) &Model::get);
  //function_model_class.def("truncate",&Model::truncate);
  function_model_class.def("evaluate",(Vector<I>(Model::*)(const Vector<I>&)const) &Model::evaluate);
  function_model_class.def("jacobian",(Matrix<I>(Model::*)(const Vector<I>&)const) &Model::jacobian);
  function_model_class.def("__add__",(Model(*)(const Model&,const Model&)) &Ariadne::operator+);
  function_model_class.def("__sub__",(Model(*)(const Model&,const Model&)) &Ariadne::operator-);
  function_model_class.def(self_ns::str(self));
 
  def("evaluate",(Vector<I>(Model::*)(const Vector<I>&)const) &Model::evaluate);
  def("project",(Model(*)(const Model&,const Slice&)) &project);
  def("join",(Model(*)(const Model&,const Model&)) &join);
  def("compose",(Model(*)(const Model&,const Model&)) &compose);
  //def("inverse",(Model(*)(const Model&)) &inverse);
  def("implicit",(Model(*)(const Model&)) &implicit);
  //def("derivative",(Model(*)(const Model&, uint)) &derivative);
  def("flow",(Model(*)(const Model&)) &flow);
  //def("integrate",(Model(*)(const Model&,const R&)) &integrate);
  def("hitting",(Model(*)(const Model&,const Model&)) &hitting);
  //def("disjoint",(Model(*)(const Model&,const Vector<I>&)) &disjoint);
  def("solve",(Vector<I>(*)(const Model&,const Vector<R>&)) &solve);
}


BOOST_PYTHON_MODULE(function)
{
  export_function_interface();
  export_python_function();
  export_model();
}
