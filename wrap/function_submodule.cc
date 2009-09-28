/***************************************************************************
 *            function_submodule.cc
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

#include <iostream>
#include <iomanip>
#include "array.h"
#include "container.h"
#include "numeric.h"
#include "real.h"
#include "vector.h"
#include "expansion.h"
#include "multi_index.h"
#include "taylor_model.h"
#include "differential.h"
#include "polynomial.h"
#include "affine.h"
#include "taylor_function.h"
#include "function.h"
#include "expression.h"
#include "formula.h"

#include "utilities.h"


#include <boost/python.hpp>
using namespace boost::python;

namespace Ariadne {

template<>
struct from_python< MultiIndex >
{
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<MultiIndex>()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::tuple tup=extract<boost::python::tuple>(obj_ptr);
        void* storage = ((converter::rvalue_from_python_storage<MultiIndex>*)   data)->storage.bytes;
        MultiIndex res(len(tup));
        for(uint i=0; i!=res.size(); ++i) { res.set(i,extract<uint>(tup[i])); }
        new (storage) MultiIndex(res);
        data->convertible = storage;
    }
};




template<>
struct from_python<VectorFunction>
{
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<VectorFunction>()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        list lst=extract<list>(obj_ptr);
        void* storage = ((converter::rvalue_from_python_storage< VectorFunction >*)   data)->storage.bytes;
        VectorFunction res(len(lst),0);
        for(uint i=0; i!=res.result_size(); ++i) { res.set(i,ScalarFunction(extract<ScalarFunction&>(lst[i]))); }
        new (storage) VectorFunction(res);
        data->convertible = storage;
    }
};




template<class X, class D>
inline Matrix<X> get_jacobian(const Vector<D>& d) {
    const uint rs=d.size(); const uint as=d[0].argument_size();
    Matrix<X> J(rs,as);
    for(uint i=0; i!=rs; ++i) {
        for(uint j=0; j!=as; ++j) {
            J[i][j]=d[i][j];
        }
    }
    return J;
}



class ScalarPythonFunction
    : public ScalarFunctionInterface
{
  public:
    ScalarPythonFunction(std::string& nm, uint as, const object& pyf) : _name(nm), _argument_size(as), _pyf(pyf) { }
    ScalarPythonFunction(uint as, const object& pyf) : _name(),  _argument_size(as), _pyf(pyf) { }
    ScalarPythonFunction(const object& pyf)
        : _name(),
          _argument_size(extract<int>(pyf.attr("argument_size"))),
          _pyf(pyf) { }

    ScalarPythonFunction* clone() const { return new ScalarPythonFunction(*this); }
    virtual uint argument_size() const { return this->_argument_size; }
    virtual ushort smoothness() const { return 255; }

    virtual Float evaluate (const Vector<Float>& x) const {
        return boost::python::extract<Float>(this->_pyf(x)); }
    virtual Interval evaluate (const Vector<Interval>& x) const {
        return boost::python::extract<Interval>(this->_pyf(x)); }
    virtual TaylorModel evaluate (const Vector<TaylorModel>& x) const {
        return boost::python::extract<TaylorModel>(this->_pyf(x)); }
    virtual Differential<Float> evaluate (const Vector< Differential<Float> >& x) const {
        return boost::python::extract< Differential<Float> >(this->_pyf(x)); }
    virtual Differential<Interval> evaluate (const Vector< Differential<Interval> >& x) const {
        return boost::python::extract< Differential<Interval> >(this->_pyf(x)); }
    virtual Differential<TaylorModel> evaluate (const Vector< Differential<TaylorModel> >& x) const {
        return boost::python::extract< Differential<TaylorModel> >(this->_pyf(x)); }

    virtual Vector<Float> gradient(const Vector<Float>& x) const {
        return this->evaluate(Differential<Float>::variables(1u,x)).gradient(); }
    virtual Vector<Interval> gradient(const Vector<Interval>& x) const {
        return this->evaluate(Differential<Interval>::variables(1u,x)).gradient(); }

    virtual ScalarPythonFunction* derivative (uint j) const {
        ARIADNE_ASSERT_MSG(false,"Cannot differentiate a Python function"); return 0; }
    virtual std::ostream& write(std::ostream& os) const {
        os << "ScalarUserFunction( ";
        if(this->_name.size()>0) { os << "name=" << this->_name << ", "; }
        os << "argument_size="<<this->_argument_size;
        return os << " )"; }
  private:
    std::string _name;
    uint _argument_size;
    boost::python::object _pyf;
};


class VectorPythonFunction
    : public VectorFunctionInterface
{
  public:
    VectorPythonFunction(std::string& nm, uint rs, uint as, const object& pyf) : _name(nm), _result_size(rs), _argument_size(as), _pyf(pyf) { }
    VectorPythonFunction(uint rs, uint as, const object& pyf) : _name(), _result_size(rs), _argument_size(as), _pyf(pyf) { }
    VectorPythonFunction(const object& pyf)
        : _name(),
          _result_size(extract<int>(pyf.attr("result_size"))),
          _argument_size(extract<int>(pyf.attr("argument_size"))),
          _pyf(pyf) { }

    VectorPythonFunction* clone() const { return new VectorPythonFunction(*this); }
    virtual uint result_size() const { return this->_result_size; }
    virtual uint argument_size() const { return this->_argument_size; }
    virtual ushort smoothness() const { return 255; }

    virtual Vector<Float> evaluate (const Vector<Float>& x) const {
        return boost::python::extract< Vector<Float> >(this->_pyf(x)); }
    virtual Vector<Interval> evaluate (const Vector<Interval>& x) const {
        return boost::python::extract< Vector<Interval> >(this->_pyf(x)); }
    virtual Vector<TaylorModel> evaluate (const Vector<TaylorModel>& x) const {
        return boost::python::extract< Vector<TaylorModel> >(this->_pyf(x)); }
    virtual Vector< Differential<Float> > evaluate (const Vector< Differential<Float> >& x) const {
        return boost::python::extract< Vector< Differential<Float> > >(this->_pyf(x)); }
    virtual Vector< Differential<Interval> > evaluate (const Vector< Differential<Interval> >& x) const {
        return boost::python::extract< Vector< Differential<Interval> > >(this->_pyf(x)); }

    virtual Matrix<Float> jacobian (const Vector<Float>& x) const {
        return this->evaluate(Differential<Float>::variables(1u,x)).jacobian(); }
    virtual Matrix<Interval> jacobian (const Vector<Interval>& x) const {
        return this->evaluate(Differential<Interval>::variables(1u,x)).jacobian(); }


    virtual std::ostream& write(std::ostream& os) const {
        os << "VectorUserFunction( ";
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



}

using namespace Ariadne;

typedef Float F;
typedef Interval I;
typedef Vector<Float> FV;
typedef Vector<Interval> IV;
typedef Matrix<Float> FMx;
typedef Matrix<Interval> IMx;
typedef Vector< Differential<Float> > FSDV;
typedef Vector< Differential<Interval> > ISDV;
typedef Vector<TaylorModel> TMV;
typedef VectorTaylorFunction TFM;
typedef TaylorModel TM;





void export_multi_index()
{
    class_< MultiIndex > multi_index_class("MultiIndex", init<uint>());
    multi_index_class.def(init<MultiIndex>());
    multi_index_class.def("__getitem__",&MultiIndex::get);
    multi_index_class.def("__setitem__",&MultiIndex::set);
    multi_index_class.def("degree",&MultiIndex::degree);
    multi_index_class.def(self_ns::str(self));

    from_python<MultiIndex>();
}

template<class X>
void export_monomial()
{
    typedef ExpansionValue<X> M;
    class_< M > monomial_class(python_name<X>("Monomial"), init<MultiIndex,X>());
    monomial_class.def("key",(const MultiIndex&(M::*)()const)&M::key,return_value_policy<copy_const_reference>());
    monomial_class.def("data",(const X&(M::*)()const) &M::data,return_value_policy<copy_const_reference>());
    monomial_class.def(self_ns::str(self));
}

template<class X>
void export_polynomial()
{
    typedef Polynomial<X> P;
    X real;

    class_< P > polynomial_class(python_name<X>("Polynomial"), init<int>());
    polynomial_class.def("constant", (P(*)(uint,double)) &P::constant);
    polynomial_class.staticmethod("constant");
    polynomial_class.def("variable", (P(*)(uint,uint)) &P::variable);
    polynomial_class.staticmethod("variable");

    polynomial_class.def("argument_size", &P::argument_size);
    polynomial_class.def("insert", &P::insert);
    polynomial_class.def(+self);
    polynomial_class.def(-self);
    polynomial_class.def(self+self);
    polynomial_class.def(self-self);
    polynomial_class.def(self*self);
    polynomial_class.def(self+real);
    polynomial_class.def(self-real);
    polynomial_class.def(self*real);
    polynomial_class.def(self/real);
    polynomial_class.def(real+self);
    polynomial_class.def(real-self);
    polynomial_class.def(real*self);
    polynomial_class.def("__iter__",boost::python::iterator<P>());
    polynomial_class.def(self_ns::str(self));
    //polynomial_class.def(self_ns::repr(self));

    to_python< Vector< Polynomial<X> > >();
}

void export_scalar_function()
{
    class_<ScalarFunction>
        scalar_function_class("ScalarFunction", init<ScalarFunction>());
    scalar_function_class.def(init<RealExpression,RealSpace>());
    scalar_function_class.def(init< Polynomial<Real> >());
    scalar_function_class.def("argument_size", &ScalarFunction::argument_size);
    scalar_function_class.def("derivative", &ScalarFunction::derivative);
    scalar_function_class.def("polynomial", &ScalarFunction::polynomial);
    scalar_function_class.def("__call__", (Interval(ScalarFunction::*)(const Vector<Interval>&)const)&ScalarFunction::operator() );
    scalar_function_class.def("__pos__", &__pos__<ScalarFunction,ScalarFunction>);
    scalar_function_class.def("__neg__", &__neg__<ScalarFunction,ScalarFunction>);
    scalar_function_class.def("__add__", &__add__<ScalarFunction,ScalarFunction,ScalarFunction>);
    scalar_function_class.def("__sub__", &__sub__<ScalarFunction,ScalarFunction,ScalarFunction>);
    scalar_function_class.def("__mul__", &__mul__<ScalarFunction,ScalarFunction,ScalarFunction>);
    scalar_function_class.def("__div__", &__div__<ScalarFunction,ScalarFunction,ScalarFunction>);
    scalar_function_class.def("__add__", &__add__<ScalarFunction,ScalarFunction,Real>);
    scalar_function_class.def("__sub__", &__sub__<ScalarFunction,ScalarFunction,Real>);
    scalar_function_class.def("__mul__", &__mul__<ScalarFunction,ScalarFunction,Real>);
    scalar_function_class.def("__div__", &__div__<ScalarFunction,ScalarFunction,Real>);
    scalar_function_class.def("__radd__", &__radd__<ScalarFunction,ScalarFunction,Real>);
    scalar_function_class.def("__rsub__", &__rsub__<ScalarFunction,ScalarFunction,Real>);
    scalar_function_class.def("__rmul__", &__rmul__<ScalarFunction,ScalarFunction,Real>);
    scalar_function_class.def("__rdiv__", &__rdiv__<ScalarFunction,ScalarFunction,Real>);
    scalar_function_class.def(self_ns::str(self));

    scalar_function_class.def("constant", (ScalarFunction(*)(uint,Real)) &ScalarFunction::constant);
    scalar_function_class.def("variable", (ScalarFunction(*)(uint,uint)) &ScalarFunction::variable);
    scalar_function_class.staticmethod("constant");
    scalar_function_class.staticmethod("variable");

    def("evaluate_approx", (Float(*)(const ScalarFunction&,const Vector<Float>&)) &evaluate_approx);
    def("evaluate", (Interval(*)(const ScalarFunction&,const Vector<Interval>&)) &evaluate);
    def("gradient_approx",(Vector<Float>(*)(const ScalarFunction&,const Vector<Float>&)) &gradient_approx);
    def("gradient",(Vector<Interval>(*)(const ScalarFunction&,const Vector<Interval>&)) &gradient);

    def("derivative", (ScalarFunction(ScalarFunction::*)(uint)const) &ScalarFunction::derivative);

}

void export_vector_function()
{
    class_<VectorFunction>
        vector_function_class("VectorFunction", init<VectorFunction>());
    vector_function_class.def("result_size", &VectorFunction::result_size);
    vector_function_class.def("argument_size", &VectorFunction::argument_size);
    vector_function_class.def("__getitem__", &VectorFunction::get);
    vector_function_class.def("__setitem__", &VectorFunction::set);
    vector_function_class.def("__call__", (Vector<Interval>(VectorFunction::*)(const Vector<Interval>&)const)&VectorFunction::operator() );
    vector_function_class.def(self_ns::str(self));

    vector_function_class.def("constant", (VectorFunction(*)(uint,Vector<Real>)) &VectorFunction::constant);
    vector_function_class.def("identity", (VectorFunction(*)(uint)) &VectorFunction::identity);
    vector_function_class.staticmethod("constant");
    vector_function_class.staticmethod("identity");

    def("evaluate_approx", (Vector<Float>(*)(const VectorFunction&,const Vector<Float>&)) &evaluate_approx);
    def("evaluate", (Vector<Interval>(*)(const VectorFunction&,const Vector<Interval>&)) &evaluate);

    def("join", (VectorFunction(*)(const ScalarFunction&, const ScalarFunction&)) &join);
    def("join", (VectorFunction(*)(const VectorFunction&, const ScalarFunction&)) &join);
    def("join", (VectorFunction(*)(const ScalarFunction&, const VectorFunction&)) &join);
    def("join", (VectorFunction(*)(const VectorFunction&, const VectorFunction&)) &join);

    def("compose", (ScalarFunction(*)(const ScalarFunction&,const VectorFunction&)) &compose);
    def("compose", (VectorFunction(*)(const VectorFunction&,const VectorFunction&)) &compose);

    from_python<VectorFunction>();
}


void export_scalar_python_function()
{
    class_<ScalarPythonFunction, bases< ScalarFunctionInterface > > scalar_python_function_class("ScalarUserFunction", init<object>());
    scalar_python_function_class.def(init<uint,object>());
}

void export_vector_python_function()
{
    class_<VectorPythonFunction, bases< VectorFunctionInterface > > vector_python_function_class("VectorUserFunction", init<object>());
    vector_python_function_class.def(init<uint,uint,object>());
}


void function_submodule() {
    to_python< array<std::string> >();
    from_python< array<std::string> >();

    export_multi_index();
    export_monomial<Real>();
    export_polynomial<Real>();

    export_scalar_function();
    export_vector_function();

    //export_scalar_python_function();
    //export_vector_python_function();
}

