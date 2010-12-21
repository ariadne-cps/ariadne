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

#include <boost/python.hpp>

#include <iostream>
#include <iomanip>

#include "array.h"
#include "container.h"
#include "numeric.h"
#include "vector.h"
#include "expansion.h"
#include "multi_index.h"
#include "taylor_model.h"
#include "differential.h"
#include "polynomial.h"
#include "affine.h"
#include "taylor_function.h"
#include "constraint.h"
#include "function.h"
#include "function_mixin.h"
#include "expression.h"
#include "space.h"
#include "assignment.h"
#include "formula.h"
#include "constraint_solver.h"

#include "../src/function_mixin.tcc"

#include "utilities.h"

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
struct from_python<RealVectorFunction>
{
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<RealVectorFunction>()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        list lst=extract<list>(obj_ptr);
        void* storage = ((converter::rvalue_from_python_storage< RealVectorFunction >*)   data)->storage.bytes;
        RealVectorFunction res(len(lst),0);
        for(uint i=0; i!=res.result_size(); ++i) { res.set(i,extract<RealScalarFunction>(lst[i])); }
        new (storage) RealVectorFunction(res);
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
    : public ScalarFunctionMixin<ScalarPythonFunction,Real>
{
    friend class ScalarFunctionMixin<ScalarPythonFunction,Real>;
    template<class T> void _compute(T& r, const Vector<T>& a) const {
        r=boost::python::extract<T>(this->_pyf(a)); }
  public:
    ScalarPythonFunction(std::string& nm, uint as, const object& pyf) : _name(nm), _argument_size(as), _pyf(pyf) { }
    ScalarPythonFunction(uint as, const object& pyf) : _name(),  _argument_size(as), _pyf(pyf) { }
    ScalarPythonFunction(const object& pyf)
        : _name(),
          _argument_size(extract<int>(pyf.attr("argument_size"))),
          _pyf(pyf) { }

    ScalarPythonFunction* clone() const { return new ScalarPythonFunction(*this); }
    virtual uint argument_size() const { return this->_argument_size; }

    virtual Vector<Float> gradient(const Vector<Float>& x) const {
        return this->evaluate(Differential<Float>::variables(1u,x)).gradient(); }
    virtual Vector<Interval> gradient(const Vector<Interval>& x) const {
        return this->evaluate(Differential<Interval>::variables(1u,x)).gradient(); }

    virtual RealScalarFunctionInterface* _derivative (uint j) const {
        ARIADNE_FAIL_MSG("Cannot get a component of a Python function"); }
    virtual std::ostream& repr(std::ostream& os) const { return os; }
    virtual std::ostream& write(std::ostream& os) const {
        os << "ScalarUserFunction( ";
        if(this->_name.size()>0) { os << "name=" << this->_name << ", "; }
        os << "argument_size="<<this->_argument_size;
        return os << " )"; }
    RealScalarFunction derivative(uint j) const { return this->_derivative(j); }
  private:
    std::string _name;
    uint _argument_size;
    boost::python::object _pyf;
};


class VectorPythonFunction
    : public VectorFunctionMixin<VectorPythonFunction,Real>
{
    friend class VectorFunctionMixin<VectorPythonFunction,Real>;
    template<class T> void _compute(Vector<T>& r, const Vector<T>& a) const {
        r=boost::python::extract< Vector<T> >(this->_pyf(a)); }
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

    virtual Matrix<Float> jacobian (const Vector<Float>& x) const {
        return this->evaluate(Differential<Float>::variables(1u,x)).jacobian(); }
    virtual Matrix<Interval> jacobian (const Vector<Interval>& x) const {
        return this->evaluate(Differential<Interval>::variables(1u,x)).jacobian(); }

    virtual RealScalarFunctionInterface* _get(uint i) const {
        ARIADNE_FAIL_MSG("Cannot get a component of a Python function"); }
    virtual RealScalarFunction operator[](uint i) const {
        ARIADNE_FAIL_MSG("Cannot get a component of a Python function"); }

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
typedef Vector<IntervalTaylorModel> TMV;
typedef VectorTaylorFunction TFM;
typedef IntervalTaylorModel TM;





void export_multi_index()
{
    class_< MultiIndex > multi_index_class("MultiIndex", init<uint>());
    multi_index_class.def(init<MultiIndex>());
    multi_index_class.def("__getitem__",&MultiIndex::get);
    multi_index_class.def("__setitem__",&MultiIndex::set);
    multi_index_class.def("degree",&MultiIndex::degree);
    multi_index_class.def(self_ns::str(self));

    from_python<MultiIndex>();
    to_python< List<MultiIndex> >();
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
    X real;

    class_< Polynomial<X> > polynomial_class(python_name<X>("Polynomial"), init< Polynomial<X> >());
    polynomial_class.def(init<uint>());
    polynomial_class.def("constant", (Polynomial<X>(*)(uint,double)) &Polynomial<X>::constant);
    polynomial_class.staticmethod("constant");
    polynomial_class.def("variable", (Polynomial<X>(*)(uint,uint)) &Polynomial<X>::variable);
    polynomial_class.staticmethod("variable");
    polynomial_class.def("variables", (Vector< Polynomial<X> >(*)(uint)) &Polynomial<X>::variables);
    polynomial_class.staticmethod("variables");

    polynomial_class.def("argument_size", &Polynomial<X>::argument_size);
    polynomial_class.def("insert", &Polynomial<X>::insert);
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
    polynomial_class.def("__iter__",boost::python::iterator< Polynomial<X> >());
    polynomial_class.def(self_ns::str(self));
    //polynomial_class.def(self_ns::repr(self));

    to_python< Vector< Polynomial<X> > >();
}

void export_scalar_function()
{
    class_<RealScalarFunction>
        scalar_function_class("ScalarFunction", init<RealScalarFunction>());
    scalar_function_class.def(init<uint>());
    scalar_function_class.def(init< Polynomial<Real> >());
    scalar_function_class.def("argument_size", &RealScalarFunction::argument_size);
    scalar_function_class.def("derivative", &RealScalarFunction::derivative);
    scalar_function_class.def("polynomial", &RealScalarFunction::polynomial);
    scalar_function_class.def("__call__", (Interval(RealScalarFunction::*)(const Vector<Interval>&)const)&RealScalarFunction::operator() );
    scalar_function_class.def("__call__", (Float(RealScalarFunction::*)(const Vector<Float>&)const)&RealScalarFunction::operator() );
    scalar_function_class.def("__call__", (IntervalDifferential(RealScalarFunction::*)(const Vector<IntervalDifferential>&)const)&RealScalarFunction::evaluate );
    scalar_function_class.def("__call__", (FloatDifferential(RealScalarFunction::*)(const Vector<FloatDifferential>&)const)&RealScalarFunction::evaluate );
    scalar_function_class.def("gradient", (Vector<Interval>(RealScalarFunction::*)(const Vector<Interval>&)const)&RealScalarFunction::gradient );
    scalar_function_class.def("gradient", (Vector<Interval>(RealScalarFunction::*)(const Vector<Interval>&)const)&RealScalarFunction::gradient );
    scalar_function_class.def("__pos__", &__pos__<RealScalarFunction,RealScalarFunction>);
    scalar_function_class.def("__neg__", &__neg__<RealScalarFunction,RealScalarFunction>);
    scalar_function_class.def("__add__", &__add__<RealScalarFunction,RealScalarFunction,RealScalarFunction>);
    scalar_function_class.def("__sub__", &__sub__<RealScalarFunction,RealScalarFunction,RealScalarFunction>);
    scalar_function_class.def("__mul__", &__mul__<RealScalarFunction,RealScalarFunction,RealScalarFunction>);
    scalar_function_class.def("__div__", &__div__<RealScalarFunction,RealScalarFunction,RealScalarFunction>);
    scalar_function_class.def("__add__", &__add__<RealScalarFunction,RealScalarFunction,Real>);
    scalar_function_class.def("__sub__", &__sub__<RealScalarFunction,RealScalarFunction,Real>);
    scalar_function_class.def("__mul__", &__mul__<RealScalarFunction,RealScalarFunction,Real>);
    scalar_function_class.def("__div__", &__div__<RealScalarFunction,RealScalarFunction,Real>);
    scalar_function_class.def("__radd__", &__radd__<RealScalarFunction,RealScalarFunction,Real>);
    scalar_function_class.def("__rsub__", &__rsub__<RealScalarFunction,RealScalarFunction,Real>);
    scalar_function_class.def("__rmul__", &__rmul__<RealScalarFunction,RealScalarFunction,Real>);
    scalar_function_class.def("__rdiv__", &__rdiv__<RealScalarFunction,RealScalarFunction,Real>);
    scalar_function_class.def("__eq__", &__eq__<NonlinearConstraint,RealScalarFunction,Interval>);
    scalar_function_class.def("__eq__", &__eq__<NonlinearConstraint,RealScalarFunction,Real>);
    scalar_function_class.def("__le__", &__le__<NonlinearConstraint,RealScalarFunction,Real>);
    scalar_function_class.def("__ge__", &__ge__<NonlinearConstraint,RealScalarFunction,Real>);
    scalar_function_class.def(self_ns::str(self));

    scalar_function_class.def("constant", (RealScalarFunction(*)(uint,Real)) &RealScalarFunction::constant);
    scalar_function_class.def("coordinate", (RealScalarFunction(*)(uint,uint)) &RealScalarFunction::coordinate);
    scalar_function_class.staticmethod("constant");
    scalar_function_class.staticmethod("coordinate");

    def("evaluate_approx", (Float(*)(const RealScalarFunction&,const Vector<Float>&)) &evaluate_approx);
    def("evaluate", (Interval(*)(const RealScalarFunction&,const Vector<Interval>&)) &evaluate);

    def("derivative", (RealScalarFunction(RealScalarFunction::*)(uint)const) &RealScalarFunction::derivative);

    def("pow", (RealScalarFunction(*)(const RealScalarFunction&,int)) &pow);
    def("rec", (RealScalarFunction(*)(const RealScalarFunction&)) &rec);
    def("sqr", (RealScalarFunction(*)(const RealScalarFunction&)) &sqr);
    def("sqrt", (RealScalarFunction(*)(const RealScalarFunction&)) &sqrt);
    def("exp", (RealScalarFunction(*)(const RealScalarFunction&)) &exp);
    def("log", (RealScalarFunction(*)(const RealScalarFunction&)) &log);
    def("sin", (RealScalarFunction(*)(const RealScalarFunction&)) &sin);
    def("cos", (RealScalarFunction(*)(const RealScalarFunction&)) &cos);
    def("tan", (RealScalarFunction(*)(const RealScalarFunction&)) &tan);

    typedef Polynomial<Real> RealPolynomial;
    implicitly_convertible<RealPolynomial,RealScalarFunction>();
}

void export_vector_function()
{

    class_<RealVectorFunction>
        vector_function_class("VectorFunction", init<RealVectorFunction>());
    vector_function_class.def(init<uint,uint>());

    vector_function_class.def("result_size", &RealVectorFunction::result_size);
    vector_function_class.def("argument_size", &RealVectorFunction::argument_size);
    vector_function_class.def("__getitem__", &RealVectorFunction::get);
    vector_function_class.def("__setitem__", &RealVectorFunction::set);
    vector_function_class.def("__call__", (Vector<Interval>(RealVectorFunction::*)(const Vector<Interval>&)const)&RealVectorFunction::operator() );
    vector_function_class.def("__call__", (Vector<Float>(RealVectorFunction::*)(const Vector<Float>&)const)&RealVectorFunction::operator() );
    vector_function_class.def("__call__", (Vector<IntervalDifferential>(RealVectorFunction::*)(const Vector<IntervalDifferential>&)const)&RealVectorFunction::evaluate );
    vector_function_class.def("__call__", (Vector<FloatDifferential>(RealVectorFunction::*)(const Vector<FloatDifferential>&)const)&RealVectorFunction::evaluate );
    vector_function_class.def("jacobian", (Matrix<Interval>(RealVectorFunction::*)(const Vector<Interval>&)const) &RealVectorFunction::jacobian);
    vector_function_class.def("jacobian", (Matrix<Float>(RealVectorFunction::*)(const Vector<Float>&)const) &RealVectorFunction::jacobian);
    vector_function_class.def(self_ns::str(self));

    vector_function_class.def("constant", (RealVectorFunction(*)(uint,Vector<Real>)) &RealVectorFunction::constant);
    vector_function_class.def("identity", (RealVectorFunction(*)(uint)) &RealVectorFunction::identity);
    vector_function_class.staticmethod("constant");
    vector_function_class.staticmethod("identity");

    def("evaluate_approx", (Vector<Float>(*)(const RealVectorFunction&,const Vector<Float>&)) &evaluate_approx);
    def("evaluate", (Vector<Interval>(*)(const RealVectorFunction&,const Vector<Interval>&)) &evaluate);

    def("join", (RealVectorFunction(*)(const RealScalarFunction&, const RealScalarFunction&)) &join);
    def("join", (RealVectorFunction(*)(const RealVectorFunction&, const RealScalarFunction&)) &join);
    def("join", (RealVectorFunction(*)(const RealScalarFunction&, const RealVectorFunction&)) &join);
    def("join", (RealVectorFunction(*)(const RealVectorFunction&, const RealVectorFunction&)) &join);

    def("compose", (RealScalarFunction(*)(const RealScalarFunction&,const RealVectorFunction&)) &compose);
    def("compose", (RealVectorFunction(*)(const RealVectorFunction&,const RealVectorFunction&)) &compose);

    from_python<RealVectorFunction>();
}


void export_scalar_python_function()
{
    class_<ScalarPythonFunction, bases< RealScalarFunctionInterface > > scalar_python_function_class("ScalarUserFunction", init<object>());
    scalar_python_function_class.def(init<uint,object>());
}

void export_vector_python_function()
{
    class_<VectorPythonFunction, bases< RealVectorFunctionInterface > > vector_python_function_class("VectorUserFunction", init<object>());
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

