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
struct from_python< Vector<ScalarFunctionInterface> >
{
    from_python() { converter::registry::push_back(&convertible,&construct,type_id< Vector<ScalarFunctionInterface> >()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        list lst=extract<list>(obj_ptr);
        void* storage = ((converter::rvalue_from_python_storage< Vector<ScalarFunctionInterface> >*)   data)->storage.bytes;
        Vector<ScalarFunctionInterface> res(len(lst));
        for(uint i=0; i!=res.size(); ++i) { res.set(i,extract<ScalarFunctionInterface&>(lst[i])); }
        new (storage) Vector<ScalarFunctionInterface>(res);
        data->convertible = storage;
    }
};


template<>
struct from_python<VectorPolynomialFunction>
{
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<VectorPolynomialFunction>()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        list lst=extract<list>(obj_ptr);
        void* storage = ((converter::rvalue_from_python_storage< Vector<ScalarFunctionInterface> >*)   data)->storage.bytes;
        VectorPolynomialFunction res(Vector< Polynomial<Interval> >(len(lst)));
        for(uint i=0; i!=res.size(); ++i) { res.set(i,ScalarPolynomialFunction(extract<ScalarPolynomialFunction&>(lst[i]))); }
        new (storage) Vector<ScalarFunctionInterface>(res);
        data->convertible = storage;
    }
};


// Define to remove ambiguity
ScalarPolynomialFunction operator*(const ScalarPolynomialFunction& p, const Float& x) {
    return ScalarPolynomialFunction(static_cast<const Polynomial<Interval>&>(p)*Interval(x));
}

// Define explicitly since cannot automatically wrap Vector<Polynomial<Interval> > to ScalarPolynomialFunction.
VectorPolynomialFunction antiderivative(const VectorPolynomialFunction& p, uint k) {
    VectorPolynomialFunction res(p.result_size(),p.argument_size());
    for(uint i=0; i!=res.result_size(); ++i) {
        res.set(i,antiderivative(p[i],k));
    }
    return res;
}

// Define explicitly since cannot automatically wrap Vector<Polynomial<Interval> > to ScalarPolynomialFunction.
VectorPolynomialFunction derivative(const VectorPolynomialFunction& p, uint k) {
    VectorPolynomialFunction res(p.result_size(),p.argument_size());
    for(uint i=0; i!=res.result_size(); ++i) {
        res.set(i,derivative(p[i],k));
    }
    return res;
}

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


template<class T> T evaluate(const ScalarFunctionInterface& f, const Vector<T>& x) { return f.evaluate(x); }
template<class T> Vector<T> evaluate(const VectorFunctionInterface& f, const Vector<T>& x) { return f.evaluate(x); }


// RATIONALE:
// Don't use standard method for wrapping abstract base classes, as this doesn't work well with
// shared_ptr<>. Instead, write static functions calling the overridden functions without using inheritence
class ScalarFunctionPyWrap {
  public:
    static ScalarFunctionInterface* clone(const ScalarFunctionInterface& f) { return f.clone(); };
    static uint argument_size(const ScalarFunctionInterface& f) { return f.argument_size(); }
    static ushort smoothness(const ScalarFunctionInterface& f) { return f.smoothness(); }
    template<class T> static T evaluate(const ScalarFunctionInterface& f, const Vector<T>& x) { return f.evaluate(x); }
};

/*
class ScalarFunctionPyWrap
    : public ScalarFunctionInterface
    , public wrapper< ScalarFunctionInterface >
{
    virtual ScalarFunctionInterface* clone() const { return this->get_override("clone")(); };
    virtual uint argument_size() const { return this->get_override("argument_size")(); }
    virtual ushort smoothness() const { return this->get_override("smoothness")(); }
    virtual Float evaluate(const Vector<Float>&) const { return this->get_override("__call__")(); }
    virtual Interval evaluate(const Vector<Interval>&) const { return this->get_override("__call__")(); }
    virtual TaylorModel evaluate(const Vector<TaylorModel>&) const { return this->get_override("__call__")(); }
    virtual Differential<Float> evaluate(const Vector< Differential<Float> >&) const { return this->get_override("__call__")(); }
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >&) const { return this->get_override("__call__")(); }
    virtual ScalarFunctionInterface* derivative(uint j) const { return this->get_override("derivative")(); }
    virtual std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

*/

class VectorFunctionPyWrap
    : public VectorFunctionInterface
    , public wrapper< VectorFunctionInterface >
{
    virtual VectorFunctionInterface* clone() const { return this->get_override("clone")(); };
    virtual uint result_size() const { return this->get_override("result_size")(); }
    virtual uint argument_size() const { return this->get_override("argument_size")(); }
    virtual ushort smoothness() const { return this->get_override("smoothness")(); }
    virtual Vector<Float> evaluate(const Vector<Float>&) const { return this->get_override("evaluate")(); }
    virtual Vector<Interval> evaluate(const Vector<Interval>&) const { return this->get_override("evaluate")(); }
    virtual Vector<TaylorModel> evaluate(const Vector<TaylorModel>&) const { return this->get_override("evaluate")(); }
    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >&) const { return this->get_override("evaluate")(); }
    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >&) const { return this->get_override("evaluate")(); }
    virtual Matrix<Float> jacobian(const Vector<Float>&) const { return this->get_override("jacobian")(); }
    virtual Matrix<Interval> jacobian(const Vector<Interval>&) const { return this->get_override("jacobian")(); }
    virtual Vector< Differential<Float> > expansion(const Vector<Float>&, const ushort&) const { return this->get_override("expansion")(); }
    virtual Vector< Differential<Interval> > expansion(const Vector<Interval>&, const ushort&) const { return this->get_override("expansion")(); }
    virtual std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};




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
    //monomial_class.def("key",(MultiIndex&(M::*)())&M::key,return_value_policy<reference_existing_object>());
    monomial_class.def("key",(const MultiIndex&(M::*)()const)&M::key,return_value_policy<copy_const_reference>());
    monomial_class.def("data",(const X&(M::*)()const) &M::data,return_value_policy<copy_const_reference>());
    monomial_class.def(self_ns::str(self));
}

template<class X>
void export_polynomial()
{
    typedef Polynomial<X> P;

    //double flt;
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

void export_scalar_function_interface()
{
    class_<ScalarFunctionInterface, shared_ptr<ScalarFunctionInterface>, boost::noncopyable>
        scalar_function_interface_class("ScalarFunctionInterface",no_init);

    // Don't use following standard wrapping technique due to clash with shared_ptr
    //expression_interface_class.def("argument_size", pure_virtual(&ScalarFunctionInterface::argument_size));

    scalar_function_interface_class.def("argument_size", &ScalarFunctionPyWrap::argument_size);
    scalar_function_interface_class.def("smoothness", &ScalarFunctionPyWrap::smoothness);
    scalar_function_interface_class.def("__call__", &ScalarFunctionPyWrap::evaluate<Float>);
    scalar_function_interface_class.def("__call__", &ScalarFunctionPyWrap::evaluate<Interval>);
    scalar_function_interface_class.def("__call__", &ScalarFunctionPyWrap::evaluate<TaylorModel>);
    scalar_function_interface_class.def("__call__", &ScalarFunctionPyWrap::evaluate< Differential<Interval> >);
    scalar_function_interface_class.def(self_ns::str(self));
    //expression_interface_class.def(self_ns::repr(self));

    def("evaluate",&ScalarFunctionPyWrap::evaluate<Interval>);
    def("evaluate",&ScalarFunctionPyWrap::evaluate<TaylorModel>);
}


void export_scalar_python_function()
{
    class_<ScalarPythonFunction, bases< ScalarFunctionInterface > > scalar_python_function_class("ScalarUserFunction", init<object>());
    scalar_python_function_class.def(init<uint,object>());
}





void export_scalar_affine_function()
{
    typedef ScalarAffineFunction Self;

    class_<ScalarAffineFunction, bases< ScalarFunctionInterface > > scalar_affine_function_class("ScalarAffineFunction", init<Vector<Interval>, Interval> ());
    scalar_affine_function_class.def("value",&ScalarAffineFunction::value,return_value_policy<copy_const_reference>());
    scalar_affine_function_class.def("gradient",(Vector<Interval>(ScalarAffineFunction::*)(const Vector<Interval>&)const)&ScalarAffineFunction::gradient);
    scalar_affine_function_class.def("__pos__", &__pos__<Self,Self>);
    scalar_affine_function_class.def("__neg__", &__neg__<Self,Self>);
    scalar_affine_function_class.def("__add__", &__add__<Self,Self,Self>);
    scalar_affine_function_class.def("__sub__", &__sub__<Self,Self,Self>);
    scalar_affine_function_class.def("__mul__", &__mul__<Self,Self,Interval>);
    scalar_affine_function_class.def("__rmul__", &__rmul__<Self,Self,Interval>);
    scalar_affine_function_class.def("__div__", &__div__<Self,Self,Interval>);

}



void export_scalar_polynomial_function()
{
    typedef ScalarPolynomialFunction Self;

    Float flt;
    Interval ivl;

    to_python< array<ScalarPolynomialFunction> >();

    implicitly_convertible< Polynomial<Interval>, ScalarPolynomialFunction>();

    class_< ScalarPolynomialFunction, bases<ScalarFunctionInterface> > scalar_polynomial_function_class("ScalarPolynomialFunction", init<int>());
    scalar_polynomial_function_class.def(init<Expression<Real>,Space<Real> >());
    scalar_polynomial_function_class.def(init<ScalarAffineFunction>());
    scalar_polynomial_function_class.def(init<ScalarPolynomialFunction>());
    scalar_polynomial_function_class.def(init< Polynomial<Interval> >());
    scalar_polynomial_function_class.def("constant", (ScalarPolynomialFunction(*)(uint,double)) &ScalarPolynomialFunction::constant);
    scalar_polynomial_function_class.def("constant", (ScalarPolynomialFunction(*)(uint,Interval)) &ScalarPolynomialFunction::constant);
    scalar_polynomial_function_class.def("variable", (ScalarPolynomialFunction(*)(uint,uint)) &ScalarPolynomialFunction::variable);
    scalar_polynomial_function_class.def("variables", (array<ScalarPolynomialFunction>(*)(uint)) &ScalarPolynomialFunction::variables);
    scalar_polynomial_function_class.staticmethod("constant");
    scalar_polynomial_function_class.staticmethod("variable");
    scalar_polynomial_function_class.staticmethod("variables");

    scalar_polynomial_function_class.def("argument_size", &ScalarPolynomialFunction::argument_size);
    scalar_polynomial_function_class.def("insert", &ScalarPolynomialFunction::insert);
    scalar_polynomial_function_class.def("__iter__",boost::python::iterator<ScalarPolynomialFunction>());

    //implicitly_convertible<RealExpressionPointer,PolynomialExpression>();

    scalar_polynomial_function_class.def("__pos__", &__pos__<Self,Self>);
    scalar_polynomial_function_class.def("__neg__", &__neg__<Self,Self>);
    scalar_polynomial_function_class.def("__add__", &__add__<Self,Self,Self>);
    scalar_polynomial_function_class.def("__sub__", &__sub__<Self,Self,Self>);
    scalar_polynomial_function_class.def("__mul__", &__mul__<Self,Self,Self>);

    scalar_polynomial_function_class.def("__add__", &__add__<Self,Self,Float>);
    scalar_polynomial_function_class.def("__sub__", &__sub__<Self,Self,Float>);
    scalar_polynomial_function_class.def("__mul__", &__mul__<Self,Self,Float>);
    scalar_polynomial_function_class.def("__div__", &__div__<Self,Self,Float>);
    scalar_polynomial_function_class.def("__radd__", &__add__<Self,Self,Float>);
    scalar_polynomial_function_class.def("__rsub__", &__rsub__<Self,Self,Float>);
    scalar_polynomial_function_class.def("__rmul__", &__mul__<Self,Self,Float>);

    scalar_polynomial_function_class.def("__add__", &__add__<Self,Self,Interval>);
    scalar_polynomial_function_class.def("__sub__", &__sub__<Self,Self,Interval>);
    scalar_polynomial_function_class.def("__mul__", &__mul__<Self,Self,Interval>);
    scalar_polynomial_function_class.def("__div__", &__div__<Self,Self,Interval>);
    scalar_polynomial_function_class.def("__radd__", &__add__<Self,Self,Interval>);
    scalar_polynomial_function_class.def("__rsub__", &__rsub__<Self,Self,Interval>);
    scalar_polynomial_function_class.def("__rmul__", &__mul__<Self,Self,Interval>);

    scalar_polynomial_function_class.def(self+=self);
    scalar_polynomial_function_class.def(self-=self);
    scalar_polynomial_function_class.def(self+=flt);
    scalar_polynomial_function_class.def(self-=flt);
    scalar_polynomial_function_class.def(self*=flt);
    scalar_polynomial_function_class.def(self/=flt);
    scalar_polynomial_function_class.def(self+=ivl);
    scalar_polynomial_function_class.def(self-=ivl);
    scalar_polynomial_function_class.def(self*=ivl);
    scalar_polynomial_function_class.def(self/=ivl);

    scalar_polynomial_function_class.def("derivative",(ScalarPolynomialFunction*(ScalarPolynomialFunction::*)(uint)const)&ScalarPolynomialFunction::derivative,return_value_policy<manage_new_object>());
    scalar_polynomial_function_class.def("antiderivative",(ScalarPolynomialFunction*(ScalarPolynomialFunction::*)(uint)const)&ScalarPolynomialFunction::antiderivative,return_value_policy<manage_new_object>());


    def("derivative",(ScalarPolynomialFunction*(ScalarPolynomialFunction::*)(uint)const)&ScalarPolynomialFunction::derivative,return_value_policy<manage_new_object>());
    def("antiderivative",(ScalarPolynomialFunction*(ScalarPolynomialFunction::*)(uint)const)&ScalarPolynomialFunction::antiderivative,return_value_policy<manage_new_object>());

    //implicitly_convertible<Polynomial<Float>,PolynomialExpression>();
    //implicitly_convertible<Polynomial<Interval>,PolynomialExpression>();
}




void export_vector_function_interface()
{
    class_<VectorFunctionPyWrap, boost::noncopyable> vector_function_interface_class("VectorFunctionInterface");
    vector_function_interface_class.def("argument_size", pure_virtual(&VectorFunctionInterface::argument_size));
    vector_function_interface_class.def("result_size", pure_virtual(&VectorFunctionInterface::result_size));
    vector_function_interface_class.def("smoothness", pure_virtual(&VectorFunctionInterface::smoothness));
    vector_function_interface_class.def("__len__", pure_virtual(&VectorFunctionInterface::result_size));
    vector_function_interface_class.def("__call__",pure_virtual((IV(VectorFunctionInterface::*)(const IV&)const)&VectorFunctionInterface::evaluate));
    vector_function_interface_class.def("__call__",pure_virtual((TMV(VectorFunctionInterface::*)(const TMV&)const)&VectorFunctionInterface::evaluate));
    vector_function_interface_class.def("jacobian",(FMx(VectorFunctionInterface::*)(const FV&)const)&VectorFunctionInterface::jacobian);
    vector_function_interface_class.def("jacobian",(IMx(VectorFunctionInterface::*)(const IV&)const)&VectorFunctionInterface::jacobian);
    vector_function_interface_class.def(self_ns::str(self));
    //vector_function_interface_class.def(self_ns::repr(self));

    def("evaluate",(Vector<Float>(*)(const VectorFunctionInterface&,const Vector<Float>&))&evaluate);
    def("evaluate",(Vector<Interval>(*)(const VectorFunctionInterface&,const Vector<Interval>&))&evaluate);
    def("evaluate",(Vector<TaylorModel>(*)(const VectorFunctionInterface&,const Vector<TaylorModel>&))&evaluate);
}

void export_vector_function()
{
    typedef Vector<ScalarFunctionInterface> VectorFunction;

    class_<VectorFunction, bases< VectorFunctionInterface > > vector_function_class("VectorFunction",init<VectorFunction>());
    vector_function_class.def("__getitem__",(ScalarFunctionInterface const&(VectorFunction::*)(uint)const)&VectorFunction::operator[],return_value_policy<copy_const_reference>());

    from_python<VectorFunction>();
}

void export_vector_python_function()
{
    class_<VectorPythonFunction, bases< VectorFunctionInterface > > vector_python_function_class("VectorUserFunction", init<object>());
    vector_python_function_class.def(init<uint,uint,object>());
}

void export_vector_affine_function()
{
    class_<VectorAffineFunction, bases< VectorFunctionInterface > >
        vector_affine_function_class("VectorAffineFunction", init<VectorAffineFunction>());
    vector_affine_function_class.def(init< Matrix<Interval>, Vector<Interval> >());
    vector_affine_function_class.def("__getitem__",(VectorAffineFunction(VectorAffineFunction::*)(uint)const)&VectorAffineFunction::operator[]);

}

void export_vector_polynomial_function()
{
    from_python< array<String> >();
    //array_from_python_list<RealVariable>();

    class_< VectorPolynomialFunction, bases<VectorFunctionInterface> >
        vector_polynomial_function_class("VectorPolynomialFunction", init<VectorPolynomialFunction>());
    vector_polynomial_function_class.def(init<uint,uint>());
    vector_polynomial_function_class.def("__getitem__",(Interval const&(VectorPolynomialFunction::*)(uint)const)&VectorPolynomialFunction::operator[],return_value_policy<copy_const_reference>());
    vector_polynomial_function_class.def("__setitem__",(void(VectorPolynomialFunction::*)(uint,const VectorPolynomialFunction&)const)&VectorPolynomialFunction::set);

    def("derivative", (VectorPolynomialFunction(*)(const VectorPolynomialFunction&,uint)) &derivative);
    def("antiderivative", (VectorPolynomialFunction(*)(const VectorPolynomialFunction&,uint)) &antiderivative);

    def("join", (VectorPolynomialFunction(*)(const VectorPolynomialFunction&, const VectorPolynomialFunction&)) &join);
    def("join", (VectorPolynomialFunction(*)(const VectorPolynomialFunction&, const ScalarPolynomialFunction&)) &join);
    def("join", (VectorPolynomialFunction(*)(const ScalarPolynomialFunction&, const ScalarPolynomialFunction&)) &join);

    from_python<VectorPolynomialFunction>();

}




void function_submodule() {
    to_python< array<std::string> >();
    from_python< array<std::string> >();

    export_multi_index();


    export_scalar_function_interface();
    export_vector_function_interface();

    //export_polynomial<Float>();
    //export_polynomial<Interval>();
    export_monomial<Interval>();

    export_scalar_affine_function();
    export_scalar_polynomial_function();
    export_scalar_python_function();

    export_vector_function();
    export_vector_affine_function();
    export_vector_polynomial_function();
    export_vector_python_function();
}

