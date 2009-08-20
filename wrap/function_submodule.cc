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
#include "taylor_expression.h"
#include "taylor_function.h"
#include "function.h"
#include "expression.h"
#include "formula.h"

#include "utilities.h"


#include <boost/python.hpp>
using namespace boost::python;

namespace Ariadne {

struct vector_function_from_python_list
{
    vector_function_from_python_list()     {
        converter::registry::push_back(&convertible,&construct,type_id< Vector<ScalarFunctionInterface> >());
    }

    static void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr)) { return 0; }
        return obj_ptr;
    }

    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data)
    {
        list lst=extract<list>(obj_ptr);
        //if (value == 0) boost::python::throw_error_already_set();
        void* storage = ((converter::rvalue_from_python_storage< Vector<ScalarFunctionInterface> >*)   data)->storage.bytes;
        Vector<ScalarFunctionInterface> res(len(lst));
        for(uint i=0; i!=res.size(); ++i) { res.set(i,extract<ScalarFunctionInterface&>(lst[i])); }
        new (storage) Vector<ScalarFunctionInterface>(res);
        data->convertible = storage;
    }
};


// Define to remove ambiguity
ScalarPolynomialFunction operator*(const ScalarPolynomialFunction& p, const Float& x) {
    return ScalarPolynomialFunction(static_cast<const Polynomial<Interval>&>(p)*Interval(x));
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

template<class X> X* make(const boost::python::object& obj);

template<>
MultiIndex*
make(const boost::python::object& obj) {
    boost::python::tuple tup=boost::python::extract<boost::python::tuple>(obj);
    MultiIndex* a=new MultiIndex(len(tup));
    for(uint i=0; i!=a->size(); ++i) {
        //shared_ptr<ScalarFunctionInterface> expression_ptr=boost::python::extract<shared_ptr<ScalarFunctionInterface> >(elements[i]);
        uint ai=boost::python::extract<uint>(tup[i]);
        a->set(i,ai);
    }
    return a;
}

template<>
Vector<ScalarFunctionInterface>*
make(const boost::python::object& obj) {
    boost::python::list elements=boost::python::extract<boost::python::list>(obj);
    Vector<ScalarFunctionInterface>* r=new Vector<ScalarFunctionInterface>(len(elements));
    for(uint i=0; i!=r->size(); ++i) {
        //shared_ptr<ScalarFunctionInterface> expression_ptr=boost::python::extract<shared_ptr<ScalarFunctionInterface> >(elements[i]);
        const ScalarFunctionInterface& expression_ref=boost::python::extract<const ScalarFunctionInterface&>(elements[i]);
        shared_ptr<ScalarFunctionInterface> expression_ptr(expression_ref.clone());
        r->set(i,expression_ptr);
    }
    return r;
}

template<>
PolynomialFunction*
make<PolynomialFunction>(const boost::python::object& obj)
{
    boost::python::list lst=boost::python::extract<boost::python::list>(obj);
    Vector<ScalarPolynomialFunction> expressions(boost::python::len(lst));
    for(uint i=0; i!=expressions.size(); ++i) {
        expressions[i]=boost::python::extract<ScalarPolynomialFunction>(lst[i]);
    }
    return new PolynomialFunction(expressions);
}

template<class T> static T evaluate(const ScalarFunctionInterface& f, const Vector<T>& x) { return f.evaluate(x); }
template<class T> static Vector<T> evaluate(const FunctionInterface& f, const Vector<T>& x) { return f.evaluate(x); }

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

class FunctionPyWrap
    : public FunctionInterface
    , public wrapper< FunctionInterface >
{
    virtual FunctionInterface* clone() const { return this->get_override("clone")(); };
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
        os << "UserExpression( ";
        if(this->_name.size()>0) { os << "name=" << this->_name << ", "; }
        os << "argument_size="<<this->_argument_size;
        return os << " )"; }
  private:
    std::string _name;
    uint _argument_size;
    boost::python::object _pyf;
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
        os << "UserFunction( ";
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



/*
template<class F> TaylorFunction call_evaluate(const F& f, const TaylorFunction& tf) {
    return TaylorFunction(tf.domain(),f.evaluate(tf.models())); }
template<class E> TaylorExpression call_evaluate_expression(const E& f, const TaylorFunction& tf) {
    return TaylorExpression(tf.domain(),f.evaluate(tf.models())); }
*/

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
typedef TaylorFunction TFM;
typedef TaylorModel TM;

template<class T> T evaluate(const ScalarFunctionInterface& expr, const Vector<T>& ix) { return expr.evaluate(ix); }




void export_multi_index()
{
    class_< MultiIndex > multi_index_class("MultiIndex", init<uint>());
    multi_index_class.def(init<MultiIndex>());
    multi_index_class.def("__init__",make_constructor(&make<MultiIndex>));
    multi_index_class.def("__getitem__",&MultiIndex::get);
    multi_index_class.def("__setitem__",&MultiIndex::set);
    multi_index_class.def("degree",&MultiIndex::degree);
    multi_index_class.def(self_ns::str(self));
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

    to_python_converter< Vector< Polynomial<X> >, vector_to_python_list< Polynomial<X> > >();
}

void export_scalar_function_interface()
{
    class_<ScalarFunctionInterface, shared_ptr<ScalarFunctionInterface>, boost::noncopyable>
        scalar_function_interface_class("ScalarFunctionInterface",no_init);

    // Don't use following standard wrapping technique due to clash with shared_ptr
    //expression_interface_class.def("argument_size", pure_virtual(&ScalarFunctionInterface::argument_size));

    scalar_function_interface_class.def("argument_size", &ScalarFunctionPyWrap::argument_size);
    scalar_function_interface_class.def("smoothness", &ScalarFunctionPyWrap::smoothness);
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
   class_<ScalarAffineFunction, bases< ScalarFunctionInterface > > scalar_affine_function_class("Affine", init<Vector<Interval>, Interval> ());
    //scalar_affine_function_class.def(init<ExpressionPointer>());
    scalar_affine_function_class.def("value",&ScalarAffineFunction::value);
    scalar_affine_function_class.def("gradient",(Vector<Interval>(ScalarAffineFunction::*)(const Vector<Interval>&)const)&ScalarAffineFunction::gradient);
    scalar_affine_function_class.def("__add__",(ScalarAffineFunction(*)(const ScalarAffineFunction&,const ScalarAffineFunction&))&operator+);
    scalar_affine_function_class.def("__mul__",(ScalarAffineFunction(*)(const ScalarAffineFunction&,const Interval&))&operator*);
    scalar_affine_function_class.def("__rmul__",&__rmul__<ScalarAffineFunction,ScalarAffineFunction,Interval>);

    def("derivative",(Interval(*)(const ScalarAffineFunction&,uint)) &derivative);

}



void export_scalar_polynomial_function()
{
    typedef ScalarPolynomialFunction PE;

    Float flt;
    Interval ivl;

    array_to_python_list<ScalarPolynomialFunction>();

    implicitly_convertible< Polynomial<Interval>, ScalarPolynomialFunction>();

    class_< PE, bases<ScalarFunctionInterface> > scalar_polynomial_function_class("Polynomial", init<int>());
    scalar_polynomial_function_class.def(init<Expression<Real>,Space<Real> >());
    scalar_polynomial_function_class.def(init<ScalarAffineFunction>());
    scalar_polynomial_function_class.def(init<ScalarPolynomialFunction>());
    scalar_polynomial_function_class.def(init< Polynomial<Interval> >());
    scalar_polynomial_function_class.def("constant", (PE(*)(uint,double)) &PE::constant);
    scalar_polynomial_function_class.def("constant", (PE(*)(uint,Interval)) &PE::constant);
    scalar_polynomial_function_class.def("variable", (PE(*)(uint,uint)) &PE::variable);
    scalar_polynomial_function_class.def("variables", (array<PE>(*)(uint)) &PE::variables);
    scalar_polynomial_function_class.staticmethod("constant");
    scalar_polynomial_function_class.staticmethod("variable");
    scalar_polynomial_function_class.staticmethod("variables");

    scalar_polynomial_function_class.def("argument_size", &PE::argument_size);
    scalar_polynomial_function_class.def("insert", &PE::insert);
    scalar_polynomial_function_class.def("__iter__",boost::python::iterator<PE>());

    //implicitly_convertible<RealExpressionPointer,PolynomialExpression>();

    scalar_polynomial_function_class.def("__pos__", &__pos__<PE,PE>);
    scalar_polynomial_function_class.def("__neg__", &__neg__<PE,PE>);
    scalar_polynomial_function_class.def("__add__", &__add__<PE,PE,PE>);
    scalar_polynomial_function_class.def("__sub__", &__sub__<PE,PE,PE>);
    scalar_polynomial_function_class.def("__mul__", &__mul__<PE,PE,PE>);

    scalar_polynomial_function_class.def("__add__", &__add__<PE,PE,Float>);
    scalar_polynomial_function_class.def("__sub__", &__sub__<PE,PE,Float>);
    scalar_polynomial_function_class.def("__mul__", &__mul__<PE,PE,Float>);
    scalar_polynomial_function_class.def("__div__", &__div__<PE,PE,Float>);
    scalar_polynomial_function_class.def("__radd__", &__add__<PE,PE,Float>);
    scalar_polynomial_function_class.def("__rsub__", &__rsub__<PE,PE,Float>);
    scalar_polynomial_function_class.def("__rmul__", &__mul__<PE,PE,Float>);

    scalar_polynomial_function_class.def("__add__", &__add__<PE,PE,Interval>);
    scalar_polynomial_function_class.def("__sub__", &__sub__<PE,PE,Interval>);
    scalar_polynomial_function_class.def("__mul__", &__mul__<PE,PE,Interval>);
    scalar_polynomial_function_class.def("__div__", &__div__<PE,PE,Interval>);
    scalar_polynomial_function_class.def("__radd__", &__add__<PE,PE,Interval>);
    scalar_polynomial_function_class.def("__rsub__", &__rsub__<PE,PE,Interval>);
    scalar_polynomial_function_class.def("__rmul__", &__mul__<PE,PE,Interval>);

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

    def("derivative",(PE*(PE::*)(uint)const)&PE::derivative,return_value_policy<manage_new_object>());
    def("antiderivative",(PE*(PE::*)(uint)const)&PE::antiderivative,return_value_policy<manage_new_object>());

    //implicitly_convertible<Polynomial<Float>,PolynomialExpression>();
    //implicitly_convertible<Polynomial<Interval>,PolynomialExpression>();
}




void export_function_interface()
{
    class_<FunctionPyWrap, boost::noncopyable> function_interface_class("FunctionInterface");
    function_interface_class.def("argument_size", pure_virtual(&FunctionInterface::argument_size));
    function_interface_class.def("result_size", pure_virtual(&FunctionInterface::result_size));
    function_interface_class.def("smoothness", pure_virtual(&FunctionInterface::smoothness));
    function_interface_class.def("__len__", pure_virtual(&FunctionInterface::result_size));
    function_interface_class.def("__call__",pure_virtual((IV(FunctionInterface::*)(const IV&)const)&FunctionInterface::evaluate));
    function_interface_class.def("__call__",pure_virtual((TMV(FunctionInterface::*)(const TMV&)const)&FunctionInterface::evaluate));
    function_interface_class.def("jacobian",(FMx(FunctionInterface::*)(const FV&)const)&FunctionInterface::jacobian);
    function_interface_class.def("jacobian",(IMx(FunctionInterface::*)(const IV&)const)&FunctionInterface::jacobian);
    function_interface_class.def(self_ns::str(self));
    //function_interface_class.def(self_ns::repr(self));

    def("evaluate",(Vector<Float>(*)(const FunctionInterface&,const Vector<Float>&))&evaluate);
    def("evaluate",(Vector<Interval>(*)(const FunctionInterface&,const Vector<Interval>&))&evaluate);
    def("evaluate",(Vector<TaylorModel>(*)(const FunctionInterface&,const Vector<TaylorModel>&))&evaluate);
}

void export_python_function()
{
    class_<PythonFunction, bases< FunctionInterface > > python_function_class("UserFunction", init<object>());
    python_function_class.def(init<uint,uint,object>());
}

void export_vector_function()
{
    typedef Vector<ScalarFunctionInterface> VectorFunction;

    class_<VectorFunction, bases< FunctionInterface > > vector_function_class("VectorFunction",no_init);
    vector_function_class.def("__init__", make_constructor(&make<VectorFunction>) );
    vector_function_class.def("__getitem__",(ScalarFunctionInterface const&(VectorFunction::*)(uint)const)&VectorFunction::operator[],return_value_policy<copy_const_reference>());

    vector_function_from_python_list();
}

void export_affine_function()
{
    class_<AffineFunction, bases< FunctionInterface > >
        affine_function_class("AffineFunction", init<Matrix<Float>, Vector<Float> >());
    affine_function_class.def(init<AffineFunction>());
    affine_function_class.def("__getitem__",(ScalarAffineFunction(AffineFunction::*)(uint)const)&AffineFunction::operator[]);

}

void export_polynomial_function()
{
    array_from_python_list<String>();
    //array_from_python_list<RealVariable>();

    class_< PolynomialFunction, bases<FunctionInterface> >
        polynomial_function_class("PolynomialFunction", no_init);
    polynomial_function_class.def("__init__",make_constructor(&make<PolynomialFunction>));
    polynomial_function_class.def(init<PolynomialFunction>());
    polynomial_function_class.def("__getitem__",(ScalarPolynomialFunction const&(PolynomialFunction::*)(uint)const)&PolynomialFunction::operator[],return_value_policy<copy_const_reference>());
    polynomial_function_class.def("__setitem__",(void(PolynomialFunction::*)(uint,const ScalarPolynomialFunction&)const)&PolynomialFunction::set);


}




void function_submodule() {
    to_python_converter< array<ScalarPolynomialFunction>, array_to_python_list<ScalarPolynomialFunction> >();
    to_python_converter< array<std::string>, array_to_python_list<std::string> >();

    export_multi_index();


    export_scalar_function_interface();
    export_function_interface();

    //export_polynomial<Float>();
    //export_polynomial<Interval>();
    export_monomial<Interval>();

    export_scalar_affine_function();
    export_scalar_polynomial_function();
    export_scalar_python_function();

    export_vector_function();
    export_python_function();
    export_affine_function();
    export_polynomial_function();
}

