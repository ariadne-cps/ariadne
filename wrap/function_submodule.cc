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
#include "numeric.h"
#include "vector.h"
#include "differential.h"
#include "polynomial.h"
#include "taylor_model.h"
#include "taylor_function.h"
#include "function.h"


#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;


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
void read(Vector< Differential<X> >& td, const object& obj) {
    list elements=extract<list>(obj);
    ARIADNE_ASSERT(td.result_size()==uint(len(elements)));
    for(uint i=0; i!=td.size(); ++i) {
        boost::python::extract< Differential<X> > etv(elements[i]);
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
        Vector<Float> r(this->_result_size);
        read(r,this->_pyf(x));
        return r; }
    virtual Vector<Interval> evaluate (const Vector<Interval>& x) const {
        Vector<Interval> r(this->_result_size);
        read(r,this->_pyf(x));
        return r; }
    virtual Vector<TaylorModel> evaluate (const Vector<TaylorModel>& x) const {
        Vector<TaylorModel> r(this->_result_size);
        read(r,this->_pyf(x));
        return r; }
    virtual Vector< Differential<Float> > evaluate (const Vector< Differential<Float> >& x) const {
        Vector< Differential<Float> > r(this->_result_size);
        read(r,this->_pyf(x));
        return r; }
    virtual Vector< Differential<Interval> > evaluate (const Vector< Differential<Interval> >& x) const {
        Vector< Differential<Interval> > r(this->_result_size);
        read(r,this->_pyf(x));
        return r; }
    virtual Matrix<Float> jacobian (const Vector<Float>& x) const {
        Vector< Differential<Float> > rj(this->_result_size,this->_argument_size,1u);
        Vector< Differential<Float> > aj=Differential<Float>::variables(x.size(),x.size(),1u,x);
        read(rj,this->_pyf(aj));
        return get_jacobian<Float>(rj); }
    virtual Matrix<Interval> jacobian (const Vector<Interval>& x) const {
        Vector< Differential<Interval> > rj(this->_result_size,this->_argument_size,1u);
        Vector< Differential<Interval> > aj=Differential<Interval>::variables(x.size(),x.size(),1u,x);
        read(rj,this->_pyf(aj));
        return get_jacobian<Interval>(rj); }
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
typedef Matrix<Float> FMx;
typedef Matrix<Interval> IMx;
typedef Vector< Differential<Float> > FSDV;
typedef Vector< Differential<Interval> > ISDV;
typedef TaylorFunction TFM;

template<class F> TaylorFunction call_evaluate(const F& f, const TaylorFunction& tf) {
    return TaylorFunction(tf.domain(),f.evaluate(tf.models())); }

void export_function_interface()
{
    class_<FunctionPyWrap, boost::noncopyable> function_interface_class("FunctionInterface");
    function_interface_class.def("argument_size", pure_virtual(&FunctionInterface::argument_size));
    function_interface_class.def("result_size", pure_virtual(&FunctionInterface::result_size));
    function_interface_class.def("smoothness", pure_virtual(&FunctionInterface::smoothness));
    function_interface_class.def("evaluate",pure_virtual((IV(FunctionInterface::*)(const IV&)const)&FunctionInterface::evaluate));
    function_interface_class.def("jacobian",pure_virtual((IMx(FunctionInterface::*)(const IV&)const)&FunctionInterface::jacobian));
}


void export_python_function()
{
    class_<PythonFunction, bases< FunctionInterface > > python_function_class("AriadneFunction", init<object>());
    python_function_class.def(init<uint,uint,object>());
    python_function_class.def("argument_size", &PythonFunction::argument_size);
    python_function_class.def("result_size", &PythonFunction::result_size);
    python_function_class.def("smoothness", &PythonFunction::smoothness);
    python_function_class.def("__call__",pure_virtual((FV(FunctionInterface::*)(const FV&)const)&FunctionInterface::evaluate));
    python_function_class.def("__call__",pure_virtual((IV(FunctionInterface::*)(const IV&)const)&FunctionInterface::evaluate));
    python_function_class.def("evaluate",pure_virtual((FV(FunctionInterface::*)(const FV&)const)&FunctionInterface::evaluate));
    python_function_class.def("evaluate",pure_virtual((IV(FunctionInterface::*)(const IV&)const)&FunctionInterface::evaluate));
    python_function_class.def("jacobian",pure_virtual((FMx(FunctionInterface::*)(const FV&)const)&FunctionInterface::jacobian));
    python_function_class.def("jacobian",pure_virtual((IMx(FunctionInterface::*)(const IV&)const)&FunctionInterface::jacobian));
    python_function_class.def(self_ns::str(self));
}

std::string __str__(const AffineFunction& f) {
    std::stringstream ss;
    for(uint i=0; i!=f.result_size(); ++i) {
        ss<<(i==0?"AffineFunction( ":"; ");
        if(f.b()[i]!=0) { ss<<f.b()[i]; }
        for(uint j=0; j!=f.argument_size(); ++j) {
            if(f.A()[i][j]!=0) {
                if(f.A()[i][j]>0) { ss<<"+"; } else { ss<<"-"; }
                if(abs(f.A()[i][j])!=1) { ss<<abs(f.A()[i][j]); }
                //ss<<char('x'+j);
                ss<<"x"<<j;
            }
        }
    }
    ss<<" )";
    return ss.str();
}

void export_affine_function()
{
    typedef AffineFunction AF;

    class_<AffineFunction, bases< FunctionInterface > > affine_function_class("AffineFunction", init<Matrix<Float>, Vector<Float> >());
    affine_function_class.def("argument_size", &AffineFunction::argument_size);
    affine_function_class.def("result_size", &AffineFunction::result_size);
    affine_function_class.def("smoothness", &AffineFunction::smoothness);
    affine_function_class.def("__call__",(FV(AffineFunction::*)(const FV&)const)&AffineFunction::evaluate);
    affine_function_class.def("__call__",(IV(AffineFunction::*)(const IV&)const)&AffineFunction::evaluate);
    affine_function_class.def("__call__",(TFM(*)(const AffineFunction&, const TFM&))&call_evaluate<AffineFunction>);
    affine_function_class.def("evaluate",(FV(AffineFunction::*)(const FV&)const)&AffineFunction::evaluate);
    affine_function_class.def("evaluate",(IV(AffineFunction::*)(const IV&)const)&AffineFunction::evaluate);
    affine_function_class.def("evaluate",(TFM(*)(const AffineFunction&, const TFM&))&call_evaluate<AffineFunction>);
    affine_function_class.def("jacobian",(FMx(AffineFunction::*)(const FV&)const)&AffineFunction::jacobian);
    affine_function_class.def("jacobian",(IMx(AffineFunction::*)(const IV&)const)&AffineFunction::jacobian);
    //affine_function_class.def(self_ns::str(self));
    affine_function_class.def("__str__",(std::string(*)(const AffineFunction&))&__str__);
}

void export_polynomial()
{
    typedef Polynomial<Float> FP;

    Float real;
    class_< FP > polynomial_class("Polynomial", init<int>());
    polynomial_class.def("constant", (FP(*)(uint,double)) &FP::constant);
    polynomial_class.staticmethod("constant");
    polynomial_class.def("variable", (FP(*)(uint,uint)) &FP::variable);
    polynomial_class.staticmethod("variable");

    polynomial_class.def("argument_size", &FP::argument_size);
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
    polynomial_class.def(self_ns::str(self));
}



void function_submodule() {
    export_function_interface();
    export_affine_function();
    export_polynomial();
    export_python_function();
}

