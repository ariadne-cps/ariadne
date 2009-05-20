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
#include "expansion.h"
#include "multi_index.h"
#include "differential.h"
#include "polynomial.h"
#include "taylor_model.h"
#include "taylor_expression.h"
#include "taylor_function.h"
#include "function.h"

#include "utilities.h"


#include <boost/python.hpp>
using namespace boost::python;

namespace Ariadne {

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
void read(Differential<X>& df, const object& obj) {
    boost::python::extract< Differential<X> > edf(obj);
    boost::python::extract<double> eflt(obj);
    if(edf.check()) {
        df=edf();
    } else if(eflt.check()) {
        df=eflt();
    } else {
        df=edf();
    }
}

template<class X>
void read(Vector< Differential<X> >& dfv, const object& obj) {
    list elements=extract<list>(obj);
    ARIADNE_ASSERT(dfv.result_size()==uint(len(elements)));
    for(uint i=0; i!=dfv.size(); ++i) {
        read(dfv[i],elements[i]);
    }
}

void read(TaylorModel& tm, const boost::python::object& obj);

class ExpressionPyWrap
    : public ExpressionInterface
    , public wrapper< ExpressionInterface >
{
    virtual ExpressionInterface* clone() const { return this->get_override("clone")(); };
    virtual uint argument_size() const { return this->get_override("argument_size")(); }
    virtual ushort smoothness() const { return this->get_override("smoothness")(); }
    virtual Float evaluate(const Vector<Float>&) const { return this->get_override("evaluate")(); }
    virtual Interval evaluate(const Vector<Interval>&) const { return this->get_override("evaluate")(); }
    virtual TaylorModel evaluate(const Vector<TaylorModel>&) const { return this->get_override("evaluate")(); }
    virtual Differential<Float> evaluate(const Vector< Differential<Float> >&) const { return this->get_override("evaluate")(); }
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >&) const { return this->get_override("evaluate")(); }
    virtual std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

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




class PythonExpression
    : public ExpressionInterface
{
  public:
    PythonExpression(std::string& nm, uint as, const object& pyf) : _name(nm), _argument_size(as), _pyf(pyf) { }
    PythonExpression(uint as, const object& pyf) : _name(),  _argument_size(as), _pyf(pyf) { }
    PythonExpression(const object& pyf)
        : _name(),
          _argument_size(extract<int>(pyf.attr("argument_size"))),
          _pyf(pyf) { }

    PythonExpression* clone() const { return new PythonExpression(*this); }
    virtual uint argument_size() const { return this->_argument_size; }
    virtual ushort smoothness() const { return 255; }

    virtual Float evaluate (const Vector<Float>& x) const {
        Float r;
        read(r,this->_pyf(x));
        return r; }
    virtual Interval evaluate (const Vector<Interval>& x) const {
        Interval r;
        read(r,this->_pyf(x));
        return r; }
    virtual TaylorModel evaluate (const Vector<TaylorModel>& x) const {
        TaylorModel r;
        read(r,this->_pyf(x));
        return r; }
    virtual Differential<Float> evaluate (const Vector< Differential<Float> >& x) const {
        Differential<Float> r;
        read(r,this->_pyf(x));
        return r; }
    virtual Differential<Interval> evaluate (const Vector< Differential<Interval> >& x) const {
        Differential<Interval> r;
        read(r,this->_pyf(x));
        return r; }
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


template<class F> TaylorFunction call_evaluate(const F& f, const TaylorFunction& tf) {
    return TaylorFunction(tf.domain(),f.evaluate(tf.models())); }
template<class E> TaylorExpression call_evaluate_expression(const E& f, const TaylorFunction& tf) {
    return TaylorExpression(tf.domain(),f.evaluate(tf.models())); }

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
typedef TaylorFunction TFM;


void export_expression_interface()
{
    class_<ExpressionPyWrap, boost::noncopyable> expression_interface_class("ExpressionInterface");
    expression_interface_class.def("argument_size", pure_virtual(&ExpressionInterface::argument_size));
    expression_interface_class.def("smoothness", pure_virtual(&ExpressionInterface::smoothness));
    expression_interface_class.def("evaluate",pure_virtual((Interval(ExpressionInterface::*)(const IV&)const)&ExpressionInterface::evaluate));
}

void export_function_interface()
{
    class_<FunctionPyWrap, boost::noncopyable> function_interface_class("FunctionInterface");
    function_interface_class.def("argument_size", pure_virtual(&FunctionInterface::argument_size));
    function_interface_class.def("result_size", pure_virtual(&FunctionInterface::result_size));
    function_interface_class.def("smoothness", pure_virtual(&FunctionInterface::smoothness));
    function_interface_class.def("evaluate",pure_virtual((IV(FunctionInterface::*)(const IV&)const)&FunctionInterface::evaluate));
    function_interface_class.def("jacobian",pure_virtual((IMx(FunctionInterface::*)(const IV&)const)&FunctionInterface::jacobian));
}


void export_python_expression()
{
    class_<PythonExpression, bases< ExpressionInterface > > python_expression_class("UserExpression", init<object>());
    python_expression_class.def(init<uint,object>());
    python_expression_class.def("argument_size", &PythonFunction::argument_size);
    python_expression_class.def("smoothness", &PythonFunction::smoothness);
    python_expression_class.def("__call__",(F(ExpressionInterface::*)(const FV&)const)&ExpressionInterface::evaluate);
    python_expression_class.def("__call__",(I(ExpressionInterface::*)(const IV&)const)&ExpressionInterface::evaluate);
    python_expression_class.def("evaluate",(F(ExpressionInterface::*)(const FV&)const)&ExpressionInterface::evaluate);
    python_expression_class.def("evaluate",(I(ExpressionInterface::*)(const IV&)const)&ExpressionInterface::evaluate);
    python_expression_class.def(self_ns::str(self));
}

void export_python_function()
{
    class_<PythonFunction, bases< FunctionInterface > > python_function_class("UserFunction", init<object>());
    python_function_class.def(init<uint,uint,object>());
    python_function_class.def("argument_size", &PythonFunction::argument_size);
    python_function_class.def("result_size", &PythonFunction::result_size);
    python_function_class.def("smoothness", &PythonFunction::smoothness);
    python_function_class.def("__call__",(FV(FunctionInterface::*)(const FV&)const)&FunctionInterface::evaluate);
    python_function_class.def("__call__",(IV(FunctionInterface::*)(const IV&)const)&FunctionInterface::evaluate);
    python_function_class.def("evaluate",(FV(FunctionInterface::*)(const FV&)const)&FunctionInterface::evaluate);
    python_function_class.def("evaluate",(IV(FunctionInterface::*)(const IV&)const)&FunctionInterface::evaluate);
    python_function_class.def("jacobian",(FMx(FunctionInterface::*)(const FV&)const)&FunctionInterface::jacobian);
    python_function_class.def("jacobian",(IMx(FunctionInterface::*)(const IV&)const)&FunctionInterface::jacobian);
    python_function_class.def(self_ns::str(self));
}

std::string __str__(const AffineExpression& f) {
    std::stringstream ss;
    ss<<"AffineExpression( ";
    if(f.b()!=0) { ss<<f.b(); }
    for(uint j=0; j!=f.argument_size(); ++j) {
        if(f.a()[j]!=0) {
            if(f.a()[j]>0) { ss<<"+"; } else { ss<<"-"; }
            if(abs(f.a()[j])!=1) { ss<<abs(f.a()[j]); }
            //ss<<char('x'+j);
            ss<<"x"<<j;
        }
    }
    ss<<" )";
    return ss.str();
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

void export_affine_expression()
{
    typedef Float F;
    typedef Interval I;
    typedef AffineExpression AE;
    typedef TaylorExpression TE;

    class_<AffineExpression, bases< ExpressionInterface > > affine_expression_class("AffineExpression", init<Vector<Float>, Float> ());
    affine_expression_class.def("argument_size", &AffineExpression::argument_size);
    affine_expression_class.def("smoothness", &AffineExpression::smoothness);
    affine_expression_class.def("__call__",(F(AffineExpression::*)(const FV&)const)&AffineExpression::evaluate);
    affine_expression_class.def("__call__",(I(AffineExpression::*)(const IV&)const)&AffineExpression::evaluate);
    affine_expression_class.def("__call__",(TE(*)(const AffineExpression&,const TFM&))&call_evaluate_expression<AffineExpression>);
    affine_expression_class.def("__add__",(AffineExpression(*)(const AffineExpression&,const AffineExpression&))&operator+);
    affine_expression_class.def("__mul__",(AffineExpression(*)(const AffineExpression&,const Interval&))&operator*);
    affine_expression_class.def("evaluate",(F(AffineExpression::*)(const FV&)const)&AffineExpression::evaluate);
    affine_expression_class.def("evaluate",(I(AffineExpression::*)(const IV&)const)&AffineExpression::evaluate);
    affine_expression_class.def("evaluate",(TE(*)(const AffineExpression&,const TFM&))&call_evaluate_expression<AffineExpression>);
    //affine_expression_class.def(self_ns::str(self));
    affine_expression_class.def("__str__",(std::string(*)(const AffineExpression&))&__str__);
    affine_expression_class.def("__repr__",(std::string(*)(const AffineExpression&))&__str__);

    def("derivative",(Interval(*)(const AffineExpression&,uint)) &derivative);
}

void export_affine_function()
{
    typedef AffineFunction AF;

    class_<AffineFunction, bases< FunctionInterface > > affine_function_class("AffineFunction", init<Matrix<Float>, Vector<Float> >());
    affine_function_class.def("argument_size", &AffineFunction::argument_size);
    affine_function_class.def("result_size", &AffineFunction::result_size);
    affine_function_class.def("smoothness", &AffineFunction::smoothness);
    affine_function_class.def("__getitem__",(AffineExpression(AffineFunction::*)(uint)const)&AffineFunction::operator[]);
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
    affine_function_class.def("__repr__",(std::string(*)(const AffineFunction&))&__str__);
}

void export_multi_index()
{
    class_< MultiIndex > multi_index_class("MultiIndex", no_init);
    //    multi_index_class.def("__getitem__",(MultiIndex::value_type(MultiIndex::*)(const MultiIndex::size_type&)const) &MultiIndex::operator[]);
    multi_index_class.def("degree",(MultiIndex::value_type(MultiIndex::*)()const) &MultiIndex::degree);
    multi_index_class.def(self_ns::str(self));
}


template<class X>
void export_polynomial()
{
    typedef ExpansionValue<X> M;
    typedef Polynomial<X> P;

    class_< M > monomial_class(python_name<X>("Monomial"), no_init);
    monomial_class.def("key",(MultiIndex&(M::*)())&M::key,return_value_policy<reference_existing_object>());
    //monomial_class.def("data",&monomial_data<X>,return_value_policy<reference_existing_object>());
    //monomial_class.def("data",(X&(M::*)()) &M::data,return_value_policy<reference_existing_object>());
    monomial_class.def("data",(const X&(M::*)()const) &M::data,return_value_policy<copy_const_reference>());
    monomial_class.def(self_ns::str(self));

    //double flt;
    X real;
    class_< P > polynomial_class(python_name<X>("Polynomial"), init<int>());
    polynomial_class.def("constant", (P(*)(uint,double)) &P::constant);
    polynomial_class.staticmethod("constant");
    polynomial_class.def("variable", (P(*)(uint,uint)) &P::variable);
    polynomial_class.staticmethod("variable");

    polynomial_class.def("argument_size", &P::argument_size);
    polynomial_class.def(+self);
    polynomial_class.def(-self);
    polynomial_class.def(self+self);
    polynomial_class.def(self-self);
    polynomial_class.def(self*self);
    //polynomial_class.def(self+flt);
    //polynomial_class.def(self-flt);
    //polynomial_class.def(self*flt);
    //polynomial_class.def(self/flt);
    //polynomial_class.def(flt+self);
    //polynomial_class.def(flt-self);
    //polynomial_class.def(flt*self);
    polynomial_class.def(self+real);
    polynomial_class.def(self-real);
    polynomial_class.def(self*real);
    polynomial_class.def(self/real);
    polynomial_class.def(real+self);
    polynomial_class.def(real-self);
    polynomial_class.def(real*self);
    polynomial_class.def("__iter__",boost::python::iterator<P>());
    polynomial_class.def(self_ns::str(self));
    polynomial_class.def("__mul__",&__mul__< Vector<Polynomial<X> >, Polynomial<X>, Vector<Float> >);

    class_< Vector<P> > polynomial_function_class(python_name<X>("PolynomialVector"), no_init);
    polynomial_function_class.def("__getitem__",(const P&(Vector<P>::*)(size_t)const)&Vector<P>::operator[],return_value_policy<copy_const_reference>());
    polynomial_function_class.def(self_ns::str(self));
}


void export_polynomial_expression()
{
    typedef PolynomialExpression PE;
    Float flt;
    Interval ivl;

    class_< PE, bases<ExpressionInterface> > polynomial_expression_class("PolynomialExpression", init<int,int>());
    polynomial_expression_class.def("constant", (PE(*)(uint,double)) &PE::constant);
    polynomial_expression_class.staticmethod("constant");
    polynomial_expression_class.def("variable", (PE(*)(uint,uint)) &PE::variable);
    polynomial_expression_class.staticmethod("variable");

    polynomial_expression_class.def("argument_size", &PE::argument_size);

    //polynomial_expression_class.def(+self);
    //polynomial_expression_class.def(-self);
    //polynomial_expression_class.def(self+self);
    //polynomial_expression_class.def(self-self);
    //polynomial_expression_class.def(self*self);

    //polynomial_expression_class.def(self+flt);
    //polynomial_expression_class.def(self-flt);
    //polynomial_expression_class.def(self*flt);
    //polynomial_expression_class.def(self/flt);
    //polynomial_expression_class.def(flt+self);
    //polynomial_expression_class.def(flt-self);
    //polynomial_expression_class.def(flt*self);

    polynomial_expression_class.def("__pos__", &__pos__<PE,PE>);
    polynomial_expression_class.def("__neg__", &__neg__<PE,PE>);
    polynomial_expression_class.def("__add__", &__add__<PE,PE,PE>);
    polynomial_expression_class.def("__sub__", &__sub__<PE,PE,PE>);
    polynomial_expression_class.def("__mul__", &__mul__<PE,PE,PE>);

    polynomial_expression_class.def("__add__", &__add__<PE,PE,Float>);
    polynomial_expression_class.def("__sub__", &__sub__<PE,PE,Float>);
    polynomial_expression_class.def("__mul__", &__mul__<PE,PE,Float>);
    polynomial_expression_class.def("__div__", &__div__<PE,PE,Float>);
    polynomial_expression_class.def("__radd__", &__add__<PE,PE,Float>);
    polynomial_expression_class.def("__rsub__", &__rsub__<PE,PE,Float>);
    polynomial_expression_class.def("__rmul__", &__mul__<PE,PE,Float>);

    polynomial_expression_class.def("__add__", &__add__<PE,PE,Interval>);
    polynomial_expression_class.def("__sub__", &__sub__<PE,PE,Interval>);
    polynomial_expression_class.def("__mul__", &__mul__<PE,PE,Interval>);
    polynomial_expression_class.def("__div__", &__div__<PE,PE,Interval>);
    polynomial_expression_class.def("__radd__", &__add__<PE,PE,Interval>);
    polynomial_expression_class.def("__rsub__", &__rsub__<PE,PE,Interval>);
    polynomial_expression_class.def("__rmul__", &__mul__<PE,PE,Interval>);

    polynomial_expression_class.def(self+=self);
    polynomial_expression_class.def(self-=self);
    polynomial_expression_class.def(self+=flt);
    polynomial_expression_class.def(self-=flt);
    polynomial_expression_class.def(self*=flt);
    polynomial_expression_class.def(self/=flt);
    polynomial_expression_class.def(self+=ivl);
    polynomial_expression_class.def(self-=ivl);
    polynomial_expression_class.def(self*=ivl);
    polynomial_expression_class.def(self/=ivl);

    polynomial_expression_class.def(self_ns::str(self));

}

void function_submodule() {
    export_multi_index();

    //export_polynomial<Float>();
    //export_polynomial<Interval>();

    export_expression_interface();
    export_affine_expression();
    export_polynomial_expression();
    export_python_expression();

    export_function_interface();
    export_affine_function();
    export_python_function();
}

