/***************************************************************************
 *            function.h
 *
 *  Copyright 2008-9  Pieter Collins
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

/*! \file function.h
 *  \brief Built-in and user functions and expressions
 */

#ifndef ARIADNE_FUNCTION_H
#define ARIADNE_FUNCTION_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "function_interface.h"

#include "macros.h"
#include "pointer.h"

#include "vector.h"
#include "matrix.h"
#include "polynomial.h"
#include "taylor_model.h"
#include "differential.h"
#include "operators.h"
#include "expression.h"
#include "affine.h"

namespace Ariadne {


template<uint AS, uint PS=0u, uint SM=255u>
struct ScalarFunctionData
{
    const uint argument_size() const { return AS; }
    const uint parameter_size() const { return PS; }
    const uint smoothness() const { return SM; }
};


//! \brief A wrapper for converting templated C++ functions to %Ariadne functions.
//!
//! Given a C++ class T with a (static) template method
//!   <code>template<class R, class A, class P> compute(R& r, const A& a, const P& p);</code>
//! the type <code>UserExpression<T></code> is an Ariadne expression defined by \f$r=f(a)\f$.
//! The constructor for UserExpression<T> takes a Vector<Float> or Vector<Interval> argument which is used for \a p.
//!
//! The class T must also define meta-data <c>argument_size(), parameter_size()
//! and smoothness()</c>. These are most easily defined by inheriting from the
//! <tt>ExpressionData<AS,PS,SM=SMOOTH></tt> class.
//!
//! The constant \a SMOOTH is used for an arbitrarily-differentiable function.
template<class T> class UserScalarFunction
    : public ScalarFunctionInterface
{
  public:
    UserScalarFunction() : _p(this->parameter_size()) { }
    UserScalarFunction(const Vector<Float>& p) : _p(p) { }
    UserScalarFunction(const Vector<Interval>& p) : _p(p) { }

    const Vector<Interval> parameters() const { return _p; }

    virtual SizeType argument_size() const { return T::argument_size(); }
    virtual SizeType parameter_size() const { return T::parameter_size(); }
    virtual SmoothnessType smoothness() const { return T::smoothness(); }

    virtual Float evaluate(const Vector<Float>& x) const {
        Float r=0; T::compute(r,x,midpoint(_p)); return r; }
    virtual Interval evaluate(const Vector<Interval>& x) const {
        Interval r=0; T::compute(r,x,_p); return r; }

    virtual TaylorModel evaluate(const Vector<TaylorModel>& x) const {
        TaylorModel r(x[0].argument_size(),x[0].accuracy_ptr()); T::compute(r,x,_p); return r; }

    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const {
        Differential<Float> r(x[0].argument_size(),x[0].degree()); T::compute(r,x,midpoint(_p)); return r; }
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const {
        Differential<Interval> r(x[0].argument_size(),x[0].degree()); T::compute(r,x,_p); return r; }

    virtual UserScalarFunction<T>* derivative(uint j) { ARIADNE_NOT_IMPLEMENTED; }

    virtual Matrix<Float> gradient(const Vector<Float>& x) const {
        return this->evaluate(Differential<Float>::variables(1u,x)).gradient(); }
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const {
        return this->evaluate(Differential<Interval>::variables(1u,x)).gradient(); }

    virtual std::ostream& write(std::ostream& os) const  {
        return os << "UserExpression( argument_size="<<this->argument_size()<<" )"; }
  private:
    Vector<Interval> _p;
};




//! A polynomial expression \f$p:\R^n\rightarrow\R\f$.
class ScalarPolynomialFunction
    : public ScalarFunctionInterface
    , public Polynomial<Interval>
{
  public:
    ScalarPolynomialFunction() : Polynomial<Interval>() { }
    explicit ScalarPolynomialFunction(uint n) : Polynomial<Interval>(n) { }
    ScalarPolynomialFunction(const Polynomial<Float>& p) : Polynomial<Interval>(p) { }
    ScalarPolynomialFunction(const Polynomial<Interval>& p) : Polynomial<Interval>(p) { }
    ScalarPolynomialFunction& operator=(const Polynomial<Interval>& p) { this->Polynomial<Interval>::operator=(p); return *this; }

    ScalarPolynomialFunction(RealExpressionPointer e, const Array<String>& s);

    ScalarPolynomialFunction(const ScalarAffineFunction& a) : Polynomial<Interval>(a.argument_size()) {
        uint n=this->argument_size(); (*this)[MultiIndex::zero(n)]=a.b();
        for(uint i=0; i!=n; ++i) { (*this)[MultiIndex::unit(n,i)]=a.a()[i]; } }

    static ScalarPolynomialFunction constant(uint n, Float c) {
        return Polynomial<Interval>::constant(n,c); }
    static ScalarPolynomialFunction constant(uint n, Interval c) {
        return Polynomial<Interval>::constant(n,c); }
    static ScalarPolynomialFunction variable(uint n, uint j) {
        return Polynomial<Interval>::variable(n,j); }
    static array<ScalarPolynomialFunction> variables(uint n) {
        array<ScalarPolynomialFunction> r(n); for(uint i=0; i!=n; ++i) { r[i]=variable(n,i); } return r; }

    virtual ScalarPolynomialFunction* clone() const { return new ScalarPolynomialFunction(*this); }

    virtual ScalarFunctionInterface::SizeType argument_size() const {
        return static_cast<const Polynomial<Interval>&>(*this).argument_size(); }
    virtual ScalarFunctionInterface::SmoothnessType smoothness() const {
        return SMOOTH; }
    virtual Float evaluate(const Vector<Float>& x) const {
        return Ariadne::evaluate(midpoint(static_cast<const Polynomial<Interval>&>(*this)),x); }
    virtual Interval evaluate(const Vector<Interval>& x) const {
        return Ariadne::evaluate(static_cast<const Polynomial<Interval>&>(*this),x); }
    virtual TaylorModel evaluate(const Vector<TaylorModel>& x) const {
        return Ariadne::evaluate(static_cast<const Polynomial<Interval>&>(*this),x); }
    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const {
        return Ariadne::evaluate(midpoint(static_cast<const Polynomial<Interval>&>(*this)),x); }
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const {
        return Ariadne::evaluate(static_cast<const Polynomial<Interval>&>(*this),x); }
    virtual Differential<TaylorModel> evaluate(const Vector< Differential<TaylorModel> >& x) const {
        return Ariadne::evaluate(static_cast<const Polynomial<Interval>&>(*this),x); }

    virtual ScalarPolynomialFunction* derivative(uint j) const {
        return new ScalarPolynomialFunction(Ariadne::derivative(*this,j)); }
    virtual ScalarPolynomialFunction* antiderivative(uint j) const {
        return new ScalarPolynomialFunction(Ariadne::antiderivative(*this,j)); }

    virtual Vector<ScalarPolynomialFunction> gradient() const {
        Vector<ScalarPolynomialFunction> g(this->argument_size());
        for(uint i=0; i!=g.size(); ++i) { g[i]=Ariadne::derivative(*this,i); }
        return g; }

    virtual std::ostream& write(std::ostream& os) const;
};

inline std::ostream& ScalarPolynomialFunction::write(std::ostream& os) const {
    const ScalarPolynomialFunction& p=*this;
    bool first_term=true;
    //bool first_factor=true;
    os<<"P(";
    Float r;
    if(p.expansion().size()==0) {
        r=0.0;
        os << "0";
    } else {
        r=radius(p.begin()->data());
        for(Polynomial<Interval>::const_iterator iter=p.begin(); iter!=p.end(); ++iter) {
            MultiIndex a=iter->key();
            Float v=midpoint(iter->data());
            if(abs(v)<1e-15) { r+=abs(v); v=0; }
            if(v!=0) {
                if(v>0 && !first_term) { os<<"+"; }
                first_term=false;
                bool first_factor=true;
                if(v<0) { os<<"-"; }
                if(abs(v)!=1 || a.degree()==0) { os<<abs(v); first_factor=false; }
                for(uint j=0; j!=a.size(); ++j) {
                    if(a[j]!=0) {
                        if(first_factor) { first_factor=false; } else { os <<"*"; }
                        os<<"x"<<j; if(a[j]!=1) { os<<"^"<<int(a[j]); } }
                }
            }
        }
    }
    if(r>0) { os<<"+/-"<<r; }
    os<<")";
    return os;
}

/*

inline ScalarPolynomialFunction::ScalarPolynomialFunction(ExpressionPointer e, const Array<String>& s) {
    ScalarPolynomialFunction& r=*this;
    const ExpressionInterface* eptr=e.operator->();
    if(dynamic_cast<const ConstantExpression*>(eptr)) {
        const ConstantExpression* cptr=dynamic_cast<const ConstantExpression*>(eptr);
        r =  ScalarPolynomialFunction::constant(s.size(),cptr->value());
    } else if(dynamic_cast<const CoordinateExpression*>(eptr)) {
        const CoordinateExpression* pptr=dynamic_cast<const CoordinateExpression*>(eptr);
        r = ScalarPolynomialFunction::variable(s.size(),pptr->coordinate());
        //return PolynomialExpression::variable(pptr->argument_size(),pptr->coordinate());
    } else if(dynamic_cast<const UnaryExpression<Neg>*>(eptr)) {
        const UnaryExpression<Neg>* uptr=dynamic_cast<const UnaryExpression<Neg>*>(eptr);
        r = -ScalarPolynomialFunction(uptr->_arg_ptr,s);
    } else if(dynamic_cast<const BinaryExpression<Add>*>(eptr)) {
        const BinaryExpression<Add>* bptr=dynamic_cast<const BinaryExpression<Add>*>(eptr);
        r =  ScalarPolynomialFunction(bptr->_arg1_ptr,s)+ScalarPolynomialFunction(bptr->_arg2_ptr,s);
    } else if(dynamic_cast<const BinaryExpression<Sub>*>(eptr)) {
        const BinaryExpression<Sub>* bptr=dynamic_cast<const BinaryExpression<Sub>*>(eptr);
        r =  ScalarPolynomialFunction(bptr->_arg1_ptr,s)-ScalarPolynomialFunction(bptr->_arg2_ptr,s);
    } else if(dynamic_cast<const BinaryExpression<Mul>*>(eptr)) {
        const BinaryExpression<Mul>* bptr=dynamic_cast<const BinaryExpression<Mul>*>(eptr);
        r =  ScalarPolynomialFunction(bptr->_arg1_ptr,s)*ScalarPolynomialFunction(bptr->_arg2_ptr,s);
    } else if(dynamic_cast<const BinaryExpression<Div>*>(eptr)) {
        const BinaryExpression<Div>* bptr=dynamic_cast<const BinaryExpression<Div>*>(eptr);
        const ConstantExpression* cptr=dynamic_cast<const ConstantExpression*>(bptr->_arg2_ptr.operator->());
        ARIADNE_ASSERT_MSG(cptr,"Cannot convert expression "<<*e<<" to polynomial.");
        r =  ScalarPolynomialFunction(bptr->_arg1_ptr,s)*(1/cptr->value());
    } else {
        ARIADNE_ASSERT_MSG(false,"Cannot convert expression "<<*e<<" to polynomial.");
    }
}

*/

inline ScalarPolynomialFunction::ScalarPolynomialFunction(RealExpressionPointer e, const Array<String>& s) {
    ScalarPolynomialFunction& r=*this;
    const ExpressionInterface<Real>* eptr=e.operator->();
    if(dynamic_cast<const ConstantExpression<Real>*>(eptr)) {
        const ConstantExpression<Real>* cptr=dynamic_cast<const ConstantExpression<Real>*>(eptr);
        r =  ScalarPolynomialFunction::constant(s.size(),cptr->value());
    } else if(dynamic_cast<const CoordinateExpression<Real>*>(eptr)) {
        const CoordinateExpression<Real>* pptr=dynamic_cast<const CoordinateExpression<Real>*>(eptr);
        r = ScalarPolynomialFunction::variable(s.size(),pptr->coordinate());
        //return PolynomialExpression::variable(pptr->argument_size(),pptr->coordinate());
    } else if(dynamic_cast<const UnaryExpression<Real,Neg>*>(eptr)) {
        const UnaryExpression<Real,Neg>* uptr=dynamic_cast<const UnaryExpression<Real,Neg>*>(eptr);
        r = -ScalarPolynomialFunction(uptr->_arg_ptr,s);
    } else if(dynamic_cast<const BinaryExpression<Real,Add>*>(eptr)) {
        const BinaryExpression<Real,Add>* bptr=dynamic_cast<const BinaryExpression<Real,Add>*>(eptr);
        r =  ScalarPolynomialFunction(bptr->_arg1_ptr,s)+ScalarPolynomialFunction(bptr->_arg2_ptr,s);
    } else if(dynamic_cast<const BinaryExpression<Real,Sub>*>(eptr)) {
        const BinaryExpression<Real,Sub>* bptr=dynamic_cast<const BinaryExpression<Real,Sub>*>(eptr);
        r =  ScalarPolynomialFunction(bptr->_arg1_ptr,s)-ScalarPolynomialFunction(bptr->_arg2_ptr,s);
    } else if(dynamic_cast<const BinaryExpression<Real,Mul>*>(eptr)) {
        const BinaryExpression<Real,Mul>* bptr=dynamic_cast<const BinaryExpression<Real,Mul>*>(eptr);
        r =  ScalarPolynomialFunction(bptr->_arg1_ptr,s)*ScalarPolynomialFunction(bptr->_arg2_ptr,s);
    } else if(dynamic_cast<const BinaryExpression<Real,Div>*>(eptr)) {
        const BinaryExpression<Real,Div>* bptr=dynamic_cast<const BinaryExpression<Real,Div>*>(eptr);
        const ConstantExpression<Real>* cptr=dynamic_cast<const ConstantExpression<Real>*>(bptr->_arg2_ptr.operator->());
        ARIADNE_ASSERT_MSG(cptr,"Cannot convert expression "<<*e<<" to polynomial.");
        r =  ScalarPolynomialFunction(bptr->_arg1_ptr,s)*(1/cptr->value());
    } else {
        ARIADNE_ASSERT_MSG(false,"Cannot convert expression "<<*e<<" to polynomial.");
    }
}




// A wrapper for classes with non-static _compute and _compute_approx methods
template<class F>
class FunctionTemplate
    : public FunctionInterface
{
  private:
    template<class R, class A> void _base_compute(R& r, const A& a) const {
        static_cast<const F*>(this)->_compute(r,a); }
    template<class R, class A> void _base_compute_approx(R& r, const A& a) const {
        static_cast<const F*>(this)->_compute_approx(r,a); }
  protected:
    FunctionTemplate() { }
  public:
    virtual SmoothnessType smoothness() const { return SMOOTH; }

    virtual Vector<Float> evaluate(const Vector<Float>& x) const {
        Vector<Float> r(this->result_size()); _base_compute_approx(r,x); return r; }
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const {
        Vector<Interval> r(this->result_size()); _base_compute(r,x); return r; }

    virtual Vector<TaylorModel> evaluate(const Vector<TaylorModel>& x) const {
        Vector<TaylorModel> r(this->result_size(),TaylorModel(x[0].argument_size(),x[0].accuracy_ptr()));
        _base_compute(r,x); return r; }

    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const {
        Vector< Differential<Float> > r(this->result_size(),Differential<Float>(x[0].argument_size(),x[0].degree()));
        _base_compute_approx(r,x); return r; }
    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const {
        Vector< Differential<Interval> > r(this->result_size(),Differential<Interval>(x[0].argument_size(),x[0].degree()));
        _base_compute(r,x); return r; }

    virtual Matrix<Float> jacobian(const Vector<Float>& x) const {
        return Ariadne::jacobian(this->evaluate(Differential<Float>::variables(1u,x))); }
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const {
        return Ariadne::jacobian(this->evaluate(Differential<Interval>::variables(1u,x))); }

    // TODO: Find a better way for writing functions which can handle transformations which may not have a
    // write() method or operator<<.
    virtual std::ostream& write(std::ostream& os) const  {
        return os << "Function( result_size="<<this->result_size()<<", argument_size="<<this->argument_size()<<" )"; }
};


//! \brief A convenience class defining the meta-data of an %Ariadne function.
template<uint RS, uint AS, uint PS=0u, uint SM=255u>
class FunctionData
{
  public:
    //!
    static const uint result_size() { return RS; }
    //!
    static const uint argument_size() { return AS; }
    //!
    static const uint parameter_size() { return PS; }
    //!
    static const uint smoothness() { return SM; }
};


//! \brief A wrapper for converting templated C++ functions to %Ariadne functions.
//!
//! Given a C++ class T with a (static) template method
//!   <code>template<class R, class A, class P> compute(R& r, const A& a, const P& p);</code>
//! the type <code>Function<T></code> is an Ariadne function defined by \f$r=f(a)\f$.
//! The constructor for Function<T> takes a Vector<Float> argument which is used for \a p.
//!
//! The class T must also define meta-data <c>result_size(), argument_size(), parameter_size()
//! and smoothness()</c> as static functions. These are most easily defined by inheriting from the
//! <tt>FunctionData<RS,AS,PS,SM=SMOOTH></tt> class.
//!
//! The constant \a SMOOTH is used for an arbitrarily-differentiable function.
template<class T> class UserFunction
    : public FunctionInterface
{
  public:
    UserFunction() : _p(this->parameter_size()) { }
    UserFunction(const Vector<Float>& p) : _p(p) { }
    UserFunction(const Vector<Interval>& p) : _p(p) { }

    const Vector<Interval> parameters() const { return _p; }

    virtual UserFunction<T>* clone() const { return new UserFunction<T>(*this); }

    virtual SizeType result_size() const { return T::result_size(); }
    virtual SizeType argument_size() const { return T::argument_size(); }
    virtual SizeType parameter_size() const { return T::parameter_size(); }
    virtual SmoothnessType smoothness() const { return T::smoothness(); }

    virtual Vector<Float> evaluate(const Vector<Float>& x) const {
        Vector<Float> r(this->result_size(),0.0); T::compute(r,x,midpoint(_p)); return r; }
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const {
        Vector<Interval> r(this->result_size(),0.0); T::compute(r,x,_p); return r; }

    virtual Vector<TaylorModel> evaluate(const Vector<TaylorModel>& x) const {
        Vector<TaylorModel> r(this->result_size(),TaylorModel(x[0].argument_size(),x[0].accuracy_ptr()));
        T::compute(r,x,_p); return r; }

    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const {
        Vector< Differential<Float> > r(this->result_size(),Differential<Float>(x[0].argument_size(),x[0].degree()));
        T::compute(r,x,midpoint(_p)); return r; }
    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const {
        Vector< Differential<Interval> > r(this->result_size(),Differential<Interval>(x[0].argument_size(),x[0].degree()));
        T::compute(r,x,_p); return r; }

    virtual Matrix<Float> jacobian(const Vector<Float>& x) const {
        return Ariadne::jacobian(this->evaluate(Differential<Float>::variables(1u,x))); }
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const {
        return Ariadne::jacobian(this->evaluate(Differential<Interval>::variables(1u,x))); }

    // TODO: Find a better way for writing functions which can handle transformations which may not have a
    // write() method or operator<<.
    virtual std::ostream& write(std::ostream& os) const  {
        return os << "UserFunction( result_size="<<this->result_size()<<", argument_size="<<this->argument_size()<<" )"; }
  private:
    Vector<Interval> _p;
};


//! \brief A wrapper for converting templated C++ functions to %Ariadne functions.
//! \deprecated Use UserFunction instead.
template<class T> class Function
    : public UserFunction<T>
{
  public:
    Function() : UserFunction<T>() { _deprecated_warning(); }
    Function(const Vector<Float>& p) : UserFunction<T>(p) { _deprecated_warning(); }
    Function(const Vector<Interval>& p) : UserFunction<T>(p) { _deprecated_warning(); }
  private:
    void _deprecated_warning() {
        static bool firsttime=true;
        if(firsttime) { firsttime=false;
        std::cerr<<"Deprecated: template<class T> Function is deprecated; use template<class T> UserFunction instead.\n"; }
    }
};

//! \brief A vector-valued function defined by its scalar-valued components.
template<>
class Vector<ScalarFunctionInterface>
    : public FunctionInterface
{
    friend class FunctionTemplate< Vector<ScalarFunctionInterface> >;
  public:
    typedef unsigned int SizeType;
    typedef unsigned short int SmoothnessType;

    Vector(SizeType n) : _expressions(n) { }
    SizeType size() const { return _expressions.size(); }
    const ScalarFunctionInterface& operator[](SizeType i) const { return *_expressions[i]; }
    void set(SizeType i, const ScalarFunctionInterface& f) {
        ARIADNE_ASSERT(i<this->result_size());
        ARIADNE_ASSERT(this->size()==0 || !_expressions[0] || f.argument_size()==this->argument_size());
        _expressions[i]=shared_ptr<const ScalarFunctionInterface>(f.clone()); }
    void set(SizeType i, shared_ptr<const ScalarFunctionInterface> p) {
        ARIADNE_ASSERT(i<this->result_size());
        ARIADNE_ASSERT(this->size()==0 || !_expressions[0] || p->argument_size()==this->argument_size()); _expressions[i]=p; }
    void set(SizeType i, const ScalarFunctionInterface* p) {
        ARIADNE_ASSERT(i<this->result_size());
        ARIADNE_ASSERT(this->size()==0 || !_expressions[0] || p->argument_size()==this->argument_size());
        _expressions[i]=shared_ptr<const ScalarFunctionInterface>(p); }
    const ScalarFunctionInterface& get(uint i) const {
        ARIADNE_ASSERT(i<this->result_size());
        return *this->_expressions[i]; }

    virtual Vector<ScalarFunctionInterface>* clone() const { return new Vector<ScalarFunctionInterface>(*this); }

    virtual SizeType result_size() const { return _expressions.size(); }
    virtual SizeType argument_size() const { if(_expressions[0]) { return _expressions[0]->argument_size(); } else { return 0; } }
    virtual SmoothnessType smoothness() const {
        SmoothnessType res=_expressions[0]->smoothness();
        for(uint i=1; i!=this->size(); ++i) {
            res=max(res,_expressions[i]->smoothness()); } return res; }

    virtual Vector<Float> evaluate(const Vector<Float>& x) const {
        Vector<Float> r(result_size()); _compute_approx(r,x); return r; }
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const {
        Vector<Interval> r(result_size()); _compute(r,x); return r; }
    virtual Vector<TaylorModel> evaluate(const Vector<TaylorModel>& x) const {
        Vector<TaylorModel> r(result_size()); _compute(r,x); return r; }
    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const {
        Vector< Differential<Float> > r(result_size()); _compute_approx(r,x); return r; }
    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const {
        Vector< Differential<Interval> > r(result_size()); _compute_approx(r,x); return r; }

    virtual Matrix<Float> jacobian(const Vector<Float>& x) const {
        return Ariadne::jacobian(this->evaluate(Differential<Float>::variables(1u,x))); }
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const {
        return Ariadne::jacobian(this->evaluate(Differential<Interval>::variables(1u,x))); }

    virtual std::ostream& write(std::ostream& os) const {
        os<<"["; for(uint i=0; i!=this->result_size(); ++i) { if(i!=0) { os<<","; } os<<(*this)[i]; } return os<<"]"; }
  private:
    template<class Y> inline void _compute(Vector<Y>& r, const Vector<Y>& x) const {
        for(SizeType i=0; i!=this->result_size(); ++i) {
            r[i]=this->_expressions[i]->evaluate(x); } }
    template<class Y> inline void _compute_approx(Vector<Y>& r, const Vector<Y>& x) const {
        _compute(r,x); }
  private:
    array< boost::shared_ptr<const ScalarFunctionInterface> > _expressions;
};



//! A polynomial function.
class PolynomialFunction
    : public Vector<ScalarFunctionInterface>
{
    typedef Vector<ScalarFunctionInterface> Base;
  public:
    explicit PolynomialFunction(uint rs, uint as) : Base(rs) {
        for(uint i=0; i!=rs; ++i) { Base::set(i,ScalarPolynomialFunction(as)); } }
    explicit PolynomialFunction(uint rs, const ScalarPolynomialFunction& p) : Base(rs) {
        for(uint i=0; i!=rs; ++i) { Base::set(i,p); } }

    template<class E> PolynomialFunction(const ublas::vector_expression<E>& p) : Base(p().size()) {
        for(uint i=0; i!=p().size(); ++i) { Base::set(i,PolynomialExpression(p()[i])); } }

    PolynomialFunction(const Vector< Polynomial<Float> >& p) : Base(p.size()) {
        for(uint i=0; i!=p.size(); ++i) { Base::set(i,ScalarPolynomialFunction(p[i])); } }
    PolynomialFunction(const Vector< Polynomial<Interval> >& p) : Base(p.size()) {
        for(uint i=0; i!=p.size(); ++i) { Base::set(i,ScalarPolynomialFunction(p[i])); } }
    PolynomialFunction(const Vector<ScalarPolynomialFunction>& p) : Base(p.size()) {
        for(uint i=0; i!=p.size(); ++i) { Base::set(i,p[i]); } }

    const ScalarPolynomialFunction& operator[](uint i) const {
        return dynamic_cast<const ScalarPolynomialFunction&>(Base::get(i)); }

    const ScalarPolynomialFunction& get(uint i) const {
        return dynamic_cast<const ScalarPolynomialFunction&>(Base::get(i)); }
    void set(uint i,const ScalarPolynomialFunction& p) {
        return this->Base::set(i,p); }

    virtual std::ostream& write(std::ostream& os) const {
        return os << "PolynomialFunction("<<static_cast<const Vector<ScalarFunctionInterface>&>(*this)<<")"; }
  private:
    friend PolynomialFunction operator+(const PolynomialFunction&, const PolynomialFunction&);
    friend PolynomialFunction flip(const PolynomialFunction&, uint);
};

inline PolynomialFunction join(const ScalarPolynomialFunction& p1, const ScalarPolynomialFunction& p2) {
    return join(static_cast<const Polynomial<Interval>&>(p1),static_cast<const Polynomial<Interval>&>(p2));
}

inline PolynomialFunction operator*(const ScalarPolynomialFunction& p, const Vector<Float>& e) {
    PolynomialFunction r(e.size(),p.argument_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,static_cast<ScalarPolynomialFunction>(e[i]*p)); }
    return r;
}

inline PolynomialFunction operator+(const PolynomialFunction& p1, const PolynomialFunction& p2) {
    ARIADNE_ASSERT(p1.result_size()==p2.result_size());
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    PolynomialFunction r(p1.result_size(),p1.argument_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,static_cast<ScalarPolynomialFunction>(p1[i]+p2[i])); }
    return r;
}

inline PolynomialFunction operator-(const PolynomialFunction& p1, const PolynomialFunction& p2) {
    ARIADNE_ASSERT(p1.result_size()==p2.result_size());
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    PolynomialFunction r(p1.result_size(),p1.argument_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,static_cast<ScalarPolynomialFunction>(p1[i]-p2[i])); }
    return r;
}







//! A constant function \f$ x\mapsto c\f$ from \f$R^m\f$ to \f$R^n\f$.
class ConstantFunction
    : public FunctionTemplate<ConstantFunction>
{
  public:
    ConstantFunction(uint rs, double c0, ...) : _as(), _c(rs) {
        ARIADNE_ASSERT(rs>0); va_list args; va_start(args,c0);
        _c[0]=c0; for(size_t i=1; i!=rs; ++i) { _c[i]=va_arg(args,double); } _as=va_arg(args,int);
        va_end(args);
    }

    ConstantFunction(uint rs, Interval c0, ...) : _as(), _c(rs) {
        ARIADNE_ASSERT(rs>0); double l,u; va_list args; va_start(args,c0);
        _c[0]=c0; for(size_t i=1; i!=rs; ++i) { l=va_arg(args,double); u=va_arg(args,double); _c[i]=Interval(l,u); } _as=va_arg(args,int);
        va_end(args);
    }

    ConstantFunction(const Vector<Float>& c, uint as) : _as(as), _c(c) { }
    ConstantFunction(const Vector<Interval>& c, uint as) : _as(as), _c(c) { }
//    ConstantExpression<Real,Real> operator[](uint i) const { return ConstantExpression<Real,Real>(_as,this->_c[i]); }
    Interval operator[](uint i) const { return this->_c[i]; }
    const Vector<Interval>& c() const { return _c; }

    ConstantFunction* clone() const { return new ConstantFunction(*this); }

    SizeType result_size() const { return _c.size(); }
    SizeType argument_size() const { return _as; }

    std::ostream& write(std::ostream& os) const {
        return os << "ConstantFunction( argument_size=" << this->argument_size()
                  << ", c=" << this->c() << " )"; }
  private:
    friend class FunctionTemplate<ConstantFunction>;
    template<class R, class A> void _compute(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) { r[i]=_c[i]; } }
    template<class R, class A> void _compute_approx(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) { r[i]=midpoint(_c[i]); } }
  private:
    uint _as;
    Vector<Interval> _c;
};


//! \brief A projection function \f$ x'_i= x_{p(i)}\f$.
class ProjectionFunction
    : public FunctionTemplate<ProjectionFunction>
{
  public:
    //! \brief Construct the identity function in dimension \a n.
    ProjectionFunction(uint n) : _as(n), _p(n) {
        for(uint i=0; i!=n; ++i) { _p[i]=i; } }
    //! \brief Construct the projection functions \f$f_i(x)=x_{i+k}\f$ for \f$i=0,\ldots,m-1\f$. Precondition: \f$m+k\leq n\f$.
    ProjectionFunction(uint m, uint n, uint k) : _as(n), _p(m) {
        ARIADNE_ASSERT(m+k<=n); for(uint j=0; j!=m; ++j) { _p[j]=k+j; } }
    //! \brief Construct the projection function  with \f$f_i(x)=x_{p_i}\f$ for \f$i=0,\ldots,m-1\f$.
    ProjectionFunction(uint m, uint n, const array<uint>& p) : _as(n), _p(p) {
        ARIADNE_ASSERT(p.size()==m); for(uint i=0; i!=_p.size(); ++i) { ARIADNE_ASSERT(p[i]<n); } }
    //! \brief Construct the projection function with \f$f_i(x)=x_{p_i}\f$ for \f$i=0,\ldots,|p|-1\f$.
    ProjectionFunction(const array<uint>& p, uint n) : _as(n), _p(p) {
        for(uint i=0; i!=_p.size(); ++i) { ARIADNE_ASSERT(p[i]<n); } }
    ProjectionFunction(const Range& rng, uint as) : _as(as), _p(rng.size()) {
        ARIADNE_ASSERT(rng.start()+rng.size()<=as);
        for(uint i=0; i!=_p.size(); ++i) { _p[i]=rng.start()+i; } }
    std::ostream& write(std::ostream& os) const {
        return os << "ProjectionFunction( p=" << this->_p << " )"; }

    ProjectionFunction* clone() const { return new ProjectionFunction(*this); }

    SizeType result_size() const { return this->_p.size(); }
    SizeType argument_size() const { return this->_as; }
    SizeType p(unsigned int j) const { return this->_p[j]; }

    CoordinateExpression<Real> operator[](unsigned int i) const { return CoordinateExpression<Real>(_as,_p[i]); }

  private:
    friend class FunctionTemplate<ProjectionFunction>;
    template<class R, class A> void _compute(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) { r[i]=x[_p[i]]; } }
    template<class R, class A> void _compute_approx(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) { r[i]=x[_p[i]]; } }
  private:
    uint _as;
    array<uint> _p;
};


//! \brief The identity function \f$ x\mapsto x\f$ in \f$\R^n\f$.
class IdentityFunction
    : public FunctionTemplate<IdentityFunction>
{
  public:
    //! \brief Construct the identity function in dimension \a n.
    IdentityFunction(uint n) : _n(n) { }

    IdentityFunction* clone() const { return new IdentityFunction(*this); }

    CoordinateExpression<Real> operator[](unsigned int i) const { return CoordinateExpression<Real>(_n,i); }

    SizeType result_size() const { return this->_n; }
    SizeType argument_size() const { return this->_n; }

    std::ostream& write(std::ostream& os) const {
        return os << "IdentityFunction( size=" << this->result_size() << " )"; }
  private:
    friend class FunctionTemplate<IdentityFunction>;
    template<class R, class A> void _compute(R& r, const A& x) const { r=x; }
    template<class R, class A> void _compute_approx(R& r, const A& x) const { r=x; }
  private:
    uint _n;
};


//! \brief A scaling function \f$x_i' = o_i+l_ix_i\f$.
class ScalingFunction
    : public FunctionTemplate<ScalingFunction>
{
  public:
    //! \brief The scaling function \f$x_i' = o_i+l_ix_i\f$.
    explicit ScalingFunction(const Vector<Float>& origin,
                             const Vector<Float>& lengths)
        : _o(origin), _l(lengths) { ARIADNE_ASSERT(origin.size()==lengths.size()); }
    //! \brief The scaling function which takes the unit interval \f$[-1,+1]^n\f$ into \a range.
    explicit ScalingFunction(const Vector<Interval>& range)
        : _o(midpoint(range)), _l(range.size()) { for(uint i=0; i!=_l.size(); ++i) { _l[i]=range[i].radius(); } }
    const Vector<Float>& origin() const { return _o; }
    const Vector<Float>& lengths() const { return _l; }
    ScalingFunction* clone() const { return new ScalingFunction(*this); }
    SizeType result_size() const { return _l.size(); }
    SizeType argument_size() const { return _l.size(); }
    std::ostream& write(std::ostream& os) const {
        return os << "ScalingFunction( o=" << this->origin() << ", l=" << this->lengths() << " )"; }
  private:
    friend class FunctionTemplate<ScalingFunction>;
    template<class R, class A> void _compute(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) { r[i]=_o[i]+_l[i]*x[i]; } }
    template<class R, class A> void _compute_approx(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) { r[i]=_o[i]+_l[i]*x[i]; } }
  private:
    Vector<Float> _o;
    Vector<Float> _l;
};






class FunctionElement
    : public ScalarFunctionInterface
{
  public:
    typedef unsigned int SizeType;
    typedef unsigned short SmoothnessType;

    FunctionElement(const FunctionInterface& f, SizeType i)
        : _fptr(f.clone()), _i(i) { ARIADNE_ASSERT(i<f.result_size()); }
    FunctionElement(shared_ptr<const FunctionInterface> fptr, SizeType i)
        : _fptr(fptr), _i(i) { ARIADNE_ASSERT(i<fptr->result_size()); }
    FunctionElement* clone() const { return new FunctionElement(*this); }

    SizeType argument_size() const { return _fptr->argument_size(); }
    SmoothnessType smoothness() const { return _fptr->smoothness(); }

    virtual Float evaluate(const Vector<Float>& x) const {
        return this->_evaluate(x); }
    virtual Interval evaluate(const Vector<Interval>& x) const {
        return this->_evaluate(x); }
    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const {
        return this->_evaluate(x); }
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const {
        return this->_evaluate(x); }
    virtual TaylorModel evaluate(const Vector<TaylorModel>& x) const {
        return this->_evaluate(x); }
    virtual Differential<TaylorModel> evaluate(const Vector< Differential<TaylorModel> >& x) const {
        // FIXME: Incorrect
        return x[0]; }

    virtual FunctionElement* derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual std::ostream& write(std::ostream& os) const { ARIADNE_NOT_IMPLEMENTED; }
  private:
    template<class Y> inline Y _evaluate(const Vector<Y>& x) const {
        return this->_fptr->evaluate(x)[this->_i]; }
  private:
    shared_ptr<const FunctionInterface> _fptr;
    SizeType _i;
};



inline Vector<ScalarFunctionInterface>
join(uint n, const ScalarFunctionInterface* e0, const ScalarFunctionInterface* e1, ...) {
    Vector<ScalarFunctionInterface> res(n);
    res.set(0,e0->clone());
    res.set(1,e1->clone());
    va_list args; va_start(args,e1);
    for(uint i=2; i!=n; ++i) {
        const ScalarFunctionInterface* ei=va_arg(args,const ScalarFunctionInterface*);
        res.set(i,ei->clone());
    }
    return res;
}


class ComposedFunction
    : public FunctionTemplate<ComposedFunction>
{
    friend class FunctionTemplate<ComposedFunction>;
  public:
    ComposedFunction(shared_ptr<const FunctionInterface> g, shared_ptr<const FunctionInterface> f)
        : _f(f), _g(g) { ARIADNE_ASSERT(g->argument_size()==f->result_size()); }
    ComposedFunction* clone() const { return new ComposedFunction(*this); }

    virtual SizeType result_size() const { return this->_g->result_size(); }
    virtual SizeType argument_size() const { return this->_f->argument_size(); }
    virtual SmoothnessType smoothness() const { return min(_f->smoothness(),_g->smoothness()); }
  private:
    template<class X> inline void _compute(Vector<X>& r, const Vector<X>& x) const {
        r=this->_g->evaluate(this->_f->evaluate(x)); }
    template<class X> inline void _compute_approx(Vector<X>& r, const Vector<X>& x) const {
        _compute(r,x); }
  private:
    shared_ptr<const FunctionInterface> _f;
    shared_ptr<const FunctionInterface> _g;
};

class JoinedFunction
    : public FunctionTemplate<JoinedFunction>
{
    friend class FunctionTemplate<JoinedFunction>;
  public:
    JoinedFunction(shared_ptr<const FunctionInterface> f1, shared_ptr<const FunctionInterface> f2)
        : _f1(f1), _f2(f2) { ARIADNE_ASSERT(f1->argument_size()==f2->argument_size()); }
    JoinedFunction* clone() const { return new JoinedFunction(*this); }

    virtual SizeType result_size() const { return _f1->result_size()+_f2->result_size(); }
    virtual SizeType argument_size() const { return _f1->argument_size(); }
    virtual SmoothnessType smoothness() const { return min(_f1->smoothness(),_f2->smoothness()); }
  private:
    template<class X> inline void _compute(Vector<X>& r, const Vector<X>& x) const {
        r=join(_f1->evaluate(x),_f2->evaluate(x)); }
    template<class X> inline void _compute_approx(Vector<X>& r, const Vector<X>& x) const {
        _compute(r,x); }
  private:
    shared_ptr<const FunctionInterface> _f1;
    shared_ptr<const FunctionInterface> _f2;
};

class CombinedFunction
    : public FunctionTemplate<JoinedFunction>
{
  public:
    CombinedFunction(shared_ptr<const FunctionInterface> f1, shared_ptr<const FunctionInterface> f2)
        : _f1(f1), _f2(f2) { }
    CombinedFunction* clone() const { return new CombinedFunction(*this); }

    virtual SizeType result_size() const { return _f1->result_size()+_f2->result_size(); }
    virtual SizeType argument_size() const { return _f1->argument_size()+_f2->argument_size(); }
    virtual SmoothnessType smoothness() const { return min(_f1->smoothness(),_f2->smoothness()); }
  private:
    template<class X> inline void _compute(Vector<X>& r, const Vector<X>& x) const {
        return r=combine(_f1->evaluate(project(x,range(0,_f1->argument_size()))),
                         _f2->evaluate(project(x,range(_f1->argument_size(),this->argument_size())))); }
    template<class X> inline void _compute_approx(Vector<X>& r, const Vector<X>& x) const  {
        _compute(r,x); }
  private:
    shared_ptr<const FunctionInterface> _f1;
    shared_ptr<const FunctionInterface> _f2;
};



} // namespace Ariadne

#endif
