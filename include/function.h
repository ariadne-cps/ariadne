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
#include "container.h"

#include "vector.h"
#include "matrix.h"
#include "polynomial.h"
#include "affine.h"
#include "taylor_model.h"
#include "differential.h"
#include "operators.h"
#include "expression.h"
#include "formula.h"

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
template<class T> class ScalarUserFunction
    : public ScalarFunctionInterface
{
  public:
    ScalarUserFunction() : _p(this->parameter_size()) { }
    ScalarUserFunction(const Vector<Float>& p) : _p(p) { }
    ScalarUserFunction(const Vector<Interval>& p) : _p(p) { }

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

    virtual ScalarUserFunction<T>* derivative(uint j) { ARIADNE_NOT_IMPLEMENTED; }

    virtual Vector<Float> gradient(const Vector<Float>& x) const {
        return this->evaluate(Differential<Float>::variables(1u,x)).gradient(); }
    virtual Vector<Interval> gradient(const Vector<Interval>& x) const {
        return this->evaluate(Differential<Interval>::variables(1u,x)).gradient(); }

    virtual std::ostream& write(std::ostream& os) const  {
        return os << "ScalarUserFunction( argument_size="<<this->argument_size()<<" )"; }
  private:
    Vector<Interval> _p;
};

template<class X> X evaluate(const Expression<Real>& e, const Vector<X>& x);

class ScalarExpressionFunction
//    : public ScalarFunctionInterface
{
    Array<String> _space;
    Expression<Real> _expression;
  public:
    typedef unsigned int SizeType;

    ScalarExpressionFunction(const Expression<Real>& e, const Array<String>& s)
        : _space(s.begin(),s.end()), _expression(_convert(e,_variables(s))) { }

    virtual uint argument_size() const { return _space.size(); }
    virtual Float evaluate(const Vector<Float>& x) { return this->_evaluate(x); }
    virtual Interval evaluate(const Vector<Interval>& x) { return this->_evaluate(x); }
    virtual TaylorModel evaluate(const Vector<TaylorModel>& x) { return this->_evaluate(x); }
    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) { return this->_evaluate(x); }
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) { return this->_evaluate(x); }
    virtual std::ostream& write(std::ostream& os) const {
        return os << "Function( space="<<this->_space<<", expression="<<this->_expression<<" )"; }

  protected:
    template<class X> X _evaluate(const Vector<X>& x) { return Ariadne::evaluate(_expression,x); }

    static Map<String,SizeType> _variables(const Array<String>& s) {
        Map<String,SizeType> variables;
        for(uint i=0; i!=s.size(); ++i) {
            ARIADNE_ASSERT_MSG(!variables.has_key(s[i]),"Variable name "<<s[i]<<" is duplicated in space "<<s);
            variables.insert(s[i],i); }
        return variables; }

    static Expression<Real> _convert(const Expression<Real> e, const Map<String,SizeType>& s);

/*
    static Expression<Real> _convert(shared_ptr<const ExpressionInterface<Real> > e, const Map<String,SizeType>& s) {
        const ExpressionInterface<Real>* eptr=e.operator->();
        const BinaryExpression<Real,Operator>* bptr=dynamic_cast<const BinaryExpression<Real,Operator>*>(eptr);
        if(bptr) { return make_expression<Real>(bptr->_op,_convert(bptr->_arg1_ptr,s),_convert(bptr->_arg2_ptr,s)); }
        const UnaryExpression<Real,Operator>* uptr=dynamic_cast<const UnaryExpression<Real,Operator>*>(eptr);
        if(uptr) { return make_expression<Real>(uptr->_op,_convert(uptr->_arg_ptr,s)); }
        const ConstantExpression<Real>* cptr=dynamic_cast<const ConstantExpression<Real>*>(eptr);
        if(cptr) { return Expression<Real>(*cptr); }
        const VariableExpression<Real>* vptr=dynamic_cast<const VariableExpression<Real>*>(eptr);
        if(vptr) { return Expression<Real>(new CoordinateExpression<Real>(s[vptr->name()])); }
        const CoordinateExpression<Real>* iptr=dynamic_cast<const CoordinateExpression<Real>*>(eptr);
        ARIADNE_ASSERT_MSG(!iptr,"Cannot convert numbered variable");
    }
*/

};



//! An affine expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b\f$.
class ScalarConstantFunction
    : public ScalarFunctionInterface
{
  public:
    operator Interval() const;
};


//! An affine expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b\f$.
class ScalarAffineFunction
    : public ScalarFunctionInterface,
      public Affine<Interval>
{
  public:
    explicit ScalarAffineFunction(uint n) : Affine<Interval>(n) { }
    ScalarAffineFunction(const Affine<Interval>& aff) : Affine<Interval>(aff) { }
    ScalarAffineFunction(const Vector<Float>& g, const Float& c) : Affine<Interval>(g,c) { }
    ScalarAffineFunction(const Vector<Interval>& g, const Interval& c) : Affine<Interval>(g,c) { }
    ScalarAffineFunction(uint as, double c, double g0, ...) {
        Float _c; Vector<Float> _g(as);
        _g[0]=g0; va_list args; va_start(args,g0);
        for(uint i=1; i!=as; ++i) { _g[i]=va_arg(args,double); }
        *this=ScalarAffineFunction(_g,_c); }
//    ScalarAffineFunction(const ConstantScalarFunction& c)
//        : _c(c), _g(c.argument_size()) { }

    ScalarAffineFunction(Expression<Real> e, const Array<String>& s);
//    ScalarAffineFunction(Expression<Real> e, const Map<String,SizeType>& s);

    static ScalarAffineFunction constant(uint n, Float c) {
        return ScalarAffineFunction(Vector<Float>(n),c); }
    static ScalarAffineFunction constant(uint n, Interval c) {
        return ScalarAffineFunction(Vector<Interval>(n),c); }
    static ScalarAffineFunction variable(uint n, uint j) {
        return ScalarAffineFunction(Vector<Interval>::unit(n,j),Interval(0.0)); }

    virtual ScalarAffineFunction* clone() const { return new ScalarAffineFunction(*this); }

    virtual SizeType argument_size() const { return static_cast<Affine<Interval>const&>(*this).argument_size(); }
    virtual SmoothnessType smoothness() const { return SMOOTH; }
    virtual Float evaluate(const Vector<Float>& x) const { return this->_evaluate_approx<Float>(x); }
    virtual Interval evaluate(const Vector<Interval>& x) const { return this->_evaluate<Interval>(x); }
    virtual TaylorModel evaluate(const Vector<TaylorModel>& x) const { return this->_evaluate<TaylorModel>(x); }
    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const {
        return this->_evaluate_approx< Differential<Float> >(x); }
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const {
        return this->_evaluate< Differential<Interval> >(x); }
    virtual Differential<TaylorModel> evaluate(const Vector< Differential<TaylorModel> >& x) const {
        return this->_evaluate< Differential<TaylorModel> >(x); }

    virtual Vector<Float> gradient(const Vector<Float>& x) const {
        return this->evaluate(Differential<Float>::variables(1u,x)).gradient(); }
    virtual Vector<Interval> gradient(const Vector<Interval>& x) const {
        return this->evaluate(Differential<Interval>::variables(1u,x)).gradient(); }

//    virtual ConstantScalarFunction* derivative(uint j) const { return new ConstantScalarFunction(_g.size(),_g[j]); }
    virtual ScalarFunctionInterface* derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }


    virtual std::ostream& write(std::ostream& os) const { return os << static_cast<Affine<Interval>const&>(*this); }
  private:
    template<class R, class A>
    R _evaluate_approx(const A& x) const {
        R r=x[0]*0+midpoint(b()); for(uint j=0; j!=argument_size(); ++j) { r+=midpoint(a()[j])*x[j]; } return r; }
    template<class R, class A>
    R _evaluate(const A& x) const {
        R r=x[0]*0+b(); for(uint j=0; j!=argument_size(); ++j) { r+=a()[j]*x[j]; } return r; }

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

    ScalarPolynomialFunction(const Expression<Real>& e, const Space<Real>& s) { *this=ScalarPolynomialFunction(polynomial<Interval>(e,s)); }

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

    virtual Vector<Float> gradient(const Vector<Float>& x) const {
        return this->evaluate(Differential<Float>::variables(1u,x)).gradient(); }
    virtual Vector<Interval> gradient(const Vector<Interval>& x) const {
        return this->evaluate(Differential<Interval>::variables(1u,x)).gradient(); }

    virtual ScalarPolynomialFunction* derivative(uint j) const {
        return new ScalarPolynomialFunction(Ariadne::derivative(*this,j)); }
    virtual ScalarPolynomialFunction* antiderivative(uint j) const {
        return new ScalarPolynomialFunction(Ariadne::antiderivative(*this,j)); }

    virtual Vector<ScalarPolynomialFunction> derivative() const {
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


// A wrapper for classes with non-static _compute and _compute_approx methods
template<class F>
class VectorFunctionTemplate
    : public VectorFunctionInterface
{
  private:
    template<class R, class A> void _base_compute(R& r, const A& a) const {
        static_cast<const F*>(this)->_compute(r,a); }
    template<class R, class A> void _base_compute_approx(R& r, const A& a) const {
        static_cast<const F*>(this)->_compute_approx(r,a); }
  protected:
    VectorFunctionTemplate() { }
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
class VectorFunctionData
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
//! <tt>ScalarFunctionData<RS,AS,PS,SM=SMOOTH></tt> class.
//!
//! The constant \a SMOOTH is used for an arbitrarily-differentiable function.
template<class T> class VectorUserFunction
    : public VectorFunctionInterface
{
  public:
    VectorUserFunction() : _p(this->parameter_size()) { }
    VectorUserFunction(const Vector<Float>& p) : _p(p) { }
    VectorUserFunction(const Vector<Interval>& p) : _p(p) { }

    const Vector<Interval> parameters() const { return _p; }

    virtual VectorUserFunction<T>* clone() const { return new VectorUserFunction<T>(*this); }

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
        return os << "VectorUserFunction( result_size="<<this->result_size()<<", argument_size="<<this->argument_size()<<" )"; }
  private:
    Vector<Interval> _p;
};




//! \brief A vector-valued function defined by its scalar-valued components.
template<>
class Vector<ScalarFunctionInterface>
    : public VectorFunctionInterface
{
    friend class VectorFunctionTemplate< Vector<ScalarFunctionInterface> >;
  public:
    typedef unsigned int SizeType;
    typedef unsigned short int SmoothnessType;

    Vector(SizeType n) : _expressions(n) { }
    SizeType size() const { return _expressions.size(); }
    const ScalarFunctionInterface& operator[](SizeType i) const { return *_expressions[i]; }
    void set(SizeType i, const ScalarFunctionInterface& f) {
        ARIADNE_ASSERT(i<this->result_size());
        ARIADNE_ASSERT(this->size()==0 || i==0 || f.argument_size()==this->argument_size());
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



//! A constant function \f$ x\mapsto c\f$ from \f$R^m\f$ to \f$R^n\f$.
class VectorConstantFunction
    : public VectorFunctionTemplate<VectorConstantFunction>
{
  public:
    VectorConstantFunction(uint rs, double c0, ...) : _as(), _c(rs) {
        ARIADNE_ASSERT(rs>0); va_list args; va_start(args,c0);
        _c[0]=c0; for(size_t i=1; i!=rs; ++i) { _c[i]=va_arg(args,double); } _as=va_arg(args,int);
        va_end(args);
    }

    VectorConstantFunction(uint rs, Interval c0, ...) : _as(), _c(rs) {
        ARIADNE_ASSERT(rs>0); double l,u; va_list args; va_start(args,c0);
        _c[0]=c0; for(size_t i=1; i!=rs; ++i) { l=va_arg(args,double); u=va_arg(args,double); _c[i]=Interval(l,u); } _as=va_arg(args,int);
        va_end(args);
    }

    VectorConstantFunction(const Vector<Float>& c, uint as) : _as(as), _c(c) { }
    VectorConstantFunction(const Vector<Interval>& c, uint as) : _as(as), _c(c) { }
//    ConstantExpression<Real,Real> operator[](uint i) const { return ConstantExpression<Real,Real>(_as,this->_c[i]); }
    Interval operator[](uint i) const { return this->_c[i]; }
    const Vector<Interval>& c() const { return _c; }

    VectorConstantFunction* clone() const { return new VectorConstantFunction(*this); }

    SizeType result_size() const { return _c.size(); }
    SizeType argument_size() const { return _as; }

    std::ostream& write(std::ostream& os) const {
        return os << "VectorConstantFunction( argument_size=" << this->argument_size()
                  << ", c=" << this->c() << " )"; }
  private:
    friend class VectorFunctionTemplate<VectorConstantFunction>;
    template<class R, class A> void _compute(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) { r[i]=_c[i]; } }
    template<class R, class A> void _compute_approx(R& r, const A& x) const {
        for(uint i=0; i!=result_size(); ++i) { r[i]=midpoint(_c[i]); } }
  private:
    uint _as;
    Vector<Interval> _c;
};


//! An affine function \f$x\mapsto Ax+b\f$.
class VectorAffineFunction
    : public VectorFunctionInterface
{
  public:

    //! Construct an affine function from the matrix \a A and vector \a b.
    VectorAffineFunction(const Matrix<Interval>& A, const Vector<Interval>& b)
        : _fA(midpoint(A)), _fb(midpoint(b)), _iA(A), _ib(b) { ARIADNE_ASSERT(A.row_size()==b.size()); }

    VectorAffineFunction(uint rs, double b0, ...) : _fA(), _fb(rs) {
        ARIADNE_ASSERT(rs>0);
        va_list args; va_start(args,b0);
        _fb[0]=b0; for(size_t i=1; i!=rs; ++i) { _fb[i]=va_arg(args,double); }
        uint as=va_arg(args,int); ARIADNE_ASSERT(as>0); _fA=Matrix<Float>(rs,as);
        for(size_t i=0; i!=rs; ++i) { for(size_t j=0; i!=as; ++j) { _fA[i][j]=va_arg(args,double); } }
        va_end(args);
        _iA=Matrix<Interval>(_fA); _ib=Vector<Interval>(_fb);
    }

    const Matrix<Float>& A() const { return _fA; }
    const Vector<Float>& b() const { return _fb; }

    virtual VectorAffineFunction* clone() const { return new VectorAffineFunction(*this); }

    virtual SizeType result_size() const { return _fA.row_size(); }
    virtual SizeType argument_size() const { return _fA.column_size(); }
    virtual SmoothnessType smoothness() const { return SMOOTH; }

    virtual Vector<Float> evaluate(const Vector<Float>& x) const {
        return prod(_fA,x)+_fb; }
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const {
        return prod(_iA,x)+_ib; }
    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const {
        return prod(_fA,x)+_fb; }
    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const {
        return prod(_iA,x)+_ib; }
    virtual Vector<TaylorModel> evaluate(const Vector<TaylorModel>& x) const {
        return prod(_iA,x)+_ib;
        ARIADNE_ASSERT(x.size()==this->argument_size());
        for(uint i=1; i<x.size(); ++i) { ARIADNE_ASSERT(x[i].argument_size()==x[0].argument_size()); }
        Vector<TaylorModel> r(this->result_size(),x[0]*0.0);
        for(uint i=0; i!=r.size(); ++i) { r[i]=_ib[i]; for(uint j=0; j!=x.size(); ++j) { r[i]+=_iA[i][j]*x[j]; } }
        return r;
    }
    virtual Matrix<Float> jacobian(const Vector<Float>& x) const {
        return this->_fA; }
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const {
        return this->_iA; }

    virtual ScalarAffineFunction operator[](unsigned int i) const  {
        ARIADNE_ASSERT(i<this->result_size());
        Vector<Interval> g(this->argument_size()); for(unsigned int j=0; j!=g.size(); ++j) { g[j]=this->_iA[i][j]; }
        return ScalarAffineFunction(g,_ib[i]); }

    virtual std::ostream& write(std::ostream& os) const;
  private:
    Matrix<Float> _fA; Vector<Float> _fb;
    Matrix<Interval> _iA; Vector<Interval> _ib;
};

inline std::ostream& VectorAffineFunction::write(std::ostream& os) const {
    //return os << "VectorAffineFunction( A="<<midpoint(_iA)<<", b="<<midpoint(_ib)<<" )";
    const VectorAffineFunction& f=*this;
    for(uint i=0; i!=f.result_size(); ++i) {
        os<<(i==0?"VectorAffineFunction( ":"; ");
        if(f.b()[i]!=0) { os<<f.b()[i]; }
        for(uint j=0; j!=f.argument_size(); ++j) {
            if(f.A()[i][j]!=0) {
                if(f.A()[i][j]>0) { os<<"+"; } else { os<<"-"; }
                if(abs(f.A()[i][j])!=1) { os<<abs(f.A()[i][j]); }
                //ss<<char('x'+j);
                os<<"x"<<j;
            }
        }
    }
    os<<" )";
    return os;
}



//! A polynomial function.
class VectorPolynomialFunction
    : public Vector<ScalarFunctionInterface>
{
    typedef Vector<ScalarFunctionInterface> Base;
  public:
    explicit VectorPolynomialFunction(uint rs, uint as) : Base(rs) {
        for(uint i=0; i!=rs; ++i) { Base::set(i,ScalarPolynomialFunction(as)); } }
    explicit VectorPolynomialFunction(uint rs, const ScalarPolynomialFunction& p) : Base(rs) {
        for(uint i=0; i!=rs; ++i) { Base::set(i,p); } }

    template<class E> VectorPolynomialFunction(const ublas::vector_expression<E>& p) : Base(p().size()) {
        for(uint i=0; i!=p().size(); ++i) { Base::set(i,ScalarPolynomialFunction(p()[i])); } }

    VectorPolynomialFunction(const Vector< Polynomial<Float> >& p) : Base(p.size()) {
        for(uint i=0; i!=p.size(); ++i) { Base::set(i,ScalarPolynomialFunction(p[i])); } }
    VectorPolynomialFunction(const Vector< Polynomial<Interval> >& p) : Base(p.size()) {
        for(uint i=0; i!=p.size(); ++i) { Base::set(i,ScalarPolynomialFunction(p[i])); } }
    VectorPolynomialFunction(const Vector<ScalarPolynomialFunction>& p) : Base(p.size()) {
        for(uint i=0; i!=p.size(); ++i) { Base::set(i,p[i]); } }
    VectorPolynomialFunction(const array< Polynomial<Interval> >& p) : Base(p.size()) {
        for(uint i=0; i!=p.size(); ++i) { Base::set(i,ScalarPolynomialFunction(p[i])); } }

    static VectorPolynomialFunction identity(uint n) {
        return VectorPolynomialFunction(ScalarPolynomialFunction::variables(n)); }

    const ScalarPolynomialFunction& operator[](uint i) const {
        return dynamic_cast<const ScalarPolynomialFunction&>(Base::get(i)); }

    const ScalarPolynomialFunction& get(uint i) const {
        return dynamic_cast<const ScalarPolynomialFunction&>(Base::get(i)); }
    void set(uint i,const ScalarPolynomialFunction& p) {
        return this->Base::set(i,p); }

    virtual std::ostream& write(std::ostream& os) const {
        return os << "VectorPolynomialFunction("<<static_cast<const Vector<ScalarFunctionInterface>&>(*this)<<")"; }
  private:
    friend VectorPolynomialFunction operator+(const VectorPolynomialFunction&, const VectorPolynomialFunction&);
    friend VectorPolynomialFunction flip(const VectorPolynomialFunction&, uint);
};

inline VectorPolynomialFunction join(const ScalarPolynomialFunction& p1, const ScalarPolynomialFunction& p2) {
    return join(static_cast<const Polynomial<Interval>&>(p1),static_cast<const Polynomial<Interval>&>(p2));
}

inline VectorPolynomialFunction join(const VectorPolynomialFunction& p1, const ScalarPolynomialFunction& p2) {
    // Need to implement this explicitly since VectorPolynomialFunction does not inherit from Vector<Polynomial>
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    VectorPolynomialFunction r(p1.result_size()+1u,p1.argument_size());
    for(uint i=0; i!=p1.result_size(); ++i) {
        r.set(i,p1[i]); }
    r.set(p1.result_size(),p2);
    return r;
}

inline VectorPolynomialFunction join(const ScalarPolynomialFunction& p1, const VectorPolynomialFunction& p2) {
    // Need to implement this explicitly since VectorPolynomialFunction does not inherit from Vector<Polynomial>
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    VectorPolynomialFunction r(p2.result_size()+1u,p1.argument_size());
    r.set(0u,p1);
    for(uint i=0; i!=p2.result_size(); ++i) {
        r.set(i+1,p2[i]); }
    return r;
}

inline VectorPolynomialFunction join(const VectorPolynomialFunction& p1, const VectorPolynomialFunction& p2) {
    // Need to implement this explicitly since VectorPolynomialFunction does not inherit from Vector<Polynomial>
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    VectorPolynomialFunction r(p1.result_size()+p2.result_size(),p1.argument_size());
    for(uint i=0; i!=p1.result_size(); ++i) {
        r.set(i,p1[i]); }
    for(uint i=0; i!=p2.result_size(); ++i) {
        r.set(p1.result_size()+i,p2[i]); }
    return r;
}

inline VectorPolynomialFunction operator*(const ScalarPolynomialFunction& p, const Vector<Float>& e) {
    VectorPolynomialFunction r(e.size(),p.argument_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,static_cast<ScalarPolynomialFunction>(e[i]*p)); }
    return r;
}

inline VectorPolynomialFunction operator+(const VectorPolynomialFunction& p1, const VectorPolynomialFunction& p2) {
    ARIADNE_ASSERT(p1.result_size()==p2.result_size());
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    VectorPolynomialFunction r(p1.result_size(),p1.argument_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,static_cast<ScalarPolynomialFunction>(p1[i]+p2[i])); }
    return r;
}

inline VectorPolynomialFunction operator-(const VectorPolynomialFunction& p1, const VectorPolynomialFunction& p2) {
    ARIADNE_ASSERT(p1.result_size()==p2.result_size());
    ARIADNE_ASSERT(p1.argument_size()==p2.argument_size());
    VectorPolynomialFunction r(p1.result_size(),p1.argument_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,static_cast<ScalarPolynomialFunction>(p1[i]-p2[i])); }
    return r;
}








//! \brief A projection function \f$ x'_i= x_{p(i)}\f$.
class ProjectionFunction
    : public VectorFunctionTemplate<ProjectionFunction>
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

  private:
    friend class VectorFunctionTemplate<ProjectionFunction>;
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
    : public VectorFunctionTemplate<IdentityFunction>
{
  public:
    //! \brief Construct the identity function in dimension \a n.
    IdentityFunction(uint n) : _n(n) { }

    IdentityFunction* clone() const { return new IdentityFunction(*this); }

    //CoordinateExpression<Real> operator[](unsigned int i) const { return CoordinateExpression<Real>(_n,i); }

    SizeType result_size() const { return this->_n; }
    SizeType argument_size() const { return this->_n; }

    std::ostream& write(std::ostream& os) const {
        return os << "VectorIdentityFunction( size=" << this->result_size() << " )"; }
  private:
    friend class VectorFunctionTemplate<IdentityFunction>;
    template<class R, class A> void _compute(R& r, const A& x) const { r=x; }
    template<class R, class A> void _compute_approx(R& r, const A& x) const { r=x; }
  private:
    uint _n;
};


//! \brief A scaling function \f$x_i' = o_i+l_ix_i\f$.
class VectorScalingFunction
    : public VectorFunctionTemplate<VectorScalingFunction>
{
  public:
    //! \brief The scaling function \f$x_i' = o_i+l_ix_i\f$.
    explicit VectorScalingFunction(const Vector<Float>& origin,
                             const Vector<Float>& lengths)
        : _o(origin), _l(lengths) { ARIADNE_ASSERT(origin.size()==lengths.size()); }
    //! \brief The scaling function which takes the unit interval \f$[-1,+1]^n\f$ into \a range.
    explicit VectorScalingFunction(const Vector<Interval>& range)
        : _o(midpoint(range)), _l(range.size()) { for(uint i=0; i!=_l.size(); ++i) { _l[i]=range[i].radius(); } }
    const Vector<Float>& origin() const { return _o; }
    const Vector<Float>& lengths() const { return _l; }
    VectorScalingFunction* clone() const { return new VectorScalingFunction(*this); }
    SizeType result_size() const { return _l.size(); }
    SizeType argument_size() const { return _l.size(); }
    std::ostream& write(std::ostream& os) const {
        return os << "VectorScalingFunction( o=" << this->origin() << ", l=" << this->lengths() << " )"; }
  private:
    friend class VectorFunctionTemplate<VectorScalingFunction>;
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

    FunctionElement(const VectorFunctionInterface& f, SizeType i)
        : _fptr(f.clone()), _i(i) { ARIADNE_ASSERT(i<f.result_size()); }
    FunctionElement(shared_ptr<const VectorFunctionInterface> fptr, SizeType i)
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

    virtual Vector<Float> gradient(const Vector<Float>& x) const {
        return this->evaluate(Differential<Float>::variables(1u,x)).gradient(); }
    virtual Vector<Interval> gradient(const Vector<Interval>& x) const {
        return this->evaluate(Differential<Interval>::variables(1u,x)).gradient(); }

    virtual FunctionElement* derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual std::ostream& write(std::ostream& os) const { ARIADNE_NOT_IMPLEMENTED; }
  private:
    template<class Y> inline Y _evaluate(const Vector<Y>& x) const {
        return this->_fptr->evaluate(x)[this->_i]; }
  private:
    shared_ptr<const VectorFunctionInterface> _fptr;
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
    : public VectorFunctionTemplate<ComposedFunction>
{
    friend class VectorFunctionTemplate<ComposedFunction>;
  public:
    ComposedFunction(shared_ptr<const VectorFunctionInterface> g, shared_ptr<const VectorFunctionInterface> f)
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
    shared_ptr<const VectorFunctionInterface> _f;
    shared_ptr<const VectorFunctionInterface> _g;
};

class JoinedFunction
    : public VectorFunctionTemplate<JoinedFunction>
{
    friend class VectorFunctionTemplate<JoinedFunction>;
  public:
    JoinedFunction(shared_ptr<const VectorFunctionInterface> f1, shared_ptr<const VectorFunctionInterface> f2)
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
    shared_ptr<const VectorFunctionInterface> _f1;
    shared_ptr<const VectorFunctionInterface> _f2;
};

class CombinedFunction
    : public VectorFunctionTemplate<JoinedFunction>
{
  public:
    CombinedFunction(shared_ptr<const VectorFunctionInterface> f1, shared_ptr<const VectorFunctionInterface> f2)
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
    shared_ptr<const VectorFunctionInterface> _f1;
    shared_ptr<const VectorFunctionInterface> _f2;
};



} // namespace Ariadne

#endif
