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
#include "expression_interface.h"
#include "function_interface.h"

#include "macros.h"
#include "pointer.h"

#include "vector.h"
#include "matrix.h"
#include "polynomial.h"
#include "differential.h"
#include "taylor_model.h"

namespace Ariadne {

static const int SMOOTH=255;

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
    virtual smoothness_type smoothness() const { return SMOOTH; }

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

    virtual size_type result_size() const { return T::result_size(); }
    virtual size_type argument_size() const { return T::argument_size(); }
    virtual size_type parameter_size() const { return T::parameter_size(); }
    virtual smoothness_type smoothness() const { return T::smoothness(); }

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




template<uint AS, uint PS=0u, uint SM=255u>
struct ExpressionData
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
template<class T> class UserExpression
    : public ExpressionInterface
{
  public:
    UserExpression() : _p(this->parameter_size()) { }
    UserExpression(const Vector<Float>& p) : _p(p) { }
    UserExpression(const Vector<Interval>& p) : _p(p) { }

    const Vector<Interval> parameters() const { return _p; }

    virtual size_type argument_size() const { return T::argument_size(); }
    virtual size_type parameter_size() const { return T::parameter_size(); }
    virtual smoothness_type smoothness() const { return T::smoothness(); }

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

    virtual UserExpression<T>* derivative(uint j) { ARIADNE_ASSERT(false); }

    virtual Matrix<Float> gradient(const Vector<Float>& x) const {
        return this->evaluate(Differential<Float>::variables(1u,x)).gradient(); }
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const {
        return this->evaluate(Differential<Interval>::variables(1u,x)).gradient(); }

    virtual std::ostream& write(std::ostream& os) const  {
        return os << "UserExpression( argument_size="<<this->argument_size()<<" )"; }
  private:
    Vector<Interval> _p;
};




class ConstantExpression
    : public ExpressionInterface
{
  public:
    ConstantExpression(unsigned int as, const Float& c) : _as(as), _c(c) { }
    ConstantExpression(unsigned int as, const Interval& c) : _as(as), _c(c) { }
    ConstantExpression* clone() const { return new ConstantExpression(*this); }
    virtual size_type argument_size() const { return _as; }
    virtual smoothness_type smoothness() const { return SMOOTH; }
    virtual Float evaluate(const Vector<Float>& x) const { return midpoint(_c); }
    virtual Interval evaluate(const Vector<Interval>& x) const { return _c; }
    virtual TaylorModel evaluate(const Vector<TaylorModel>& x) const { return x[0]*0+_c; }
    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const { return x[0]*0+midpoint(_c); }
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const { return x[0]*0+midpoint(_c); }
    virtual ConstantExpression* derivative(uint j) const { return new ConstantExpression(_as,0.0); }
    virtual Vector<Float> gradient() const { return Vector<Float>::zero(_as); }
    virtual std::ostream& write(std::ostream& os) const {
        if(_c.lower()==_c.upper()) { os<<midpoint(_c); } else { os<<_c; } return os; }
  private:
    unsigned int _as;
    Interval _c;
};


//! A coordinate projection \f$\pi:\R^n\rightarrow\R\f$ given by \f$\pi(x)=x_j\f$.
class ProjectionExpression
    : public ExpressionInterface
{
  public:
    ProjectionExpression(unsigned int as, unsigned int j) : _as(as), _j(j) { }
    ProjectionExpression* clone() const { return new ProjectionExpression(*this); }
    virtual size_type argument_size() const { return _as; }
    virtual smoothness_type smoothness() const { return SMOOTH; }
    virtual Float evaluate(const Vector<Float>& x) const { return x[_j]; }
    virtual Interval evaluate(const Vector<Interval>& x) const { return x[_j]; }
    virtual TaylorModel evaluate(const Vector<TaylorModel>& x) const { return x[_j]; }
    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const { return x[_j]; }
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const { return x[_j]; }
    virtual ConstantExpression* derivative(uint j) const { return new ConstantExpression(_as,0.0); }
    virtual Vector<Float> gradient() const { return Vector<Float>::unit(_as,_j); }
    virtual std::ostream& write(std::ostream& os) const { return os << "x["<<_j<<"]"; }
  private:
    unsigned int _as;
    unsigned int _j;
};




//! An affine expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b\f$.
class AffineExpression
    : public ExpressionInterface
{
  public:
    AffineExpression(const Vector<Float>& g, const Float& c) : _c(c), _g(g) { }
    AffineExpression(const Vector<Interval>& g, const Interval& c) : _c(c), _g(g) { }
    AffineExpression(uint as, double c, double g0, ...) : _c(c), _g(as) {
        _g[0]=g0; va_list args; va_start(args,g0);
        for(uint i=1; i!=as; ++i) { _g[i]=va_arg(args,double); } }

    static AffineExpression constant(uint n, Float c) {
        return AffineExpression(Vector<Float>(n),c); }
    static AffineExpression constant(uint n, Interval c) {
        return AffineExpression(Vector<Interval>(n),c); }
    static AffineExpression variable(uint n, uint j) {
        return AffineExpression(Vector<Interval>::unit(n,j),Interval(0.0)); }

    const Vector<Interval>& a() const { return _g; }
    const Interval& b() const { return _c; }

    virtual AffineExpression* clone() const { return new AffineExpression(*this); }

    virtual size_type argument_size() const { return _g.size(); }
    virtual smoothness_type smoothness() const { return SMOOTH; }
    virtual Float evaluate(const Vector<Float>& x) const {
        Float r=midpoint(_c); for(uint j=0; j!=_g.size(); ++j) { r+=midpoint(_g[j])*x[j]; } return r; }
    virtual Interval evaluate(const Vector<Interval>& x) const {
        Interval r=_c; for(uint j=0; j!=_g.size(); ++j) { r+=_g[j]*x[j]; } return r; }
    virtual TaylorModel evaluate(const Vector<TaylorModel>& x) const {
        TaylorModel r=x[0]*0+_c; for(uint j=0; j!=_g.size(); ++j) { r+=_g[j]*x[j]; } return r; }
    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const {
        Differential<Float> r=x[0]*0+midpoint(_c); for(uint j=0; j!=_g.size(); ++j) { r+=midpoint(_g[j])*x[j]; } return r; };
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const {
        Differential<Interval> r=x[0]*0+_c; for(uint j=0; j!=_g.size(); ++j) { r+=_g[j]*x[j]; } return r; };

    virtual ConstantExpression* derivative(uint j) const { return new ConstantExpression(_g.size(),_g[j]); }
    virtual Vector<Interval> gradient() const { return this->_g; }

    virtual std::ostream& write(std::ostream& os) const { return os<<"AffineExpression( a="<<this->_g<<", b="<<this->_c<<" )"; }
  private:
    friend AffineExpression operator+(const AffineExpression&, const AffineExpression&);
    friend AffineExpression operator*(const Interval&, const AffineExpression&);
    friend Interval derivative(const AffineExpression&, uint);
    template<class R, class A>
    void compute(R& r, const A& x) const {
        r=x[0]*0+_c; for(uint j=0; j!=argument_size(); ++j) { r+=_g[j]*x[j]; } }
  private:
    Interval _c;
    Vector<Interval> _g;
};

inline AffineExpression operator+(const AffineExpression& f1, const AffineExpression& f2) {
    return AffineExpression(Vector<Interval>(f1._g+f2._g),f1._c+f2._c); }
inline AffineExpression operator*(const Interval& c, const AffineExpression& f) {
    return AffineExpression(Vector<Interval>(c*f._g),c*f._c); }
inline AffineExpression operator*(const AffineExpression& f, const Interval& c) { return c*f; }
inline Interval derivative(const AffineExpression& f, uint k) { return f._g[k]; }

//! A polynomial expression \f$p:\R^n\rightarrow\R\f$.
class PolynomialExpression
    : public ExpressionInterface,
      public Polynomial<Interval>
{
  public:
    PolynomialExpression() : Polynomial<Interval>() { }
    explicit PolynomialExpression(uint n) : Polynomial<Interval>(Polynomial<Interval>::constant(n,Interval(0))) { }
    PolynomialExpression(uint n, uint j) : Polynomial<Interval>(Polynomial<Interval>::variable(n,j)) { }
    PolynomialExpression(const Polynomial<Float>& p) : Polynomial<Interval>(p) { }
    PolynomialExpression(const Polynomial<Interval>& p) : Polynomial<Interval>(p) { }
    PolynomialExpression& operator=(const Polynomial<Interval>& p) { this->Polynomial<Interval>::operator=(p); return *this; }

    static PolynomialExpression constant(uint n, Float c) {
        return Polynomial<Interval>::constant(n,c); }
    static PolynomialExpression constant(uint n, Interval c) {
        return Polynomial<Interval>::constant(n,c); }
    static PolynomialExpression variable(uint n, uint j) {
        return Polynomial<Interval>::variable(n,j); }

    virtual PolynomialExpression* clone() const { return new PolynomialExpression(*this); }

    virtual ExpressionInterface::size_type argument_size() const {
        return static_cast<const Polynomial<Interval>&>(*this).argument_size(); }
    virtual ExpressionInterface::smoothness_type smoothness() const {
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

    virtual PolynomialExpression* derivative(uint j) const {
        return new PolynomialExpression(Ariadne::derivative(*this,j)); }

    virtual Vector<PolynomialExpression> gradient() const {
        Vector<PolynomialExpression> g(this->argument_size());
        for(uint i=0; i!=g.size(); ++i) { g[i]=Ariadne::derivative(*this,i); }
        return g; }

    virtual std::ostream& write(std::ostream& os) const { return os << static_cast<Polynomial<Interval>const&>(*this); }
};


typedef shared_ptr<const ExpressionInterface> ExpressionPointer;



template<class Op> class UnaryExpression
    : public ExpressionInterface
{
  public:
    UnaryExpression(Op op, const ExpressionInterface* expr) : _op(op), _subexpr(expr) { }
    UnaryExpression(Op op, ExpressionPointer expr) : _op(op), _subexpr(expr) { }
    virtual UnaryExpression<Op>* clone() const { return new UnaryExpression<Op>(_op,_subexpr); }
    virtual size_type argument_size() const { return _subexpr->argument_size(); }
    virtual smoothness_type smoothness() const { return 255; }
    virtual Float evaluate(const Vector<Float>& x) const { return _op(_subexpr->evaluate(x)); }
    virtual Interval evaluate(const Vector<Interval>& x) const { return _op(_subexpr->evaluate(x)); }
    virtual TaylorModel evaluate(const Vector<TaylorModel>& x) const { return _op(_subexpr->evaluate(x)); }
    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const { return _op(_subexpr->evaluate(x)); }
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const { return _op(_subexpr->evaluate(x)); }
    virtual ExpressionInterface* derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual Vector<Float> gradient() const { ARIADNE_NOT_IMPLEMENTED; }
    virtual std::ostream& write(std::ostream& os) const { return os << _op << "(" << *_subexpr << ")"; }
  private:
    template<class R, class A> inline
    void compute(R& r, const A& a) { Op _op; _op(r,_subexpr->evaluate(a)); }
  private:
    Op _op;
    shared_ptr<const ExpressionInterface> _subexpr;
};


template<class Op> class BinaryExpression
    : public ExpressionInterface
{
  public:
    BinaryExpression(Op op, const ExpressionInterface* expr1, const ExpressionInterface* expr2)
        : _op(op), _subexpr1(expr1), _subexpr2(expr2) { }
    BinaryExpression(Op op, ExpressionPointer expr1, ExpressionPointer expr2)
        : _op(op), _subexpr1(expr1), _subexpr2(expr2)  { }
    virtual BinaryExpression<Op>* clone() const { return new BinaryExpression<Op>(_op,_subexpr1,_subexpr2); }
    virtual size_type argument_size() const { return _subexpr1->argument_size(); }
    virtual smoothness_type smoothness() const { return 255; }
    virtual Float evaluate(const Vector<Float>& x) const {
        return _op(_subexpr1->evaluate(x),_subexpr2->evaluate(x)); }
    virtual Interval evaluate(const Vector<Interval>& x) const {
        return _op(_subexpr1->evaluate(x),_subexpr2->evaluate(x)); }
    virtual TaylorModel evaluate(const Vector<TaylorModel>& x) const {
        return _op(_subexpr1->evaluate(x),_subexpr2->evaluate(x)); }
    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const {
        return _op(_subexpr1->evaluate(x),_subexpr2->evaluate(x)); }
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const {
        return _op(_subexpr1->evaluate(x),_subexpr2->evaluate(x)); }
    virtual ExpressionInterface* derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual std::ostream& write(std::ostream& os) const {
        return os << *_subexpr1 << _op << *_subexpr2; }
  private:
    template<class R, class A> inline
    void compute(R& r, const A& a) { r=Op()(_subexpr1->evaluate(a),_subexpr2->evaluate(a)); }
  private:
    Op _op;
    shared_ptr<const ExpressionInterface> _subexpr1;
    shared_ptr<const ExpressionInterface> _subexpr2;
};



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
    ConstantExpression operator[](uint i) const { return ConstantExpression(_as,this->_c[i]); }
    const Vector<Interval>& c() const { return _c; }

    ConstantFunction* clone() const { return new ConstantFunction(*this); }

    size_type result_size() const { return _c.size(); }
    size_type argument_size() const { return _as; }

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

    size_type result_size() const { return this->_p.size(); }
    size_type argument_size() const { return this->_as; }
    size_type p(unsigned int j) const { return this->_p[j]; }

    ProjectionExpression operator[](unsigned int i) const { return ProjectionExpression(_as,_p[i]); }

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


//! The identity function \f$ x\mapsto x\f$ in \f$\R^n\f$.
class IdentityFunction
    : public FunctionTemplate<IdentityFunction>
{
  public:
    //! Construct the identity function in dimension \a n.
    IdentityFunction(uint n) : _n(n) { }

    IdentityFunction* clone() const { return new IdentityFunction(*this); }

    ProjectionExpression operator[](unsigned int i) const { return ProjectionExpression(_n,i); }

    size_type result_size() const { return this->_n; }
    size_type argument_size() const { return this->_n; }

    std::ostream& write(std::ostream& os) const {
        return os << "IdentityFunction( size=" << this->result_size() << " )"; }
  private:
    friend class FunctionTemplate<IdentityFunction>;
    template<class R, class A> void _compute(R& r, const A& x) const { r=x; }
    template<class R, class A> void _compute_approx(R& r, const A& x) const { r=x; }
  private:
    uint _n;
};


//! An scaling function \f$x_i' = o_i+l_ix_i\f$.
class ScalingFunction
    : public FunctionTemplate<ScalingFunction>
{
  public:
    explicit ScalingFunction(const Vector<Float>& origin,
                             const Vector<Float>& lengths)
        : _o(origin), _l(lengths) { ARIADNE_ASSERT(origin.size()==lengths.size()); }
    explicit ScalingFunction(const Vector<Interval>& range)
        : _o(midpoint(range)), _l(range.size()) { for(uint i=0; i!=_l.size(); ++i) { _l[i]=range[i].radius(); } }
    const Vector<Float>& origin() const { return _o; }
    const Vector<Float>& lengths() const { return _l; }
    ScalingFunction* clone() const { return new ScalingFunction(*this); }
    size_type result_size() const { return _l.size(); }
    size_type argument_size() const { return _l.size(); }
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



//! An affine function \f$x\mapsto Ax+b\f$.
class AffineFunction
    : public FunctionInterface
{
  public:

    //! Construct an affine function from the matrix \a A and vector \a b.
    AffineFunction(const Matrix<Interval>& A, const Vector<Interval>& b)
        : _fA(midpoint(A)), _fb(midpoint(b)), _iA(A), _ib(b) { ARIADNE_ASSERT(A.row_size()==b.size()); }

    AffineFunction(uint rs, double b0, ...) : _fA(), _fb(rs) {
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

    virtual AffineFunction* clone() const { return new AffineFunction(*this); }

    virtual size_type result_size() const { return _fA.row_size(); }
    virtual size_type argument_size() const { return _fA.column_size(); }
    virtual smoothness_type smoothness() const { return SMOOTH; }

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

    virtual AffineExpression operator[](unsigned int i) const  {
        ARIADNE_ASSERT(i<this->result_size());
        Vector<Interval> g(this->argument_size()); for(unsigned int j=0; j!=g.size(); ++j) { g[j]=this->_iA[i][j]; }
        return AffineExpression(g,_ib[i]); }

    virtual std::ostream& write(std::ostream& os) const {
        return os << "AffineFunction( A="<<midpoint(_iA)<<", b="<<midpoint(_ib)<<" )"; }
  private:
    Matrix<Float> _fA; Vector<Float> _fb;
    Matrix<Interval> _iA; Vector<Interval> _ib;
};

//! A polynomial function.
class PolynomialFunction
    : public FunctionInterface
{
  public:
    explicit PolynomialFunction(uint rs, uint as) : _p(rs,PolynomialExpression(as)) { }
    explicit PolynomialFunction(const PolynomialExpression& p) : _p(1,p) { }

    PolynomialFunction(const Vector< Polynomial<Float> >& p) : _p(p) { }
    PolynomialFunction(const Vector< Polynomial<Interval> >& p) : _p(p) { }

    template<class E> PolynomialFunction(const ublas::vector_expression<E>& p)
        : _p(static_cast< Vector< Polynomial<Interval> > >(p)) { }

    virtual PolynomialFunction* clone() const { return new PolynomialFunction(*this); }

    virtual size_type result_size() const { return _p.size(); }
    virtual size_type argument_size() const { return _p[0].argument_size(); }
    virtual smoothness_type smoothness() const { return SMOOTH; }

    virtual Vector<Float> evaluate(const Vector<Float>& x) const {
        return Ariadne::evaluate(midpoint(_p),x); }
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const {
        return Ariadne::evaluate(_p,x); }
    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const {
        return Ariadne::evaluate(midpoint(_p),x); }
    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const {
        return Ariadne::evaluate(_p,x); }
    virtual Vector<TaylorModel> evaluate(const Vector<TaylorModel>& x) const {
        return Ariadne::evaluate(_p,x); }

    virtual Matrix<Float> jacobian(const Vector<Float>& x) const {
        return Ariadne::jacobian(this->evaluate(Differential<Float>::variables(1u,x))); }
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const {
        return Ariadne::jacobian(this->evaluate(Differential<Interval>::variables(1u,x))); }

    virtual PolynomialExpression operator[](unsigned int i) { return _p[i]; }

    virtual std::ostream& write(std::ostream& os) const {
        return os << "PolynomialFunction"<<_p; }
  private:
    Vector< Polynomial<Interval> > _p;
};

inline PolynomialFunction join(const PolynomialExpression& p1, const PolynomialExpression& p2) {
    return join(static_cast<const Polynomial<Interval>&>(p1),static_cast<const Polynomial<Interval>&>(p2));
}

class FunctionElement
    : public ExpressionInterface
{
  public:
    typedef unsigned int size_type;
    typedef unsigned short smoothness_type;

    FunctionElement(const FunctionInterface& f, size_type i)
        : _fptr(f.clone()), _i(i) { ARIADNE_ASSERT(i<f.result_size()); }
    FunctionElement(shared_ptr<const FunctionInterface> fptr, size_type i)
        : _fptr(fptr), _i(i) { ARIADNE_ASSERT(i<fptr->result_size()); }
    FunctionElement* clone() const { return new FunctionElement(*this); }

    size_type argument_size() const { return _fptr->argument_size(); }
    smoothness_type smoothness() const { return _fptr->smoothness(); }

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

    virtual FunctionElement* derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }
    virtual std::ostream& write(std::ostream& os) const { ARIADNE_NOT_IMPLEMENTED; }
  private:
    template<class Y> inline Y _evaluate(const Vector<Y>& x) const {
        return this->_fptr->evaluate(x)[this->_i]; }
  private:
    shared_ptr<const FunctionInterface> _fptr;
    size_type _i;
};


template<>
class Vector<ExpressionInterface>
    : public FunctionTemplate< Vector<ExpressionInterface> >
{
    friend class FunctionTemplate< Vector<ExpressionInterface> >;
  public:
    typedef unsigned int size_type;
    typedef unsigned short int smoothness_type;

    Vector(size_type n) : _expressions(n) { }
    size_type size() const { return _expressions.size(); }
    const ExpressionInterface& operator[](size_type i) const { return *_expressions[i]; }
    void set(size_type i, shared_ptr<const ExpressionInterface> p) {
        ARIADNE_ASSERT(this->size()==0 || p->argument_size()==this->argument_size()); _expressions[i]=p; }
    void set(size_type i, const ExpressionInterface* p) {
        ARIADNE_ASSERT(this->size()==0 || p->argument_size()==this->argument_size());
        _expressions[i]=shared_ptr<const ExpressionInterface>(p); }

    virtual Vector<ExpressionInterface>* clone() const { return new Vector<ExpressionInterface>(*this); }

    virtual size_type result_size() const { return _expressions.size(); }
    virtual size_type argument_size() const { if(_expressions[0]) { return _expressions[0]->argument_size(); } else { return 0; } }
    virtual smoothness_type smoothness() const {
        smoothness_type res=_expressions[0]->smoothness();
        for(uint i=1; i!=this->size(); ++i) {
            res=max(res,_expressions[i]->smoothness()); } return res; }

    virtual std::ostream& write(std::ostream& os) const {
        for(uint i=0; i!=this->result_size(); ++i) { os<<(*this)[i]; } return os; }
  private:
    template<class Y> inline void _compute(Vector<Y>& r, const Vector<Y>& x) const {
        for(size_type i=0; i!=this->result_size(); ++i) {
            r[i]=this->_expressions[i]->evaluate(x); } }
    template<class Y> inline void _compute_approx(Vector<Y>& r, const Vector<Y>& x) const {
        _compute(r,x); }
  private:
    array< boost::shared_ptr<const ExpressionInterface> > _expressions;
};


inline Vector<ExpressionInterface>
join(uint n, const ExpressionInterface* e0, const ExpressionInterface* e1, ...) {
    Vector<ExpressionInterface> res(n);
    res.set(0,e0->clone());
    res.set(1,e1->clone());
    va_list args; va_start(args,e1);
    for(uint i=2; i!=n; ++i) {
        const ExpressionInterface* ei=va_arg(args,const ExpressionInterface*);
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

    virtual size_type result_size() const { return this->_g->result_size(); }
    virtual size_type argument_size() const { return this->_f->argument_size(); }
    virtual smoothness_type smoothness() const { return min(_f->smoothness(),_g->smoothness()); }
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

    virtual size_type result_size() const { return _f1->result_size()+_f2->result_size(); }
    virtual size_type argument_size() const { return _f1->argument_size(); }
    virtual smoothness_type smoothness() const { return min(_f1->smoothness(),_f2->smoothness()); }
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

    virtual size_type result_size() const { return _f1->result_size()+_f2->result_size(); }
    virtual size_type argument_size() const { return _f1->argument_size()+_f2->argument_size(); }
    virtual smoothness_type smoothness() const { return min(_f1->smoothness(),_f2->smoothness()); }
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
