/***************************************************************************
 *            function.h
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

/*! \file function.h
 *  \brief Concrete function types.
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
#include "transformation.h"
#include "differential.h"
#include "taylor_model.h"

namespace Ariadne {

// A wrapper for transformations
// This class is for internal use only; we can easily specify parameter types.
template<class T>
class FunctionBase
    : public FunctionInterface,
      public T
{
  protected:
    FunctionBase() : T() { }
    template<class S> FunctionBase(const S& s) : T(s) { }
    template<class S1, class S2> FunctionBase(const S1& s1, const S2& s2) : T(s1,s2) { }
    template<class S1, class S2, class S3> FunctionBase(const S1& s1, const S2& s2, const S3& s3) : T(s1,s2,s3) { }
  public:
    virtual FunctionBase<T>* clone() const { return new FunctionBase<T>(*this); }
    virtual uint result_size() const { return this->T::result_size(); }
    virtual uint argument_size() const { return this->T::argument_size(); }
    virtual ushort smoothness() const { return this->T::smoothness(); }
    virtual Vector<Float> evaluate(const Vector<Float>& x) const {
        Vector<Float> r(this->T::result_size()); this->T::compute(r,x); return r; }
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const {
        Vector<Interval> r(this->T::result_size()); this->T::compute(r,x); return r; }
    virtual Vector<TaylorModel> evaluate(const Vector<TaylorModel>& x) const {
        Vector<TaylorModel> r(this->result_size(),TaylorModel(x[0].argument_size()));
        this->T::compute(r,x); return r; }
    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const {
        Vector< Differential<Float> > r(this->T::result_size()); this->T::compute(r,x); return r; }
    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const {
        Vector< Differential<Interval> > r(this->T::result_size()); this->T::compute(r,x); return r; }

    virtual Matrix<Float> jacobian(const Vector<Float>& x) const {
        return Ariadne::jacobian(this->_expansion(x,1u)); }
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const {
        return Ariadne::jacobian(this->_expansion(x,1u)); }
    virtual Vector< Differential<Float> > expansion(const Vector<Float>& x, const ushort& s) const {
        return this->_expansion(x,s); }
    virtual Vector< Differential<Interval> > expansion(const Vector<Interval>& x, const ushort& s) const {
        return this->_expansion(x,s); }

    // TODO: Find a better way for writing functions which can handle transformations which may not have a
    // write() method or operator<<.
    virtual std::ostream& write(std::ostream& os) const  {
        return os << "Function( result_size="<<this->result_size()<<", argument_size="<<this->argument_size()<<" )"; }
  private:
    template<class X> Vector< Differential<X> > _expansion(const Vector<X>& x, const ushort& s) const {
        const uint rs=this->T::result_size();
        const uint as=this->T::argument_size();
        Vector< Differential<X> > dx(as,Differential<X>(as,s));
        Vector< Differential<X> > dr(rs,Differential<X>(as,s));
        for(uint i=0; i!=as; ++i) { dx[i]=x[i]; }
        for(uint i=0; i!=as; ++i) { dx[i][i]=1; }
        this->T::compute(dr,dx);
        return dr;
    }
};


/*

// A wrapper for transformations
// This class is for internal use only; we can easily specify parameter types.
template<class F>
class FunctionTemplate
    : public FunctionInterface
{
  private:
    Vector<Float> p;
  private:
    template<class R,class A, class P>
    void _compute(R& r, const A& x, const P& p) const {
        static_cast<const F&>(*this)._compute(r,x,p); }
    template<class X> Vector<X> _evaluate(const Vector<X>& x) const {
        return static_cast<const F&>(*this)._evaluate(x); }
  protected:
    FunctionTemplate() { }
  public:
    virtual Vector<Float> evaluate(const Vector<Float>& x) const {
        Vector<Float> r(this->result_size()); this->_compute(r,x,p); return r; }
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const {
        Vector<Interval> r(this->result_size()); this->_compute(r,x,p); return r; }
    virtual Vector<TaylorModel> evaluate(const Vector<TaylorModel>& x) const {
        Vector<TaylorModel> r=TaylorModel::zeros(this->result_size(),this->argument_size());
        this->_compute(r,x,p); return r; }
    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const {
        Vector< Differential<Float> > r(this->result_size()); this->_compute(r,x,p); return r; }
    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const {
        Vector< Differential<Interval> > r(this->result_size()); this->_compute(r,x,p); return r; }

    virtual Matrix<Float> jacobian(const Vector<Float>& x) const {
        return this->_expansion(x,1u).jacobian(); }
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const {
        return this->_expansion(x,1u).jacobian(); }
    virtual Vector< Differential<Float> > expansion(const Vector<Float>& x, const ushort& s) const {
        return this->_expansion(x,s); }
    virtual Vector< Differential<Interval> > expansion(const Vector<Interval>& x, const ushort& s) const {
        return this->_expansion(x,s); }

    // TODO: Find a better way for writing functions which can handle transformations which may not have a
    // write() method or operator<<.
    virtual std::ostream& write(std::ostream& os) const  {
        return os << "Function( result_size="<<this->result_size()<<", argument_size="<<this->argument_size()<<" )"; }

  private:
    template<class X> Vector< Differential<X> > _expansion(const Vector<X>& x, const ushort& s) const {
        const uint rs=this->result_size();
        const uint as=this->argument_size();
        Vector< Differential<X> > dx(as,Differential<X>(as,s));
        Vector< Differential<X> > dr(rs,Differential<X>(as,s));
        for(uint i=0; i!=as; ++i) { dx[i]=x[i]; }
        for(uint i=0; i!=as; ++i) { dx[i][i]=1; }
        this->_compute(dr,dx,p);
        return dr;
    }
};

*/

// A wrapper for transformations
// This class is for internal use only; we can easily specify parameter types.
template<class F>
class FunctionTemplate
    : public FunctionInterface
{
  private:
    template<class X> Vector<X> _evaluate(const Vector<X>& x) const {
        return static_cast<const F&>(*this)._evaluate(x); }
  protected:
    FunctionTemplate() { }
  public:
    virtual Vector<Float> evaluate(const Vector<Float>& x) const {
        return this->_evaluate(x); }
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const {
        return this->_evaluate(x); }

    virtual Vector<TaylorModel> evaluate(const Vector<TaylorModel>& x) const {
        return this->_evaluate(x); }

    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const {
        return this->_evaluate(x); }
    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const {
        return this->_evaluate(x); }

    virtual Matrix<Float> jacobian(const Vector<Float>& x) const {
        return Ariadne::jacobian(this->_evaluate(Differential<Float>::variables(1u,x))); }
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const {
        return Ariadne::jacobian(this->_evaluate(Differential<Interval>::variables(1u,x))); }

    // TODO: Find a better way for writing functions which can handle transformations which may not have a
    // write() method or operator<<.
    virtual std::ostream& write(std::ostream& os) const  {
        return os << "Function( result_size="<<this->result_size()<<", argument_size="<<this->argument_size()<<" )"; }
};


template<uint RS, uint AS, uint PS=0u, uint SM=255u>
struct FunctionData
{
    const uint result_size() const { return RS; }
    const uint argument_size() const { return AS; }
    const uint parameter_size() const { return PS; }
    const uint smoothness() const { return SM; }
};

// FIXME: Make interval computations use interval version of parameters!
template<class T> struct FunctionWrapper
    : public T
{
    FunctionWrapper() : T(), _p(Vector<Float>(T::parameter_size())) { }
    FunctionWrapper(const Vector<Float>& p) : T(), _p(p) { ARIADNE_ASSERT(p.size()==T::parameter_size()); }
    FunctionWrapper(const Vector<Interval>& p) : T(), _p(midpoint(p)) { ARIADNE_ASSERT(p.size()==T::parameter_size()); }
    const Vector<Float>& parameters() const { return this->_p; }
    template<class R, class A> void compute(R& r, const A& x) const { this->T::compute(r,x,_p); }
  private:
    Vector<Float> _p;
};

//! \brief A wrapper for converting templated C++ functions to %Ariadne functions.
//!
//! Given a C++ class T with a (static) template method
//!   <code>template<class R, class A, class P> compute(R& r, const A& a, const P& p);</code>
//! the type <code>Function<T></code> is an Ariadne function defined by \f$r=f(a)\f$.
//! The constructor for Function<T> takes a Vector<Float> argument which is used for \a p.
//!
//! The class T must also define meta-data <c>result_size(), argument_size(), parameter_size()
//! and smoothness()</c>. These are most easily defined by inheriting from the
//! <tt>FunctionData<RS,AS,PS,SM=SMOOTH></tt> class.
//!
//! The constant \a SMOOTH is used for an arbitrarily-differentiable function.
template<class T> class Function
    : public FunctionBase< FunctionWrapper<T> >
{
  public:
    Function() : FunctionBase< FunctionWrapper<T> >() { }
    Function(const Vector<Float>& p) : FunctionBase< FunctionWrapper<T> >(p) { }
    Function(const Vector<Interval>& p) : FunctionBase< FunctionWrapper<T> >(p) { }
};




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

    virtual Vector<ExpressionInterface>* clone() const {
        return new Vector<ExpressionInterface>(*this); }

    virtual size_type result_size() const { return _expressions.size(); }
    virtual size_type argument_size() const { return _expressions[0]->argument_size(); }
    virtual smoothness_type smoothness() const {
        smoothness_type res=_expressions[0]->smoothness();
        for(uint i=1; i!=this->size(); ++i) {
            res=max(res,_expressions[i]->smoothness()); } return res; }

    virtual std::ostream& write(std::ostream& os) const {
        for(uint i=0; i!=this->result_size(); ++i) { os<<(*this)[i]; } return os; }
  private:
    template<class Y> inline Vector<Y> _evaluate(const Vector<Y>& x) const {
        Vector<Y> res(this->result_size());
        for(size_type i=0; i!=this->result_size(); ++i) {
            res[i]=this->_expressions[i]->evaluate(x); }
        return res;
    }
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
    template<class X> inline Vector<X> _evaluate(const Vector<X>& x) const {
        return this->_g->evaluate(this->_f->evaluate(x)); }
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
    template<class X> inline Vector<X> _evaluate(const Vector<X>& x) const {
        return join(_f1->evaluate(x),_f2->evaluate(x)); }
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
    template<class X> inline Vector<X> _evaluate(const Vector<X>& x) const {
        return combine(_f1->evaluate(project(x,range(0,_f1->argument_size()))),
                       _f2->evaluate(project(x,range(_f1->argument_size(),this->argument_size())))); }
  private:
    shared_ptr<const FunctionInterface> _f1;
    shared_ptr<const FunctionInterface> _f2;
};



//! A constant function \f$ x\mapsto c\f$ from \f$R^m\f$ to \f$R^n\f$.
class ConstantFunction
    : public FunctionBase<ConstantTransformation>
{
  public:
    //! A constant function with value \f$c\in\R^n\f$ and domain \f$\R^m\f$.
    ConstantFunction(const Vector<Float>& c, uint m)
        : FunctionBase<ConstantTransformation>(c,m) { }
    std::ostream& write(std::ostream& os) const {
        return os << "ConstantFunction( argument_size=" << this->argument_size()
                  << ", c=" << this->c() << " )"; }

};


//! A projection function \f$ x'_i= x_{p(i)}\f$.
class ProjectionFunction
    : public FunctionBase<ProjectionTransformation>
{
  public:
    //! Construct the identity function in dimension \a n.
    ProjectionFunction(uint n)
        : FunctionBase<ProjectionTransformation>(Range(0,n),n) { }
    //! Construct the projection functions which maps variables \f$x_i,\ldots,x_{i+m-1}\f$ in \f$\R^n\f$ to \f$x'_0\ldots x'_{m-1}\f$.
    ProjectionFunction(uint m, uint n, uint i)
        : FunctionBase<ProjectionTransformation>(m,n,i) { }
    //! Construct a projection function based on the array \a p.
    ProjectionFunction(const array<uint>& p, uint as)
        : FunctionBase<ProjectionTransformation>(p,as) { }
    std::ostream& write(std::ostream& os) const {
        return os << "ProjectionFunction( p=" << this->p() << " )"; }
};


//! The identity function \f$ x\mapsto x\f$ in \f$\R^n\f$.
class IdentityFunction
    : public FunctionBase<IdentityTransformation>
{
  public:
    //! Construct the identity function in dimension \a n.
    IdentityFunction(uint n)
        : FunctionBase<IdentityTransformation>(n) { }
    std::ostream& write(std::ostream& os) const {
        return os << "IdentityFunction( size=" << this->result_size() << " )"; }
};


//! An scaling function \f$x_i' = o_i+l_ix_i\f$.
class ScalingFunction
    : public FunctionBase<ScalingTransformation>
{
  public:
    //! Construct an affine function from the matrix \a A and vector \a b.
    ScalingFunction(const Vector<Float>& o, const Vector<Float>& l)
        : FunctionBase<ScalingTransformation>(o,l) { }
    explicit ScalingFunction(const Vector<Interval>& bx)
        : FunctionBase<ScalingTransformation>(bx) { }
    std::ostream& write(std::ostream& os) const {
        return os << "ScalingFunction( o=" << this->origin() << ", l=" << this->lengths() << " )"; }
};


/*
//! An affine function \f$x\mapsto Ax+b\f$.
class AffineFunction
    : public FunctionBase<AffineTransformation>
{
  public:
    //! Construct an affine function from the matrix \a A and vector \a b.
    AffineFunction(const Matrix<Float>& A, const Vector<Float>& b)
        : FunctionBase<AffineTransformation>(A,b) { }
    std::ostream& write(std::ostream& os) const {
        return os << "AffineFunction( A=" << this->A() << ", b=" << this->b() << " )"; }
};
*/

//! An affine function \f$x\mapsto Ax+b\f$.
class AffineFunction
    : public FunctionInterface
{
  public:
    //! Construct an affine function from the matrix \a A and vector \a b.
    AffineFunction(const Matrix<Interval>& A, const Vector<Interval>& b)
        : _fA(midpoint(A)), _fb(midpoint(b)), _iA(A), _ib(b) { ARIADNE_ASSERT(A.row_size()==b.size()); }

    const Matrix<Float>& A() const { return _fA; }
    const Vector<Float>& b() const { return _fb; }

    virtual AffineFunction* clone() const { return new AffineFunction(*this); }

    virtual size_type result_size() const { return _fA.row_size(); }
    virtual size_type argument_size() const { return _fA.column_size(); }
    virtual smoothness_type smoothness() const { return 255u; }

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
        return Ariadne::jacobian(this->evaluate(Differential<Float>::variables(1u,x))); }
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const {
        return Ariadne::jacobian(this->evaluate(Differential<Interval>::variables(1u,x))); }

    virtual std::ostream& write(std::ostream& os) const {
        return os << "AffineFunction( A="<<_iA<<", b="<<_ib<<" )"; }
  private:
    Matrix<Float> _fA; Vector<Float> _fb;
    Matrix<Interval> _iA; Vector<Interval> _ib;
};

//! A polynomial function.
class PolynomialFunction
    : public FunctionInterface
{
  public:
    PolynomialFunction(const Vector< Polynomial<Float> >& p) : _fp(p), _ip(p) { }
    PolynomialFunction(const Vector< Polynomial<Interval> >& p) : _fp(midpoint(p)), _ip(p) { }

    template<class E> PolynomialFunction(const ublas::vector_expression<E>& p)
        : _ip(static_cast< Vector< Polynomial<Interval> > >(p)) { _fp=midpoint(_ip); }

    virtual PolynomialFunction* clone() const { return new PolynomialFunction(*this); }

    virtual size_type result_size() const { return _fp.size(); }
    virtual size_type argument_size() const { return _fp[0].argument_size(); }
    virtual smoothness_type smoothness() const { return 255u; }

    virtual Vector<Float> evaluate(const Vector<Float>& x) const {
        return Ariadne::evaluate(_fp,x); }
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const {
        return Ariadne::evaluate(_ip,x); }
    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const {
        return Ariadne::evaluate(_fp,x); }

    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const {
        return Ariadne::evaluate(_ip,x); }
    virtual Vector<TaylorModel> evaluate(const Vector<TaylorModel>& x) const {
        return Ariadne::evaluate(_ip,x); }

    virtual Matrix<Float> jacobian(const Vector<Float>& x) const {
        return Ariadne::jacobian(this->evaluate(Differential<Float>::variables(1u,x))); }
    virtual Matrix<Interval> jacobian(const Vector<Interval>& x) const {
        return Ariadne::jacobian(this->evaluate(Differential<Interval>::variables(1u,x))); }

    virtual std::ostream& write(std::ostream& os) const {
        return os << "PolynomialFunction"<<_ip; }
  private:
    Vector< Polynomial<Float> > _fp;
    Vector< Polynomial<Interval> > _ip;
};


//! A polynomial expression.
class PolynomialExpression
    : public ExpressionInterface
{
  public:
    PolynomialExpression(const Polynomial<Float>& p) : _fp(p), _ip(p) { }
    PolynomialExpression(const Polynomial<Interval>& p) : _fp(midpoint(p)), _ip(p) { }

    virtual PolynomialExpression* clone() const { return new PolynomialExpression(*this); }

    virtual size_type argument_size() const { return _fp.argument_size(); }
    virtual smoothness_type smoothness() const { return 255u; }

    virtual Float evaluate(const Vector<Float>& x) const {
        return Ariadne::evaluate(_fp,x); }
    virtual Interval evaluate(const Vector<Interval>& x) const {
        return Ariadne::evaluate(_ip,x); }
    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const {
        return Ariadne::evaluate(_fp,x); }

    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const {
        return Ariadne::evaluate(_ip,x); }
    virtual TaylorModel evaluate(const Vector<TaylorModel>& x) const {
        return Ariadne::evaluate(_ip,x); }

    virtual std::ostream& write(std::ostream& os) const {
        return os << "PolynomialExpression"<<_ip; }
  private:
    Polynomial<Float> _fp;
    Polynomial<Interval> _ip;
};


} // namespace Ariadne

#endif
