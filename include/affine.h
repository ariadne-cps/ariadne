/***************************************************************************
 *            affine.h
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

/*! \file affine.h
 *  \brief Affine scalar and vector functions
 */

#ifndef ARIADNE_AFFINE_H
#define ARIADNE_AFFINE_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "function_interface.h"

#include "macros.h"
#include "pointer.h"

#include "vector.h"
#include "matrix.h"

namespace Ariadne {

typedef std::string String;
template<class T> class Array;
template<class K, class V> class Map;

class Real;
template<class T> class Variable;
template<class R> class Expression;

//! An affine expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b\f$.
class ConstantScalarFunction
    : public ScalarFunctionInterface
{
  public:
    operator Interval() const;
};


//! An affine expression \f$f:\R^n\rightarrow\R\f$ given by \f$f(x)=\sum_{i=0}^{n-1} a_i x_i + b\f$.
class ScalarAffineFunction
    : public ScalarFunctionInterface
{
  public:
    explicit ScalarAffineFunction(uint n) : _c(0), _g(n) { }
    ScalarAffineFunction(const Vector<Float>& g, const Float& c) : _c(c), _g(g) { }
    ScalarAffineFunction(const Vector<Interval>& g, const Interval& c) : _c(c), _g(g) { }
    ScalarAffineFunction(uint as, double c, double g0, ...) : _c(c), _g(as) {
        _g[0]=g0; va_list args; va_start(args,g0);
        for(uint i=1; i!=as; ++i) { _g[i]=va_arg(args,double); } }
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

    const Vector<Interval>& a() const { return _g; }
    const Interval& b() const { return _c; }

    Vector<Interval> gradient() const { return this->_g; }
    Interval value() const { return _c; }

    virtual ScalarAffineFunction* clone() const { return new ScalarAffineFunction(*this); }

    virtual SizeType argument_size() const { return _g.size(); }
    virtual SmoothnessType smoothness() const { return SMOOTH; }
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
    virtual Differential<TaylorModel> evaluate(const Vector< Differential<TaylorModel> >& x) const {
        Differential<TaylorModel> r=x[0]*0+_c; for(uint j=0; j!=_g.size(); ++j) { r+=_g[j]*x[j]; } return r; };

//    virtual ConstantScalarFunction* derivative(uint j) const { return new ConstantScalarFunction(_g.size(),_g[j]); }
    virtual ScalarFunctionInterface* derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }


    virtual std::ostream& write(std::ostream& os) const;
  private:
    friend ScalarAffineFunction operator-(const ScalarAffineFunction&);
    friend ScalarAffineFunction operator+(const ScalarAffineFunction&, const ScalarAffineFunction&);
    friend ScalarAffineFunction operator-(const ScalarAffineFunction&, const ScalarAffineFunction&);
    friend ScalarAffineFunction operator*(const Interval&, const ScalarAffineFunction&);
    friend ScalarAffineFunction operator*(const ScalarAffineFunction&, const Interval&);
    friend ScalarAffineFunction operator/(const ScalarAffineFunction&, const Interval&);
    friend Interval derivative(const ScalarAffineFunction&, uint);
    template<class R, class A>
    void compute(R& r, const A& x) const {
        r=x[0]*0+_c; for(uint j=0; j!=argument_size(); ++j) { r+=_g[j]*x[j]; } }
  private:
    Interval _c;
    Vector<Interval> _g;
};

//! \relates ScalarAffineFunction
//! \brief Negation of an ScalarAffineFunction.
inline ScalarAffineFunction operator-(const ScalarAffineFunction& f) {
    return ScalarAffineFunction(Vector<Interval>(-f._g),-f._c); }
//! \relates ScalarAffineFunction
//! \brief Addition of two AffineExpressions.
inline ScalarAffineFunction operator+(const ScalarAffineFunction& f1, const ScalarAffineFunction& f2) {
    return ScalarAffineFunction(Vector<Interval>(f1._g+f2._g),f1._c+f2._c); }
//! \relates ScalarAffineFunction
//! \brief Subtraction of two AffineExpressions.
inline ScalarAffineFunction operator-(const ScalarAffineFunction& f1, const ScalarAffineFunction& f2) {
    return ScalarAffineFunction(Vector<Interval>(f1._g-f2._g),f1._c-f2._c); }
//! \relates ScalarAffineFunction
//! \brief Scalar multiplication of an ScalarAffineFunction.
inline ScalarAffineFunction operator*(const Interval& c, const ScalarAffineFunction& f) {
    return ScalarAffineFunction(Vector<Interval>(c*f._g),c*f._c); }
//! \relates ScalarAffineFunction
//! \brief Scalar multiplication of an ScalarAffineFunction.
inline ScalarAffineFunction operator*(const ScalarAffineFunction& f, const Interval& c) { return c*f; }
//! \relates ScalarAffineFunction
//! \brief Scalar division of an ScalarAffineFunction.
inline ScalarAffineFunction operator/(const ScalarAffineFunction& f, const Interval& c) { return (1/c)*f; }
//! \relates ScalarAffineFunction
//! \brief The derivative of an ScalarAffineFunction gives a constant.
inline Interval derivative(const ScalarAffineFunction& f, uint k) { return f._g[k]; }

inline std::ostream& ScalarAffineFunction::write(std::ostream& os) const {
    const ScalarAffineFunction& f=*this;
    os<<"A(";
    if(f.b()!=0) { os<<f.b(); }
    for(uint j=0; j!=f.argument_size(); ++j) {
        if(f.a()[j]!=0) {
            if(f.a()[j]>0) { os<<"+"; } else { os<<"-"; }
            if(abs(f.a()[j])!=1) { os<<abs(f.a()[j])<<"*"; }
            //ss<<char('x'+j);
            os<<"x"<<j;
        }
    }
    os<<")";
    return os;
}

/*
inline ScalarAffineFunction::ScalarAffineFunction(RealExpressionPointer e, const Map<String,SizeType>& s) {
    ScalarAffineFunction& r=*this;
    const ExpressionInterface<Real>* eptr=e.operator->();
    Operator op=e->type();
    const ConstantExpression<Real>* cptr;
    const VariableExpression<Real>* vptr;
    const CoordinateExpression<Real>* iptr;
    const UnaryExpression<Real,Operator>* uptr;
    const BinaryExpression<Real,Operator>* bptr;
    const ConstantExpression<Real>* cptr1;
    const ConstantExpression<Real>* cptr2;

    switch(op) {
        case CNST:
            cptr=dynamic_cast<const ConstantExpression<Real>*>(eptr);
            r = ScalarAffineFunction::constant(s.size(),cptr->value()); break;
        case VAR:
            vptr=dynamic_cast<const VariableExpression<Real>*>(eptr);
            r = ScalarAffineFunction::variable(s.size(),s[vptr->name()]); break;
        case IND:
            iptr=dynamic_cast<const CoordinateExpression<Real>*>(eptr);
            r = ScalarAffineFunction::variable(s.size(),iptr->coordinate()); break;
        case NEG:
            uptr=static_cast<const UnaryExpression<Real,Operator>*>(uptr);
            r = -ScalarAffineFunction(uptr->_arg_ptr,s); break;
        case ADD:
            bptr=static_cast<const BinaryExpression<Real,Operator>*>(eptr);
            r = ScalarAffineFunction(bptr->_arg1_ptr,s)+ScalarAffineFunction(bptr->_arg2_ptr,s); break;
        case SUB:
            bptr=static_cast<const BinaryExpression<Real,Operator>*>(eptr);
            r = ScalarAffineFunction(bptr->_arg1_ptr,s)-ScalarAffineFunction(bptr->_arg2_ptr,s); break;
        case DIV:
            bptr=dynamic_cast<const BinaryExpression<Real,Operator>*>(eptr);
            cptr=dynamic_cast<const ConstantExpression<Real>*>(bptr->_arg2_ptr.operator->());
            ARIADNE_ASSERT_MSG(cptr,"Cannot convert expression "<<*e<<" to affine form.");
            r = ScalarAffineFunction(bptr->_arg1_ptr,s)/cptr->value(); break;
        case MUL:
            bptr=static_cast<const BinaryExpression<Real,Operator>*>(eptr);
            cptr1=dynamic_cast<const ConstantExpression<Real>*>(bptr->_arg1_ptr.operator->());
            cptr2=dynamic_cast<const ConstantExpression<Real>*>(bptr->_arg2_ptr.operator->());
            ARIADNE_ASSERT_MSG(cptr1 || cptr2,"Cannot convert expression "<<*e<<" to affine form.");
            if(cptr1) { r =  cptr1->value() * ScalarAffineFunction(bptr->_arg2_ptr,s); }
            else { r =  ScalarAffineFunction(bptr->_arg1_ptr,s) * cptr2->value(); }
            break;
        default:
            ARIADNE_ASSERT_MSG(false,"Cannot convert expression "<<*e<<" to affine form.");
    }
}
*/


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

inline std::ostream& AffineFunction::write(std::ostream& os) const {
    //return os << "AffineFunction( A="<<midpoint(_iA)<<", b="<<midpoint(_ib)<<" )";
    const AffineFunction& f=*this;
    for(uint i=0; i!=f.result_size(); ++i) {
        os<<(i==0?"AffineFunction( ":"; ");
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


} // namespace Ariadne

#endif /* ARIADNE_AFFINE_H */
