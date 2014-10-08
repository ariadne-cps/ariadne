/***************************************************************************
 *            formula.h
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


/*! \file formula.h
 *  \brief Formulae over variables
 */

#ifndef ARIADNE_FORMULA_H
#define ARIADNE_FORMULA_H

#include <cstdarg>
#include <iostream>
#include <iomanip>
#include <string>


#include "utility/macros.h"
#include "utility/pointer.h"
#include "utility/container.h"
#include "utility/stlio.h"

#include "expression/operators.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/expansion.h"
#include <boost/concept_check.hpp>

namespace Ariadne {

template<class X> class Formula;
typedef Formula<ApproximateNumber> ApproximateFormula;
typedef Formula<ValidatedNumber> ValidatedFormula;

struct Index {
    Nat _i;
  public:
    explicit Index(Nat i) : _i(i) { }
    operator Nat () const { return _i; }
};
inline OutputStream& operator<<(OutputStream& os, const Index& ind) {
    return os << Nat(ind); }

template<class X> class FormulaNode;

//! \brief A formula defining a real function.
//!
//! The Formula class is implemented as a directed acyclic graph, with
//! each node being an atomic operation.
template<class X>
class Formula {
    typedef Index I;
  public:
    typedef X NumericType;
    typedef X ConstantType;
    typedef I IndexType;
  public:
    Formula(const FormulaNode<X>* fptr) : _root(fptr) { }
    Formula(counted_pointer< const FormulaNode<X> > fptr) : _root(fptr) { }
  public:
    //! \brief Construct the constant expression with the default value of \a X.
    Formula();
    //! \brief Construct an expression from a value.
    Formula(const X& c);
    //! \brief Construct an expression from an index.
    Formula(const Index& i);

    //! \brief Set equal to a constant.
    Formula<X>& operator=(const X& c);
    Formula<X>& operator=(double c);
  public:
    static Formula<X> zero();
    template<class R> static Formula<X> constant(const R& c);
    static Formula<X> coordinate(Nat i);
    static Vector< Formula<X> > identity(Nat n);
  public:
    const Operator& op() const;
    OperatorCode code() const;
    OperatorKind kind() const;
    const X& val() const;
    const I& ind() const;
    const Formula<X>& arg() const;
    const Int& num() const;
    const Formula<X>& arg1() const;
    const Formula<X>& arg2() const;
  public:
    const FormulaNode<X>* node_ptr() const { return _root.operator->(); }
  private:
    counted_pointer< const FormulaNode<X> > _root;
};
template<class X> OutputStream& operator<<(OutputStream& os, const Formula<X>& f);

template<class X>
struct FormulaNode {
    mutable uint count;
    Operator op;
    virtual ~FormulaNode();
    explicit FormulaNode(const Operator& o) : count(0u), op(o) { }
    explicit FormulaNode(OperatorCode cd, OperatorKind knd) : count(0u), op(cd,knd) { }
};

template<class X> struct ConstantFormulaNode : public FormulaNode<X> {
    X val;
    ConstantFormulaNode(const X& v) : FormulaNode<X>(CNST,NULLARY), val(v) { }
};
template<class X> struct IndexFormulaNode : public FormulaNode<X> {
    Index ind;
    IndexFormulaNode(Nat i) : FormulaNode<X>(IND,COORDINATE), ind(i) { }
    IndexFormulaNode(const Index& i) : FormulaNode<X>(IND,COORDINATE), ind(i) { }
};
template<class X, class A=X> struct UnaryFormulaNode : public FormulaNode<X> {
    Formula<X> arg;
    UnaryFormulaNode(const Operator& op, Formula<X> const& a) : FormulaNode<X>(op), arg(a) { }
};
template<class X, class A1=X, class A2=A1> struct BinaryFormulaNode {
    Formula<X> arg1; Formula<X> arg2;
};
template<class X> struct BinaryFormulaNode<X> : public FormulaNode<X> {
    Formula<X> arg1; Formula<X> arg2;
    BinaryFormulaNode(const Operator& op, Formula<X> const& a1, Formula<X> const& a2)
        : FormulaNode<X>(op), arg1(a1), arg2(a2) { }
};
template<class X> struct ScalarFormulaNode : public UnaryFormulaNode<X> {
    Int num;
    ScalarFormulaNode(const Operator& op, Formula<X> const& a, Int n)
        : UnaryFormulaNode<X>(op,a), num(n) { }
};

template<class X> FormulaNode<X>::~FormulaNode() { }

template<class X> inline const Operator& Formula<X>::op() const {
    return node_ptr()->op; }
template<class X> inline OperatorCode Formula<X>::code() const {
    return node_ptr()->op.code(); }
template<class X> inline OperatorKind Formula<X>::kind() const {
    return node_ptr()->op.kind(); }
template<class X> inline const X& Formula<X>::val() const {
    return static_cast<const ConstantFormulaNode<X>*>(node_ptr())->val; }
template<class X> inline const Index& Formula<X>::ind() const {
    return static_cast<const IndexFormulaNode<X>*>(node_ptr())->ind; }
template<class X> inline const Formula<X>& Formula<X>::arg() const {
    return static_cast<const UnaryFormulaNode<X>*>(node_ptr())->arg; }
template<class X> inline const Int& Formula<X>::num() const {
    return static_cast<const ScalarFormulaNode<X>*>(node_ptr())->num; }
template<class X> inline const Formula<X>& Formula<X>::arg1() const {
    return static_cast<const BinaryFormulaNode<X>*>(node_ptr())->arg1; }
template<class X> inline const Formula<X>& Formula<X>::arg2() const {
    return static_cast<const BinaryFormulaNode<X>*>(node_ptr())->arg2; }


template<class X> inline Formula<X>::Formula() : _root(new ConstantFormulaNode<X>(X())) { }
template<class X> inline Formula<X>::Formula(const X& c) : _root(new ConstantFormulaNode<X>(c)) { }
template<class X> inline Formula<X>::Formula(const Index& i) : _root(new IndexFormulaNode<X>(i)) { }
template<class X> inline Formula<X>& Formula<X>::operator=(const X& c) { return *this=Formula(c); }
template<class X> inline Formula<X>& Formula<X>::operator=(double c) { return *this=Formula<X>::constant(static_cast<X>(c)); }

template<class X> inline Formula<X> Formula<X>::zero() {
    return new ConstantFormulaNode<X>(numeric_cast<X>(0)); }
template<class X> template<class R> inline Formula<X> Formula<X>::constant(const R& c) {
    return new ConstantFormulaNode<X>(numeric_cast<X>(c)); }
template<class X> inline Formula<X> Formula<X>::coordinate(Nat j) {
    return new IndexFormulaNode<X>(Index(j)); }
template<class X> inline Vector< Formula<X> > Formula<X>::identity(Nat n) {
    Vector< Formula<X> > r(n); for(uint i=0; i!=n; ++i) { r[i]=Formula<X>::coordinate(i); } return r; }

template<class X, class R> inline Formula<X> make_formula(const R& c) {
    return new ConstantFormulaNode<X>(static_cast<X>(c)); }
template<class X, class R> inline Formula<X> make_formula(Cnst op, const R& c) {
    return new ConstantFormulaNode<X>(static_cast<X>(c)); }
template<class X> inline Formula<X> make_formula(Ind op, Nat j) {
    return new IndexFormulaNode<X>(j); }
template<class X> inline Formula<X> make_formula(const Operator& op, const Formula<X>& arg) {
    return new UnaryFormulaNode<X>(op,arg); }
template<class X> inline Formula<X> make_formula(const Operator& op, const Formula<X>& arg1, const Formula<X>& arg2) {
    return new BinaryFormulaNode<X>(op,arg1,arg2); }
template<class X> inline Formula<X> make_formula(const Operator& op, const Formula<X>& arg, Int num) {
    return new ScalarFormulaNode<X>(op,arg,num); }

template<class X> inline Formula<X>& operator+=(Formula<X>& f1, const Formula<X>& f2) { Formula<X> r=f1+f2; return f1=r; }
template<class X> inline Formula<X>& operator*=(Formula<X>& f1, const Formula<X>& f2) { Formula<X> r=f1*f2; return f1=r; }

template<class X> inline Formula<X> operator+(const Formula<X>& f) { return make_formula(Pos(),f); }
template<class X> inline Formula<X> operator-(const Formula<X>& f) { return make_formula(Neg(),f); }
template<class X> inline Formula<X> operator+(const Formula<X>& f1, const Formula<X>& f2) { return make_formula(Add(),f1,f2); }
template<class X> inline Formula<X> operator-(const Formula<X>& f1, const Formula<X>& f2) { return make_formula(Sub(),f1,f2); }
template<class X> inline Formula<X> operator*(const Formula<X>& f1, const Formula<X>& f2) { return make_formula(Mul(),f1,f2); }
template<class X> inline Formula<X> operator/(const Formula<X>& f1, const Formula<X>& f2) { return make_formula(Div(),f1,f2); }

template<class X> inline Formula<X> neg(const Formula<X>& f) { return make_formula(Neg(),f); }
template<class X> inline Formula<X> rec(const Formula<X>& f) { return make_formula(Rec(),f); }
template<class X> inline Formula<X> sqr(const Formula<X>& f) { return make_formula(Sqr(),f); }
template<class X> inline Formula<X> pow(const Formula<X>& f, int n) { return make_formula(Pow(),f,n); }
template<class X> inline Formula<X> sqrt(const Formula<X>& f) { return make_formula(Sqrt(),f); }
template<class X> inline Formula<X> exp(const Formula<X>& f) { return make_formula(Exp(),f); }
template<class X> inline Formula<X> log(const Formula<X>& f) { return make_formula(Log(),f); }
template<class X> inline Formula<X> sin(const Formula<X>& f) { return make_formula(Sin(),f); }
template<class X> inline Formula<X> cos(const Formula<X>& f) { return make_formula(Cos(),f); }
template<class X> inline Formula<X> tan(const Formula<X>& f) { return make_formula(Tan(),f); }
template<class X> inline Formula<X> atan(const Formula<X>& f) { return make_formula(Atan(),f); }

template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator+(Formula<X> f, R c) { return f + make_formula<X>(c); }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator-(Formula<X> f, R c) { return f - make_formula<X>(c); }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator*(Formula<X> f, R c) { return f * make_formula<X>(c); }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator/(Formula<X> f, R c) { return f / make_formula<X>(c); }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator+(R c, Formula<X> f) { return make_formula<X>(c) + f; }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator-(R c, Formula<X> f) { return make_formula<X>(c) - f; }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator*(R c, Formula<X> f) { return make_formula<X>(c) * f; }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type operator/(R c, Formula<X> f) { return make_formula<X>(c) / f; }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type& operator+=(Formula<X>& f, const R& c) { return f+=make_formula<X>(c); }
template<class X, class R> inline typename EnableIfNumeric<R,Formula<X> >::Type& operator*=(Formula<X>& f, const R& c) { return f*=make_formula<X>(c); }


// Make a constant of type T with value c based on a prototype vector v
template<class X, class T> inline T make_constant(const X& c, const Vector<T>& v, typename EnableIf< Not< IsNumeric<T> >, Void >::Type* =0 ) {
    return v.zero_element()+numeric_cast<typename T::NumericType>(c);
}

// Make a constant of type T with value c based on a prototype vector v
template<class X, class T> inline T make_constant(const X& c, const Vector<T>& v, typename EnableIf<IsNumeric<T>,Void>::Type* = 0) {
    return v.zero_element()+numeric_cast<T>(c);
}

// Make a constant of type T with value c based on a prototype vector v
template<class X, class T> inline Formula<T> make_constant(const X& c, const Vector< Formula<T> >& v) {
    return Formula<T>::constant(c);
}

template<class X, class T> T evaluate(const Formula<X>& f, const Vector<T>& x) {
    switch(f.kind()) {
        case COORDINATE: return x[f.ind()];
        case NULLARY: return make_constant(f.val(),x);
        case UNARY: return compute(f.op(),evaluate(f.arg(),x));
        case BINARY: return compute(f.op(),evaluate(f.arg1(),x),evaluate(f.arg2(),x));
        case SCALAR: return compute(f.op(),evaluate(f.arg(),x),f.num());
        default: ARIADNE_FAIL_MSG("Cannot evaluate formula "<<f<<" on "<<x<<"; unknown operator "<<f.op()<<" of kind "<<f.kind()<<"\n");
    }
}


//! \brief Convert the expression with index type \c I to one with variables indexed by \a J.
template<class X, class T> const T& cached_evaluate(const Formula<X>& f, const Vector<T>& x, Map<const void*,T>& cache) {
    const FormulaNode<X>* fptr=f.node_ptr();
    if(cache.has_key(fptr)) { return cache.get(fptr); }
    switch(f.kind()) {
        case VARIABLE: return insert( cache, fptr, x[f.ind()] );
        case NULLARY: return insert( cache, fptr, make_constant(f.val(),x) );
        case UNARY: return insert( cache, fptr, compute(f.op(),cached_evaluate(f.arg(),x,cache)) );
        case BINARY: return insert( cache, fptr, compute(f.op(),cached_evaluate(f.arg1(),x,cache),cached_evaluate(f.arg2(),x,cache)) );
        case SCALAR: return insert( cache, fptr, compute(f.op(),cached_evaluate(f.arg(),x,cache),f.num()) );
        default: ARIADNE_FAIL_MSG("Cannot evaluate formula "<<f<<" on "<<x<<"; unknown operator "<<f.op()<<" of kind "<<f.kind()<<"\n");
    }
}

template<class X, class T> inline T cached_evaluate(const Formula<X>& f, const Vector<T>& v) {
    Map<const void*,T> cache;
    return cached_evaluate(f.handle(),v,cache);
}

template<class X, class T> Vector<T> cached_evaluate(const Vector< Formula<X> >& f, const Vector<T>& v) {
    assert(v.size()!=0);
    Vector<T> r(f.size());
    Map<const void*,T> cache;
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=cached_evaluate(f[i],v,cache);
    }
    return r;
}

template<class X> Formula<X> derivative(const Formula<X>& f, uint j)
{
    switch(f.op()) {
        case CNST:
            return Formula<X>::constant(0);
        case IND:
            if(f.ind()==j) { return Formula<X>::constant(1); }
            else { return Formula<X>::constant(0); }
        case ADD:
            return derivative(f.arg1(),j)+derivative(f.arg2(),j);
        case SUB:
            return derivative(f.arg1(),j)-derivative(f.arg2(),j);
        case MUL:
            return f.arg1()*derivative(f.arg2(),j)+derivative(f.arg1(),j)*f.arg2();
        case DIV:
            return derivative(f.arg1() * rec(f.arg2()),j);
        case NEG:
            return  - derivative(f.arg(),j);
        case REC:
            return  - derivative(f.arg(),j) * rec(sqr(f.arg()));
        case SQR:
            return static_cast<X>(2) * derivative(f.arg(),j) * f.arg();
        case EXP:
            return derivative(f.arg(),j) * f.arg();
        case LOG:
            return derivative(f.arg(),j) * rec(f.arg());
        case SIN:
            return derivative(f.arg(),j) * cos(f.arg());
        case COS:
            return -derivative(f.arg(),j) * sin(f.arg());
        case TAN:
            return derivative(f.arg(),j) * (static_cast<X>(1)-sqr(f.arg()));
        default:
            ARIADNE_THROW(std::runtime_error,"derivative(Formual<X>)",
                          "Cannot compute derivative of "<<f<<"\n");
    }
}

//! \brief Write to an output stream
template<class X> OutputStream& operator<<(std::ostream& os, const Formula<X>& f) {
    switch(f.op()) {
        //case CNST: return os << std::fixed << std::setprecision(4) << fptr->val;
        case CNST:
            os << f.val(); return os;
            //if(f.val()==0.0) { return os << 0.0; } if(abs(f.val())<1e-4) { os << std::fixed << f.val(); } else { os << f.val(); } return os;
        case IND:
            return os << "x" << f.ind();
        case ADD:
            return os << f.arg1() << '+' << f.arg2();
        case SUB:
            os << f.arg1() << '-';
            switch(f.arg2().op()) { case ADD: case SUB: os << '(' << f.arg2() << ')'; break; default: os << f.arg2(); }
            return os;
        case MUL:
            switch(f.arg1().op()) { case ADD: case SUB: case DIV: os << '(' << f.arg1() << ')'; break; default: os << f.arg1(); }
            os << '*';
            switch(f.arg2().op()) { case ADD: case SUB: os << '(' << f.arg2() << ')'; break; default: os << f.arg2(); }
            return os;
        case DIV:
            switch(f.arg1().op()) { case ADD: case SUB: case DIV: os << '(' << f.arg1() << ')'; break; default: os << f.arg1(); }
            os << '/';
            switch(f.arg2().op()) { case ADD: case SUB: case MUL: case DIV: os << '(' << f.arg2() << ')'; break; default: os << f.arg2(); }
            return os;
        case POW:
            return os << "pow" << '(' << f.arg() << ',' << f.num() << ')';
        default:
            switch(f.kind()) {
                case UNARY: return os << f.op() << "(" << f.arg() << ")";
                case BINARY: return os << f.op() << "(" << f.arg1() << "," << f.arg2() << ")";
                case COMPARISON: return os << "(" << f.arg1() << symbol(f.op()) << f.arg2() << ")";
                default: ARIADNE_FAIL_MSG("Cannot output formula with operator "<<f.op()<<" of kind "<<f.kind()<<"\n");
            }
    }
}

// Declare conversion operators from an expression
template<class X> class Expression;
template<class X> class Space;
Formula<Real> formula(const Expression<Real>& e, const Space<Real>& spc);

//! \ingroup FunctionModule
//! \brief Convert a power-series expansion into a formula using a version of Horner's rule.
//! See J. M. Pena and T. Sauer, "On the multivariate Horner scheme", SIAM J. Numer. Anal. 37(4) 1186-1197, 2000.
//!
template<class X> Formula<X> formula(const Expansion<X>& e)
{
    Vector< Formula<X> > identity(e.argument_size());
    for(uint i=0; i!=identity.size(); ++i) { identity[i]=Formula<X>::coordinate(i); }
    return horner_evaluate(e,identity);
}

// Class for which an object x produces coordinate \f$x_j\f$ when calling \c x[j].
struct Coordinate { Formula<Real> operator[](uint j) { return Formula<Real>::coordinate(j); } };

} // namespace Ariadne


#endif // ARIADNE_FORMULA_H
