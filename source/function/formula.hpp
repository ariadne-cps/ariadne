/***************************************************************************
 *            formula.hpp
 *
 *  Copyright 2008-17 Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file formula.hpp
 *  \brief Formulae over variables
 */

#ifndef ARIADNE_FORMULA_HPP
#define ARIADNE_FORMULA_HPP

#include <cstdarg>
#include <iostream>
#include <iomanip>
#include <string>


#include "../utility/macros.hpp"
#include "../utility/pointer.hpp"
#include "../utility/container.hpp"
#include "../utility/stlio.hpp"

#include "../numeric/operators.hpp"
#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/expansion.hpp"
#include "../algebra/operations.hpp"

namespace Ariadne {

template<class A, class X >
struct DispatchAlgebraOperations;

template<class Y> class Formula;
typedef Formula<ApproximateNumber> ApproximateFormula;
typedef Formula<ValidatedNumber> ValidatedFormula;
typedef Formula<EffectiveNumber> EffectiveFormula;

template<> inline Real compute(OperatorCode op, const EffectiveNumber& x1, const Real& x2) {
    return compute(op,Real(x1),x2);
}
template<template<class>class A> inline A<Real> compute(OperatorCode op, const EffectiveNumber& x1, const A<Real>& x2) {
    return compute(op,Real(x1),x2);
}

struct Index {
    Nat _i;
  public:
    explicit Index(Nat i) : _i(i) { }
    operator Nat () const { return _i; }
    friend OutputStream& operator<<(OutputStream& os, const Index& ind) {
        return os << Nat(ind); }
};

template<class Y> class FormulaNode;

template<class Y> inline Formula<Y> make_formula(const Y& c) {
    return Formula<Y>::constant(c); }
template<class Y> inline Formula<Y> make_formula(Cnst op, const Y& c) {
    return Formula<Y>::constant(c); }
template<class Y> inline Formula<Y> make_formula(Ind op, Nat j) {
    return Formula<Y>::index(j); }
template<class Y> inline Formula<Y> make_formula(const Operator& op, const Formula<Y>& arg) {
    return Formula<Y>::unary(op,arg); }
template<class Y> inline Formula<Y> make_formula(const Operator& op, const Formula<Y>& arg1, const Formula<Y>& arg2) {
    return Formula<Y>::binary(op,arg1,arg2); }
template<class Y> inline Formula<Y> make_formula(const Operator& op, const Formula<Y>& arg1, const Y& arg2) {
    return Formula<Y>::binary(op,arg1,make_formula(arg2)); }
template<class Y> inline Formula<Y> make_formula(const Operator& op, const Formula<Y>& arg, Int num) {
    return Formula<Y>::graded(op,arg,num); }
template<class Y, class OP> inline Formula<Y> make_formula(const OP& op, Y const& cnst, const Formula<Y>& arg) {
    OperatorCode op_code=static_cast<OperatorCode>((char)op.code()+((char)OperatorCode::SADD-(char)OperatorCode::ADD));
    return Formula<Y>::scalar(Operator(op_code,OperatorKind::SCALAR),cnst,arg); }

class FormulaOperations {
    template<class Y, class R, EnableIf<IsConstructible<Y,R>> =dummy> friend Formula<Y> operator+(Formula<Y> f, R c) { return f + Y(c); }
    template<class Y, class R, EnableIf<IsConstructible<Y,R>> =dummy> friend Formula<Y> operator-(Formula<Y> f, R c) { return f - Y(c); }
    template<class Y, class R, EnableIf<IsConstructible<Y,R>> =dummy> friend Formula<Y> operator*(Formula<Y> f, R c) { return f * Y(c); }
    template<class Y, class R, EnableIf<IsConstructible<Y,R>> =dummy> friend Formula<Y> operator/(Formula<Y> f, R c) { return f / Y(c); }
    template<class Y, class R, EnableIf<IsConstructible<Y,R>> =dummy> friend Formula<Y> operator+(R c, Formula<Y> f) { return Y(c)+f; }
    template<class Y, class R, EnableIf<IsConstructible<Y,R>> =dummy> friend Formula<Y> operator-(R c, Formula<Y> f) { return Y(c)-f; }
    template<class Y, class R, EnableIf<IsConstructible<Y,R>> =dummy> friend Formula<Y> operator*(R c, Formula<Y> f) { return Y(c)*f; }
    template<class Y, class R, EnableIf<IsConstructible<Y,R>> =dummy> friend Formula<Y> operator/(R c, Formula<Y> f) { return Y(c)/f; }

    template<class Y, class R, EnableIf<IsConstructible<Y,R>> =dummy> friend Formula<Y>& operator+=(Formula<Y>& f, const R& c) { return f+=Y(c); }
    template<class Y, class R, EnableIf<IsConstructible<Y,R>> =dummy> friend Formula<Y>& operator-=(Formula<Y>& f, const R& c) { return f-=Y(c); }
    template<class Y, class R, EnableIf<IsConstructible<Y,R>> =dummy> friend Formula<Y>& operator*=(Formula<Y>& f, const R& c) { return f*=Y(c); }
    template<class Y, class R, EnableIf<IsConstructible<Y,R>> =dummy> friend Formula<Y>& operator/=(Formula<Y>& f, const R& c) { return f/=Y(c); }
};

template<class Y> struct AlgebraOperations<Formula<Y>,Y> {
  public:
    template<class OP> static Formula<Y> apply(OP op, Formula<Y> const& f) {
        return Formula<Y>::unary(op,f); }
    template<class OP> static Formula<Y> apply(OP op, Formula<Y> const& f1, Formula<Y> const& f2) {
        return Formula<Y>::binary(op,f1,f2); }
    template<class OP> static Formula<Y> apply(OP op, Formula<Y> const& f1, Y const& c2) {
        return Formula<Y>::binary(op,f1,Formula<Y>::constant(c2)); }
    template<class OP> static Formula<Y> apply(OP op, Y const& c1, Formula<Y> const& f2) {
        return Formula<Y>::binary(op,Formula<Y>::constant(c1),f2); }
    static Formula<Y> apply(Pow op, Formula<Y> const& f, Int n) {
        return Formula<Y>::graded(op,f,n); }
};

//! \brief A formula defining a real function.
//!
//! The Formula class is implemented as a directed acyclic graph, with
//! each node being an atomic operation.
template<class Y>
class Formula
    : public DispatchTranscendentalAlgebraOperations<Formula<Y>,Y>
    , public FormulaOperations
{
    typedef Index I;
  public:
    typedef typename Y::Paradigm Paradigm;
    typedef Y NumericType;
    typedef Y ConstantType;
    typedef I IndexType;
  private:
    explicit Formula(const FormulaNode<Y>* fptr, PointerTag) : _root(fptr) { }
    explicit Formula(counted_pointer<const FormulaNode<Y>> fptr, PointerTag) : _root(fptr) { }
  public:
    //! \brief Construct the constant expression with the default value of \a Y.
    Formula();
    //! \brief Construct the constant expression with value \a c.
    Formula(Y const& c);
    //! \brief Set equal to a constant.
    Formula<Y>& operator=(const Y& c);
    template<class X, EnableIf<IsConstructible<Y,X>> =dummy> Formula<Y>& operator=(const X& c);
  public:
    //! \brief Return the constant formula zero.
    Formula<Y> create_zero() const;
    //! \brief Return the constant formula with value \a c.
    Formula<Y> create_constant(const Y& c) const;
  public:
    static Formula<Y> zero();
    static Formula<Y> constant(const Y& c);
    static Formula<Y> constant(Int c);
    static Formula<Y> coordinate(Nat i);
    static Vector<Formula<Y>> coordinates(Nat n);
    static Formula<Y> unary(const Operator& op, Formula<Y> const& a);
    static Formula<Y> binary(const Operator& op, Formula<Y> const& a1, Formula<Y> const& a2);
    static Formula<Y> graded(const Operator& op, Formula<Y> const& a1, Int n2);
    static Formula<Y> scalar(const Operator& op, Y const& c1, Formula<Y> const& a2);
    static Vector<Formula<Y>> identity(Nat n);
  public:
    const Operator& op() const;
    OperatorCode code() const;
    OperatorKind kind() const;
    const Y& val() const;
    const I& ind() const;
    const Formula<Y>& arg() const;
    const Y& cnst() const;
    const Int& num() const;
    const Formula<Y>& arg1() const;
    const Formula<Y>& arg2() const;
  public:
    Formula<Y> _derivative(SizeType j) const;
    OutputStream& _write(OutputStream& os) const;
  public:
    const FormulaNode<Y>* node_ptr() const { return _root.operator->(); }
  private:
    counted_pointer<const FormulaNode<Y>> _root;
};
template<class Y> Formula<Y> derivative(const Formula<Y>& f, SizeType j) {
    return f._derivative(j); }
template<class Y> OutputStream& operator<<(OutputStream& os, const Formula<Y>& f) {
    return f._write(os); }

template<class Y>
class FormulaNode {
  public:
    mutable Nat count;
    Operator op;
    virtual ~FormulaNode() = default;
    explicit FormulaNode(const Operator& o) : count(0u), op(o) { }
    explicit FormulaNode(OperatorCode cd, OperatorKind knd) : count(0u), op(cd,knd) { }
};

template<class Y> struct ConstantFormulaNode : public FormulaNode<Y> {
    Y val;
    ConstantFormulaNode(const Y& v) : FormulaNode<Y>(OperatorCode::CNST,OperatorKind::NULLARY), val(v) { }
};
template<class Y> struct IndexFormulaNode : public FormulaNode<Y> {
    Index ind;
    IndexFormulaNode(Nat i) : FormulaNode<Y>(OperatorCode::IND,OperatorKind::COORDINATE), ind(i) { }
    IndexFormulaNode(const Index& i) : FormulaNode<Y>(OperatorCode::IND,OperatorKind::COORDINATE), ind(i) { }
};
template<class Y, class A=Y> struct UnaryFormulaNode : public FormulaNode<Y> {
    Formula<Y> arg;
    UnaryFormulaNode(const Operator& oper, Formula<Y> const& a)
        : FormulaNode<Y>(oper), arg(a) { }
};
template<class Y, class A1=Y, class A2=A1> struct BinaryFormulaNode {
    Formula<Y> arg1; Formula<Y> arg2;
};
template<class Y> struct BinaryFormulaNode<Y> : public FormulaNode<Y> {
    Formula<Y> arg1; Formula<Y> arg2;
    BinaryFormulaNode(const Operator& oper, Formula<Y> const& a1, Formula<Y> const& a2)
        : FormulaNode<Y>(oper), arg1(a1), arg2(a2) { }
};
template<class Y> struct GradedFormulaNode : public UnaryFormulaNode<Y> {
    Int num;
    GradedFormulaNode(const Operator& oper, Formula<Y> const& a, Int n)
        : UnaryFormulaNode<Y>(oper,a), num(n) { }
};
template<class Y> struct ScalarFormulaNode : public UnaryFormulaNode<Y> {
    Y cnst;
    ScalarFormulaNode(const Operator& oper, Y const& c, Formula<Y> const& a)
        : UnaryFormulaNode<Y>(oper,a), cnst(c) { }
};

template<class Y> inline const Operator& Formula<Y>::op() const {
    return node_ptr()->op; }
template<class Y> inline OperatorCode Formula<Y>::code() const {
    return node_ptr()->op.code(); }
template<class Y> inline OperatorKind Formula<Y>::kind() const {
    return node_ptr()->op.kind(); }
template<class Y> inline const Y& Formula<Y>::val() const {
    return static_cast<const ConstantFormulaNode<Y>*>(node_ptr())->val; }
template<class Y> inline const Index& Formula<Y>::ind() const {
    return static_cast<const IndexFormulaNode<Y>*>(node_ptr())->ind; }
template<class Y> inline const Formula<Y>& Formula<Y>::arg() const {
    return static_cast<const UnaryFormulaNode<Y>*>(node_ptr())->arg; }
template<class Y> inline const Int& Formula<Y>::num() const {
    return static_cast<const GradedFormulaNode<Y>*>(node_ptr())->num; }
template<class Y> inline const Y& Formula<Y>::cnst() const {
    return static_cast<const ScalarFormulaNode<Y>*>(node_ptr())->cnst; }
template<class Y> inline const Formula<Y>& Formula<Y>::arg1() const {
    return static_cast<const BinaryFormulaNode<Y>*>(node_ptr())->arg1; }
template<class Y> inline const Formula<Y>& Formula<Y>::arg2() const {
    return static_cast<const BinaryFormulaNode<Y>*>(node_ptr())->arg2; }


template<class Y> inline Formula<Y>::Formula() : _root(new ConstantFormulaNode<Y>(Y())) { }
template<class Y> inline Formula<Y>::Formula(const Y& c) : _root(new ConstantFormulaNode<Y>(c)) { }
template<class Y> inline Formula<Y>& Formula<Y>::operator=(const Y& c) { return *this=Formula<Y>::constant(c); }
template<class Y> template<class X, EnableIf<IsConstructible<Y,X>>> inline Formula<Y>& Formula<Y>::operator=(const X& c) { return *this=Y(c); }

template<class Y> inline Formula<Y> Formula<Y>::create_zero() const { return Formula<Y>::constant(0); }
template<class Y> inline Formula<Y> Formula<Y>::create_constant(const Y& c) const { return Formula<Y>::constant(c); }

template<class Y> inline Formula<Y> Formula<Y>::zero() {
    return Formula<Y>(new ConstantFormulaNode<Y>(static_cast<Y>(0)),PointerTag()); }
template<class Y> inline Formula<Y> Formula<Y>::constant(const Y& c) {
    return Formula<Y>(new ConstantFormulaNode<Y>(c),PointerTag()); }
template<class Y> inline Formula<Y> Formula<Y>::coordinate(Nat j) {
    return Formula<Y>(new IndexFormulaNode<Y>(Index(j)),PointerTag()); }
template<class Y> inline Vector<Formula<Y>> Formula<Y>::coordinates(Nat n) {
    Vector<Formula<Y>> r(n); for(SizeType i=0; i!=n; ++i) { r[i]=Formula<Y>::coordinate(i); } return r; }
template<class Y> inline Formula<Y> Formula<Y>::unary(const Operator& op, Formula<Y> const& a) {
    return Formula<Y>(new UnaryFormulaNode<Y>(op,a),PointerTag()); }
template<class Y> inline Formula<Y> Formula<Y>::binary(const Operator& op, Formula<Y> const& a1, Formula<Y> const& a2) {
    return Formula<Y>(new BinaryFormulaNode<Y>(op,a1,a2),PointerTag()); }
template<class Y> inline Formula<Y> Formula<Y>::graded(const Operator& op, Formula<Y> const& a1, Int n2) {
    return Formula<Y>(new GradedFormulaNode<Y>(op,a1,n2),PointerTag()); }
template<class Y> inline Formula<Y> Formula<Y>::scalar(const Operator& op, Y const& c1, Formula<Y> const& a2) {
    return Formula<Y>(new ScalarFormulaNode<Y>(op,c1,a2),PointerTag()); }
template<class Y> inline Vector<Formula<Y>> Formula<Y>::identity(Nat n) {
    Vector<Formula<Y>> r(n); for(Nat i=0; i!=n; ++i) { r[i]=Formula<Y>::coordinate(i); } return r; }

// DEPRECATED
template<class Y> inline Formula<Y> Formula<Y>::constant(Int c) {
    return Formula<Y>::constant(Y(c)); }


template<class X> using NumericType = typename X::NumericType;
template<class X> using GenericType = typename X::GenericType;
template<class X> using GenericNumericType = GenericType<NumericType<X>>;

template<class X, class Y> inline X make_constant(const Y& c, X r) {
    r=c; return std::move(r); }

inline Real make_constant(const EffectiveNumber& c, const Real& x) {
    return Real(c); }
inline Formula<Real> make_constant(const EffectiveNumber& c, const Formula<Real>& x) {
    return Formula<Real>::constant(Real(c)); }
template<class X, EnableIf<IsSame<X,Real>> =dummy> Algebra<X> make_constant(const EffectiveNumber& c, const Algebra<X>& x) {
    return make_constant(Real(c),x); }

// Make a constant of type Y with value c based on a prototype vector v
template<class X, class Y> inline X make_constant(const Y& c, const Vector<X>& v) {
    return make_constant(c,v.zero_element());
}

template<class X, class Y> X direct_evaluate(const Formula<Y>& f, const Vector<X>& x) {
    switch(f.kind()) {
        case OperatorKind::COORDINATE: return x[f.ind()];
        case OperatorKind::NULLARY: return make_constant(f.val(),x);
        case OperatorKind::UNARY: return compute(f.op(),evaluate(f.arg(),x));
        case OperatorKind::BINARY: return compute(f.op(),evaluate(f.arg1(),x),evaluate(f.arg2(),x));
        case OperatorKind::SCALAR: return compute(f.op(),f.cnst(),evaluate(f.arg(),x));
        case OperatorKind::GRADED: return compute(f.op(),evaluate(f.arg(),x),f.num());
        default: ARIADNE_FAIL_MSG("Cannot evaluate formula "<<f<<" on "<<x<<"; unknown operator "<<f.op()<<" of kind "<<f.kind()<<"\n");
    }
}


//! \brief Convert the expression with index type \c I to one with variables indexed by \a J.
template<class X, class Y> const X& cached_evaluate(const Formula<Y>& f, const Vector<X>& x, Map<const Void*,X>& cache) {
    const FormulaNode<Y>* fptr=f.node_ptr();
    if(cache.has_key(fptr)) { return cache.get(fptr); }
    switch(f.kind()) {
        case OperatorKind::COORDINATE: return insert( cache, fptr, x[f.ind()] );
        case OperatorKind::NULLARY: return insert( cache, fptr, make_constant(f.val(),x) );
        case OperatorKind::UNARY: return insert( cache, fptr, compute(f.op(),cached_evaluate(f.arg(),x,cache)) );
        case OperatorKind::BINARY: return insert( cache, fptr, compute(f.op(),cached_evaluate(f.arg1(),x,cache),cached_evaluate(f.arg2(),x,cache)) );
        case OperatorKind::SCALAR: return insert( cache, fptr, compute(f.op(),f.cnst(),cached_evaluate(f.arg(),x,cache)) );
        case OperatorKind::GRADED: return insert( cache, fptr, compute(f.op(),cached_evaluate(f.arg(),x,cache),f.num()) );
        default: ARIADNE_FAIL_MSG("Cannot evaluate formula "<<f<<" on "<<x<<"; unknown operator "<<f.op()<<" of kind "<<f.kind()<<"\n");
    }
}

template<class X, class Y> inline X cached_evaluate(const Formula<Y>& f, const Vector<X>& v) {
    Map<const Void*,X> cache;
    return cached_evaluate(f,v,cache);
}

template<class X, class Y> Vector<X> cached_evaluate(const Vector<Formula<Y>>& f, const Vector<X>& v) {
    assert(v.size()!=0);
    Vector<X> r(f.size(),zero_element(v));
    Map<const Void*,X> cache;
    for(Nat i=0; i!=r.size(); ++i) {
        r[i]=cached_evaluate(f[i],v,cache);
    }
    return r;
}

template<class X, class Y> X evaluate(const Formula<Y>& f, const Vector<X>& x) {
    return cached_evaluate(f,x);
}
template<class X, class Y> Vector<X> evaluate(const Vector<Formula<Y>>& f, const Vector<X>& x) {
    return cached_evaluate(f,x);
}

//! \ingroup FunctionModule
//! \brief Convert a power-series expansion into a formula using a version of Horner's rule.
//! See J. M. Pena and T. Sauer, "On the multivariate Horner scheme", SIAM J. Numer. Anal. 37(4) 1186-1197, 2000.
//!
template<class X> Formula<X> formula(const Expansion<MultiIndex,X>& e)
{
    Vector<Formula<X>> identity(e.argument_size());
    for(Nat i=0; i!=identity.size(); ++i) { identity[i]=Formula<X>::coordinate(i); }
    return horner_evaluate(e,identity);
}

} // namespace Ariadne


#endif // ARIADNE_FORMULA_HPP
