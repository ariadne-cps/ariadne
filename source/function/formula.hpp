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
typedef Formula<ExactNumber> ExactFormula;

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

//! \brief Substitute all occurrences of index \a i with constant value \a c.
template<class Y> Formula<Y> substitute(const Formula<Y>& a, const Nat& i, const Y& c);
//! \brief Substitute all occurrences of indices with formulae \a as into \a a.
template<class Y> Formula<Y> substitute(const Formula<Y>& a, const List<Pair<Nat,Formula<Y>>>& as);
//! \brief Substitute all occurrences of indices with formulae \a as into \a av.
template<class Y> Vector<Formula<Y>> substitute(const Vector<Formula<Y>>& av, const List<Pair<Nat,Formula<Y>>>& as);

//! \brief Simplify the formula to reduce the number of nodes.
template<class Y> Formula<Y> simplify(const Formula<Y>& a);
//! \brief Simplify the formulae to reduce the number of nodes.
template<class Y> Vector<Formula<Y>> simplify(const Vector<Formula<Y>>& a);

//! \brief Tests whether two formulas are identical.
template<class Y> Bool identical(const Formula<Y>& a1, const Formula<Y>& a2);
//! \brief Tests whether two vector formulas are identical.
template<class Y> Bool identical(const Vector<Formula<Y>>& a1, const Vector<Formula<Y>>& a2);

//! \brief Returns \a true if the formula \a a is syntactically constant in the indices \a is.
template<class Y> Bool is_constant_in(const Formula<Y>& a, const Set<Nat>& is);
//! \brief Returns \a true if the formula \a a is syntactically affine in the indices \a is.
template<class Y> Bool is_affine_in(const Formula<Y>& a, const Set<Nat>& is);
//! \brief Returns \a true if the vector formula \a e is syntactically affine in the indices \a is.
template<class Y> Bool is_affine_in(const Vector<Formula<Y>>& a, const Set<Nat>& is);
//! \brief Returns \a true if the vector formula \a a is syntactically additive (possibly with multipliers) in the indices \a is.
template<class Y> Bool is_additive_in(const Vector<Formula<Y>>& a, const Set<Nat>& is);

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

namespace {
template<class X, class Y> inline const Formula<Y>& _substitute_coordinate(const Nat& ie, const Nat& is, const Formula<Y>& a, const Formula<X>& s) {
    ARIADNE_ASSERT_MSG(ie!=is,"Cannot substitute formula "<<s<<" for coordinate "<<ie<<"\n");
    return a; }
template<class Y> inline const Formula<Y>& _substitute_coordinate(const Nat& ie, const Nat& is, const Formula<Y>& e, const Formula<Y>& s) {
    return ie==is ? s : e; }
} // namespace

template<class Y> Formula<Y> substitute(const Formula<Y>& a, const Nat& i, const Formula<Y>& is) {
    switch(a.kind()) {
        case OperatorKind::BINARY: return make_formula<Y>(a.op(),substitute(a.arg1(),i,is),substitute(a.arg2(),i,is));
        case OperatorKind::UNARY: return make_formula<Y>(a.op(),substitute(a.arg(),i,is));
        case OperatorKind::GRADED: return make_formula<Y>(a.op(),substitute(a.arg(),i,is),a.num());
        case OperatorKind::NULLARY: return make_formula<Y>(a.val());
        case OperatorKind::COORDINATE: return _substitute_coordinate(a.ind(),i,a,is);
        default: ARIADNE_FAIL_MSG("Cannot substitute "<<is<<" for index "<<i<<" in an unknown formula "<<a<<"\n");
    }
}

template<class Y> Formula<Y> substitute(const Formula<Y>& a, const Nat& i, const Y& c) {
    return substitute(a,i,Formula<Y>::constant(c));
}

template<class Y> Formula<Y> substitute(const Formula<Y>& a, const Pair<Nat,Formula<Y>>& as) {
    return substitute(a,as.first,as.second);
}

template<class Y> Formula<Y> substitute(const Formula<Y>& a, const List<Pair<Nat,Formula<Y>>>& as) {
    Formula<Y> r=a;
    for(SizeType i=0; i!=as.size(); ++i) {
        r=substitute(r,as[i]);
    }
    return r;
}

template<class Y> Vector<Formula<Y>> substitute(const Vector<Formula<Y>>& av, const List<Pair<Nat,Formula<Y>>>& as) {
    Vector<Formula<Y>> r(av.size());
    for(SizeType i=0; i!=av.size(); ++i) {
        r[i]=substitute(av[i],as);
    }
    return r;
}

template<class Y> inline Formula<Y> simplify(const Formula<Y>& a) {

    if(a.kind() == OperatorKind::UNARY) {
        Formula<Y> sarg=simplify(a.arg());
        if(sarg.op()==OperatorCode::CNST) {
            return Formula<Y>(compute(a.op(),sarg.val()));
        } else {
            return make_formula<Y>(a.op(),sarg);
        }
    }

    if(a.kind() == OperatorKind::GRADED) {
        Formula<Y> sarg=simplify(a.arg());
        Formula<Y> one(static_cast<Y>(1));
        switch(a.op()) {
            case OperatorCode::POW:
                switch (a.num()) {
                case 0: return one;
                case 1: return sarg;
                default: return make_formula<Y>(OperatorCode::POW,sarg,a.num());
                }
            default:
                return make_formula<Y>(a.op(),sarg,a.num());
        }
    }

    if(a.kind() != OperatorKind::BINARY) { return a; }

    Formula<Y> sarg1=simplify(a.arg1());
    Formula<Y> sarg2=simplify(a.arg2());
    Formula<Y> zero(Formula<Y>::zero());
    Formula<Y> one(static_cast<Y>(1));
    switch(a.op()) {
        case OperatorCode::ADD:
            if(identical(sarg2,zero)) { return sarg1; }
            if(identical(sarg1,zero)) { return sarg2; }
            break;
        case OperatorCode::SUB:
            if(identical(sarg2,zero)) { return sarg1; }
            if(identical(sarg1,zero)) { return -sarg2; }
            break;
        case OperatorCode::MUL:
            if(identical(sarg1,zero)) { return zero; }
            if(identical(sarg2,zero)) { return zero; }
            if(identical(sarg1,one)) { return sarg2; }
            if(identical(sarg2,one)) { return sarg1; }
            break;
        case OperatorCode::DIV:
            if(identical(sarg1,zero)) { return sarg1; }
            if(identical(sarg1,one)) { return rec(sarg2); }
            if(identical(sarg2,one)) { return sarg1; }
        default:
            break;
    }
    return make_formula<Y>(a.op(),sarg1,sarg2);
}

template<class Y> inline Vector<Formula<Y>> simplify(const Vector<Formula<Y>>& a) {

    Vector<Formula<Y>> r(a.size());
    for(SizeType i=0; i!=a.size(); ++i) {
        r[i]=simplify(a[i]);
    }
    return r;
}

inline Bool same(EffectiveNumber const& v1, EffectiveNumber const& v2) {
    // FIXME: Use symbolic approach
    DoublePrecision pr;
    FloatDPBounds x1(v1,pr);
    FloatDPBounds x2(v2,pr);
    return x1.lower_raw()==x2.upper_raw() && x1.upper_raw() == x2.lower_raw();
}

template<class Y> Bool identical(const Formula<Y>& a1, const Formula<Y>& a2)
{
    if(a1.node_ptr()==a2.node_ptr()) { return true; }
    if(a1.op()!=a2.op()) { return false; }
    switch(a1.kind()) {
        case OperatorKind::COORDINATE:
            return a1.ind() == a2.ind();
        case OperatorKind::NULLARY:
            return same(a1.val(),a2.val());
        case OperatorKind::UNARY:
            return identical(a1.arg(),a2.arg());
        case OperatorKind::GRADED:
            return identical(a1.arg(),a2.arg()) && a1.num() == a2.num();
        case OperatorKind::BINARY:
            switch(a1.op()) {
            case OperatorCode::MUL: case OperatorCode::ADD:
                return (identical(a1.arg1(),a2.arg1()) && identical(a1.arg2(),a2.arg2())) ||
                       (identical(a1.arg1(),a2.arg2()) && identical(a1.arg2(),a2.arg1()));
            default:
                return identical(a1.arg1(),a2.arg1()) && identical(a1.arg2(),a2.arg2());
            }
        default:
            return false;
    }
}

template<class Y> Bool identical(const Vector<Formula<Y>>& a1, const Vector<Formula<Y>>& a2) {
    if (a1.size() != a2.size()) return false;

    for (auto i : range(a1.size())) {
        if (not identical(a1[i],a2[i])) return false;
    }
    return true;
}

template<class Y> Bool is_constant_in(const Formula<Y>& a, const Set<Nat>& is) {
    switch(a.kind()) {
        case OperatorKind::COORDINATE: return not is.contains(a.ind());
        case OperatorKind::NULLARY: return true;
        case OperatorKind::UNARY: case OperatorKind::SCALAR: case OperatorKind::GRADED: return is_constant_in(a.arg(),is);
        case OperatorKind::BINARY: return is_constant_in(a.arg1(),is) and is_constant_in(a.arg2(),is);
        default: ARIADNE_FAIL_MSG("Cannot evaluate if formula "<<a<<" is constant in "<<is<<"\n");
    }
}


template<class Y> Bool is_affine_in(const Formula<Y>& a, const Set<Nat>& is) {
    switch(a.op()) {
        case OperatorCode::CNST: return true;
        case OperatorCode::IND: return true;
        case OperatorCode::ADD: case OperatorCode::SUB: return is_affine_in(a.arg1(),is) and is_affine_in(a.arg2(),is);
        case OperatorCode::MUL: return (is_affine_in(a.arg1(),is) and is_constant_in(a.arg2(),is)) or (is_constant_in(a.arg1(),is) and is_affine_in(a.arg2(),is));
        case OperatorCode::DIV: return (is_affine_in(a.arg1(),is) and is_constant_in(a.arg2(),is));
        case OperatorCode::POS: case OperatorCode::NEG: return is_affine_in(a.arg(),is);
        case OperatorCode::POW: case OperatorCode::SQR: case OperatorCode::COS: case OperatorCode::SIN: case OperatorCode::TAN: return is_constant_in(a.arg(),is);
        default: ARIADNE_FAIL_MSG("Not currently supporting code '"<<a.op()<<"' for evaluation of affinity in given indices\n");
    }
}

template<class Y> Bool is_affine_in(const Vector<Formula<Y>>& as, const Set<Nat>& is) {
    for (auto idx : range(as.size()))
        if (not is_affine_in(as[idx],is)) return false;
    return true;
}


template<class Y> Bool is_additive_in(const Vector<Formula<Y>>& as, const Set<Nat>& is) {
    // We treat the vector of formulas as additive in is if each variable in is appears at most once in all expressions,
    // with a constant multiplier
    // (FIXME: this simplifies the case of a diagonalisable matrix of constant multipliers)

    for (auto i : is) {
        Bool already_found = false;
        for (auto idx : range(as.size())) {
            const Formula<Y>& a = as[idx];
            auto der = simplify(derivative(a, i));
            if (not identical(der,Formula<Y>::zero())) {
                if (already_found) {
                    return false;
                } else {
                    already_found = true;
                    if (der.op() != OperatorCode::CNST) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}


} // namespace Ariadne


#endif // ARIADNE_FORMULA_HPP
