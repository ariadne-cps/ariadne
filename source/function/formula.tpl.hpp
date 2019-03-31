/***************************************************************************
 *            formula.tpl.hpp
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

namespace Ariadne {


template<class OP, class... AS> using Symbolic = ExpressionTemplate<OP,AS...>;

template<class Y> using ConstantFormulaNode = Symbolic<Cnst,Y>;

template<class Y> class IndexFormulaNode : public Symbolic<Var,Index> { public: using ExpressionTemplate<Var,Index>::ExpressionTemplate; };

template<class Y> using UnaryFormulaNode = ExpressionTemplate<UnaryElementaryOperator,Formula<Y>>;
template<class Y> using BinaryFormulaNode = ExpressionTemplate<BinaryElementaryOperator,Formula<Y>,Formula<Y>>;
template<class Y> using GradedFormulaNode = ExpressionTemplate<GradedElementaryOperator,Formula<Y>,Int>;
template<class Y> using ScalarFormulaNode = ExpressionTemplate<BinaryElementaryOperator,Y,Formula<Y>>;

template<class O, class A, template<class>class E> struct ExpressionTemplate<O,A,E<A>> {
    O _op; A _cnst; E<A> _arg;
    ExpressionTemplate(O o, A c, E<A> a) : _op(o), _cnst(c), _arg(a) { }
    template<class T> operator T() const { return _op(static_cast<T>(_cnst),static_cast<T>(_arg)); }
    template<class... AS> auto operator() (AS... vals) const -> decltype(_op(_cnst,_arg(vals...))) {
        return _op(_cnst(vals...),_arg(vals...)); }
    friend OutputStream& operator<<(OutputStream& os, ExpressionTemplate expr) {
        return os << expr._op.code() << "(" << expr._arg1 << "," << expr._arg2 << ")"; }
};

template<class Y> using FormulaNodeVariantType = Variant<ConstantFormulaNode<Y>,IndexFormulaNode<Y>, UnaryFormulaNode<Y>,BinaryFormulaNode<Y>,GradedFormulaNode<Y>,ScalarFormulaNode<Y>>;

template<class Y>
class FormulaNode
    : public FormulaNodeVariantType<Y>
{
  public:
    template<class FN> FormulaNode(FN&& fn) : FormulaNodeVariantType<Y>(std::forward<FN>(fn)) { }
    template<class VIS> decltype(auto) visit(VIS&& vis) const {
        return std::visit(std::forward<VIS>(vis),static_cast<FormulaNodeVariantType<Y> const&>(*this)); }

};



template<class Y> inline Formula<Y> make_formula(const Y& c) {
    return Formula<Y>::constant(c); }
template<class Y> inline Formula<Y> make_formula(Cnst op, const Y& c) {
    return Formula<Y>::constant(c); }
template<class Y> inline Formula<Y> make_formula(Ind op, SizeType j) {
    return Formula<Y>::index(j); }
template<class Y> inline Formula<Y> make_formula(const UnaryElementaryOperator& op, const Formula<Y>& arg) {
    return Formula<Y>::unary(op,arg); }
template<class Y> inline Formula<Y> make_formula(const BinaryElementaryOperator& op, const Formula<Y>& arg1, const Formula<Y>& arg2) {
    return Formula<Y>::binary(op,arg1,arg2); }
template<class Y> inline Formula<Y> make_formula(const BinaryElementaryOperator& op, const Formula<Y>& arg1, const Y& arg2) {
    return Formula<Y>::binary(op,arg1,make_formula(arg2)); }
template<class Y> inline Formula<Y> make_formula(const GradedElementaryOperator& op, const Formula<Y>& arg, Int num) {
    return Formula<Y>::graded(op,arg,num); }
template<class Y, class OP> inline Formula<Y> make_formula(const OP& op, Y const& cnst, const Formula<Y>& arg) {
    OperatorCode op_code=static_cast<OperatorCode>((char)op.code()+((char)OperatorCode::SADD-(char)OperatorCode::ADD));
    return Formula<Y>::scalar(Operator(op_code,OperatorKind::SCALAR),cnst,arg); }

template<class Y> Formula<Y> make_formula(const Operator& op, const Formula<Y>& arg) {
    return make_formula(UnaryElementaryOperator(op.code()),arg); }
template<class Y> Formula<Y> make_formula(const Operator& op, const Formula<Y>& arg1, const Formula<Y>& arg2) {
    return make_formula(BinaryElementaryOperator(op.code()),arg1,arg2); }
template<class Y> Formula<Y> make_formula(const Operator& op, const Formula<Y>& arg1, Int n2) {
    return make_formula(GradedElementaryOperator(op.code()),arg1,n2); }
template<class Y> Formula<Y> make_formula(const Operator& op, const Y& c1, const Formula<Y>& arg2) {
    return make_formula(BinaryElementaryOperator(op.code()),c1,arg2); }


namespace {
template<class X, class Y> decltype(auto) evaluate(Symbolic<Cnst,Y> c, Vector<X> const& v) {
    return make_constant(c._val,v); }
template<class X, class I> decltype(auto) evaluate(Symbolic<Var,I> i, Vector<X> const& v) {
    return v[i._ind]; }
template<class X, class OP, class A> decltype(auto) evaluate(Symbolic<OP,A> s, Vector<X> const& v) {
    return s._op(evaluate(s._arg,v)); }
template<class X, class OP, class A1, class A2> decltype(auto) evaluate(Symbolic<OP,A1,A2> s, Vector<X> const& v) {
    return s._op(evaluate(s._arg1,v),evaluate(s._arg2,v)); }
template<class X, class OP, class A1, class A2, class A3> decltype(auto) evaluate(Symbolic<OP,A1,A2,A3> s, Vector<X> const& v) {
    return s._op(evaluate(s._arg1,v),evaluate(s._arg2,v),evaluate(s._arg3,v)); }
template<class X, class OP, class A> decltype(auto) evaluate(Symbolic<OP,A,Formula<A>> s, Vector<X> const& v) {
    return s._op(evaluate(Symbolic<Cnst,A>(s._cnst),v),evaluate(s._arg,v)); }
template<class X, class OP, class A> decltype(auto) evaluate(Symbolic<OP,A,Int> s, Vector<X> const& v) {
    return s._op(evaluate(s._arg,v),s._num); }

template<class X, class Y> decltype(auto) _cached_evaluate(Symbolic<Cnst,Y> const& c, Vector<X> const& v, Map<const Void*,X>& cache) {
    return make_constant(c._val,v); }
template<class X, class I> decltype(auto) _cached_evaluate(Symbolic<Var,I> const& i, Vector<X> const& v, Map<const Void*,X>& cache) {
    return v[i._ind]; }
template<class X, class OP, class A> decltype(auto) _cached_evaluate(Symbolic<OP,A> const& s, Vector<X> const& v, Map<const Void*,X>& cache) {
    return s._op(cached_evaluate(s._arg,v)); }
template<class X, class OP, class A> decltype(auto) _cached_evaluate(Symbolic<OP,A,Formula<A>> const& s, Vector<X> const& v, Map<const Void*,X>& cache) {
    return s._op(evaluate(ConstantFormulaNode<A>(s._cnst),v),cached_evaluate(s._arg,v)); }
template<class X, class OP, class A> decltype(auto) _cached_evaluate(Symbolic<OP,Formula<A>,A> const& s, Vector<X> const& v, Map<const Void*,X>& cache) {
    return s._op(s._cnst,cached_evaluate(s._arg,v)); }
template<class X, class OP, class A> decltype(auto) _cached_evaluate(Symbolic<OP,A,Int> const& s, Vector<X> const& v, Map<const Void*,X>& cache) {
    return s._op(cached_evaluate(s._arg,v),s._num); }
template<class X, class OP, class A1, class A2> decltype(auto) _cached_evaluate(Symbolic<OP,A1,A2> const& s, Vector<X> const& v, Map<const Void*,X>& cache) {
    return s._op(cached_evaluate(s._arg1,v),cached_evaluate(s._arg2,v)); }
template<class X, class OP, class A1, class A2, class A3> decltype(auto) _cached_evaluate(Symbolic<OP,A1,A2,A3> const& s, Vector<X> const& v, Map<const Void*,X>& cache) {
    return s._op(cached_evaluate(s._arg1,v),cached_evaluate(s._arg2,v),cached_evaluate(s._arg3,v)); }

template<class X, class OP, class... AS> decltype(auto) _cached_evaluate(Symbolic<OP,AS...> const& s, Vector<X> const& v, Map<const Void*,X>& cache) {
    return v[0]; }

} // namespace


template<class X, class Y> X direct_evaluate(const Formula<Y>& f, const Vector<X>& x) {
    return f.node_ref().visit([&x](auto s){return evaluate(s,x);});
}


//! \brief Convert the expression with index type \c I to one with variables indexed by \a J.
template<class X, class Y> const X& cached_evaluate(const Formula<Y>& f, const Vector<X>& x, Map<const Void*,X>& cache) {
    const FormulaNode<Y>* fptr=f.node_ptr();
    if(cache.has_key(fptr)) { return cache.get(fptr); }
    else { X r=f.node_ref().visit([&x,&cache](auto s){return _cached_evaluate(s,x,cache);}); return insert(cache,fptr,r); }
}

template<class X, class Y> inline X cached_evaluate(const Formula<Y>& f, const Vector<X>& v) {
    Map<const Void*,X> cache;
    return cached_evaluate(f,v,cache);
}

template<class X, class Y> Vector<X> cached_evaluate(const Vector<Formula<Y>>& f, const Vector<X>& v) {
    assert(v.size()!=0);
    Vector<X> r(f.size(),zero_element(v));
    Map<const Void*,X> cache;
    for(SizeType i=0; i!=r.size(); ++i) {
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




namespace {
template<class Y, class J> decltype(auto) derivative(IndexFormulaNode<Y> const& s, J j) { return ConstantFormulaNode<Y>(Y(s._ind==j?1:0)); }
} // namespace



//! \ingroup FunctionModule
//! \brief Convert a power-series expansion into a formula using a version of Horner's rule.
//! See J. M. Pena and T. Sauer, "On the multivariate Horner scheme", SIAM J. Numer. Anal. 37(4) 1186-1197, 2000.
//!
template<class X> Formula<X> formula(const Expansion<MultiIndex,X>& e)
{
    Vector<Formula<X>> identity(e.argument_size());
    for(SizeType i=0; i!=identity.size(); ++i) { identity[i]=Formula<X>::coordinate(i); }
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
        case OperatorKind::COORDINATE: case OperatorKind::VARIABLE: return _substitute_coordinate(a.ind(),i,a,is);
        default: ARIADNE_FAIL_MSG("Cannot substitute "<<is<<" for index "<<i<<" in an unknown formula "<<a<<"\n");
    }
}



template<class Y> inline Formula<Y> simplify(const Formula<Y>& a) {
    if(a.kind() == OperatorKind::UNARY) {
        Formula<Y> sarg=simplify(a.arg());
        if(sarg.op()==OperatorCode::CNST) {
            UnaryElementaryOperator uop(a.op());
            return Formula<Y>(uop(sarg.val()));
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
        case OperatorKind::COORDINATE: case OperatorKind::VARIABLE:
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

template<class Y> Vector<Formula<Y>> simplify(const Vector<Formula<Y>>& a) {
    return Vector<Formula<Y>>(a.size(),[&a](SizeType i){return simplify(a[i]);});
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
        case OperatorKind::COORDINATE: case OperatorKind::VARIABLE: return not is.contains(a.ind());
        case OperatorKind::NULLARY: return true;
        case OperatorKind::UNARY: case OperatorKind::SCALAR: case OperatorKind::GRADED: return is_constant_in(a.arg(),is);
        case OperatorKind::BINARY: return is_constant_in(a.arg1(),is) and is_constant_in(a.arg2(),is);
        default: ARIADNE_FAIL_MSG("Cannot evaluate if formula "<<a<<" is constant in "<<is<<"\n");
    }
}


template<class Y> Bool is_affine_in(const Formula<Y>& a, const Set<Nat>& is) {
    switch(a.op()) {
        case OperatorCode::CNST: return true;
        case OperatorCode::IND: case OperatorCode::VAR: return true;
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
