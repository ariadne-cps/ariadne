/***************************************************************************
 *            expression.cpp
 *
 *  Copyright 2009--17  Pieter Collins
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

#include "../utility/standard.hpp"
#include "../config.hpp"


#include "../algebra/algebra.hpp"
#include "../algebra/algebra_wrapper.hpp"

#include "../function/formula.hpp"
#include "../function/formula.tpl.hpp"

#include "../symbolic/constant.hpp"
#include "../symbolic/variables.hpp"
#include "../symbolic/expression.hpp"
#include "../symbolic/assignment.hpp"
#include "../symbolic/space.hpp"
#include "../symbolic/valuation.hpp"

#include "expression.tpl.hpp"

namespace Ariadne {

template<> template<>
Expression<Real>::operator ElementaryAlgebra<Real>() const {
    return ElementaryAlgebra<Real>(new ElementaryAlgebraWrapper<Expression<Real>,Real>(*this));
}


template class Expression<Boolean>;
template class Expression<Kleenean>;
template class Expression<String>;
template class Expression<Integer>;
template class Expression<Real>;

template Bool before<Real>(Expression<Real> const& e1, Expression<Real> const& e2);
template Nat count_nodes<Real>(const Expression<Real>& e);
template Nat count_distinct_nodes<Real>(const Expression<Real>& e);
template Nat count_distinct_node_ptrs<Real>(const Expression<Real>& e);

template Expression<Real> ElementaryAlgebra<Real>::extract() const;

Expression<Boolean> operator&&(Expression<Boolean> const& e1, Expression<Boolean> const& e2) {
    return make_expression<Boolean>(AndOp(),e1,e2); }
Expression<Boolean> operator||(Expression<Boolean> const& e1, Expression<Boolean> const& e2) {
    return make_expression<Boolean>(OrOp(),e1,e2); }
Expression<Boolean> operator!(Expression<Boolean> const& e) {
    return make_expression<Boolean>(NotOp(),e); }


Expression<Kleenean> operator&&(Expression<Kleenean> const& e1, Expression<Kleenean> const& e2) {
    return make_expression<Kleenean>(AndOp(),e1,e2); }
Expression<Kleenean> operator||(Expression<Kleenean> const& e1, Expression<Kleenean> const& e2) {
    return make_expression<Kleenean>(OrOp(),e1,e2); }
Expression<Kleenean> operator!(Expression<Kleenean> const& e) {
    return make_expression<Kleenean>(NotOp(),e); }


Expression<Boolean> operator==(Variable<String> v1, const String& s2) {
    return make_expression<Boolean>(Equal(),Expression<String>(v1),Expression<String>::constant(s2)); }
Expression<Boolean> operator!=(Variable<String> v1, const String& s2) {
    return make_expression<Boolean>(Unequal(),Expression<String>(v1),Expression<String>::constant(s2)); }


Expression<Boolean> operator==(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Boolean>(Equal(),e1,e2); }
Expression<Boolean> operator!=(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Boolean>(Unequal(),e1,e2); }
Expression<Boolean> operator>=(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Boolean>(Geq(),e1,e2); }
Expression<Boolean> operator<=(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Boolean>(Leq(),e1,e2); }
Expression<Boolean> operator>(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Boolean>(Gtr(),e1,e2); }
Expression<Boolean> operator<(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Boolean>(Less(),e1,e2); }



Expression<Integer> operator+(Expression<Integer> const& e) {
    return make_expression<Integer>(Pos(),e); }
Expression<Integer> operator-(Expression<Integer> const& e) {
    return make_expression<Integer>(Neg(),e); }
Expression<Integer> operator+(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Integer>(Add(),e1,e2); }
Expression<Integer> operator-(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Integer>(Sub(),e1,e2); }
Expression<Integer> operator*(Expression<Integer> const& e1, Expression<Integer> const& e2) {
    return make_expression<Integer>(Mul(),e1,e2); }

Expression<Integer>& operator+=(Expression<Integer>& e1, Expression<Integer> const& e2) {
    return e1=e1+e2; }
Expression<Integer>& operator-=(Expression<Integer>& e1, Expression<Integer> const& e2) {
    return e1=e1-e2; }
Expression<Integer>& operator*=(Expression<Integer>& e1, Expression<Integer> const& e2) {
    return e1=e1*e2; }

Expression<Kleenean> sgn(Expression<Real> const& e) {
    return make_expression<Kleenean>(Sgn(),e); }

Expression<Kleenean> operator==(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Kleenean>(Equal(),e1,e2); }
Expression<Kleenean> operator!=(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Kleenean>(Unequal(),e1,e2); }
Expression<Kleenean> operator>=(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Kleenean>(Geq(),e1,e2); }
Expression<Kleenean> operator<=(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Kleenean>(Leq(),e1,e2); }
Expression<Kleenean> operator>(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Kleenean>(Gtr(),e1,e2); }
Expression<Kleenean> operator<(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Kleenean>(Less(),e1,e2); }


Expression<Real> operator+(Expression<Real> const& e) {
    return make_expression<Real>(Pos(),e); }
Expression<Real> operator-(Expression<Real> const& e) {
    return make_expression<Real>(Neg(),e); }
Expression<Real> operator+(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Add(),e1,e2); }
Expression<Real> operator-(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Sub(),e1,e2); }
Expression<Real> operator*(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Mul(),e1,e2); }
Expression<Real> operator/(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Div(),e1,e2); }

Expression<Real>& operator+=(Expression<Real>& e1, Expression<Real> const& e2) {
    return e1=e1+e2; }
Expression<Real>& operator-=(Expression<Real>& e1, Expression<Real> const& e2) {
    return e1=e1-e2; }
Expression<Real>& operator*=(Expression<Real>& e1, Expression<Real> const& e2) {
    return e1=e1*e2; }
Expression<Real>& operator/=(Expression<Real>& e1, Expression<Real> const& e2) {
    return e1=e1/e2; }

Expression<Real> add(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Add(),e1,e2); }
Expression<Real> sub(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Sub(),e1,e2); }
Expression<Real> mul(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Mul(),e1,e2); }
Expression<Real> div(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Div(),e1,e2); }
Expression<Real> pow(Expression<Real> const& e, Int n) {
    return make_expression<Real>(Pow(),e,n); }

Expression<Real> nul(Expression<Real> const& e) {
    return make_expression<Real>(Real(0)); }
Expression<Real> pos(Expression<Real> const& e) {
    return make_expression<Real>(Pos(),e); }
Expression<Real> neg(Expression<Real> const& e) {
    return make_expression<Real>(Neg(),e); }
Expression<Real> rec(Expression<Real> const& e) {
    return make_expression<Real>(Rec(),e); }
Expression<Real> sqr(Expression<Real> const& e) {
    return make_expression<Real>(Sqr(),e); }
Expression<Real> hlf(Expression<Real> const& e) {
    return make_expression<Real>(Hlf(),e); }
Expression<Real> sqrt(Expression<Real> const& e) {
    return make_expression<Real>(Sqrt(),e); }
Expression<Real> exp(Expression<Real> const& e) {
    return make_expression<Real>(Exp(),e); }
Expression<Real> log(Expression<Real> const& e) {
    return make_expression<Real>(Log(),e); }
Expression<Real> sin(Expression<Real> const& e) {
    return make_expression<Real>(Sin(),e); }
Expression<Real> cos(Expression<Real> const& e) {
    return make_expression<Real>(Cos(),e); }
Expression<Real> tan(Expression<Real> const& e) {
    return make_expression<Real>(Tan(),e); }
Expression<Real> asin(Expression<Real> const& e) {
    return make_expression<Real>(Asin(),e); }
Expression<Real> acos(Expression<Real> const& e) {
    return make_expression<Real>(Acos(),e); }
Expression<Real> atan(Expression<Real> const& e) {
    return make_expression<Real>(Atan(),e); }

Expression<Real> max(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Max(),e1,e2); }
Expression<Real> min(Expression<Real> const& e1, Expression<Real> const& e2) {
    return make_expression<Real>(Min(),e1,e2); }
Expression<Real> abs(Expression<Real> const& e) {
    return make_expression<Real>(Abs(),e); }

template String evaluate(const Expression<String>& e, const Valuation<String>& x);
template Integer evaluate(const Expression<Integer>& e, const Valuation<Integer>& x);
template Real evaluate(const Expression<Real>& e, const Valuation<Real>& x);
template Boolean evaluate(const Expression<Boolean>& e, const Valuation<String>& x);
template Boolean evaluate(const Expression<Boolean>& e, const Valuation<Integer>& x);
template Kleenean evaluate(const Expression<Kleenean>& e, const Valuation<Real>& x);

template Real evaluate(Expression<Real> const&, Map<Identifier, Real> const&);


template Set<Identifier> arguments(const Expression<Boolean>& e);
template Set<Identifier> arguments(const Expression<Kleenean>& e);
template Set<Identifier> arguments(const Expression<Real>& e);


template Expression<Kleenean> substitute(const Expression<Kleenean>& e, const Variable<Kleenean>& v, const Kleenean& c);
template Expression<Kleenean> substitute(const Expression<Kleenean>& e, const Variable<Real>& v, const Real& c);
template Expression<Real> substitute(const Expression<Real>& e, const Variable<Real>& v, const Real& c);
template Expression<Real> substitute(const Expression<Real>& e, const Variable<Real>& v, const Expression<Real>& c);
template Expression<Real> substitute(const Expression<Real>& e, const List< Assignment< Variable<Real>, Expression<Real> > >& c);
template Vector<Expression<Real>> substitute(const Vector<Expression<Real>>& e, const List< Assignment< Variable<Real>, Expression<Real> > >& c);
template Expression<Kleenean> substitute(const Expression<Kleenean>& e, const List< Assignment< Variable<Real>, Expression<Real> > >& c);

template Expression<Real> simplify(const Expression<Real>& e);
template Expression<Kleenean> simplify(const Expression<Kleenean>& e);

template Void eliminate_common_subexpressions(Expression<Real>&);
template Void eliminate_common_subexpressions(Vector<Expression<Real>>&);




Expression<Real> indicator(Expression<Kleenean> e, Sign sign) {
    switch(e.op()) {
        case OperatorCode::CNST: {
            Kleenean value=( sign==Sign::POSITIVE ? e.val() : !e.val() );
            ValidatedKleenean checked_value = value.check(Effort::get_default());
            if(definitely(checked_value)) { return Expression<Real>::constant(+1); }
            else if(not possibly(checked_value)) {  return Expression<Real>::constant(-1); }
            else { return Expression<Real>::constant(0); }
        }
        case OperatorCode::VAR:
            return Expression<Real>(Variable<Real>(e.var()));
        case OperatorCode::GEQ: case OperatorCode::GT:
            if(sign==Sign::POSITIVE) { return e.cmp1<Real>()-e.cmp2<Real>(); }
            else { return e.cmp2<Real>()-e.cmp1<Real>(); }
        case OperatorCode::LEQ: case OperatorCode::LT:
            if(sign==Sign::POSITIVE) { return e.cmp2<Real>()-e.cmp1<Real>(); }
            else { return e.cmp1<Real>()-e.cmp2<Real>(); }
        case OperatorCode::AND:
            return min(indicator(e.arg1(),sign),indicator(e.arg2(),sign));
        case OperatorCode::OR:
            return max(indicator(e.arg1(),sign),indicator(e.arg2(),sign));
        case OperatorCode::NOT:
            return neg(indicator(e.arg(),sign));
        default:
            ARIADNE_FAIL_MSG("Cannot compute indicator function of expression " << e);
    }
}



template Bool is_constant(const Expression<Real>&, const Real&);
template Bool is_constant(const Expression<Kleenean>&, const Kleenean&);

template Bool is_variable(const Expression<Real>&, const Variable<Real>&);
template Bool identical(const Expression<Real>&, const Expression<Real>&);

template Bool is_constant_in<Real>(const Expression<Real>& e, const Set<Variable<Real>>& vs);

Bool is_affine_in(const Expression<Real>& e, const Set<Variable<Real>>& vs) {
    switch(e.op()) {
        case OperatorCode::CNST: return true;
        case OperatorCode::VAR: return true;
        case OperatorCode::ADD: case OperatorCode::SUB: return is_affine_in(e.arg1(),vs) and is_affine_in(e.arg2(),vs);
        case OperatorCode::MUL: return (is_affine_in(e.arg1(),vs) and is_constant_in(e.arg2(),vs)) or (is_constant_in(e.arg1(),vs) and is_affine_in(e.arg2(),vs));
        case OperatorCode::DIV: return (is_affine_in(e.arg1(),vs) and is_constant_in(e.arg2(),vs));
        case OperatorCode::POS: case OperatorCode::NEG: return is_affine_in(e.arg(),vs);
        case OperatorCode::POW: case OperatorCode::SQR: case OperatorCode::COS: case OperatorCode::SIN: case OperatorCode::TAN: return is_constant_in(e.arg(),vs);
        default: ARIADNE_FAIL_MSG("Not currently supporting code '"<<e.op()<<"' for evaluation of affinity in given variables\n");
    }
}

Bool is_affine_in(const Vector<Expression<Real>>& e, const Set<Variable<Real>>& vs) {
    for (auto i : range(e.size()))
        if (not is_affine_in(e[i],vs)) return false;
    return true;
}

Bool is_additive_in(const Vector<Expression<Real>>& ev, const Set<Variable<Real>>& vs) {
    // We treat the vector of expressions as additive in vs if each variable in vs appears at most once in all expressions,
    // with a constant value of 1
    // (FIXME: this simplifies the case of a constant multiplier, for which would need to rescale the variable)
    // (FIXME: more generally, this simplifies the case of a diagonalisable matrix of constant multipliers)

    auto one = Expression<Real>::constant(1);
    auto zero = Expression<Real>::constant(0);

    for (auto v : vs) {
        Bool already_found = false;
        Bool already_found_one = false;
        for (auto i : range(ev.size())) {
            const Expression<Real>& e = ev[i];
            auto der = simplify(derivative(e, v));
            if (not identical(der,zero)) {
                if (already_found) {
                    return false;
                } else {
                    already_found = true;
                    if (identical(der,one)) {
                        if (already_found_one) {
                            return false;
                        } else {
                            already_found_one = true;
                        }
                    } else {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}


Bool opposite(Expression<Kleenean> e1, Expression<Kleenean> e2) {

    OperatorCode e1op;
    OperatorCode e2op;
    switch(e1.op()) {
        case OperatorCode::GEQ: case OperatorCode::GT: e1op=OperatorCode::GEQ; break;
        case OperatorCode::LEQ: case OperatorCode::LT: e1op=OperatorCode::LEQ; break;
        default: return false;
    }
    switch(e2.op()) {
        case OperatorCode::GEQ: case OperatorCode::GT: e2op=OperatorCode::GEQ; break;
        case OperatorCode::LEQ: case OperatorCode::LT: e2op=OperatorCode::LEQ; break;
        default: return false;
    }

    // Both expressions are <=,<,>=,> comparisons
    Expression<Real> const& e1arg1=e1.cmp1<Real>();
    Expression<Real> const& e1arg2=e2.cmp2<Real>();
    Expression<Real> const& e2arg1=e1.cmp1<Real>();
    Expression<Real> const& e2arg2=e2.cmp2<Real>();

    // Test if the expressions are of the form a1<=a2; a1>=a2 or a1<=a2; a2<=a1
    if(e1op==e2op) {
        if(identical(e1arg1,e2arg2) && identical(e1arg2,e2arg1)) { return true; }
        else { return false; }
    } else {
        if(identical(e1arg1,e2arg1) && identical(e1arg2,e2arg2)) { return true; }
        else { return false; }
    }

}


Expression<Real> derivative(const Expression<Real>& e, Variable<Real> v)
{
    /*
    switch(e.op()) {
        case OperatorCode::CNST:
            return Expression<Real>::constant(0);
        case OperatorCode::VAR:
            if(e.var()==v.name()) { return Expression<Real>::constant(1); }
            else { return Expression<Real>::constant(0); }
        case OperatorCode::ADD:
            return derivative(e.arg1(),v)+derivative(e.arg2(),v);
        case OperatorCode::SUB:
            return derivative(e.arg1(),v)-derivative(e.arg2(),v);
        case OperatorCode::MUL:
            return e.arg1()*derivative(e.arg2(),v)+derivative(e.arg1(),v)*e.arg2();
        case OperatorCode::DIV:
            return derivative(e.arg1() * rec(e.arg2()),v);
        case OperatorCode::NEG:
            return  - derivative(e.arg(),v);
        case OperatorCode::REC:
            return  - derivative(e.arg(),v) * rec(sqr(e.arg()));
        case OperatorCode::SQR:
            return static_cast<Real>(2) * derivative(e.arg(),v) * e.arg();
        case OperatorCode::POW:
            return Expression<Real>::constant(e.num())*make_expression<Real>(OperatorCode::POW,e.arg(),e.num()-1)*derivative(e.arg(),v);
        case OperatorCode::EXP:
            return derivative(e.arg(),v) * e.arg();
        case OperatorCode::LOG:
            return derivative(e.arg(),v) * rec(e.arg());
        case OperatorCode::SIN:
            return derivative(e.arg(),v) * cos(e.arg());
        case OperatorCode::COS:
            return -derivative(e.arg(),v) * sin(e.arg());
        case OperatorCode::TAN:
            return derivative(e.arg(),v) * (static_cast<Real>(1)-sqr(e.arg()));*/
    switch(e.kind()) {
        case OperatorKind::NULLARY: return Expression<Real>::constant(0);
        case OperatorKind::VARIABLE: return Expression<Real>::constant(e.var()==v.name()?1:0);
        case OperatorKind::UNARY: return compute_derivative(e.op().code(),e.arg(),derivative(e.arg(),v));
        case OperatorKind::BINARY: return compute_derivative(e.op().code(),e.arg1(),derivative(e.arg1(),v),e.arg2(),derivative(e.arg2(),v));
        case OperatorKind::GRADED: return compute_derivative(e.op().code(), e.arg(), derivative(e.arg(),v), e.num());
        default:
            ARIADNE_THROW(std::runtime_error,"derivative(Expression<Real> e, Variable<Real> v)",
                          "Cannot compute derivative of "<<e<<"\n");
    }
}



const Formula<EffectiveNumber>& cached_make_formula(const Expression<Real>& e, const Map<Identifier,Nat>& v, Map< const Void*, Formula<EffectiveNumber> >& cache)
{
    typedef EffectiveNumber Y;
    const ExpressionNode<Real>* eptr=e.node_ptr().operator->();
    if(cache.has_key(eptr)) { return cache.get(eptr); }
    switch(e.kind()) {
        case OperatorKind::VARIABLE: return insert( cache, eptr, make_formula<Y>(v[e.var()]) );
        case OperatorKind::NULLARY: return insert( cache, eptr, make_formula<Y>(e.val()) );
        case OperatorKind::UNARY: return insert( cache, eptr, make_formula<Y>(e.op(),cached_make_formula(e.arg(),v,cache)));
        case OperatorKind::BINARY: return insert( cache, eptr, make_formula<Y>(e.op(),cached_make_formula(e.arg1(),v,cache),cached_make_formula(e.arg2(),v,cache)) );
        case OperatorKind::GRADED: return insert( cache, eptr, make_formula<Y>(e.op(),cached_make_formula(e.arg(),v,cache),e.num()) );
        default: ARIADNE_FAIL_MSG("Cannot convert expression "<<e<<" to use variables "<<v<<"\n");
    }
}

Formula<EffectiveNumber> make_formula(const Expression<Real>& e, const Map<Identifier,Nat>& v)
{
    Map< const Void*, Formula<EffectiveNumber> > cache;
    return cached_make_formula(e,v,cache);
}

Formula<EffectiveNumber> make_formula(const Expression<Real>& e, const Space<Real>& spc)
{
    typedef EffectiveNumber Y;
    switch(e.kind()) {
        case OperatorKind::GRADED: return make_formula<Y>(e.op(),make_formula(e.arg(),spc),e.num());
        case OperatorKind::BINARY: return make_formula<Y>(e.op(),make_formula(e.arg1(),spc),make_formula(e.arg2(),spc));
        case OperatorKind::UNARY: return make_formula<Y>(e.op(),make_formula(e.arg(),spc));
        case OperatorKind::NULLARY: return Formula<Y>::constant(e.val());
        case OperatorKind::VARIABLE: return Formula<Y>::coordinate(spc.index(e.var()));
        default: ARIADNE_FAIL_MSG("Cannot compute formula for expression "<<e.op()<<"of kind "<<e.kind()<<" in space "<<spc);
    }
}

Formula<EffectiveNumber> make_formula(const Expression<Real>& e, const Variable<Real>& var)
{
    return make_formula(e,Space<Real>({var}));
}

Vector<Formula<EffectiveNumber>> make_formula(const Vector<Expression<Real>>& e, const Space<Real>& spc)
{
    Vector<Formula<EffectiveNumber>> res(e.size());
    for(SizeType i=0; i!=e.size(); ++i) {
        res[i]=make_formula(e[i],spc);
    }
    return res;
}

Formula<EffectiveNumber> make_formula(const Expression<Real>& e, const List<Variable<Real>>& vars)
{
    return make_formula(e,Space<Real>(vars));
}

Formula<EffectiveNumber> make_formula(const Expression<Real>& out, const List<Assignment<Variable<Real>,Expression<Real>>>& aux, const Space<Real> spc)
{
    return make_formula(substitute(out,aux),spc);
}

Vector<Formula<EffectiveNumber>> make_formula(const Vector<Expression<Real>>& out, const List<Assignment<Variable<Real>,Expression<Real>>>& aux, const Space<Real> spc)
{
    Vector<Formula<EffectiveNumber>> res(out.size());
    for(SizeType i=0; i!=out.size(); ++i) {
        res[i]=make_formula(out[i],aux,spc);
    }
    return res;
}


Expression<Real> make_expression(const Formula<Real>& f, const Space<Real>& s) {
    const List<RealVariable>& vars=s.variables();
    typedef ElementaryAlgebra<Real> RealElementaryAlgebra;
    RealElementaryAlgebra az(RealExpression::constant(0));
    Vector<RealElementaryAlgebra> va(vars.size(),az);
    for(SizeType i=0; i!=va.size(); ++i) { va[i]=RealElementaryAlgebra(RealExpression(vars[i])); }
    RealElementaryAlgebra fa=evaluate(f,va);
    return fa.template extract<RealExpression>();
}

Formula<Real> make_formula(const EffectiveScalarMultivariateFunction& f);
Expression<Real> make_expression(const ScalarMultivariateFunction<EffectiveTag>& f, const Space<Real>& s) {
    return make_expression(make_formula(f),s); }

} // namespace Ariadne
