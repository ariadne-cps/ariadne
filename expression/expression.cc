/***************************************************************************
 *            expression.cc
 *
 *  Copyright 2009  Pieter Collins
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
 *  MERCHANTABILITY or FITNESS FOperatorCode::OR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "utility/standard.h"
#include "config.h"

#include "expression/expression.h"
#include "expression/assignment.h"
#include "expression/space.h"
#include "expression/valuation.h"

#include "expression/formula.h"
#include "expression/operators.tcc"

namespace Ariadne {

inline Real operator*(const Real& x1, Int x2) { return x1*Real(x2); }
inline Real operator/(Int x1, const Real& x2) { return Real(x1)/x2; }




template class Expression<Boolean>;
template class Expression<Tribool>;
template class Expression<String>;
template class Expression<Integer>;
template class Expression<Real>;

const Formula<Real>& cached_formula(const Expression<Real>& e, const Map<Identifier,Nat>& v, Map< const Void*, Formula<Real> >& cache)
{
    const ExpressionNode<Real>* eptr=e.node_ptr();
    if(cache.has_key(eptr)) { return cache.get(eptr); }
    switch(e.kind()) {
        case OperatorKind::VARIABLE: return insert( cache, eptr, make_formula<Real>(v[e.var()]) );
        case OperatorKind::NULLARY: return insert( cache, eptr, make_formula<Real>(e.val()) );
        case OperatorKind::UNARY: return insert( cache, eptr, make_formula<Real>(e.op(),cached_formula(e.arg(),v,cache)));
        case OperatorKind::BINARY: return insert( cache, eptr, make_formula<Real>(e.op(),cached_formula(e.arg1(),v,cache),cached_formula(e.arg2(),v,cache)) );
        case OperatorKind::SCALAR: return insert( cache, eptr, make_formula<Real>(e.op(),cached_formula(e.arg(),v,cache),e.num()) );
        default: ARIADNE_FAIL_MSG("Cannot convert expression "<<e<<" to use variables "<<v<<"\n");
    }
}

//! \brief Convert the expression with index type \c I to one with variables indexed by \a J.
Formula<Real> formula(const Expression<Real>& e, const Map<Identifier,Nat>& v)
{
    Map< const Void*, Formula<Real> > cache;
    return cached_formula(e,v,cache);
}



inline Void _set_constant(Real& r, const Real& c) { r=c; }
inline Void _set_constant(Formula<Real>& r, const Real& c) { r=c; }

/*
template<class R, class A> inline
Expression<R> make_expression(OperatorCode op, Expression<A> a);
template<class R, class A1, class A2> inline
Expression<R> make_expression(OperatorCode op, Expression<A1> a1, Expression<A2> a2);
*/

Expression<Boolean> operator&&(Expression<Boolean> e1, Expression<Boolean> e2) {
    return make_expression<Boolean>(AndOp(),e1,e2); }
Expression<Boolean> operator||(Expression<Boolean> e1, Expression<Boolean> e2) {
    return make_expression<Boolean>(OrOp(),e1,e2); }
Expression<Boolean> operator!(Expression<Boolean> e) {
    return make_expression<Boolean>(NotOp(),e); }


Expression<Tribool> operator&&(Expression<Tribool> e1, Expression<Tribool> e2) {
    return make_expression<Tribool>(AndOp(),e1,e2); }
Expression<Tribool> operator||(Expression<Tribool> e1, Expression<Tribool> e2) {
    return make_expression<Tribool>(OrOp(),e1,e2); }
Expression<Tribool> operator!(Expression<Tribool> e) {
    return make_expression<Tribool>(NotOp(),e); }


Expression<Boolean> operator==(Variable<String> v1, const String& s2) {
    return make_expression<Boolean>(Equal(),Expression<String>(v1),Expression<String>(s2)); }
Expression<Boolean> operator!=(Variable<String> v1, const String& s2) {
    return make_expression<Boolean>(Unequal(),Expression<String>(v1),Expression<String>(s2)); }


Expression<Boolean> operator==(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(Equal(),e1,e2); }
Expression<Boolean> operator!=(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(Unequal(),e1,e2); }
Expression<Boolean> operator>=(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(Geq(),e1,e2); }
Expression<Boolean> operator<=(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(Leq(),e1,e2); }
Expression<Boolean> operator>(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(Gtr(),e1,e2); }
Expression<Boolean> operator<(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Boolean>(Less(),e1,e2); }



Expression<Integer> operator+(Expression<Integer> e) {
    return make_expression<Integer>(Pos(),e); }
Expression<Integer> operator-(Expression<Integer> e) {
    return make_expression<Integer>(Neg(),e); }
Expression<Integer> operator+(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Integer>(Add(),e1,e2); }
Expression<Integer> operator-(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Integer>(Sub(),e1,e2); }
Expression<Integer> operator*(Expression<Integer> e1, Expression<Integer> e2) {
    return make_expression<Integer>(Mul(),e1,e2); }



Expression<Tribool> sgn(Expression<Real> e) {
    return make_expression<Tribool>(Sgn(),e); }

Expression<Tribool> operator==(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(Equal(),e1,e2); }
Expression<Tribool> operator!=(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(Unequal(),e1,e2); }
Expression<Tribool> operator>=(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(Geq(),e1,e2); }
Expression<Tribool> operator<=(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(Leq(),e1,e2); }
Expression<Tribool> operator>(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(Gtr(),e1,e2); }
Expression<Tribool> operator<(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Tribool>(Less(),e1,e2); }


Expression<Real> operator+(Expression<Real> e) {
    return make_expression<Real>(Pos(),e); }
Expression<Real> operator-(Expression<Real> e) {
    return make_expression<Real>(Neg(),e); }
Expression<Real> operator+(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(Add(),e1,e2); }
Expression<Real> operator-(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(Sub(),e1,e2); }
Expression<Real> operator*(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(Mul(),e1,e2); }
Expression<Real> operator/(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(Div(),e1,e2); }

Expression<Real> pow(Expression<Real> e, Int n) {
    ARIADNE_NOT_IMPLEMENTED;
    //return make_expression(POW,e,n);
}

Expression<Real> neg(Expression<Real> e) {
    return make_expression<Real>(Neg(),e); }
Expression<Real> rec(Expression<Real> e) {
    return make_expression<Real>(Rec(),e); }
Expression<Real> sqr(Expression<Real> e) {
    return make_expression<Real>(Sqr(),e); }
Expression<Real> sqrt(Expression<Real> e) {
    return make_expression<Real>(Sqrt(),e); }
Expression<Real> exp(Expression<Real> e) {
    return make_expression<Real>(Exp(),e); }
Expression<Real> log(Expression<Real> e) {
    return make_expression<Real>(Log(),e); }
Expression<Real> sin(Expression<Real> e) {
    return make_expression<Real>(Sin(),e); }
Expression<Real> cos(Expression<Real> e) {
    return make_expression<Real>(Cos(),e); }
Expression<Real> tan(Expression<Real> e) {
    return make_expression<Real>(Tan(),e); }

Expression<Real> max(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(Max(),e1,e2); }
Expression<Real> min(Expression<Real> e1, Expression<Real> e2) {
    return make_expression<Real>(Min(),e1,e2); }
Expression<Real> abs(Expression<Real> e) {
    return make_expression<Real>(Abs(),e); }


Expression<Real> operator+(Expression<Real> e, Real c) { return e + Expression<Real>::constant(c); }
Expression<Real> operator-(Expression<Real> e, Real c) { return e - Expression<Real>::constant(c); }
Expression<Real> operator*(Expression<Real> e, Real c) { return e * Expression<Real>::constant(c); }
Expression<Real> operator/(Expression<Real> e, Real c) { return e / Expression<Real>::constant(c); }
Expression<Real> operator+(Real c, Expression<Real> e) { return Expression<Real>::constant(c) + e; }
Expression<Real> operator-(Real c, Expression<Real> e) { return Expression<Real>::constant(c) - e; }
Expression<Real> operator*(Real c, Expression<Real> e) { return Expression<Real>::constant(c) * e; }
Expression<Real> operator/(Real c, Expression<Real> e) { return Expression<Real>::constant(c) / e; }

Expression<Tribool> operator<=(Expression<Real> e, Real c) { return e <= Expression<Real>::constant(c); }
Expression<Tribool> operator< (Expression<Real> e, Real c) { return e <  Expression<Real>::constant(c); }
Expression<Tribool> operator>=(Expression<Real> e, Real c) { return e >= Expression<Real>::constant(c); }
Expression<Tribool> operator> (Expression<Real> e, Real c) { return e >  Expression<Real>::constant(c); }
Expression<Tribool> operator<=(Real c, Expression<Real> e) { return Expression<Real>::constant(c) <= e; }
Expression<Tribool> operator< (Real c, Expression<Real> e) { return Expression<Real>::constant(c) <  e; }
Expression<Tribool> operator>=(Real c, Expression<Real> e) { return Expression<Real>::constant(c) >= e; }
Expression<Tribool> operator> (Real c, Expression<Real> e) { return Expression<Real>::constant(c) >  e; }




template<class A>
typename Logic<A>::Type
evaluate(const Expression<typename Logic<A>::Type>& e, const Map<Identifier,A>& x)
{
    typedef typename Logic<A>::Type R;
    A* aptr=0;
    switch(e.kind()) {
        case OperatorKind::NULLARY: return static_cast<R>(e.val());
        case OperatorKind::UNARY: return compute(e.op(),evaluate(e.arg(),x));
        case OperatorKind::BINARY: return compute(e.op(),evaluate(e.arg1(),x),evaluate(e.arg2(),x));
        case OperatorKind::COMPARISON: return compare(e.op(),evaluate(e.cmp1(aptr),x),evaluate(e.cmp2(aptr),x));
        default: ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" on "<<x<<"\n");
    }
}

//! \brief Evaluate expression \a e on argument \a x which is a map of variable identifiers to values of type \c A.
template<class T>
T
evaluate(const Expression<T>& e, const Map<Identifier,T>& x)
{
    switch(e.kind()) {
        case OperatorKind::VARIABLE: return x[e.var()];
        case OperatorKind::NULLARY: return e.val();
        case OperatorKind::UNARY: return compute(e.op(),evaluate(e.arg(),x));
        case OperatorKind::BINARY: return compute(e.op(),evaluate(e.arg1(),x),evaluate(e.arg2(),x));
        case OperatorKind::SCALAR: return compute(e.op(),evaluate(e.arg(),x),e.num());
        default: ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" on "<<x<<"\n");
    }
}




/*
Boolean evaluate(const Expression<Boolean>& e, const DiscreteValuation& x) {
    const ExpressionInterface<Boolean,Identifier>* eptr=e.raw_pointer();
    const BinaryExpression<Boolean,Identifier>* bptr=dynamic_cast<const BinaryExpression<Boolean,OperatorCode>*>(eptr);
    if(bptr) { return compute(bptr->op,evaluate(bptr->arg1,x),evaluate(bptr->arg2,x)); }
    const UnaryExpression<Boolean,Identifier>* uptr=dynamic_cast<const UnaryExpression<Boolean>*>(eptr);
    if(uptr) { return compute(uptr->op,evaluate(uptr->arg,x)); }
    const ConstantExpression<Boolean,Identifier>* cptr=dynamic_cast<const ConstantExpression<Boolean>*>(eptr);
    if(cptr) { return cptr->value(); }
    const BinaryExpression<Boolean,OperatorCode,String>* bsptr=dynamic_cast<const BinaryExpression<Boolean,OperatorCode,String>*>(eptr);
    if(bsptr) { return compare(bsptr->op,evaluate(bsptr->arg1,x),evaluate(bsptr->arg2,x)); }
    const BinaryExpression<Boolean,OperatorCode,Integer>* bzptr=dynamic_cast<const BinaryExpression<Boolean,OperatorCode,Integer>*>(eptr);
    if(bzptr) { return compare(bsptr->op,evaluate(bzptr->arg1,x),evaluate(bzptr->arg2,x)); }
    ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" to a Boolean using variables "<<x);
}
*/

String evaluate(const Expression<String>& e, const StringValuation& x) {
    switch(e.op()) {
        case OperatorCode::CNST: return e.val();
        case OperatorCode::VAR: return x[e.var()];
        default: ARIADNE_FAIL_MSG("Cannot evaluate expression "<<e<<" to a String using variables "<<x);
    }
}

Integer evaluate(const Expression<Integer>& e, const IntegerValuation& x) {
    return evaluate(e,x.values());
}

Boolean evaluate(const Expression<Boolean>& e, const StringValuation& x) {
    return evaluate(e,x.values());
}

template<class X> Tribool evaluate(const Expression<Tribool>& e, const ContinuousValuation<X>& x) {
    return evaluate(e,x.values());
}

template<class X> X evaluate(const Expression<Real>& e, const ContinuousValuation<X>& x) {
    return evaluate(e,x.values());
}


template Tribool evaluate(const Expression<Tribool>& e, const ContinuousValuation<Real>& x);
template Real evaluate(const Expression<Real>& e, const ContinuousValuation<Real>& x);



template Set<Identifier> arguments(const Expression<Boolean>& e);
template Set<Identifier> arguments(const Expression<Tribool>& e);
template Set<Identifier> arguments(const Expression<Real>& e);






namespace {
template<class I, class X, class Y> inline const Expression<X>& _substitute_variable(const I& ie, const I& is, const Expression<X>& e, const Expression<Y>& s) {
    ARIADNE_ASSERT_MSG(ie!=is,"Cannot substitute expression "<<s<<" for variable "<<ie<<"\n");
    return e; }
template<class I, class X> inline const Expression<X>& _substitute_variable(const I& ie, const I& is, const Expression<X>& e, const Expression<X>& s) {
    return ie==is ? s : e; }
} // namespace

template<class X, class I, class Y> Expression<X> substitute(const Expression<X>& e, const I& v, const Expression<Y>& s) {
    switch(e.kind()) {
        case OperatorKind::COMPARISON: {
            Y* yptr=0;
            const Expression<Y>& c1=e.cmp1(yptr);
            const Expression<Y>& c2=e.cmp2(yptr);
            return make_expression<X>(e.op(),substitute(c1,v,s),substitute(c2,v,s)); }
        case OperatorKind::BINARY: return make_expression<X>(e.op(),substitute(e.arg1(),v,s),substitute(e.arg2(),v,s));
        case OperatorKind::UNARY: return make_expression<X>(e.op(),substitute(e.arg(),v,s));
        case OperatorKind::NULLARY: return make_expression<X>(e.val());
        case OperatorKind::VARIABLE: return _substitute_variable(e.var(),v,e,s);
        default: ARIADNE_FAIL_MSG("Cannot substitute "<<s<<" for a named variable "<<v<<" in an unknown expression "<<e<<"\n");
    }
}

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Variable<Y>& v, const Expression<Y>& s) {
    return substitute(e,v.name(),s);
}

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const Variable<Y>& v, const Y& c) {
    return substitute(e,v.name(),Expression<Y>(c));
}

template<class X, class Y> Expression<X> substitute(const Expression<X>& e, const List< Assignment< Variable<Y>, Expression<Y> > >& a) {
    Expression<X> r=e;
    for(Nat i=0; i!=a.size(); ++i) {
        r=substitute(r,a[i].lhs,a[i].rhs);
    }
    return r;
}
template Expression<Tribool> substitute(const Expression<Tribool>& e, const Variable<Tribool>& v, const Tribool& c);
template Expression<Tribool> substitute(const Expression<Tribool>& e, const Variable<Real>& v, const Real& c);
template Expression<Real> substitute(const Expression<Real>& e, const Variable<Real>& v, const Real& c);
template Expression<Real> substitute(const Expression<Real>& e, const Variable<Real>& v, const Expression<Real>& c);
template Expression<Real> substitute(const Expression<Real>& e, const List< Assignment< Variable<Real>, Expression<Real> > >& c);
template Expression<Tribool> substitute(const Expression<Tribool>& e, const List< Assignment< Variable<Real>, Expression<Real> > >& c);


namespace {

template<class X> inline Expression<X> _simplify(const Expression<X>& e) {
    return e;
}


template<class I> inline Expression<Real> _simplify(const Expression<Real>& e) {
    typedef Real R;

    if(e.kind() == OperatorKind::UNARY) {
        Expression<R> sarg=simplify(e.arg());
        if(sarg.op()==OperatorCode::CNST) {
            return Expression<R>(compute(e.op(),sarg.val()));
        } else {
            return make_expression<R>(e.op(),sarg);
        }
    }

    if(e.kind() != OperatorKind::BINARY) { return e; }

    Expression<R> sarg1=simplify(e.arg1());
    Expression<R> sarg2=simplify(e.arg2());
    Expression<R> zero(static_cast<R>(0));
    Expression<R> one(static_cast<R>(1));
    switch(e.op()) {
        case OperatorCode::ADD:
            if(identical(sarg2,zero)) { return sarg1; }
            if(identical(sarg1,zero)) { return sarg2; }
            break;
        case OperatorCode::SUB:
            if(identical(sarg2,zero)) { return sarg1; }
            if(identical(sarg1,zero)) { return -sarg2; }
            break;
        case OperatorCode::MUL:
            if(identical(sarg1,zero)) { return sarg1; }
            if(identical(sarg2,zero)) { return sarg2; }
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
    return e;

}

template<class I> inline Expression<Tribool> _simplify(const Expression<Tribool>& e) {
    typedef Tribool T;

    if( e.kind()==OperatorKind::UNARY ) {
        Expression<T> sarg=simplify(e.arg());
        if(e.op()==OperatorCode::NOT) {
            if( sarg.op()==OperatorCode::NOT ) {
                return sarg.arg();
            }
            if( sarg.op()==OperatorCode::CNST ) {
                return Expression<T>(compute(e.op(),sarg.val()));
            }
        }
        return make_expression<T>(e.op(),sarg);
    }

    if( e.kind()==OperatorKind::BINARY ) {
        Expression<T> sarg1=simplify(e.arg1());
        Expression<T> sarg2=simplify(e.arg2());
        if( sarg1.op()==OperatorCode::CNST && sarg2.op()==OperatorCode::CNST ) {
            if(e.op()==OperatorCode::AND) { return Expression<T>(sarg1.val() && sarg2.val()); }
            if(e.op()==OperatorCode::OR) { return Expression<T>(sarg1.val() || sarg2.val()); }
            return Expression<T>(compute(e.op(),sarg1.val(),sarg2.val()));
        } else if(sarg1.op()==OperatorCode::CNST) {
            if(e.op()==OperatorCode::AND && sarg1.val()==true) { return sarg2; }
            if(e.op()==OperatorCode::AND && sarg1.val()==false) { return sarg1; }
            if(e.op()==OperatorCode::OR && sarg1.val()==true) { return sarg1; }
            if(e.op()==OperatorCode::OR && sarg1.val()==false) { return sarg2; }
        } else if(sarg2.op()==OperatorCode::CNST) {
            if(e.op()==OperatorCode::AND && sarg2.val()==true) { return sarg1; }
            if(e.op()==OperatorCode::AND && sarg2.val()==false) { return sarg2; }
            if(e.op()==OperatorCode::OR && sarg2.val()==true) { return sarg2; }
            if(e.op()==OperatorCode::OR && sarg2.val()==false) { return sarg1; }
        } else {
            return make_expression<T>(e.op(),sarg1,sarg2);
        }
    }
    return e;
}


} // namespace

template<class X> Expression<X> simplify(const Expression<X>& e) {
    return Ariadne::_simplify(e);
}

template Expression<Real> simplify(const Expression<Real>& e);
template Expression<Tribool> simplify(const Expression<Tribool>& e);



Expression<Real> indicator(Expression<Tribool> e, Sign sign) {
    Tribool value;
    switch(e.op()) {
        case OperatorCode::CNST:
            value=( sign==POSITIVE ? e.val() : Tribool(!e.val()) );
            if(definitely(value)) { return Expression<Real>(+1); }
            else if(not possibly(value)) {  return Expression<Real>(-1); }
            else { return Expression<Real>(0); }
        case OperatorCode::VAR:
            return Expression<Real>(Variable<Real>(e.var()));
        case OperatorCode::GEQ: case OperatorCode::GT:
            if(sign==POSITIVE) { return e.cmp1<Real>()-e.cmp2<Real>(); }
            else { return e.cmp2<Real>()-e.cmp1<Real>(); }
        case OperatorCode::LEQ: case OperatorCode::LT:
            if(sign==POSITIVE) { return e.cmp2<Real>()-e.cmp1<Real>(); }
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


template<class R> Bool identical(const Expression<R>& e1, const Expression<R>& e2)
{
    if(e1.node_ptr()==e2.node_ptr()) { return true; }
    if(e1.op()!=e2.op()) { return false; }
    switch(e1.kind()) {
        case OperatorKind::VARIABLE:
            return e1.var()==e2.var();
        case OperatorKind::NULLARY:
            return same(e1.val(),e2.val());
        case OperatorKind::UNARY:
            return identical(e1.arg(),e2.arg());
        case OperatorKind::BINARY:
            return identical(e1.arg1(),e2.arg1()) && identical(e1.arg2(),e2.arg2());
        default:
            return false;
    }
}

template Bool identical(const Expression<Real>&, const Expression<Real>&);

Bool opposite(Expression<Tribool> e1, Expression<Tribool> e2) {

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




Nat dimension(const Space<Real>& spc)
{
    return spc.size();
}

Nat len(const List< Variable<Real> >& vars)
{
    return vars.size();
}


Formula<Real> formula(const Expression<Real>& e, const Space<Real>& spc)
{
    typedef Real X;
    typedef Identifier I;
    switch(e.kind()) {
        case OperatorKind::SCALAR: return make_formula(e.op(),formula(e.arg(),spc),e.num());
        case OperatorKind::BINARY: return make_formula(e.op(),formula(e.arg1(),spc),formula(e.arg2(),spc));
        case OperatorKind::UNARY: return make_formula(e.op(),formula(e.arg(),spc));
        case OperatorKind::NULLARY: return Formula<X>::constant(e.val());
        case OperatorKind::VARIABLE: return Formula<X>::coordinate(spc.index(e.var()));
        default: ARIADNE_FAIL_MSG("Cannot compute formula for expression "<<e.op()<<"of kind "<<e.kind()<<" in space "<<spc);
    }
}


Formula<Real> formula(const Expression<Real>& e, const List< Variable<Real> >& vars)
{
    return formula(e,Space<Real>(vars));
}

Formula<Real> formula(const Expression<Real>& out, const List< Assignment< Variable<Real>, Expression<Real> > >& aux, const Space<Real> spc)
{
    return formula(substitute(out,aux),spc);
}

List< Formula<Real> > formula(const List< Expression<Real> >& out, const List< Assignment< Variable<Real>, Expression<Real> > >& aux, const Space<Real> spc)
{
    List< Formula<Real> > res;
    for(Nat i=0; i!=out.size(); ++i) {
        res.append(formula(out[i],aux,spc));
    }
    return res;
}



} // namespace Ariadne
