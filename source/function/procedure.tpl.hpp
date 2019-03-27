/***************************************************************************
 *            procedure.tcc
 *
 *  Copyright 2010-17  Pieter Collins
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

#include <iostream>

#include "../utility/container.hpp"
#include "../algebra/vector.hpp"

#include "../numeric/operators.hpp"
#include "../function/formula.hpp"
#include "../algebra/expansion.hpp"
#include "../geometry/box.hpp"

#include "../algebra/evaluate.tpl.hpp"
#include "../function/formula.tpl.hpp"

// FIXME: Added to prevent compilation error in Clang++-5.0. Should not be necessary.
#include "../function/taylor_model.hpp"
#include "../algebra/algebra.hpp"

namespace Ariadne {

extern template class Procedure<ApproximateNumber>;
extern template class Procedure<ValidatedNumber>;

template<class P> Procedure<Number<P>> make_procedure(ScalarMultivariateFunction<P> const& f) {
    typedef Number<P> Y;
    Formula<Y> e=f(Formula<Y>::identity(f.argument_size()));
    return Procedure<Y>(f.argument_size(),e);
}


namespace {
typedef Map<const Void*,SizeType> CacheType;
template<class Y> Void _convert_impl(List<ProcedureInstruction>& p, List<Y>& c, ConstantFormulaNode<Y>& f, CacheType& cache) {
    c.append(f._val); p.append(ProcedureInstruction(Cnst(),c.size()-1)); }
template<class Y> Void _convert_impl(List<ProcedureInstruction>& p, List<Y>& c, IndexFormulaNode<Y>& f, CacheType& cache) {
    p.append(ProcedureInstruction(Var(),f._ind)); }
template<class Y> Void _convert_impl(List<ProcedureInstruction>& p, List<Y>& c, UnaryFormulaNode<Y>& f, CacheType& cache) {
    p.append(ProcedureInstruction(f._op,_convert(p,c,f._arg,cache))); }
template<class Y> Void _convert_impl(List<ProcedureInstruction>& p, List<Y>& c, BinaryFormulaNode<Y>& f, CacheType& cache) {
    p.append(ProcedureInstruction(f._op,_convert(p,c,f._arg1,cache),_convert(p,c,f._arg2,cache))); }
template<class Y> Void _convert_impl(List<ProcedureInstruction>& p, List<Y>& c, GradedFormulaNode<Y>& f, CacheType& cache) {
    p.append(ProcedureInstruction(f._op,_convert(p,c,f._arg,cache),f._num)); }
template<class Y> Void _convert_impl(List<ProcedureInstruction>& p, List<Y>& c, ScalarFormulaNode<Y>& f, CacheType& cache) {
    c.append(f._cnst); p.append(ProcedureInstruction(f._op,c.size()-1u,_convert(p,c,f._arg,cache))); }
} // namespace

template<class Y> SizeType _convert(List<ProcedureInstruction>& p, List<Y>& c, const Formula<Y>& f, Map<const Void*,SizeType>& cache) {
    if(cache.has_key(f.node_ptr())) { return cache[f.node_ptr()]; }
    f.node_ref().visit([&p,&c,&cache](auto s){_convert_impl(p,c,s,cache);});
    const SizeType num=p.size()-1;
    cache.insert(f.node_ptr(),num);
    return num;
}



namespace {
Void _write(OutputStream& os, ConstantProcedureInstruction const& s) { os << "c[" << s._val << "]"; }
Void _write(OutputStream& os, IndexProcedureInstruction const& s) { os << "x[" << s._ind << "]"; }
Void _write(OutputStream& os, UnaryProcedureInstruction const& s) { os << s._op << "(v[" << s._arg << "])"; }
Void _write(OutputStream& os, BinaryProcedureInstruction const& s) { os << s._op << "(v[" << s._arg1 << "],v[" << s._arg2 << "])"; }
Void _write(OutputStream& os, GradedProcedureInstruction const& s) { os << s._op << "(v[" << s._arg << "]," << s._num << ")"; }
Void _write(OutputStream& os, ScalarProcedureInstruction const& s) { os << s._op << "(c[" << s._arg1 << "],v[" << s._arg2 << "])"; }
Void _write(OutputStream& os, ProcedureInstruction const& pi) { pi.visit([&os](auto s){_write(os,s);}); }
} // namespace

template<class Y>
Void _write(OutputStream& os, const List<ProcedureInstruction>& p, const List<Y>& c)
{
    for(SizeType i=0; i!=p.size(); ++i) {
        os << "v[" << i << "]=";
        _write(os,p[i]);
        os << "; ";
    }
}


//template<class Y, class X> Formula<Y> to_formula(const Expansion<MultiIndex,X>& e) {
//    return horner_evaluate(e,Formula<Y>::identity(e.argument_size()));
//};


template<class Y>
Procedure<Y>::Procedure(SizeType as) : _argument_size(as)
{
}

template<class Y>
Procedure<Y>::Procedure(SizeType as, const Formula<Y>& f) : _argument_size(as)
{
    Map<const Void*, SizeType> ind;
    _convert(this->_instructions,this->_constants,f, ind);
}

template<class Y>
Procedure<Y>::Procedure(const ScalarMultivariateFunction<P>& f)
    : Procedure<Y>(f.argument_size(),f(Formula<Y>::identity(f.argument_size())))
{
}


template<class Y> template<class X, EnableIf<IsConvertible<X,Y>>>
Procedure<Y>::Procedure(const Expansion<MultiIndex,X>& e)
    : Procedure(e.argument_size(),horner_evaluate(e,Formula<Y>::identity(e.argument_size())))
{
}


template<class Y> OutputStream& Procedure<Y>::_write(OutputStream& os) const {
    os<<"Procedure( ";
    Ariadne::_write(os,this->_instructions,this->_constants);
    os << "r=v[" << this->_instructions.size()-1u << "] )";
    return os;
}

template<class Y> Procedure<Y>& operator+=(Procedure<Y>& f, const Y& c) {
    f._constants.append(c);
    f.new_unary_instruction(OperatorCode::CNST,f._constants.size()-1);
    f.new_binary_instruction(OperatorCode::ADD,f._instructions.size()-1,f._instructions.size()-2);
    return f;
}


template<class Y>
Vector<Procedure<Y>>::Vector(SizeType as, const Vector<Formula<Y>>& f)
    : _argument_size(as),  _results(f.size(),0u)
{
    Map<const Void*, SizeType> ind;
    for(SizeType i=0; i!=f.size(); ++i) {
        _convert(this->_instructions,this->_constants,f[i], ind);
    }
    for(Nat i=0; i!=f.size(); ++i) {
        this->_results[i]=ind[f[i].node_ptr()];
    }
}

template<class Y>
Vector<Procedure<Y>>::Vector(VectorMultivariateFunction<P> const& f)
    : Vector<Procedure<Y>>(f.argument_size(),f(Formula<Y>::identity(f.argument_size())))
{
}

template<class Y>
Vector<Procedure<Y>>::Vector(const Procedure<Y>& f)
    : _constants(f._constants)
    , _instructions(f._instructions)
    , _results(1u,f._instructions.size()-1u)
{
}

template<class Y>
OutputStream& Vector<Procedure<Y>>::_write(OutputStream& os) const {
    os<<"Procedure( ";
    Ariadne::_write(os,this->_instructions,this->_constants);
    os << "r=v" << this->_results << " )";
    return os;
}

template<class OP> struct AlgebraScalar { };

template<class V> decltype(auto) visit(OperatorCode op(), V const& v) {
    switch(op()) {
//        case OperatorCode::CNST: { X r=x.zero_element(); r=c[instruction.arg]; v[i]=r; } break;
        case OperatorCode::SADD:  return v(AlgebraScalar<Add>());
        case OperatorCode::SSUB:  return v(AlgebraScalar<Sub>());
        case OperatorCode::SMUL:  return v(AlgebraScalar<Mul>());
        case OperatorCode::SDIV:  return v(AlgebraScalar<Div>());
        case Add::code():  return v(Add());
        case Sub::code():  return v(Sub());
        case Mul::code():  return v(Mul());
        case Div::code():  return v(Div());
        case Pow::code():  return v(Pow());
        case Abs::code():  return v(Abs());
        case Pos::code():  return v(Pos());
        case Neg::code():  return v(Neg());
        case Rec::code():  return v(Rec());
        case Sqr::code():  return v(Sqr());
        case Sqrt::code(): return v(Sqrt());
        case Exp::code():  return v(Exp());
        case Log::code():  return v(Log());
        case Sin::code():  return v(Sin());
        case Cos::code():  return v(Cos());
        case Tan::code():  return v(Tan());
        case Atan::code(): return v(Atan());
        default:   ARIADNE_FAIL_MSG("Unrecognised operator "<<op());
    }
}



template<class X, class OP, class A> X evaluate(Symbolic<OP,A> const& s, X const& x);
template<class X, class OP, class A1, class A2> X evaluate(Symbolic<OP,A1,A2> const& s, X const& x1, X const& x2);
template<class X, class OP, class A1, class A2> X evaluate(Symbolic<OP,A1,A2> const& s, X const& x1, Int n2);

namespace {

template<class X, class OP> Void propagate(X& r, OP op, X const& x1, X const& x2) { r=op(x1,x2); }
template<class X, class OP> Void propagate(X& r, OP op, X const& x) { r=op(x); }

template<class X> Void propagate(X& r, BinaryElementaryOperator op, X const& x1, X const& x2) { op.visit([&](auto op){r=op(x1,x2);}); }
template<class X> Void propagate(X& r, UnaryElementaryOperator op, X const& x) { op.visit([&](auto op){r=op(x);}); }
template<class X, class N> Void propagate(X& r, GradedElementaryOperator op, X const& x, N n) { op.visit([&](auto op){r=op(x,n);}); }
template<class Y, class X> Void propagate(X& r, BinaryElementaryOperator op, Y const& y1, X const& x2) { op.visit([&](auto op){r=op(y1,x2);}); }

template<class X, class Y> Void _execute_impl(SizeType r, List<X>& v, const ConstantProcedureInstruction& pi, const List<Y>& c, const Vector<X>& x) {
    v[r]=make_constant(c[pi._val],x.zero_element()); }
template<class X, class Y> Void _execute_impl(SizeType r, List<X>& v, const IndexProcedureInstruction& pi, const List<Y>& c, const Vector<X>& x) {
    v[r]=x[pi._ind]; }
template<class X, class Y> Void _execute_impl(SizeType r, List<X>& v, const UnaryProcedureInstruction& pi, const List<Y>& c, const Vector<X>& x) {
    propagate(v[r],pi._op,v[pi._arg]); }
template<class X, class Y> Void _execute_impl(SizeType r, List<X>& v, const BinaryProcedureInstruction& pi, const List<Y>& c, const Vector<X>& x) {
    propagate(v[r],pi._op,v[pi._arg1],v[pi._arg2]); }
template<class X, class Y> Void _execute_impl(SizeType r, List<X>& v, const GradedProcedureInstruction& pi, const List<Y>& c, const Vector<X>& x) {
    propagate(v[r],pi._op,v[pi._arg],pi._num); }
template<class X, class Y> Void _execute_impl(SizeType r, List<X>& v, const ScalarProcedureInstruction& pi, const List<Y>& c, const Vector<X>& x) {
    propagate(v[r],pi._op,c[pi._arg1],v[pi._arg2]); }
}

// \related Procedure \brief Evaluate a function \a f defined by an algorithmic procedure
template<class X, class Y> Void _execute(List<X>& v, const List<ProcedureInstruction>& p, const List<Y>& c, const Vector<X>& x)
{
    ARIADNE_ASSERT(v.size()==p.size());
    for(SizeType i=0; i!=p.size(); ++i) {
        p[i].visit([i,&v,&c,&x](auto s){return _execute_impl(i,v,s,c,x);});
    }
}

template<class X> Void _backpropagate(X const& r, Add, X& a1, X& a2) { restrict(a1,r-a2); restrict(a2,r-a1); }
template<class X> Void _backpropagate(X const& r, Sub, X& a1, X& a2) { restrict(a1,a2+r); restrict(a2,a1-r); }
template<class X> Void _backpropagate(X const& r, Mul, X& a1, X& a2) { restrict(a1,r/a2); restrict(a2,r/a1); }
template<class X> Void _backpropagate(X const& r, Div, X& a1, X& a2) { restrict(a1,a2*r); restrict(a2,a1/r); }
template<class X> Void _backpropagate(X const& r, Max, X& a1, X& a2) { restrict(a1,max(r,a2)); restrict(a2,max(r,a1)); }
template<class X> Void _backpropagate(X const& r, Min, X& a1, X& a2) { restrict(a1,min(r,a2)); restrict(a2,min(r,a1)); }

template<class X, class Y> Void _backpropagate_cnst(X const& r, Add, Y const& c1, X& a2) { restrict(a2,r-c1); }
template<class X, class Y> Void _backpropagate_cnst(X const& r, Sub, Y const& c1, X& a2) { restrict(a2,c1-r); }
template<class X, class Y> Void _backpropagate_cnst(X const& r, Mul, Y const& c1, X& a2) { restrict(a2,r/c1); }
template<class X, class Y> Void _backpropagate_cnst(X const& r, Div, Y const& c1, X& a2) { restrict(a2,c1/r); }
template<class X, class Y> Void _backpropagate_cnst(X const& r, Max, Y const& c1, X& a2) { restrict(a2,max(r,c1)); }
template<class X, class Y> Void _backpropagate_cnst(X const& r, Min, Y const& c1, X& a2) { restrict(a2,min(r,c1)); }

template<class X> Void _backpropagate(X const& r, Pos, X& a) { restrict(a,r); }
template<class X> Void _backpropagate(X const& r, Neg, X& a) { restrict(a,neg(r)); }
template<class X> Void _backpropagate(X const& r, Rec, X& a) { restrict(a,rec(r)); }
template<class X> Void _backpropagate(X const& r, Sqr, X& a) { restrict(a,sqrt(r)); }
template<class X> Void _backpropagate(X const& r, Pow, X& a, Int n) { restrict(a,exp(log(r)/n)); }
template<class X> Void _backpropagate(X const& r, Sqrt, X& a) { restrict(a,sqr(r)); }
template<class X> Void _backpropagate(X const& r, Exp, X& a) { restrict(a,log(r)); }
template<class X> Void _backpropagate(X const& r, Log, X& a) { restrict(a,exp(r)); }
            // FIXME: restricting asin/acos of interval should not be dependent on any branch of asin/acos
template<class X> Void _backpropagate(X const& r, Sin, X& a) { restrict(a,asin(r)); }
template<class X> Void _backpropagate(X const& r, Cos, X& a) { restrict(a,acos(r)); }
template<class X> Void _backpropagate(X const& r, Tan, X& a) { restrict(a,atan(r)); }
template<class X> Void _backpropagate(X const& r, Atan, X& a) { restrict(a,tan(r)); }
template<class X> Void _backpropagate(X const& r, Equal, X& a1, X& a2) {  restrict(a1,r); restrict(a2,r); }
template<class X> Void _backpropagate(X const& r, Leq, X& a1, X& a2) {
    FloatDPValue inf_(inf); restrict(a1,X(-inf_,a2.upper())); restrict(a1,X(a2.lower(),+inf_)); }

template<class X> Void _backpropagate(X const& r, UnaryElementaryOperator op, X& a) {
    return op.visit([&r,&a](auto op){_backpropagate(r,op,a);}); }
template<class X> Void _backpropagate(X const& r, BinaryElementaryOperator op, X& a1, X& a2) {
    return op.visit([&r,&a1,&a2](auto op){_backpropagate(r,op,a1,a2);}); }
template<class X> Void _backpropagate(X const& r, GradedElementaryOperator op, X& a1, Int n2) {
    return op.visit([&r,&a1,n2](auto op){_backpropagate(r,op,a1,n2);}); }
template<class X, class Y> Void _backpropagate_cnst(X const& r, BinaryElementaryOperator op, Y const& c1, X& a2) {
    return op.visit([&r,&c1,&a2](auto op){_backpropagate_cnst(r,op,c1,a2);}); }

template<class X, class Y> Void _backpropagate(SizeType r, Vector<X>& x, List<X>& v, ConstantProcedureInstruction s, List<Y> const& c) { }

template<class X, class Y> Void _backpropagate(SizeType r, Vector<X>& x, List<X>& v, IndexProcedureInstruction s, List<Y> const& c) {
    restrict(x[s._ind],v[r]); }
template<class X, class Y> Void _backpropagate(SizeType r, Vector<X>& x, List<X>& v, UnaryProcedureInstruction s, List<Y> const& c) {
    _backpropagate(v[r],s._op,v[s._arg]); }
template<class X, class Y> Void _backpropagate(SizeType r, Vector<X>& x, List<X>& v, BinaryProcedureInstruction s, List<Y> const& c) {
    _backpropagate(v[r],s._op,v[s._arg1],v[s._arg2]); }
template<class X, class Y> Void _backpropagate(SizeType r, Vector<X>& x, List<X>& v, GradedProcedureInstruction s, List<Y> const& c) {
    _backpropagate(v[r],s._op,v[s._arg],s._num); }
template<class X, class Y> Void _backpropagate(SizeType r, Vector<X>& x, List<X>& v, ScalarProcedureInstruction s, List<Y> const& c) {
    _backpropagate_cnst(v[r],s._op,c[s._arg1],v[s._arg2]); }



template<class X, class Y> Void _backpropagate(Vector<X>& x, List<X>& v, const List<ProcedureInstruction>& p, const List<Y>& c)
{
    FloatDPValue inf_(inf);
    ARIADNE_ASSERT(v.size()==p.size());
    SizeType r=p.size();
    while(r!=0u) {
        --r;
        p[r].visit([r,&x,&v,&c](auto s){_backpropagate(r,x,v,s,c);});
    }
    // POSTCONDITION: No nan's get propagated to x
}




template<class X, class D> void back_gradient(D const& dr, Add, X const& a1, D& da1, X const& a2, D& da2) { da1+=dr; da2+=dr; }
template<class X, class D> void back_gradient(D const& dr, Sub, X const& a1, D& da1, X const& a2, D& da2) { da1+=dr; da2-=dr; }
template<class X, class D> void back_gradient(D const& dr, Mul, X const& a1, D& da1, X const& a2, D& da2) { da1+=dr*a2; da2+=a1*dr; }
template<class X, class D> void back_gradient(D const& dr, Div, X const& a1, D& da1, X const& a2, D& da2) { da1+=dr/a2; da2-=dr/a2*(a1/a2); }
template<class X, class D> void back_gradient(D const& dr, Max, X const& a1, D& da1, X const& a2, D& da2) { assert(false); }
template<class X, class D> void back_gradient(D const& dr, Min, X const& a1, D& da1, X const& a2, D& da2) { assert(false); }

template<class X, class D> void back_gradient(D const& dr, BinaryElementaryOperator op, X const& a1, D& da1, X const& a2, D& da2) {
    op.visit([&dr,&a1,&da1,&a2,&da2](auto op){return back_gradient(dr,op,a1,da1,a2,da2);}); }


template<class OP, class X, class D> void back_gradient(D const& dr, OP op, X const& a, D& da) { da+=op.derivative(a,dr); }
template<class X, class D> void back_gradient(D const& dr, Abs op, X const& a, D& da) { assert(false); }

template<class X, class D> void back_gradient(D const& dr, UnaryElementaryOperator op, X const& a, D& da) {
    op.visit([&dr,&a,&da](auto op){return back_gradient(dr,op,a,da);}); }

template<class Y, class X, class D> void back_gradient(D const& dr, Add op, Y const& c1, X const& a2, D& da2) { da2+=dr; }
template<class Y, class X, class D> void back_gradient(D const& dr, Sub op, Y const& c1, X const& a2, D& da2) { da2-=dr; }
template<class Y, class X, class D> void back_gradient(D const& dr, Mul op, Y const& c1, X const& a2, D& da2) { da2+=c1*dr; }
template<class Y, class X, class D> void back_gradient(D const& dr, Div op, Y const& c1, X const& a2, D& da2) { da2-=(dr/a2)*(c1/a2); }
template<class Y, class X, class D> void back_gradient(D const& dr, Max op, Y const& c1, X const& a2, D& da2) { assert(false); }
template<class Y, class X, class D> void back_gradient(D const& dr, Min op, Y const& c1, X const& a2, D& da2) { assert(false); }

template<class Y, class X, class D> void back_gradient(D const& dr, BinaryElementaryOperator op, Y const& c1, X const& a2, D& da2) {
    op.visit([&dr,&c1,&a2,&da2](auto op){return back_gradient(dr,op,c1,a2,da2);}); }


template<class N, class X, class D> void back_gradient(D const& dr, Pow op, X const& a, D& da, N n) {
    da+=op.derivative(a,dr,n); }

template<class N, class X, class D> void back_gradient(D const& dr, GradedElementaryOperator op, X const& a, D& da, N n) {
    op.visit([&dr,&a,&da,n](auto op){return back_gradient(dr,op,a,da,n);}); }



template<class X, class Y> void back_gradient_impl_(SizeType r, ConstantProcedureInstruction const& s, Vector<X> const& x, List<Y> const& c, List<X> const& v, List<X>& dv, Covector<X>& dfx) { }
template<class X, class Y> void back_gradient_impl_(SizeType r, IndexProcedureInstruction const& s, Vector<X> const& x, List<Y> const& c, List<X> const& v, List<X>& dv, Covector<X>& dfx) {
    dfx[s._ind]+=dv[r]; }
template<class X, class Y> void back_gradient_impl_(SizeType r, UnaryProcedureInstruction const& s, Vector<X> const& x, List<Y> const& c, List<X> const& v, List<X>& dv, Covector<X>& dfx) {
    back_gradient(dv[r],s._op,v[s._arg],dv[s._arg]); }
template<class X, class Y> void back_gradient_impl_(SizeType r, BinaryProcedureInstruction const& s, Vector<X> const& x, List<Y> const& c, List<X> const& v, List<X>& dv, Covector<X>& dfx) {
    back_gradient(dv[r],s._op,v[s._arg1],dv[s._arg1],v[s._arg2],dv[s._arg2]); }
template<class X, class Y> void back_gradient_impl_(SizeType r, ScalarProcedureInstruction const& s, Vector<X> const& x, List<Y> const& c, List<X> const& v, List<X>& dv, Covector<X>& dfx) {
    back_gradient(dv[r],s._op,c[s._arg1],v[s._arg2],dv[s._arg2]); }
template<class X, class Y> void back_gradient_impl_(SizeType r, GradedProcedureInstruction const& s, Vector<X> const& x, List<Y> const& c, List<X> const& v, List<X>& dv, Covector<X>& dfx) {
    back_gradient(dv[r],s._op,v[s._arg],dv[s._arg],s._num); }


template<class X, class Y> void back_gradient_impl(SizeType r, ProcedureInstruction const& pri, Vector<X> const& x, List<Y> const& c, List<X>& v, List<X>& dv, Covector<X>& dfx) {
    pri.visit([r,&x,&c,&v,&dv,&dfx](auto s){back_gradient_impl_(r,s,x,c,v,dv,dfx);}); }


template<class X, class Y> Covector<X> gradient(Procedure<Y> const& f, Vector<X> const& x) {
    X z=x.zero_element();
    auto& p=f._instructions;
    auto& c=f._constants;

    List<X> v(p.size(),z);
    List<X> dv(v.size(),z);
    Covector<X> dfx(x.size(),z);

    execute(v,f,x);

    SizeType r=p.size();
    dv[r-1]=1;
    while(r!=0u) {
        --r;
        back_gradient_impl(r,p[r],x,c,v,dv,dfx);
    }
    return dfx;
}

} // namespace Ariadne

#include "../algebra/fixed_differential.hpp"

namespace Ariadne {

template<class X, class Y> X hessian(Procedure<Y> const& f, Vector<X> const& x, Vector<X> const& s) {
    auto da=UnivariateSecondDifferential<X>::variable(s.zero_element());
    Vector<UnivariateSecondDifferential<X>> dxpas=x+da*s;
    UnivariateSecondDifferential<X> dfxpas=evaluate(f,dxpas);
    return dfxpas.hessian();
}

}

