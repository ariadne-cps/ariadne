/***************************************************************************
 *            symbolic/object.hpp
 *
 *  Copyright 2013-18  Pieter Collins
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

/*! \file symbolic/object.hpp
 *  \brief Symbolic objects
 */

#ifndef ARIADNE_OBJECT_HPP
#define ARIADNE_OBJECT_HPP

#include <variant>

#include "../utility/writable.hpp"
#include "../numeric/real.hpp"
#include "../numeric/real_interface.hpp"
#include "../numeric/float_bounds.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"

#include "templates.hpp"
#include "variables.hpp"

namespace Ariadne {

template<class SIG> using ResultOfType = typename std::result_of<SIG>::type;

struct Id { };


template<> class Vector<Real>;
template<> class Matrix<Real>;

template<class... TS> using Variant = std::variant<TS...>;

struct Get { template<class V, class I> decltype(auto) operator() (V const& v, I const& i) { return v[i]; } };

template<class... AS> struct Arguments;

template<class A> struct Argument : WritableInterface { A arg; };
template<> struct Arguments<> : WritableInterface { };
template<class A1> struct Arguments<A1> : Arguments<> { A1 arg1; };
template<class A1, class A2> struct Arguments<A1,A2> : Arguments<A1> { A2 arg2;
    Arguments(A1 a1, A2 a2) : Arguments<A1>{a1}, arg2(std::move(a2)) { } };
template<class A1, class A2, class A3> struct Arguments<A1,A2,A3> : Arguments<A1,A2> { A3 arg3;
    Arguments(A1 a1, A2 a2, A3 a3) : Arguments<A1,A2>(a1,a2), arg3(std::move(a3)) { } };

struct Write { Write(OutputStream& os) : _os(os) { } template<class T> Void operator()(T const& t) { _os <<t; } OutputStream& _os; };

template<class T, class... TS> OutputStream& operator<<(OutputStream& os, Variant<T,TS...> const& var) {
    Write writer(os); std::visit(writer,var); return os; }


template<class T> class SymbolicVariable {
    Identifier _name;
  public:
    SymbolicVariable(Identifier name) : _name(name) { }
    Identifier name() const { return _name; }
    friend OutputStream& operator<<(OutputStream& os, SymbolicVariable<T> const& v) { return os << v.name(); }
};

template<class... TS> class Valuations {
    Map<Identifier, Variant<TS...>> _vals;
  public:
    template<class T> T get(SymbolicVariable<T> v) const { return std::get<T>(_vals[v.name()]); }
    template<class T> Void set(SymbolicVariable<T> v, SelfType<T> t) { _vals[v.name()]=t; }
    template<class T> T& operator[] (SymbolicVariable<T> v) { return std::get<T>(_vals[v.name()]); }
    friend OutputStream& operator<<(OutputStream& os, Valuations<TS...> const& ts) { return os << ts._vals; }
};

struct ContinuousValuations : Valuations<Real,RealVector,RealMatrix> { };

template<class S> class Symbolic;
template<class S> class Func;

template<class R, class... AS> class Func<R(AS...)> {
  public:
    struct Interface : public WritableInterface {
        virtual R operator() (AS const& ...) const = 0;
    };
};

template<class S> class Coordinate;
template<class X> class Coordinate<X(Vector<X>)> {
    SizeType _n, _i;
  public:
    Coordinate(SizeType n, SizeType i) : _n(n), _i(i) { }
    X operator() (Vector<X> const& v) { return v[_i]; }
};


template<class S> class Coordinate { Coordinate()  { } };

template<class S> struct InterfaceTypedef { typedef typename S::Interface Type; };
template<class R, class... AS> struct InterfaceTypedef<R(AS...)> { typedef typename Func<R(AS...)>::Interface Type; };
template<class R> struct InterfaceTypedef<R()> { typedef typename R::Interface Type; };
template<> struct InterfaceTypedef<Real> { typedef Real::Interface Type; };
template<class S> using InterfaceType = typename InterfaceTypedef<S>::Type;

template<class OP, class T1, class T2, class... AS> Symbolic<ResultOfType<OP(T1,T2)>(AS...)> make_symbolic(OP op, Symbolic<T1(AS...)> t1, Symbolic<T2(AS...)> t2);
template<class OP, class T, class... AS> Symbolic<ResultOfType<OP(T)>(AS...)> make_symbolic(OP op, Symbolic<T(AS...)> t);
template<class OP, class T> Symbolic<ResultOfType<OP(T)>> make_symbolic(OP op, Symbolic<T> t);
template<class OP, class T1, class T2> Symbolic<ResultOfType<OP(T1,T2)>> make_symbolic(OP op, Symbolic<T1> t1, Symbolic<T2> t2);

class SymbolicOperators {
    template<class S1, class S2> friend decltype(auto) add(Symbolic<S1> s1, Symbolic<S2> s2) { return make_symbolic(Add(),s1,s2); }
    template<class S1, class S2> friend decltype(auto) sub(Symbolic<S1> s1, Symbolic<S2> s2) { return make_symbolic(Sub(),s1,s2); }
    template<class S1, class S2> friend decltype(auto) mul(Symbolic<S1> s1, Symbolic<S2> s2) { return make_symbolic(Mul(),s1,s2); }
    template<class S1, class S2> friend decltype(auto) div(Symbolic<S1> s1, Symbolic<S2> s2) { return make_symbolic(Div(),s1,s2); }
    template<class S> friend decltype(auto) neg(Symbolic<S> s) { return make_symbolic(Neg(),s); }
    template<class S> friend decltype(auto) rec(Symbolic<S> s) { return make_symbolic(Rec(),s); }
    template<class S> friend decltype(auto) sqrt(Symbolic<S> s) { return make_symbolic(Sqrt(),s); }
};

template<class R, class... AS> class Symbolic<R(AS...)> : public SymbolicOperators {
  public:
    typedef InterfaceType<R(AS...)> Interface;
  private:
    SharedPointer<Interface> _ptr;
  public:
    template<class E> Symbolic<E(AS...)> arg1() const;
    template<class E> Symbolic<E(AS...)> arg2() const;
  public:
    Symbolic(SharedPointer<Interface> ptr) : _ptr(ptr) { }
    template<class... ARGS> decltype(auto) operator() (ARGS... args) const { return this->_ptr->operator() (args...); }
    friend OutputStream& operator<<(OutputStream& os, Symbolic<R(AS...)> const& s) { return s._ptr->write(os); }
};

template<class R> class Symbolic : public SymbolicOperators {
  public:
    typedef InterfaceType<R> Interface;
  private:
    SharedPointer<Interface> _ptr;
  private: public:
    Symbolic<R>(SharedPointer<Interface> ptr) : _ptr(ptr) { }
    operator Interface const& () const { return *_ptr; }
    operator R const () const { return R(_ptr); }
  public:
    Symbolic<R>(R const& r) { *this = make_symbolic<R>(Cnst(),r); }
    friend OutputStream& operator<<(OutputStream& os, Symbolic<R> const& s) { return s._ptr->_write(os); }
};

template<class S, class OP, class... ES> class SymbolicExpressionTemplate;

template<class R, class... AS, class OP, class... ES>
class SymbolicExpressionTemplate<R(AS...),OP,ES...> : public ExpressionTemplate<OP,ES...>, public Symbolic<R(AS...)>::Interface {
  public:
    using ExpressionTemplate<OP,ES...>::ExpressionTemplate;
    virtual R operator() (AS const& ... as) const override final {
        ExpressionTemplate<OP,ES...> const& expr = *this; return expr(as...); }
        //return this->ExpressionTemplate<OP,ES...>::operator()(a); }
    virtual OutputStream& _write(OutputStream& os) const override final { return os << static_cast<ExpressionTemplate<OP,ES...>const&>(*this); }
};

template<class OP, class... ES>
class SymbolicExpressionTemplate<Real,OP,ES...> : public ExpressionTemplate<OP,ES...>, public Symbolic<Real>::Interface {
  public:
    using ExpressionTemplate<OP,ES...>::ExpressionTemplate;
#warning
//    virtual ValidatedReal _compute(Accuracy acc) const override final;
    virtual ValidatedReal _compute(Effort eff) const override final { assert(false); }
    virtual FloatDPBounds _compute_get(DP pr) const override final { return this->vist([&](Real const& r){return FloatDPBounds(r.compute_get(pr));}); }
    virtual FloatMPBounds _compute_get(MP pr) const override final { assert(false); }
        //return this->ExpressionTemplate<OP,ES...>::operator()(a); }
    virtual OutputStream& _write(OutputStream& os) const override final { return os << static_cast<ExpressionTemplate<OP,ES...>const&>(*this); }
};

/*
template<class OP, class... ES>
class SymbolicExpressionTemplate<Vector<Real>,Void,OP,ES...> : public ExpressionTemplate<OP,ES...>, public Symbolic<Vector<Real>>::Interface {
  public:
    using ExpressionTemplate<OP,ES...>::ExpressionTemplate;
    virtual Vector<FloatMPBounds> get(Effort eff, MP pr) const override final {
        ExpressionTemplate<OP,ES...> const& expr = *this; return expr.get(eff,pr); }
        //return this->ExpressionTemplate<OP,ES...>::operator()(a); }
    virtual OutputStream& _write(OutputStream& os) const override final { return os << static_cast<ExpressionTemplate<OP,ES...>const&>(*this); }
};
*/

template<class R, class... AS> template<class E> auto Symbolic<R(AS...)>::arg1() const -> Symbolic<E(AS...)> {
    typedef Symbolic<E(AS...)> SE; return std::dynamic_pointer_cast<Arguments<SE,SE>>(_ptr)->arg1; }
template<class R, class... AS> template<class E> auto Symbolic<R(AS...)>::arg2() const -> Symbolic<E(AS...)> {
    typedef Symbolic<E(AS...)> SE; return std::dynamic_pointer_cast<Arguments<SE,SE>>(_ptr)->arg2; }


template<class OP, class T>
Symbolic<ResultOfType<OP(T)>> make_symbolic(OP op, Symbolic<T> t) {
    typedef ResultOfType<OP(T)> R; return Symbolic<R>(std::make_shared<SymbolicExpressionTemplate<R,OP,Symbolic<T>>>(op,t));
}

template<class OP, class T1, class T2>
Symbolic<ResultOfType<OP(T1,T2)>> make_symbolic(OP op, Symbolic<T1> t1, Symbolic<T2> t2) {
    typedef ResultOfType<OP(T1,T2)> R; return Symbolic<R>(std::make_shared<SymbolicExpressionTemplate<R,OP,Symbolic<T1>,Symbolic<T2>>>(op,t1,t2));
}




template<class OP, class T, class... AS>
Symbolic<ResultOfType<OP(T)>(AS...)> make_unary_symbolic(OP op, Symbolic<T(AS...)> t) {
    typedef ResultOfType<OP(T)> R;
    return Symbolic<R(AS...)>(std::make_shared<SymbolicExpressionTemplate<R(AS...),OP,Symbolic<T(AS...)>>>(op,t));
}

template<class OP, class T1, class T2, class... AS>
Symbolic<ResultOfType<OP(T1,T2)>(AS...)> make_binary_symbolic(OP op, Symbolic<T1(AS...)> t1, Symbolic<T2(AS...)> t2) {
    typedef ResultOfType<OP(T1,T2)> R;
    return Symbolic<R(AS...)>(std::make_shared<SymbolicExpressionTemplate<R(AS...),OP,Symbolic<T1(AS...)>,Symbolic<T2(AS...)>>>(op,t1,t2));
}


template<class OP, class T1, class T2, class... AS>
Symbolic<ResultOfType<OP(T1,T2)>(AS...)> make_symbolic(OP op, Symbolic<T1(AS...)> t1, Symbolic<T2(AS...)> t2) {
    typedef ResultOfType<OP(T1,T2)> R;
    return Symbolic<R(AS...)>(std::make_shared<SymbolicExpressionTemplate<R(AS...),OP,Symbolic<T1(AS...)>,Symbolic<T2(AS...)>>>(op,t1,t2)); }

template<class OP, class T, class... AS>
Symbolic<ResultOfType<OP(T)>(AS...)> make_symbolic(OP op, Symbolic<T(AS...)> t) {
    typedef ResultOfType<OP(T)> R;
    return Symbolic<R(AS...)>(std::make_shared<SymbolicExpressionTemplate<R(AS...),OP,Symbolic<T(AS...)>>>(op,t)); }

/*
template<class OP, class A>
Symbolic<decltype(declval<OP>()(declval<A>())),A> make_symbolic(OP op) {
    typedef decltype(declval<OP>()(declval<A>())) R;
    return Symbolic<R,A>(std::make_shared<SymbolicExpressionTemplate<R,A,OP>>(op)); }
*/

template<class T>
Symbolic<T(ContinuousValuations)> make_symbolic(Var op, SymbolicVariable<T> const& v) {
    return Symbolic<T(ContinuousValuations)>(std::make_shared<SymbolicExpressionTemplate<T,ContinuousValuations,Var,SymbolicVariable<T>>>(op,v));
}

template<class R, class C> Symbolic<R> make_symbolic(Cnst op, C const& c) {
    return Symbolic<C>(std::make_shared<SymbolicExpressionTemplate<R,Cnst,C>>(op,c)); }

//template<class SIG, class C> Symbolic<SIG> make_symbolic(Cnst op, C const& c);
//    return Symbolic<SIG>(std::make_shared<SymbolicExpressionTemplate<SIG,Cnst,ConstantFunction<C,A>>>(op,c)); }

/*
template<>
Symbolic<Real,Void> inline make_symbolic(Cnst op, Real const& c) {
    return Symbolic<Real()>(std::make_shared<SymbolicExpressionTemplate<Real,Void,Cnst,Real>>(op,c)); }


template<class R, class A> inline Symbolic<R,A>::Symbolic(R const& c) {
    *this = make_symbolic<A>(Cnst(),c); }
*/


template<class R> class SymbolicExpression : public Symbolic<R(ContinuousValuations)> {
    using A = ContinuousValuations;
    using Base=Symbolic<R(ContinuousValuations)>;
  public:
    SymbolicExpression(R const& c) : Base(make_symbolic<R,A>(Cnst(),c)) { }
    SymbolicExpression(SymbolicVariable<R> const& v)  : Base(std::make_shared<SymbolicExpressionTemplate<R,ContinuousValuations,Var,SymbolicVariable<R>>>(Var(),v)) { }
};

template<class S> class SymbolicFunction;

template<class... AS> struct SymbolicFunctionConstructors;

template<class... AS> struct SymbolicFunctionConstructors {
    template<std::size_t N> static SymbolicFunction<decltype(std::get<N>(declval<Tuple<AS...>>()))(AS...)> argument() {
        typedef decltype(std::get<N>(declval<Tuple<AS...>>())) R; }
};

template<class T> struct IdentityFunction { T const& operator() (T const& t) { return t; } };


template<class X> struct SymbolicFunctionConstructors<Vector<X>> {
    static SymbolicFunction<X(Vector<X>)> coordinate(SizeType n, SizeType i);
    static SymbolicFunction<Vector<X>(Vector<X>)> identity(SizeType n);
};

template<> struct SymbolicFunctionConstructors<Real> {
    static SymbolicFunction<Real(Real)> identity();
};

template<class R, class... AS> class SymbolicFunction<R(AS...)>
    : public Symbolic<R(AS...)>, public SymbolicFunctionConstructors<AS...>
{
  public:
    using Symbolic<R(AS...)>::Symbolic;
    SymbolicFunction(R const& c) : Symbolic<R(AS...)>(std::make_shared<SymbolicExpressionTemplate<R(AS...),Cnst,R>>(Cnst(),c)) { }
};

SymbolicFunction<Real(Real)> SymbolicFunctionConstructors<Real>::identity() {
    return SymbolicFunction<Real(Real)>(std::make_shared<SymbolicExpressionTemplate<Real(Real),Var>>(Var())); }

template<class X1, class X2> decltype(auto) mul(Matrix<X1> const& A1, Vector<X2> v2) { return A1*v2; }

template<class... AS> decltype(auto) sin(Symbolic<Real(AS...)> const& sr) { return make_unary_symbolic(Sin(),sr); }

template<class R1, class R2, class... AS> decltype(auto) operator+(Symbolic<R1(AS...)> const& s1, Symbolic<R2(AS...)> const& s2) { return make_symbolic(Add(),s1,s2); }
template<class R1, class R2, class... AS> decltype(auto) operator+(Symbolic<R1(AS...)> const& s1, Constant<R2> const& c2) { return s1 + make_symbolic<R2(AS...)>(Cnst(),c2); }
template<class X1, class X2, class... AS> decltype(auto) operator*(Symbolic<X1(AS...)> const& e1, Symbolic<X2(AS...)> const& e2) { return make_symbolic(Mul(),e1,e2); }
template<class X1, class X2> decltype(auto) operator*(SymbolicVariable<X1> const& v1, SymbolicVariable<X2> const& v2) {
    SymbolicExpression<X1> e1(v1); SymbolicExpression<X2> e2(v2); return e1*e2; }



template<class... AS> decltype(auto) exp(Symbolic<Real(AS...)> const& s) { return make_symbolic(Exp(),s); }

template<class OP, class R, class... AS> auto compose(OP op, Symbolic<R(AS...)> const& sym) -> Symbolic<ResultOfType<OP(R)>(AS...)> {
    return make_symbolic(op,sym); }



template<class I, class X> class SubscriptableInterface
    : public WritableInterface
{
  public:
    typedef I IndexType;
    typedef I SizesType;
    typedef X ValueType;
  public:
    virtual ~SubscriptableInterface() = default;
    virtual SizesType size() const = 0;
    virtual ValueType operator[] (IndexType i) const = 0;
};

template<class OP, class A1, class A2> decltype(auto) compute_get(ExpressionTemplate<OP,A1,A2> const& expr, MP pr) {
    return expr.op(compute_get(expr.arg1,pr),compute_get(expr.arg2,pr)); }
FloatMPBounds compute_get(Real const& r, MP pr) {
    Effort eff(pr.bits()); return r.compute(eff).get(pr); }



template<class PR> struct ComputeGet {
    Effort _eff; PR _pr;
    FloatBounds<PR> operator() (Real const& r) const { return r.compute(_eff,_pr); }
};

template<class... AS> struct EvaluateAt {
    Tuple<AS...> _args;
    EvaluateAt(AS... as) : _args(as...) { }
    template<class F> decltype(auto) operator() (F const& f) { return std::invoke(f,_args); }
};








template<> class Vector<Real> : public VectorExpression<Vector<Real>> {
    class Interface {
      public:
        virtual ~Interface() = default;
        virtual SizeType size() const = 0;
        virtual Real get(SizeType i) const = 0;
        virtual Vector<FloatMPBounds> compute_get(MP pr) const {
            return Vector<FloatMPBounds>(this->size(),[&](SizeType i){
                Real ri=this->get(i); Effort eff(pr.bits()); ValidatedReal vri=ri.compute(eff); return vri.get(pr); }); }
    };
    SharedPointer<Interface> _ptr;
  public:
    typedef Real ScalarType;

    explicit Vector<Real>(SharedPointer<Interface> ptr) : _ptr(ptr) { }
    Real zero_element() const { return Real(0); }
    SizeType size() const { return _ptr->size(); }
    Real operator[](SizeType i) const { return _ptr->get(i); }
    Vector<FloatMPBounds> compute_get(MP pr) { return _ptr->compute_get(pr); }
    friend OutputStream& operator<<(OutputStream& os, Vector<Real> const& v) {
        os << "["; for(SizeType i=0; i!=v.size(); ++i) { if(i!=0) { os << ","; } os << v[i]; } os << "]"; return os; }

    template<class G> struct Expression : public Vector<Real>::Interface {
        SizeType _n; G _g;
        Expression(SizeType n, G const& g) : _n(n), _g(g) { }
        virtual SizeType size() const final override { return _n; }
        virtual Real get(SizeType i) const final override { return _g(i); }
    };

    template<class G> static Vector<Real> make_real_vector(SizeType n, G const& g) {
        return Vector<Real>(std::make_shared<Expression<G>>(n,g)); }

    Vector<Real>(Vector<Rational> const& v) {
        *this = make_real_vector(v.size(),[v](SizeType i){return Real(v[i]);}); }
    explicit Vector<Real>(SizeType n, Real const& z) {
        *this = make_real_vector(n,[z](SizeType i){return z;}); }
    explicit Vector<Real>(Array<Real> const& ary) {
        *this = make_real_vector(ary.size(),[ary](SizeType i){return ary[i];}); }
    Vector<Real>(InitializerList<Real> const& lst) : Vector<Real>(Array<Real>(lst)) { }
    template<class V> Vector<Real>(VectorExpression<V> const& v) : Vector<Real>(Array<Real>(v().size(),[&](SizeType i){return v()[i];})) { }
    friend Vector<Real> operator+(Vector<Real> const& v1, Vector<Real> const& v2) {
        return make_real_vector(v1.size(),[v1,v2](SizeType i){return v1[i]+v2[i];}); }
    friend Vector<Real> operator-(Vector<Real> const& v1, Vector<Real> const& v2) {
        return make_real_vector(v1.size(),[v1,v2](SizeType i){return v1[i]-v2[i];}); }
    friend Vector<Real> operator*(Real const& s1, Vector<Real> const& v2) {
        return make_real_vector(v2.size(),[s1,v2](SizeType i){return s1*v2[i];}); }
    friend Vector<Real> operator*(Vector<Real> const& v1, Real const& s2) {
        return make_real_vector(v1.size(),[v1,s2](SizeType i){return v1[i]*s2;}); }
    friend Vector<Real> operator/(Vector<Real> const& v1, Real const& s2) {
        return make_real_vector(v1.size(),[v1,s2](SizeType i){return v1[i]/s2;}); }
};


template<> class Matrix<Real> : public MatrixExpression<Matrix<Real>> {
    typedef Pair<SizeType,SizeType> IndexType;
    typedef Pair<SizeType,SizeType> SizesType;
    class Interface {
      public:
        virtual ~Interface() = default;
        virtual SizesType sizes() const = 0;
        virtual Real get(SizeType i, SizeType j) const = 0;
    };
    SharedPointer<Interface> _ptr;
  public:
    typedef Real ScalarType;

    explicit Matrix<Real>(SharedPointer<Interface> ptr) : _ptr(ptr) { }
    Real zero_element() const { return Real(0); }
    SizesType sizes() const { return _ptr->sizes(); }
    SizeType row_size() const { return std::get<0>(_ptr->sizes()); }
    SizeType column_size() const { return std::get<1>(_ptr->sizes()); }
    decltype(auto) operator[](SizeType i) const { return MatrixRow<const Matrix<Real>>(*this,i); }
    Real get(SizeType i, SizeType j) const { return _ptr->get(i,j); }
    Real at(SizeType i, SizeType j) const { return _ptr->get(i,j); }

    template<class G> struct Expression : public RealMatrix::Interface {
        SizesType _n; G _g;
        Expression(SizesType n, G const& g) : _n(n), _g(g) { }
        virtual SizesType sizes() const final override { return _n; }
        virtual Real get(SizeType i, SizeType j) const final override { return _g(i,j); }
    };

    template<class G> static Vector<Real> make_real_vector(SizesType n, G const& g) {
        return Vector<Real>(std::make_shared<RealVector::Expression<G>>(n,g)); }
    template<class G> static Matrix<Real> make_real_matrix(SizesType n, G const& g) {
        return Matrix<Real>(std::make_shared<Expression<G>>(n,g)); }

    Matrix<Real>(Matrix<Rational> const& A) {
        *this = make_real_matrix(make_pair(A.row_size(),A.column_size()),[A](SizeType i, SizeType j){return Real(A[i][j]);}); }
    explicit Matrix<Real>(SizeType m, SizeType n, Real const& z) {
        *this = make_real_matrix(make_pair(m,n),[z](SizeType i, SizeType j){return z;}); }
    explicit Matrix<Real>(Array<Array<Real>> const& ary) {
        *this = make_real_matrix(make_pair(ary.size(),ary[0].size()),[ary](SizeType i, SizeType j){return ary[i][j];}); }
    Matrix<Real>(InitializerList<InitializerList<Real>> const& lst) : Matrix<Real>(Array<Array<Real>>(lst.begin(),lst.end())) { }
    template<class V> Matrix<Real>(MatrixExpression<V> const& A); // : Matrix<Real>(Array<Real>(A().size(),[&](SizeType i){return v()[i];})) { }
    friend Matrix<Real> operator+(Matrix<Real> const& A1, Matrix<Real> const& A2) {
        return make_real_matrix(A1.sizes(),[A1,A2](SizeType i, SizeType j){return A1[i][j]+A2[i][j];}); }
    friend Matrix<Real> operator-(Matrix<Real> const& A1, Matrix<Real> const& A2) {
        return make_real_matrix(A1.sizes(),[A1,A2](SizeType i, SizeType j){return A1[i][j]-A2[i][j];}); }
    friend Matrix<Real> operator*(Matrix<Real> const& A1, Matrix<Real> const& A2) {
        ARIADNE_ASSERT(A1.column_size()==A2.row_size()); ARIADNE_ASSERT(A1.column_size()>0);
        return make_real_matrix( A1.sizes(), [A1,A2](SizeType i, SizeType j){
            Real r=A1[i][0]*A2[0][j]; for (SizeType k=1; k!=A1.column_size(); ++k) { r=r+A1[i][k]*A2[k][j];} return r;} ); }
    friend Matrix<Real> operator*(Real const& s1, Matrix<Real> const& A2) {
        return make_real_matrix(A2.sizes(),[s1,A2](SizeType i, SizeType j){return s1*A2[i][j];}); }
    friend Matrix<Real> operator*(Matrix<Real> const& A1, Real const& s2) {
        return make_real_matrix(A1.sizes(),[A1,s2](SizeType i, SizeType j){return A1[i][j]*s2;}); }
    friend Matrix<Real> operator/(Matrix<Real> const& A1, Real const& s2) {
        return make_real_matrix(A1.sizes(),[A1,s2](SizeType i, SizeType j){return A1[i][j]/s2;}); }
    friend Vector<Real> operator*(Matrix<Real> const& A1, Vector<Real> const& v2) {
        ARIADNE_ASSERT(A1.column_size()==v2.size()); ARIADNE_ASSERT(v2.size()>0);
        return RealVector::make_real_vector( A1.row_size(), [A1,v2](SizeType i){
            Real r=A1[i][0]*v2[0]; for (SizeType k=1; k!=v2.size(); ++k) { r=r+A1[i][k]*v2[k];} return r;} ); }

    friend OutputStream& operator<<(OutputStream& os, Matrix<Real> const& A) {
        for(SizeType i=0; i!=A.row_size(); ++i) {
            os << "\n  [";
            auto Ai=A[i];
            for(SizeType j=0; j!=A.column_size(); ++j) {
                if(j!=0) { os << ","; } os << Ai[j];
            }
            os << "]";
        }
        return os << "\n";
    }

};

} //namespace Ariadne


#endif /* ARIADNE_OBJECT_HPP */
