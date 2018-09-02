/***************************************************************************
 *            algebra_mixin.hpp
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

/*! \file algebra_mixin.hpp
 *  \brief Mixin class providing operations for (Banach) algebras.
 */

#ifndef ARIADNE_ALGEBRA_MIXIN_HPP
#define ARIADNE_ALGEBRA_MIXIN_HPP

#include "../algebra/algebra_interface.hpp"
#include "../numeric/operators.hpp"

namespace Ariadne {

template<class A, class X> class AlgebraMixin
    : public virtual AlgebraInterface<X>
{
    typedef X NumericType;
  private:
    template<class OP> static AlgebraInterface<X>* _eval(OP op, AlgebraMixin<A,X> const& am1, AlgebraInterface<X> const& ai2) {
        AlgebraMixin<A,X>const* amp2 = dynamic_cast<AlgebraMixin<A,X>const*>(&ai2); assert(amp2);
        A const& a1=static_cast<A const&>(am1); A const& a2=static_cast<A const&>(*amp2);
        return new A(op(a1,a2)); }
    template<class OP> static AlgebraInterface<X>* _eval(OP op, AlgebraMixin<A,X> const& am, X const& c) {
        A const& a=static_cast<A const&>(am); return new A(op(a,c)); }
    template<class OP> static AlgebraInterface<X>* _eval(OP op, AlgebraMixin<A,X> const& am) {
        A const& a=static_cast<A const&>(am); return new A(op(a)); }
  public:
    virtual AlgebraInterface<X>* _create_zero() const { return new A(static_cast<const A&>(*this).A::create_zero()); }
    virtual AlgebraInterface<X>* _create_constant(X const& c) const { return new A(static_cast<const A&>(*this).A::create_constant(c)); }
    virtual AlgebraInterface<X>* _create_copy() const { return new A(static_cast<const A&>(*this)); }
    virtual AlgebraInterface<X>* _neg() { return new AlgebraMixin<A,X>(-static_cast<A const&>(*this)); }
    virtual AlgebraInterface<X>* _add(AlgebraInterface<X> const& other) const { return _eval(Add(),*this,other); }
    virtual AlgebraInterface<X>* _sub(AlgebraInterface<X> const& other) const { return _eval(Sub(),*this,other); }
    virtual AlgebraInterface<X>* _mul(AlgebraInterface<X> const& other) const { return _eval(Mul(),*this,other); }
    virtual AlgebraInterface<X>* _add(X const& cnst) const { return _eval(Add(),*this,cnst); }
    virtual AlgebraInterface<X>* _sub(X const& cnst) const { return _eval(Sub(),*this,cnst); }
    virtual AlgebraInterface<X>* _mul(X const& cnst) const { return _eval(Mul(),*this,cnst); }
    virtual AlgebraInterface<X>* _div(X const& cnst) const { return _eval(Div(),*this,cnst); }
    virtual AlgebraInterface<X>* _radd(X const& cnst) const { return _eval(RAdd(),*this,cnst); }
    virtual AlgebraInterface<X>* _rsub(X const& cnst) const { return _eval(RSub(),*this,cnst); }
    virtual AlgebraInterface<X>* _rmul(X const& cnst) const { return _eval(RMul(),*this,cnst); }
    virtual Void _iadd(const X& c) { static_cast<A*>(this)->A::iadd(c); }
    virtual Void _imul(const X& c) { static_cast<A*>(this)->A::imul(c); }
    virtual Void _isma(const X& c, const AlgebraInterface<X>& x) {
        static_cast<A*>(this)->A::isma(c,dynamic_cast<const A&>(x)); }
    virtual Void _ifma(const AlgebraInterface<X>& x1, const AlgebraInterface<X>& x2)  {
        static_cast<A*>(this)->A::ifma(dynamic_cast<const A&>(x1),dynamic_cast<const A&>(x2)); }
    virtual OutputStream& write(OutputStream& os) const { os << static_cast<const A&>(*this); return os; }

    virtual AlgebraInterface<X>* _apply(Neg op) const { return _eval(op,*this); }
    virtual AlgebraInterface<X>* _apply(Add op, AlgebraInterface<X>const& other) const { return _eval(op,*this,other); }
    virtual AlgebraInterface<X>* _apply(Sub op, AlgebraInterface<X>const& other) const { return _eval(op,*this,other); }
    virtual AlgebraInterface<X>* _apply(Mul op, AlgebraInterface<X>const& other) const { return _eval(op,*this,other); }
    virtual AlgebraInterface<X>* _apply(Add op, X const& cnst) const { return _eval(op,*this,cnst); }
    virtual AlgebraInterface<X>* _apply(Mul op, X const& cnst) const { return _eval(op,*this,cnst); }
};

template<class A, class X> class NormedAlgebraMixin
    : public virtual NormedAlgebraInterface<X>
    , public AlgebraMixin<A,X>
{
    virtual NormedAlgebraInterface<X>* _create_ball(ErrorType r) const { return new A(static_cast<const A&>(*this).A::create_ball(r)); }
    virtual NormedAlgebraInterface<X>* _create_constant(X c) const { return new A(static_cast<const A&>(*this).A::create_constant(c)); }
    virtual NormedAlgebraInterface<X>* _create_zero() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual NormedAlgebraInterface<X>* _create_copy() const { return new A(static_cast<const A&>(*this)); }
};

template<class A, class X> class GradedAlgebraMixin
    : public virtual GradedAlgebraInterface<X>
    , public AlgebraMixin<A,X>
{
    virtual GradedAlgebraMixin<A,X>* _create_zero() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual GradedAlgebraMixin<A,X>* _create_copy() const { return new A(static_cast<const A&>(*this)); }
    virtual GradedAlgebraMixin<A,X>* _apply(const Series<X>& f) const { return new A(compose(f,static_cast<const A&>(*this))); }
};

template<class A, class X> class SymbolicAlgebraMixin
    : public virtual SymbolicAlgebraInterface<X>
    , public AlgebraMixin<A,X>
{
    virtual SymbolicAlgebraMixin<A,X>* _create_copy() const { return new A(static_cast<const A&>(*this)); }
    virtual SymbolicAlgebraMixin<A,X>* _create_zero() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual SymbolicAlgebraMixin<A,X>* _create_constant(X const& c) const { return new A(static_cast<const A&>(*this).A::create_constant(c)); }
    virtual SymbolicAlgebraMixin<A,X>* _apply(OperatorCode op) { return new A(op,static_cast<const A&>(*this)); }
};

} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_MIXIN_HPP */
