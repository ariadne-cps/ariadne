/***************************************************************************
 *            function.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file function.hpp
 *  \brief Built-in and user functions and expressions
 */

#ifndef ARIADNE_FUNCTION_HPP
#define ARIADNE_FUNCTION_HPP

#include "../config.hpp"

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "../utility/declarations.hpp"
#include "../utility/macros.hpp"
#include "../utility/pointer.hpp"
#include "../utility/container.hpp"
#include "../utility/metaprogramming.hpp"

#include "../function/function_interface.hpp"

#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/covector.hpp"
#include "../algebra/differential.hpp"
#include "../function/domain.hpp"

namespace Ariadne {

template<class P, class D> struct VectorFunctionElementReference;

template<class T> class Variable;
typedef Variable<Real> RealVariable;
template<class T> class Expression;
typedef Expression<Real> RealExpression;
template<class P, class D, class C> class FunctionExpression;

template<class P>
class FunctionConstructors {
    static_assert(Or<IsSame<P,ApproximateTag>,IsSame<P,ValidatedTag>,IsSame<P,EffectiveTag>>::value,"P must be an information level/paradigm.");
    typedef Number<P> Y;
  public:
    typedef Y NumericType;

    static ScalarFunction<P,BoxDomainType> zero(SizeType as);
    static ScalarFunction<P,BoxDomainType> constant(SizeType as, NumericType c);
    static ScalarFunction<P,BoxDomainType> coordinate(SizeType as, SizeType j);
    static List<ScalarFunction<P,BoxDomainType>> coordinates(SizeType ns);
    static VectorFunction<P,BoxDomainType> zeros(SizeType rs, SizeType as);
    static VectorFunction<P,BoxDomainType> identity(SizeType ns);

    static ScalarFunction<P,BoxDomainType> zero(BoxDomainType dom);
    static ScalarFunction<P,BoxDomainType> constant(BoxDomainType dom, NumericType c);
    static ScalarFunction<P,BoxDomainType> coordinate(BoxDomainType dom, SizeType j);
    static List<ScalarFunction<P,BoxDomainType>> coordinates(BoxDomainType dom);
    static VectorFunction<P,BoxDomainType> zeros(SizeType rs, BoxDomainType dom);
    static VectorFunction<P,BoxDomainType> identity(BoxDomainType dom);

    static ScalarFunction<P,IntervalDomainType> zero();
    static ScalarFunction<P,IntervalDomainType> constant(NumericType c);
    static ScalarFunction<P,IntervalDomainType> coordinate();
    static VectorFunction<P,IntervalDomainType> zeros(SizeType rs);
    static ScalarFunction<P,IntervalDomainType> identity();

    static ScalarFunction<P,IntervalDomainType> zero(IntervalDomainType dom);
    static ScalarFunction<P,IntervalDomainType> constant(IntervalDomainType dom, NumericType c);
    static ScalarFunction<P,IntervalDomainType> coordinate(IntervalDomainType dom, SizeOne);
    static ScalarFunction<P,IntervalDomainType> coordinate(IntervalDomainType dom);
    static VectorFunction<P,IntervalDomainType> zeros(SizeType rs, IntervalDomainType dom);
    static ScalarFunction<P,IntervalDomainType> identity(IntervalDomainType dom);

    static VectorFunction<P,BoxDomainType> constant(BoxDomainType dom, Vector<NumericType> c);
    static VectorFunction<P,IntervalDomainType> constant(IntervalDomainType dom, Vector<NumericType> c);
    static VectorFunction<P,IntervalDomainType> constant(Vector<NumericType> c);

};

template<class P, class X> using EvaluateType = decltype(declval<ScalarMultivariateFunctionInterface<P>>()._evaluate(declval<Vector<X>>()));

template<class P, class D, class C> class FunctionFacade {
};

template<class P> class FunctionFacade<P,IntervalDomainType,IntervalDomainType> {
  public:
    template<class X> Scalar<EvaluateType<P,X>> derivative(X const& x) const;
    FunctionExpression<P,IntervalDomainType,IntervalDomainType> operator() (const RealVariable& x) const;
};

template<class P> class FunctionFacade<P,IntervalDomainType,BoxDomainType> {
  public:
    template<class X> Vector<EvaluateType<P,X>> tangent(X const& x) const;
    FunctionExpression<P,IntervalDomainType,BoxDomainType> operator() (const RealVariable& x) const;
};

template<class P> class FunctionFacade<P,BoxDomainType,IntervalDomainType> {
    typedef Number<P> Y;
  public:
    template<class X> Covector<EvaluateType<P,X>> gradient(Vector<X> const& x) const;
    FunctionExpression<P,BoxDomainType,IntervalDomainType> operator() (const Vector<RealVariable>& x) const;
    //FunctionExpression<P,BoxDomainType,IntervalDomainType> operator() (const Vector<RealExpression>& x) const;
};

template<class P> class FunctionFacade<P,BoxDomainType,BoxDomainType> {
    typedef Number<P> Y;
  public:
    template<class X> Matrix<EvaluateType<P,X>> jacobian(Vector<X> const& x) const;
    FunctionExpression<P,BoxDomainType,BoxDomainType> operator() (const Vector<RealVariable>& x) const;
};

template<class P, class D, class C> class DeclareFunctionOperations;
template<class P, class D> class DeclareFunctionOperations<P,D,IntervalDomainType>
    : DeclareTranscendentalAlgebraOperations<Function<P,D,IntervalDomainType>,Number<P>> { };
template<class P, class D> class DeclareFunctionOperations<P,D,BoxDomainType>
    : DeclareVectorAlgebraOperators<Function<P,D,BoxDomainType>,Function<P,D,IntervalDomainType>,Vector<Number<P>>,Number<P>> { };

template<class P, class D, class C> class DispatchFunctionOperations;
template<class P, class D> class DispatchFunctionOperations<P,D,IntervalDomainType>
    : DispatchTranscendentalAlgebraOperations<Function<P,D,IntervalDomainType>,Number<P>> { };
template<class P, class D> class DispatchFunctionOperations<P,D,BoxDomainType>
    : DeclareVectorAlgebraOperators<Function<P,D,BoxDomainType>,Function<P,D,IntervalDomainType>,Vector<Number<P>>,Number<P>> { };

//! \ingroup FunctionModule
//! \brief A generic scalar function which can be evaluated over the number type \a X,  \f$f:\X^n\rightarrow\X\f$.
template<class P, class D, class C>
class Function
    : public FunctionConstructors<P>
    , public FunctionFacade<P,D,C>
//    , public DeclareFunctionOperations<P,D,C>
    , public DispatchFunctionOperations<P,D,C>
{
    static_assert(IsStronger<P,ApproximateTag>::value,"P must be an information level/paradigm.");
    typedef Number<P> Y;
  protected:
    SharedPointer< const FunctionInterface<P,D,C> > _ptr;
  public:
    typedef P InformationTag;
    typedef P Paradigm;
    typedef Y NumericType;
    typedef D DomainType;
    typedef C CodomainType;
    typedef decltype(declval<C>().dimension()) ResultSizeType;
    typedef decltype(declval<D>().dimension()) ArgumentSizeType;
    typedef decltype(declval<D>().dimension()) ArgumentIndexType;

    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;
    template<class Y> using Result = typename ElementTraits<C>::template Type<Y>;

    explicit Function(EuclideanDomain dom);
    explicit Function(ResultSizeType rs, EuclideanDomain dom);

    explicit Function(DomainType dom);
    explicit Function(ResultSizeType rs, DomainType dom);

    explicit Function(DomainType dom, Result<Formula<Y>>const& e);
    explicit Function(ResultSizeType rs, ScalarFunction<P,D> sf);

    Function(InitializerList<ScalarFunction<P,D>> const& lsf);
    Function(List<ScalarFunction<P,D>> const& lsf);
    Function(Vector<ScalarFunction<P,D>> const& lsf);

    ScalarFunction<P,D> create_zero() const { return ScalarFunction<P,D>::zero(this->domain()); }
    ScalarFunction<P,D> create_constant(NumericType c) const { return ScalarFunction<P,D>::constant(this->domain(),c); }
    ScalarFunction<P,D> create_coordinate(ArgumentIndexType j) const { return ScalarFunction<P,D>::coordinate(this->domain(),j); }
    VectorFunction<P,D> create_constant(Vector<NumericType> c) const { return VectorFunction<P,D>::constant(this->domain(),c); }

    Function();
    explicit Function(FunctionInterface<P,D,C>* p) : _ptr(p) { }
    explicit Function(SharedPointer<FunctionInterface<P,D,C>> p) : _ptr(p) { }
    Function(const FunctionInterface<P,D,C>& t) : _ptr(t._clone()) { }
    Function<P,D,C>& operator=(const FunctionInterface<P,D,C>& f) {
        _ptr=std::shared_ptr< FunctionInterface<P,D,C> >(f._clone()); return *this; }
    Function<P,D,C>& operator=(const Result<NumericType>& c) {
        return (*this)=this->create_constant(c); }

    template<class PP, EnableIf<IsStronger<PP,P>> =dummy>
    Function(const Function<PP,D,C>& f)
        : _ptr(std::dynamic_pointer_cast< const FunctionInterface<P,D,C> >(f.managed_pointer())) { }
    template<class PP, EnableIf<IsStronger<PP,P>> =dummy>
        Function<P,D,C>& operator=(Result<NumericType> const& c); // { return *this=this->create_constant(c); }

    SharedPointer< const FunctionInterface<P,D,C> > managed_pointer() const  { return _ptr; }
    const FunctionInterface<P,D,C>* raw_pointer() const  { return _ptr.operator->(); }
    const FunctionInterface<P,D,C>& reference() const  { return _ptr.operator*(); }
    operator const FunctionInterface<P,D,C>& () const { return _ptr.operator*(); }

    DomainType domain() const {
        return this->reference().domain(); }
    CodomainType codomain() const {
        return this->reference().codomain(); }
    ArgumentSizeType argument_size() const {
        return this->reference().argument_size(); }
    ResultSizeType result_size() const {
        return this->reference().result_size(); }

    template<class X> auto operator() (const Argument<X>& x) const -> decltype(this->reference()._evaluate(x)) {
        return this->reference()._evaluate(x); }
    template<class X> auto evaluate(const Argument<X>& x) const -> decltype(this->reference()._evaluate(x)) {
        return this->reference()._evaluate(x); }

    friend VectorFunction<P,D> operator*(ScalarFunction<P,D> const&, Vector<Y> const&);

    Function<P,D,C> derivative(ElementIndexType<D> k) const {
        return Function<P,D,C>(this->reference()._derivative(k)); }
    friend Function<P,D,C> derivative(Function<P,D,C> const& f, ElementIndexType<D> k) {
        return f.derivative(k); }

    template<class X> decltype(auto) differential(const Argument<X>& x, DegreeType d) const {
        return this->_ptr->_evaluate(Differential<EvaluateType<P,X>>::identity(d,x)); }

    Void set(SizeType i, ScalarFunction<P,D>);
    Function<P,D,IntervalDomainType> get(SizeType i) const;
    const Function<P,D,IntervalDomainType> operator[](SizeType i) const;
    const Function<P,D,BoxDomainType> operator[](Range rng) const;
    VectorFunctionElementReference<P,D> operator[](SizeType i);

    friend OutputStream& operator<<(OutputStream& os, Function<P,D,C> const& f) { f._ptr->write(os); return os; }
};

template<class A> struct AlgebraOperationsBase;

template<class P, class D> struct AlgebraOperationsBase<ScalarFunction<P,D>> {
    template<class OP> static ScalarFunction<P,D> apply(OP op, ScalarFunction<P,D> const& f);
    template<class OP> static ScalarFunction<P,D> apply(OP op, ScalarFunction<P,D> const& f1, ScalarFunction<P,D> const& f2);
    template<class OP> static ScalarFunction<P,D> apply(OP op, ScalarFunction<P,D> const& f1, Number<P> const& c2);
    template<class OP> static ScalarFunction<P,D> apply(OP op, Number<P> const& c1, ScalarFunction<P,D> const& f2);
    template<class OP> static ScalarFunction<P,D> apply(OP op, ScalarFunction<P,D> const& f, Int n);
};

template<class P, class D> struct AlgebraOperations<ScalarFunction<P,D>,Number<P>>
    : public AlgebraOperationsBase<ScalarFunction<P,D>>
{
    using F=ScalarFunction<P,D>;
    using C=Number<P>;
    using Base=AlgebraOperationsBase<ScalarFunction<P,D>>;
    static ScalarFunction<P,D> apply(Pos, ScalarFunction<P,D> const& f);
    static ScalarFunction<P,D> apply(Neg, ScalarFunction<P,D> const& f);
    static ScalarFunction<P,D> apply(Sqr, ScalarFunction<P,D> const& f);
    static ScalarFunction<P,D> apply(Rec, ScalarFunction<P,D> const& f);
    static ScalarFunction<P,D> apply(Add, ScalarFunction<P,D> const& f, Number<P> const& c);
    static ScalarFunction<P,D> apply(Mul, ScalarFunction<P,D> const& f, Number<P> const& c);
    static ScalarFunction<P,D> apply(Add, ScalarFunction<P,D> const& f1, ScalarFunction<P,D> const& f2);
    static ScalarFunction<P,D> apply(Sub, ScalarFunction<P,D> const& f1, ScalarFunction<P,D> const& f2);
    static ScalarFunction<P,D> apply(Mul, ScalarFunction<P,D> const& f1, ScalarFunction<P,D> const& f2);
    static ScalarFunction<P,D> apply(Div, ScalarFunction<P,D> const& f1, ScalarFunction<P,D> const& f2);
    static ScalarFunction<P,D> apply(Pow, ScalarFunction<P,D> const& f, Int n);
    static ScalarFunction<P,D> apply(Sqrt, ScalarFunction<P,D> const& f);
    static ScalarFunction<P,D> apply(Exp, ScalarFunction<P,D> const& f);
    static ScalarFunction<P,D> apply(Log, ScalarFunction<P,D> const& f);
    static ScalarFunction<P,D> apply(Sin, ScalarFunction<P,D> const& f);
    static ScalarFunction<P,D> apply(Cos, ScalarFunction<P,D> const& f);
    static ScalarFunction<P,D> apply(Tan, ScalarFunction<P,D> const& f);
    static ScalarFunction<P,D> apply(Asin, ScalarFunction<P,D> const& f);
    static ScalarFunction<P,D> apply(Acos, ScalarFunction<P,D> const& f);
    static ScalarFunction<P,D> apply(Atan, ScalarFunction<P,D> const& f);
    static ScalarFunction<P,D> apply(Min, ScalarFunction<P,D> const& f1, ScalarFunction<P,D> const& f2);
    static ScalarFunction<P,D> apply(Max, ScalarFunction<P,D> const& f1, ScalarFunction<P,D> const& f2);
    static ScalarFunction<P,D> apply(Abs, ScalarFunction<P,D> const& f);
};

template<class P, class D, class C> inline OutputStream&
operator<<(OutputStream& os, const Function<P,D,C>& f) {
    return f.write(os); }

template<class P, class C, class X> inline decltype(auto)
evaluate(const Function<P,IntervalDomainType,C>& f, const Scalar<X>& x) {
    return f(x); }

template<class P, class C, class X> inline decltype(auto)
evaluate(const Function<P,BoxDomainType,C>& f, const Vector<X>& x) {
    return f(x); }


template<class P, class C, class X> inline decltype(auto)
differential(const Function<P,IntervalDomainType,C>& f, const X& x, DegreeType d) {
    return f.differential(x,d); }

template<class P, class C, class X> inline decltype(auto)
differential(const Function<P,BoxDomainType,C>& f, const Vector<X>& x, DegreeType d) {
    return f.differential(x,d); }

template<class P, class D, class C, EnableIf<IsSame<ElementSizeType<D>,SizeOne>> =dummy> inline
Function<P,D,C> derivative(Function<P,D,C> const& f) {
    return f.derivative(SizeOne()); }

template<class P, class X> Scalar<EvaluateType<P,X>>
derivative(const ScalarUnivariateFunction<P>& f, const Scalar<X>& x) {
    return differential(f,x,1u).gradient(); }

template<class P, class X> Vector<EvaluateType<P,X>>
tangent(const VectorUnivariateFunction<P>& f, const Scalar<X>& x) {
    return column(differential(f,x,1u).jacobian(),0u); }

template<class P, class X> Covector<EvaluateType<P,X>>
gradient(const ScalarMultivariateFunction<P>& f, const Vector<X>& x) {
    return differential(f,x,1u).gradient(); }

template<class P, class X> Matrix<EvaluateType<P,X>>
jacobian(const VectorMultivariateFunction<P>& f, const Vector<X>& x) {
    return differential(f,x,1u).jacobian(); }

template<class P> template<class X> Covector<EvaluateType<P,X>>
FunctionFacade<P,BoxDomainType,IntervalDomainType>::gradient(Vector<X> const& x) const {
    return Ariadne::gradient(static_cast<Function<P,BoxDomainType,IntervalDomainType>const&>(*this),x);
}

template<class P> template<class X> Matrix<EvaluateType<P,X>>
FunctionFacade<P,BoxDomainType,BoxDomainType>::jacobian(Vector<X> const& x) const {
    return Ariadne::jacobian(static_cast<Function<P,BoxDomainType,BoxDomainType>const&>(*this),x);
}


ValidatedScalarMultivariateFunction& operator*=(ValidatedScalarMultivariateFunction& sf, const ExactNumber& c);
EffectiveVectorMultivariateFunction operator*(const EffectiveNumericType& c, const EffectiveVectorMultivariateFunction& vf);

EffectiveScalarMultivariateFunction embed(SizeType as1, const EffectiveScalarMultivariateFunction& f2, SizeType as3);
EffectiveVectorMultivariateFunction embed(SizeType as1, const EffectiveVectorMultivariateFunction& f2, SizeType as3);

EffectiveVectorMultivariateFunction join(const EffectiveScalarMultivariateFunction& f1, const EffectiveScalarMultivariateFunction& f2);
EffectiveVectorMultivariateFunction join(const EffectiveScalarMultivariateFunction& f1, const EffectiveVectorMultivariateFunction& f2);
EffectiveVectorMultivariateFunction join(const EffectiveVectorMultivariateFunction& f1, const EffectiveScalarMultivariateFunction& f2);
EffectiveVectorMultivariateFunction join(const EffectiveVectorMultivariateFunction& f1, const EffectiveVectorMultivariateFunction& f2);

EffectiveScalarMultivariateFunction compose(const EffectiveScalarMultivariateFunction& f, const EffectiveVectorMultivariateFunction& g);
EffectiveVectorMultivariateFunction compose(const EffectiveVectorMultivariateFunction& f, const EffectiveVectorMultivariateFunction& g);

EffectiveScalarMultivariateFunction lie_derivative(const EffectiveScalarMultivariateFunction& g, const EffectiveVectorMultivariateFunction& f);
EffectiveVectorMultivariateFunction lie_derivative(const EffectiveVectorMultivariateFunction& g, const EffectiveVectorMultivariateFunction& f);

Formula<EffectiveNumericType> make_formula(const EffectiveScalarMultivariateFunction& f);
Vector<Formula<EffectiveNumericType>> make_formula(const EffectiveVectorMultivariateFunction& f);
//RealExpression evaluate(EffectiveScalarMultivariateFunction const& f, Vector<RealVariable> const& vars);



ValidatedVectorMultivariateFunction join(const ValidatedScalarMultivariateFunction& f1, const ValidatedScalarMultivariateFunction& f2);
ValidatedVectorMultivariateFunction join(const ValidatedScalarMultivariateFunction& f1, const ValidatedVectorMultivariateFunction& f2);
ValidatedVectorMultivariateFunction join(const ValidatedVectorMultivariateFunction& f1, const ValidatedScalarMultivariateFunction& f2);
ValidatedVectorMultivariateFunction join(const ValidatedVectorMultivariateFunction& f1, const ValidatedVectorMultivariateFunction& f2);
ValidatedScalarMultivariateFunction compose(const ValidatedScalarMultivariateFunction& f, const ValidatedVectorMultivariateFunction& g);
ValidatedVectorMultivariateFunction compose(const ValidatedVectorMultivariateFunction& f, const ValidatedVectorMultivariateFunction& g);


ApproximateVectorMultivariateFunction join(const ApproximateScalarMultivariateFunction& f1, const ApproximateScalarMultivariateFunction& f2);
ApproximateVectorMultivariateFunction join(const ApproximateScalarMultivariateFunction& f1, const ApproximateVectorMultivariateFunction& f2);
ApproximateVectorMultivariateFunction join(const ApproximateVectorMultivariateFunction& f1, const ApproximateScalarMultivariateFunction& f2);
ApproximateVectorMultivariateFunction join(const ApproximateVectorMultivariateFunction& f1, const ApproximateVectorMultivariateFunction& f2);
ApproximateScalarMultivariateFunction compose(const ApproximateScalarMultivariateFunction& f, const ApproximateVectorMultivariateFunction& g);
ApproximateVectorMultivariateFunction compose(const ApproximateVectorMultivariateFunction& f, const ApproximateVectorMultivariateFunction& g);



template<class P, class D>
struct VectorFunctionElementReference
    : DispatchFunctionOperations<P,D,IntervalDomainType>
{
    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;
    typedef IntervalDomainType SC; typedef BoxDomainType VC;
    typedef VectorFunctionElementReference<P,D> SelfType;
    typedef Function<P,D,SC> ElementType;
    typedef typename ElementType::NumericType NumericType;
    Function<P,D,VC>& _vf; SizeType _i;
    VectorFunctionElementReference<P,D>(Function<P,D,VC>& vf, SizeType i) : _vf(vf), _i(i) { }
    template<class WP> operator Function<WP,D,SC> () const;
    Void operator=(const Function<P,D,SC>& sf);
    VectorFunctionElementReference<P,D>& operator=(const VectorFunctionElementReference<P,D>& sfr);
    D domain() const;
    template<class XX> XX evaluate(const Argument<XX> & x) const;
    template<class XX> XX operator()(const Argument<XX>& x) const;
};

template<class P, class D> template<class WP> inline
VectorFunctionElementReference<P,D>::operator Function<WP,D,SC> () const {
    return this->_vf.get(this->_i);
}

template<class P, class D> inline VectorFunctionElementReference<P,D> make_element_reference(ScalarFunction<P,D>& sf, SizeType i) {
    throw std::runtime_error(""); }
template<class P, class D> inline VectorFunctionElementReference<P,D> make_element_reference(VectorFunction<P,D>& vf, SizeType i) {
    return VectorFunctionElementReference<P,D>(vf,i); }

template<class P, class D, class C> inline VectorFunctionElementReference<P,D> Function<P,D,C>::operator[](SizeType i) {
    return make_element_reference(*this,i); }

template<class P, class D, class C> inline ScalarFunction<P,D> const Function<P,D,C>::operator[](SizeType i) const {
    return this->get(i); }

template<class P, class D> inline OutputStream& operator<<(OutputStream& os, const VectorFunctionElementReference<P,D>& vfe) {
    return  os << static_cast< Function<P,D,IntervalDomainType> >(vfe); }

template<class P, class D> inline Void VectorFunctionElementReference<P,D>::operator=(const Function<P,D,SC>& sf) {
    _vf.set(_i,sf); }
template<class P, class D> inline VectorFunctionElementReference<P,D>& VectorFunctionElementReference<P,D>::operator=(const VectorFunctionElementReference<P,D>& sfr) {
    _vf.set(_i,static_cast<Function<P,D,SC>>(sfr)); return *this; }
template<class P, class D> inline D VectorFunctionElementReference<P,D>::domain() const {
    return _vf.domain(); }
template<class P, class D> template<class XX> inline XX VectorFunctionElementReference<P,D>::evaluate(const Argument<XX> & x) const {
    return static_cast<Function<P,D,SC>>(*this).evaluate(x); }
template<class P, class D> template<class XX> inline XX VectorFunctionElementReference<P,D>::operator()(const Argument<XX> & x) const {
    return static_cast<Function<P,D,SC>>(*this).evaluate(x); }



UpperIntervalType evaluate_range(ScalarMultivariateFunction<ValidatedTag>const& f, const Vector<UpperIntervalType>& x);
Vector<UpperIntervalType> evaluate_range(VectorMultivariateFunction<ValidatedTag>const& f, const Vector<UpperIntervalType>& x);
Vector<Differential<UpperIntervalType>> derivative_range(VectorMultivariateFunction<ValidatedTag>const& f, const Vector<Differential<UpperIntervalType>>& x);
Covector<UpperIntervalType> gradient_range(ValidatedScalarMultivariateFunction const& f, const Vector<UpperIntervalType>& x);
Matrix<UpperIntervalType> jacobian_range(ValidatedVectorMultivariateFunction const& f, const Vector<UpperIntervalType>& x);



/*
inline Matrix<UpperIntervalType> jacobian(VectorMultivariateFunction<ValidatedTag>const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<Matrix<UpperIntervalType>>(f.jacobian(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
inline Matrix<UpperIntervalType> jacobian(VectorMultivariateFunction<ValidatedTag>const& f, const Vector<ExactInterval>& x) {
    return static_cast<Matrix<UpperIntervalType>>(f.jacobian(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
inline Matrix<UpperIntervalType> jacobian_range(VectorMultivariateFunction<ValidatedTag>const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<Matrix<UpperIntervalType>>(f.jacobian(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }

// FIXME: Needed to override templated gradient and jacobian
inline Covector<UpperIntervalType> gradient(ScalarMultivariateFunction<EffectiveTag>const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<Covector<UpperIntervalType>>(gradient(f,reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
inline Covector<UpperIntervalType> gradient(ScalarMultivariateFunction<EffectiveTag>const& f, const Vector<ExactInterval>& x) {
    return gradient(f,static_cast<Vector<UpperIntervalType>>(x)); }
inline Matrix<UpperIntervalType> jacobian(VectorMultivariateFunction<EffectiveTag>const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<Matrix<UpperIntervalType>>(f.jacobian(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
inline Matrix<UpperIntervalType> jacobian(VectorMultivariateFunction<EffectiveTag>const& f, const Vector<ExactInterval>& x) {
    return jacobian(f,static_cast<Vector<UpperIntervalType>>(x)); }
*/

template<class P> class FunctionFactory;
typedef FunctionFactory<ValidatedTag> ValidatedFunctionFactory;

template<>
class FunctionFactory<ValidatedTag>
{
    SharedPointer< const FunctionFactoryInterface<ValidatedTag> > _ptr;
  public:
    FunctionFactory(const FunctionFactoryInterface<ValidatedTag>& ref) : _ptr(ref.clone()) { }
    FunctionFactory(const FunctionFactoryInterface<ValidatedTag>* ptr) : _ptr(ptr) { }
    FunctionFactory(SharedPointer< const FunctionFactoryInterface<ValidatedTag> > ptr) : _ptr(ptr) { }
    inline ValidatedScalarMultivariateFunction create(const BoxDomainType& d, const ValidatedScalarMultivariateFunctionInterface& f) const;
    inline ValidatedVectorMultivariateFunction create(const BoxDomainType& d, const ValidatedVectorMultivariateFunctionInterface& f) const;
    inline ValidatedScalarMultivariateFunction create_zero(const BoxDomainType& d) const;
    inline ValidatedVectorMultivariateFunction create_identity(const BoxDomainType& d) const;
    friend OutputStream& operator<<(OutputStream& os, const FunctionFactory<ValidatedTag>& factory);
};

inline ValidatedScalarMultivariateFunction FunctionFactoryInterface<ValidatedTag>::create(const BoxDomainType& domain, const ValidatedScalarMultivariateFunctionInterface& function) const {
    return ValidatedScalarMultivariateFunction(SharedPointer<ValidatedScalarMultivariateFunctionInterface>(this->_create(domain,function))); }
inline ValidatedVectorMultivariateFunction FunctionFactoryInterface<ValidatedTag>::create(const BoxDomainType& domain, const ValidatedVectorMultivariateFunctionInterface& function) const {
    return ValidatedVectorMultivariateFunction(SharedPointer<ValidatedVectorMultivariateFunctionInterface>(this->_create(domain,function))); }
inline ValidatedScalarMultivariateFunction FunctionFactoryInterface<ValidatedTag>::create_zero(const BoxDomainType& domain) const {
    return this->create(domain,EffectiveScalarMultivariateFunction::zero(domain.dimension())); }
inline ValidatedVectorMultivariateFunction FunctionFactoryInterface<ValidatedTag>::create_identity(const BoxDomainType& domain) const {
    return this->create(domain,EffectiveVectorMultivariateFunction::identity(domain.dimension())); }

inline ValidatedScalarMultivariateFunction FunctionFactory<ValidatedTag>::create(const BoxDomainType& domain, const ValidatedScalarMultivariateFunctionInterface& function) const {
    return this->_ptr->create(domain,function); }
inline ValidatedVectorMultivariateFunction FunctionFactory<ValidatedTag>::create(const BoxDomainType& domain, const ValidatedVectorMultivariateFunctionInterface& function) const {
    return this->_ptr->create(domain,function); }
inline ValidatedScalarMultivariateFunction FunctionFactory<ValidatedTag>::create_zero(const BoxDomainType& domain) const {
    return this->_ptr->create_zero(domain); }
inline ValidatedVectorMultivariateFunction FunctionFactory<ValidatedTag>::create_identity(const BoxDomainType& domain) const {
    return this->_ptr->create_identity(domain); }

inline OutputStream& operator<<(OutputStream& os, const ValidatedFunctionFactory& factory) {
    factory._ptr->write(os); return os; }

} // namespace Ariadne

#endif
