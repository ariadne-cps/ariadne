/***************************************************************************
 *            function/function_model.hpp
 *
 *  Copyright  2011-20  Pieter Collins
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

/*! \file function/function_model.hpp
 *  \brief Built-in and user functions and expressions
 */

#ifndef ARIADNE_FUNCTION_MODEL_HPP
#define ARIADNE_FUNCTION_MODEL_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "../function/function.decl.hpp"
#include "../function/function_model_interface.hpp"

#include "../numeric/operators.hpp"
#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/operations.hpp"
#include "../function/domain.hpp"

#include "../function/function_interface.hpp"
#include "../function/function_mixin.hpp"
#include "../function/function.hpp"

namespace Ariadne {

template<class P, class D, class PR, class PRE> struct AlgebraOperations<ScalarFunctionModel<P,D,PR,PRE>>;

// FIXME: Extend with univariate case
template<class P, class PR, class PRE> class FunctionModelFactory {
    SharedPointer<const FunctionModelFactoryInterface<P,PR,PRE>> _ptr;
    typedef IntervalDomainType SD;
    typedef BoxDomainType VD;
  public:
    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;
    typedef SD ScalarDomainType;
    typedef VD VectorDomainType;

    operator const FunctionModelFactoryInterface<P,PR,PRE>& () const { return *_ptr; }

    explicit FunctionModelFactory(const FunctionModelFactoryInterface<P,PR,PRE>* p) : _ptr(p) { }
    explicit FunctionModelFactory(SharedPointer<const FunctionModelFactoryInterface<P,PR,PRE>> p) : _ptr(p) { }

    CanonicalNumericType<P,PR,PRE> create(Number<P> const& c) const {
        return CanonicalNumericType<P,PR,PRE>(this->_ptr->_create(c)); }
    ScalarFunctionModel<P,VD,PR,PRE> create(VectorDomainType const& dom, ScalarFunction<P,VD> const& f) const {
        return ScalarFunctionModel<P,VD,PR,PRE>(this->_ptr->_create(dom,f)); }
    VectorFunctionModel<P,VD,PR,PRE> create(VectorDomainType const& dom, VectorFunction<P,VD> const& f) const {
        return VectorFunctionModel<P,VD,PR,PRE>(this->_ptr->_create(dom,f)); }

    CanonicalNumericType<P,PR,PRE> create_number(Number<P> const& c) const {
        return CanonicalNumericType<P,PR,PRE>(this->_ptr->_create(c)); }
    ScalarFunctionModel<P,VD,PR,PRE> create_zero(VectorDomainType const& dom) const {
        return ScalarFunctionModel<P,VD,PR,PRE>(this->_ptr->_create_zero(dom)); }
    ScalarFunctionModel<P,VD,PR,PRE> create_constant(VectorDomainType const& dom, Number<P> const& c) const {
        return ScalarFunctionModel<P,VD,PR,PRE>(this->_ptr->_create_constant(dom,c)); }
    ScalarFunctionModel<P,VD,PR,PRE> create_constant(VectorDomainType const& dom, CanonicalNumericType<P,PR,PRE> const& c) const {
        return ScalarFunctionModel<P,VD,PR,PRE>(this->_ptr->_create_constant(dom,c)); }
    ScalarFunctionModel<P,VD,PR,PRE> create_coordinate(VectorDomainType const& dom, SizeType index) const {
        return ScalarFunctionModel<P,VD,PR,PRE>(this->_ptr->_create_coordinate(dom,index)); }
    VectorFunctionModel<P,VD,PR,PRE> create_zeros(SizeType rsize, VectorDomainType const& dom) const {
        return VectorFunctionModel<P,VD,PR,PRE>(this->_ptr->_create_zeros(rsize,dom)); }
    VectorFunctionModel<P,VD,PR,PRE> create_projection(VectorDomainType const& dom, Range indices) const {
        return ScalarFunctionModel<P,VD,PR,PRE>(this->_ptr->_create_projection(dom,indices)); }
    VectorFunctionModel<P,VD,PR,PRE> create_identity(VectorDomainType const& dom) const {
        return VectorFunctionModel<P,VD,PR,PRE>(this->_ptr->_create_identity(dom)); }

    // FIXME: Should return a univariate model
    ScalarFunctionModel<P,VD,PR,PRE> create_identity(ScalarDomainType const& dom) const {
        return ScalarFunctionModel<P,VD,PR,PRE>(this->_ptr->_create_coordinate(VectorDomainType(1u,dom),0u)); }
    VectorFunctionModel<P,SD,PR,PRE> create_zeros(SizeType n, ScalarDomainType const& dom) const {
        ARIADNE_NOT_IMPLEMENTED; }

    friend OutputStream& operator<<(OutputStream& os, FunctionModelFactory<P,PR,PRE> const& factory) { return factory._ptr->_write(os); }
};

template<class P, class PR, class PRE> auto inline
FunctionModelFactoryInterface<P,PR,PRE>::create_identity(IntervalDomainType const& dom) const -> ScalarFunctionModel<P,VD,PR,PRE> {
    return ScalarFunctionModel<P,VD,PR,PRE>(this->create_coordinate(BoxDomainType(1u,dom),0u)); }

template<class FCTRY, class D> class FunctionModelCreator {
    typedef typename FCTRY::Paradigm P;
    typedef typename FCTRY::PrecisionType PR;
    typedef typename FCTRY::ErrorPrecisionType PRE;
  public:
    typedef FCTRY FactoryType;
    typedef D DomainType;
    typedef P Paradigm;

    explicit FunctionModelCreator(DomainType domain, FactoryType factory) : _factory(factory), _domain(domain) { }

    decltype(auto) create(Number<P> const& c) const { return this->_factory.create(c); }
    decltype(auto) create(ScalarFunction<P,D> const& f) { return this->_factory.create(this->_domain,f); }
    decltype(auto) create(VectorFunction<P,D> const& f) { return this->_factory.create(this->_domain,f); }
    decltype(auto) create_zero() { return this->_factory.create_zero(this->_domain); }
    decltype(auto) create_zeros(SizeType n) { return this->_factory.create_zeros(n,this->_domain); }
    decltype(auto) create_constant(Number<P> const& c) const { return this->_factory.create_constant(this->_domain,c); }
    decltype(auto) create_identity() { return this->_factory.create_identity(this->_domain); }

//    ScalarFunctionModel<P,D,PR,PRE> const& create(ScalarFunctionModel<P,D,PR,PRE> const& f) const { return f; }
    CanonicalNumericType<P,PR,PRE> const& create(CanonicalNumericType<P,PR,PRE> const& c) const { return c; }
  protected:
    FactoryType _factory;
    DomainType _domain;
};

// FIXME: Merge with multivariate case
template<class FCTRY> class FunctionModelCreator<FCTRY,IntervalDomainType> {
    typedef typename FCTRY::Paradigm P;
    typedef typename FCTRY::PrecisionType PR;
    typedef typename FCTRY::ErrorPrecisionType PRE;
    typedef IntervalDomainType D;
  public:
    typedef FCTRY FactoryType;
    typedef D DomainType;
    typedef P Paradigm;

    explicit FunctionModelCreator(DomainType domain, FactoryType factory) : _factory(factory), _domain(domain) { }

    CanonicalNumericType<P,PR,PRE> create(Number<P> const& c) const { return this->_factory.create(c); }
    ScalarFunctionModel<P,D,PR,PRE> create(ScalarFunction<P,D> const& f) { ARIADNE_NOT_IMPLEMENTED; }
    VectorFunctionModel<P,D,PR,PRE> create(VectorFunction<P,D> const& f) { ARIADNE_NOT_IMPLEMENTED; }
    ScalarFunctionModel<P,D,PR,PRE> create_zero() { ARIADNE_NOT_IMPLEMENTED; }
    VectorFunctionModel<P,D,PR,PRE> create_zeros(SizeType n) { ARIADNE_NOT_IMPLEMENTED; }
    ScalarFunctionModel<P,D,PR,PRE> create_identity() { ARIADNE_NOT_IMPLEMENTED; }

    CanonicalNumericType<P,PR,PRE> const& create(CanonicalNumericType<P,PR,PRE> const& c) const { return c; }
  protected:
    FactoryType _factory;
    DomainType _domain;
};


//! \ingroup FunctionModule
//! \ingroup FunctionModelSubModule
//! \brief Generic functions on bounded domains.
template<class P, class D, class C, class PR, class PRE> class FunctionModel;

//! \ingroup FunctionModelSubModule
//! \brief Generic scalar functions on bounded domains.
template<class P, class D, class PR, class PRE> class FunctionModel<P,D,IntervalDomainType,PR,PRE>
//    : public DispatchTranscendentalAlgebraOperations<ScalarFunctionModel<P,D,PR,PRE>, CanonicalNumericType<P,PR,PRE>>
    : public DispatchElementaryAlgebraOperations<ScalarFunctionModel<P,D,PR,PRE>, CanonicalNumericType<P,PR,PRE>>
    , public ProvideConcreteGenericElementaryOperations<ScalarFunctionModel<P,D,PR,PRE>,ScalarMultivariateFunction<P>>
    , public ProvideConcreteGenericElementaryOperations<ScalarFunctionModel<P,D,PR,PRE>,Number<P>>
{
    static_assert(IsSame<D,IntervalDomainType>::value or IsSame<D,BoxDomainType>::value,"");
    using C = IntervalDomainType;
  public:
    typedef ScalarFunction<P,D> GenericType;
    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef D DomainType;
    typedef C CodomainType;
    typedef typename FunctionModelInterface<P,D,C,PR,PRE>::ValueType ValueType;
    typedef typename FunctionModelInterface<P,D,C,PR,PRE>::ErrorType ErrorType;
    typedef typename FunctionModelInterface<P,D,C,PR,PRE>::NumericType NumericType;
    typedef typename FunctionModelInterface<P,D,C,PR,PRE>::GenericNumericType GenericNumericType;
    typedef typename FunctionModelInterface<P,D,C,PR,PRE>::NormType NormType;
    typedef typename FunctionModelInterface<P,D,C,PR,PRE>::RangeType RangeType;

    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;
    template<class Y> using Result = ElementTraits<C>::template Type<Y>;
  public:
    clone_on_copy_ptr< ScalarFunctionModelInterface<P,D,PR,PRE> > _ptr;
  public:
    FunctionModel() : _ptr() { }
    explicit FunctionModel(FunctionModelInterface<P,D,C,PR,PRE>* p) : _ptr(p) { }
    FunctionModel(const SharedPointer<const FunctionModelInterface<P,D,C,PR,PRE>> p) : _ptr(p->_clone()) { }
    FunctionModel(const FunctionModel<P,D,C,PR,PRE>& f) : _ptr(f._ptr) { }
    FunctionModel& operator=(const FunctionModel<P,D,C,PR,PRE>& f) { this->_ptr=f._ptr; return *this; }
        FunctionModel(const FunctionModelInterface<P,D,C,PR,PRE>& f) : _ptr(f._clone()) { }
    FunctionModel(const Function<P,D,C>& f) : _ptr(dynamic_cast<FunctionModelInterface<P,D,C,PR,PRE>*>(f.raw_pointer()->_clone())) { }
    operator Function<P,D,C>() const { return Function<P,D,C>(this->_ptr->_clone()); }
    operator FunctionModelInterface<P,D,C,PR,PRE>& () { return *_ptr; }
    operator const FunctionModelInterface<P,D,C,PR,PRE>& () const { return *_ptr; }
    const FunctionModelInterface<P,D,C,PR,PRE>* raw_pointer() const { return _ptr.operator->(); }
    FunctionModelInterface<P,D,C,PR,PRE>& reference() { return *_ptr; }
    const FunctionModelInterface<P,D,C,PR,PRE>& reference() const { return *_ptr; }

    ScalarFunctionModel<P,D,PR,PRE>& operator=(const Number<P>& c);
    ScalarFunctionModel<P,D,PR,PRE>& operator=(const CanonicalNumericType<P,PR,PRE>& c);
    ScalarFunctionModel<P,D,PR,PRE>& operator=(const ScalarFunction<P,D>& f);
    ScalarFunctionModel<P,D,PR,PRE>& operator=(const ScalarFunctionModelInterface<P,D,PR,PRE>& f);
//    ScalarFunctionModel<P,D,PR,PRE>& operator=(const ValidatedScalarMultivariateTaylorFunctionModelDP& f);

    inline PrecisionType const precision() const { return this->value().precision(); }
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    template<class X> X operator() (const Vector<X>& x) const {
        return this->_ptr->_evaluate(x); }
    template<class X> X evaluate(const Vector<X>& x) const {
        return this->_ptr->_evaluate(x); }
    inline DomainType const domain() const { return this->_ptr->domain(); }
    inline CodomainType const codomain() const { return this->_ptr->codomain(); }
    inline RangeType const range() const { return this->_ptr->range(); }

    inline ValueType value() const { return this->_ptr->_value(); }
    inline ErrorType error() const { return this->_ptr->_error(); }

    inline Void clobber() { return this->_ptr->clobber(); }

//    inline ScalarFunctionModel<P,D,PR,PRE> apply(UnaryElementaryOperator op) const { return ScalarFunctionModel<P,D,PR,PRE>(this->_ptr->_apply(op)); }
    inline Void restrict(const DomainType& d) { *this=restriction(*this,d); }
  public:
    friend FunctionModelCreator<FunctionModelFactory<P,PR,PRE>,D> factory(ScalarFunctionModel<P,D,PR,PRE> const& f) {
        FunctionModelFactory<P,PR,PRE> factory(f._ptr->_factory()); return FunctionModelCreator<FunctionModelFactory<P,PR,PRE>,D>(f.domain(),factory); }
  public:
  public:
    friend CanonicalNumericType<P,PR,PRE> evaluate(const ScalarFunctionModel<P,D,PR,PRE>& f, const Vector<CanonicalNumericType<P,PR,PRE>>& x) {
        return f._ptr->_evaluate(x); }
    friend Number<P> evaluate(const ScalarFunctionModel<P,D,PR,PRE>& f, const Vector<Number<P>>& x) {
        return f._ptr->_evaluate(Vector<CanonicalNumericType<P,PR,PRE>>(x,f.precision())); }
    friend CanonicalNumericType<P,PR,PRE> unchecked_evaluate(const ScalarFunctionModel<P,D,PR,PRE>& f, const Vector<CanonicalNumericType<P,PR,PRE>>& x) {
        return f._ptr->_unchecked_evaluate(x); }
    friend Number<P> unchecked_evaluate(const ScalarFunctionModel<P,D,PR,PRE>& f, const Vector<Number<P>>& x) {
        return f._ptr->_unchecked_evaluate(Vector<CanonicalNumericType<P,PR,PRE>>(x,f.precision())); }

    friend ScalarFunctionModel<P,D,PR,PRE> partial_evaluate(const ScalarFunctionModel<P,D,PR,PRE>& f, SizeType j, const CanonicalNumericType<P,PR,PRE>& c) {
        return ScalarFunctionModel<P,D,PR,PRE>(f._ptr->_partial_evaluate(j,c)); }
    friend ScalarFunctionModel<P,D,PR,PRE> partial_evaluate(const ScalarFunctionModel<P,D,PR,PRE>& f, SizeType j, const Number<P>& c) {
        return partial_evaluate(f,j,CanonicalNumericType<P,PR,PRE>(c,f.precision())); }

    friend NormType norm(const ScalarFunctionModel<P,D,PR,PRE>& f) {
        return f._ptr->_norm(); }
    friend ScalarFunctionModel<P,D,PR,PRE> derivative(const ScalarFunctionModel<P,D,PR,PRE>& f, SizeType j) {
        return ScalarFunctionModel<P,D,PR,PRE>(f._ptr->_derivative(j)); }
    friend ScalarFunctionModel<P,D,PR,PRE> antiderivative(const ScalarFunctionModel<P,D,PR,PRE>& f, SizeType j) {
        return ScalarFunctionModel<P,D,PR,PRE>(f._ptr->_antiderivative(j)); }
    friend ScalarFunctionModel<P,D,PR,PRE> antiderivative(const ScalarFunctionModel<P,D,PR,PRE>& f, SizeType j, CanonicalNumericType<P,PR,PRE> c) {
        return ScalarFunctionModel<P,D,PR,PRE>(f._ptr->_antiderivative(j,c)); }
    friend ScalarFunctionModel<P,D,PR,PRE> antiderivative(const ScalarFunctionModel<P,D,PR,PRE>& f, SizeType j, const Number<P>& c) {
        return antiderivative(f,j,CanonicalNumericType<P,PR,PRE>(c,f.value().precision())); }

    friend ScalarFunctionModel<P,D,PR,PRE> embed(const DomainType& d1, const ScalarFunctionModel<P,D,PR,PRE>& f, const DomainType& d2) {
        return ScalarFunctionModel<P,D,PR,PRE>(f._ptr->_embed(d1,d2)); }
    friend ScalarFunctionModel<P,D,PR,PRE> embed(const DomainType& d, const ScalarFunctionModel<P,D,PR,PRE>& f) {
        return embed(d,f,DomainType()); }
    friend ScalarFunctionModel<P,D,PR,PRE> embed(const ScalarFunctionModel<P,D,PR,PRE>& f, const BoxDomainType& d) {
        return embed(DomainType(),f,d); }
    friend ScalarFunctionModel<P,D,PR,PRE> embed(const ScalarFunctionModel<P,D,PR,PRE>& f, const IntervalDomainType& d) {
        return embed(f,DomainType(1,d)); }
    friend ScalarFunctionModel<P,D,PR,PRE> restrict(const ScalarFunctionModel<P,D,PR,PRE>& f, const DomainType& d) {
        return ScalarFunctionModel<P,D,PR,PRE>(f._ptr->_restriction(d)); }
    friend ScalarFunctionModel<P,D,PR,PRE> restriction(const ScalarFunctionModel<P,D,PR,PRE>& f, const DomainType& d) {
        return ScalarFunctionModel<P,D,PR,PRE>(f._ptr->_restriction(d)); }

    friend VectorFunctionModel<P,D,PR,PRE> join(const ScalarFunctionModel<P,D,PR,PRE>& f1, const ScalarFunctionModel<P,D,PR,PRE>& f2) {
        return join(VectorFunctionModel<P,D,PR,PRE>(1,f1),f2); }
    friend VectorFunctionModel<P,D,PR,PRE> combine(const ScalarFunctionModel<P,D,PR,PRE>& f1, const ScalarFunctionModel<P,D,PR,PRE>& f2);
  public:
    typedef ValidatedTag VP;
    friend ScalarFunctionModel<VP,D,PR,PRE> refinement(const ScalarFunctionModel<VP,D,PR,PRE>& f1, const ScalarFunctionModel<VP,D,PR,PRE>& f2) {
        return ScalarFunctionModel<VP,D,PR,PRE>(f1._ptr->_refinement(f2)); }
    friend Boolean inconsistent(const ScalarFunctionModel<VP,D,PR,PRE>& f1, const ScalarFunctionModel<VP,D,PR,PRE>& f2) {
        return f1._ptr->_inconsistent(f2); }
    friend Boolean refines(const ScalarFunctionModel<VP,D,PR,PRE>& f1, const ScalarFunctionModel<VP,D,PR,PRE>& f2) {
        return f1._ptr->_refines(f2); }
  public:
    friend OutputStream& operator<<(OutputStream& os, const ScalarFunctionModel<P,D,PR,PRE>& f) {
        return os <<  f.operator ScalarMultivariateFunction<P>(); }
};

template<class P, class D, class PR, class PRE> struct AlgebraOperations<ScalarFunctionModel<P,D,PR,PRE>> {

/*
    static ScalarFunctionModel<P,D,PR,PRE> apply(Nul, ScalarFunctionModel<P,D,PR,PRE> f) {
        f._ptr->_imul(CanonicalNumericType<P,PR,PRE>(0)); return f; }
    static ScalarFunctionModel<P,D,PR,PRE> apply(Neg, ScalarFunctionModel<P,D,PR,PRE> f) {
        f._ptr->_imul(CanonicalNumericType<P,PR,PRE>(-1)); return f; }
    static ScalarFunctionModel<P,D,PR,PRE> apply(Add, ScalarFunctionModel<P,D,PR,PRE> f1, const ScalarFunctionModel<P,D,PR,PRE>& f2) {
        f1._ptr->_isma(CanonicalNumericType<P,PR,PRE>(+1),f2); return f1; }
    static ScalarFunctionModel<P,D,PR,PRE> apply(Sub, ScalarFunctionModel<P,D,PR,PRE> f1, const ScalarFunctionModel<P,D,PR,PRE>& f2) {
        f1._ptr->_isma(CanonicalNumericType<P,PR,PRE>(-1),f2); return f1; }
    static ScalarFunctionModel<P,D,PR,PRE> apply(Mul, const ScalarFunctionModel<P,D,PR,PRE>& f1, const ScalarFunctionModel<P,D,PR,PRE>& f2) {
        ScalarFunctionModel<P,D,PR,PRE> r=factory(f1).create_zero(); r._ptr->_ifma(f1,f2); return r; }
    static ScalarFunctionModel<P,D,PR,PRE> apply(Add, ScalarFunctionModel<P,D,PR,PRE> f1, const CanonicalNumericType<P,PR,PRE>& c2) {
        f1._ptr->_iadd(c2); return f1; }
    static ScalarFunctionModel<P,D,PR,PRE> apply(Mul, ScalarFunctionModel<P,D,PR,PRE> f1, const CanonicalNumericType<P,PR,PRE>& c2) {
        f1._ptr->_imul(c2); return f1; }
    static ScalarFunctionModel<P,D,PR,PRE> apply(Add, const CanonicalNumericType<P,PR,PRE>& c1, ScalarFunctionModel<P,D,PR,PRE> f2) {
        f2._ptr->_iadd(c1); return f2; }
    static ScalarFunctionModel<P,D,PR,PRE> apply(Mul, const CanonicalNumericType<P,PR,PRE>& c1, ScalarFunctionModel<P,D,PR,PRE> f2) {
        f2._ptr->_imul(c1); return f2; }

    static ScalarFunctionModel<P,D,PR,PRE> apply(Sub, ScalarFunctionModel<P,D,PR,PRE> f1, const CanonicalNumericType<P,PR,PRE>& c2) {
        return add(std::move(f1),neg(c2)); }
    static ScalarFunctionModel<P,D,PR,PRE> apply(Div, ScalarFunctionModel<P,D,PR,PRE> f1, const CanonicalNumericType<P,PR,PRE>& c2) {
        return mul(std::move(f1),rec(c2)); }
    static ScalarFunctionModel<P,D,PR,PRE> apply(Sub, const CanonicalNumericType<P,PR,PRE>& c1, ScalarFunctionModel<P,D,PR,PRE> f2) {
        return add(neg(std::move(f2)),c1); }
    static ScalarFunctionModel<P,D,PR,PRE> apply(Div, const CanonicalNumericType<P,PR,PRE>& c1, ScalarFunctionModel<P,D,PR,PRE> f2) {
        return mul(rec(std::move(f2)),c1); }

    static ScalarFunctionModel<P,D,PR,PRE> apply(Rec, const ScalarFunctionModel<P,D,PR,PRE>& f) {
        return f.apply(Rec()); }
    static ScalarFunctionModel<P,D,PR,PRE> apply(Div, const ScalarFunctionModel<P,D,PR,PRE>& f1, const ScalarFunctionModel<P,D,PR,PRE>& f2) {
        return mul(f1,rec(f2)); }
    static ScalarFunctionModel<P,D,PR,PRE> apply(Pow, const ScalarFunctionModel<P,D,PR,PRE>& f1, Int n2) {
        return generic_pow(f1,n2); }

    static ScalarFunctionModel<P,D,PR,PRE> apply(Abs, const ScalarFunctionModel<P,D,PR,PRE>& f) {
        return f.apply(Abs()); }
    static ScalarFunctionModel<P,D,PR,PRE> apply(Max, const ScalarFunctionModel<P,D,PR,PRE>& f1, const ScalarFunctionModel<P,D,PR,PRE>& f2) {
        return hlf(add(add(f1,f2),abs(sub(f1,f2)))); }
    static ScalarFunctionModel<P,D,PR,PRE> apply(Min, const ScalarFunctionModel<P,D,PR,PRE>& f1, const ScalarFunctionModel<P,D,PR,PRE>& f2) {
        return hlf(sub(add(f1,f2),abs(sub(f1,f2)))); }
*/

    typedef ScalarFunctionModel<P,D,PR,PRE> FM; typedef CanonicalNumericType<P,PR,PRE> X;
    typedef ScalarFunctionModelInterface<P,D,PR,PRE> FMI;

    static FM apply(BinaryElementaryOperator op, const FM& f1, const FM& f2) {
        return FM(&dynamic_cast<FMI&>(*f1._ptr->_apply(op,*f2._ptr))); }
    static FM apply(BinaryElementaryOperator op, const FM& f1, const X& c2) {
        return FM(&dynamic_cast<FMI&>(*f1._ptr->_apply(op,c2))); }
    static FM apply(BinaryElementaryOperator op, const X& c1, const FM& f2) {
        return FM(&dynamic_cast<FMI&>(*f2._ptr->_rapply(op,c1))); }
    static FM apply(UnaryElementaryOperator op, const FM& f) {
        return FM(&dynamic_cast<FMI&>(*f._ptr->_apply(op))); }
    static FM apply(GradedElementaryOperator op, const FM& f, Int n) {
        return FM(&dynamic_cast<FMI&>(*f._ptr->_apply(op,n))); }

/*
    static FM apply(BinaryElementaryOperator op, const X& c1, const FM& f2) {
        FM s1=factory(f2).create_constant(c1); return op(s1,f2); }

    static FM apply(BinaryElementaryOperator op, const FM& f1, const Number<P>& c2) {
        return op(f1,f1.create_constant(c2)); }
    static FM apply(BinaryElementaryOperator, const Number<P>& c1, const FM& f2) {
        return op(f2.create_constant(c1),f2); }

*/

};


// FIXME: Should not be needed since ScalarMultivariateFunctionModel has a representation
template<class P> inline ScalarMultivariateFunctionModel<P,DoublePrecision> embed(const ScalarMultivariateFunction<P>& f, const IntervalDomainType& d) {
    return embed(ScalarMultivariateFunctionModel<P,DoublePrecision>(f),d); }

template<class P, class D, class PR, class PRE> inline
ScalarFunctionModel<P,D,PR,PRE>& ScalarFunctionModel<P,D,PR,PRE>::operator=(const CanonicalNumericType<P,PR,PRE>& c) {
    (*this)*=CanonicalNumericType<P,PR,PRE>(0); (*this)+=c; return *this; }
template<class P, class D, class PR, class PRE> inline
ScalarFunctionModel<P,D,PR,PRE>& ScalarFunctionModel<P,D,PR,PRE>::operator=(const Number<P>& c) {
    return (*this)=factory(*this).create(c); }
template<class P, class D, class PR, class PRE> inline
ScalarFunctionModel<P,D,PR,PRE>& ScalarFunctionModel<P,D,PR,PRE>::operator=(const ScalarFunction<P,D>& f) {
    return (*this)=factory(*this).create(f); }
template<class P, class D, class PR, class PRE> inline
ScalarFunctionModel<P,D,PR,PRE>& ScalarFunctionModel<P,D,PR,PRE>::operator=(const ScalarFunctionModelInterface<P,D,PR,PRE>& f) {
    return (*this)=ScalarFunctionModel<P,D,PR,PRE>(f._clone()); }



template<class V> struct Element;

template<class M> class ScaledFunctionPatch;
template<class M> class VectorScaledFunctionPatch;
template<class M> struct Element<VectorScaledFunctionPatch<M>> { typedef ScaledFunctionPatch<M> Type; };

typedef ScaledFunctionPatch<ValidatedTaylorModelDP> ValidatedScalarMultivariateTaylorFunctionModelDP;

template<class P, class D, class PR, class PRE> class VectorFunctionModelElement
    : public DispatchTranscendentalAlgebraOperations<ScalarFunctionModel<P,D,PR,PRE>, CanonicalNumericType<P,PR,PRE>>
    , public ProvideConcreteGenericArithmeticOperators<ScalarFunctionModel<P,D,PR,PRE>>
{
    VectorFunctionModel<P,D,PR,PRE>* _p; SizeType _i;
  public:
    typedef typename ScalarFunctionModel<P,D,PR,PRE>::GenericType GenericType;
    operator const ScalarFunctionModel<P,D,PR,PRE> () const;
    VectorFunctionModelElement(VectorFunctionModel<P,D,PR,PRE>* p, SizeType i) : _p(p), _i(i) { }
    VectorFunctionModelElement<P,D,PR,PRE>& operator=(const ScalarFunctionModel<P,D,PR,PRE>& sf) {
        _p->set(_i,sf); return *this; }
    VectorFunctionModelElement<P,D,PR,PRE>& operator=(const VectorFunctionModelElement<P,D,PR,PRE>& sf) {
        return this->operator=(static_cast<ScalarFunctionModel<P,D,PR,PRE>const>(sf)); }
    Void clobber() { ScalarFunctionModel<P,D,PR,PRE> sf=_p->get(_i); sf.clobber(); _p->set(_i,sf); }
    const ScalarFunctionModel<P,D,PR,PRE> model() const { ScalarFunctionModel<P,D,PR,PRE> sf=_p->get(_i); return sf.model(); }
    const CanonicalErrorType<P,PRE> error() const { ScalarFunctionModel<P,D,PR,PRE> sf=_p->get(_i); return sf.error(); }
    Void set_error(CanonicalErrorType<P,PRE> e) const { ScalarFunctionModel<P,D,PR,PRE> sf=_p->get(_i); sf.set_error(e); _p->set(_i,sf); }
    Void set_error(Nat e) const { ScalarFunctionModel<P,D,PR,PRE> sf=_p->get(_i); sf.set_error(e); _p->set(_i,sf); }
    friend Boolean refines(const ScalarFunctionModel<ValidatedTag,D,PR,PRE>& f1, const ScalarFunctionModel<ValidatedTag,D,PR,PRE>& f2);
    friend Boolean inconsistent(const ScalarFunctionModel<ValidatedTag,D,PR,PRE>& f1, const ScalarFunctionModel<ValidatedTag,D,PR,PRE>& f2);
    friend ScalarFunctionModel<P,D,PR,PRE> antiderivative(VectorFunctionModelElement<P,D,PR,PRE> const& f, SizeType k) {
        return antiderivative(ScalarFunctionModel<P,D,PR,PRE>(f),k); }
    friend inline OutputStream& operator<<(OutputStream& os, const VectorFunctionModelElement<P,D,PR,PRE>& function) {
        return os << static_cast< const ScalarFunctionModel<P,D,PR,PRE> >(function); }
};

//! \ingroup FunctionModelSubModule
//! \brief Generic vector functions on bounded domains.
template<class P, class D, class PR, class PRE> class FunctionModel<P,D,BoxDomainType,PR,PRE>
{
    static_assert(IsSame<D,IntervalDomainType>::value or IsSame<D,BoxDomainType>::value,"");
    static_assert(IsSame<PRE,DoublePrecision>::value or IsSame<PRE,MultiplePrecision>::value,"");
    typedef BoxDomainType C;
  public:
    clone_on_copy_ptr< VectorFunctionModelInterface<P,D,PR,PRE> > _ptr;
  public:
    typedef VectorFunction<P,D> GenericType;
    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef D DomainType;
    typedef C CodomainType;

    typedef typename FunctionModelInterface<P,D,C,PR,PRE>::RangeType RangeType;
    typedef typename FunctionModelInterface<P,D,C,PR,PRE>::NormType NormType;
    typedef typename FunctionModelInterface<P,D,C,PR,PRE>::ValueType ValueType;
    typedef typename FunctionModelInterface<P,D,C,PR,PRE>::ErrorType ErrorType;
    typedef typename FunctionModelInterface<P,D,C,PR,PRE>::NumericType NumericType;
    typedef typename FunctionModelInterface<P,D,C,PR,PRE>::GenericNumericType GenericNumericType;

    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;
    template<class Y> using Result = ElementTraits<C>::template Type<Y>;
  public:
    inline FunctionModel() : _ptr() { }
    inline FunctionModel(SharedPointer<const FunctionModelInterface<P,D,C,PR,PRE>> vfp)
        : _ptr(vfp->_clone()) { }
    inline FunctionModel(SizeType n, const ScalarFunctionModelInterface<P,D,PR,PRE>& sf) {
        FunctionModelFactory<P,PR,PRE> factory(sf._factory()); *this=factory.create_zeros(n,sf.domain());
        for(SizeType i=0; i!=n; ++i) { (*this)[i]=sf; } }
    inline FunctionModel(Array<ScalarFunctionModel<P,D,PR,PRE>> const& asf)
        : FunctionModel(asf.size(),asf[0]) { for(SizeType i=0; i!=asf.size(); ++i) { (*this)[i]=asf[i]; } }
    inline FunctionModel(List<ScalarFunctionModel<P,D,PR,PRE>> const& lsf)
        : FunctionModel(lsf.size(),lsf[0]) { for(SizeType i=0; i!=lsf.size(); ++i) { (*this)[i]=lsf[i]; } }
    inline explicit FunctionModel(FunctionModelInterface<P,D,C,PR,PRE>* p) : _ptr(p) { }
    inline FunctionModel(const FunctionModelInterface<P,D,C,PR,PRE>& f) : _ptr(f._clone()) { }
    inline FunctionModel(const FunctionModel<P,D,C,PR,PRE>& f) : _ptr(f._ptr) { }
    inline FunctionModel& operator=(const FunctionModel<P,D,C,PR,PRE>& f) { this->_ptr=f._ptr; return *this; }
    inline operator const FunctionModelInterface<P,D,C,PR,PRE>& () const { return *_ptr; }
    inline operator Function<P,D,C> () const { return Function<P,D,C>(*_ptr); }
    inline const FunctionModelInterface<P,D,C,PR,PRE>* raw_pointer() const { return _ptr.operator->(); }
    inline const FunctionModelInterface<P,D,C,PR,PRE>& reference() const { return *_ptr; }
    inline FunctionModelInterface<P,D,C,PR,PRE>& reference() { return *_ptr; }

    inline SizeType result_size() const { return this->_ptr->result_size(); }
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    inline SizeType size() const { return this->_ptr->result_size(); }
    template<class XX> inline Vector<XX> operator()(const Vector<XX>& v) const { return this->_ptr->_evaluate(v); }
    template<class XX> inline Vector<XX> evaluate(const Vector<XX>& v) const { return this->_ptr->_evaluate(v); }
    inline ScalarFunctionModel<P,D,PR,PRE> const get(SizeType i) const { return ScalarFunctionModel<P,D,PR,PRE>(this->_ptr->_get(i)); }
    inline Void set(SizeType i, ScalarFunctionModel<P,D,PR,PRE> const& sf) { this->_ptr->_set(i,sf); }
    inline ScalarFunctionModel<P,D,PR,PRE> const operator[](SizeType i) const { return this->get(i); }
    inline VectorFunctionModelElement<P,D,PR,PRE> operator[](SizeType i) { return VectorFunctionModelElement<P,D,PR,PRE>(this,i); }
    inline VectorFunctionModel<P,D,PR,PRE> operator[](Range rng) { VectorFunctionModel<P,D,PR,PRE> r=factory(*this).create_zeros(rng.size());
        for(SizeType i=0; i!=rng.size(); ++i) { r[i]=this->operator[](rng[i]); } return r; }
    inline PrecisionType const precision() const { return this->get(0).precision(); }
    inline DomainType const domain() const { return this->_ptr->domain(); }
    inline CodomainType const codomain() const { return this->_ptr->codomain(); }
    inline RangeType const range() const { return this->_ptr->range(); }
    inline Vector<ValueType> const values() const { return this->_ptr->_values(); }
    inline Vector<ErrorType> const errors() const { return this->_ptr->_errors(); }
    inline ErrorType const error() const { return this->_ptr->_error(); }
    inline Void clobber() { this->_ptr->clobber(); }
    inline Matrix<NumericType> const jacobian(const Vector<NumericType>& x) const;
//        Vector<Differential<NumericType>> dfx=this->_ptr->_evaluate(Differential<NumericType>::variables(1u,x));
//        return dfx.jacobian(); }

    inline Void restrict(const DomainType& d) { *this=restriction(*this,d); }
  public:
    friend FunctionModelCreator<FunctionModelFactory<P,PR,PRE>,D> factory(VectorFunctionModel<P,D,PR,PRE> const& f) {
        FunctionModelFactory<P,PR,PRE> factory(f._ptr->_factory()); return FunctionModelCreator<FunctionModelFactory<P,PR,PRE>,D>(f.domain(),factory); }
  public:
    friend inline ScalarFunctionModel<P,D,PR,PRE> compose(const ScalarMultivariateFunction<P>& f, const VectorFunctionModel<P,D,PR,PRE>& g) {
        return ScalarFunctionModel<P,D,PR,PRE>(g._ptr->_compose(f)); }
    friend inline ScalarFunctionModel<P,D,PR,PRE> compose(const ScalarFunctionModel<P,D,PR,PRE>& f, const VectorFunctionModel<P,D,PR,PRE>& g) {
        return ScalarFunctionModel<P,D,PR,PRE>(g._ptr->_compose(f)); }
    friend inline VectorFunctionModel<P,D,PR,PRE> compose(const VectorMultivariateFunction<P>& f, const VectorFunctionModel<P,D,PR,PRE>& g) {
        return VectorFunctionModel<P,D,PR,PRE>(g._ptr->_compose(f)); }
    friend inline VectorFunctionModel<P,D,PR,PRE> compose(const VectorFunctionModel<P,D,PR,PRE>& f, const VectorFunctionModel<P,D,PR,PRE>& g) {
        return VectorFunctionModel<P,D,PR,PRE>(g._ptr->_compose(f)); }

    friend inline ScalarFunctionModel<P,D,PR,PRE> unchecked_compose(const ScalarFunctionModel<P,D,PR,PRE>& f, const VectorFunctionModel<P,D,PR,PRE>& g) {
        return ScalarFunctionModel<P,D,PR,PRE>(g._ptr->_unchecked_compose(f)); }
    friend inline VectorFunctionModel<P,D,PR,PRE> unchecked_compose(const VectorFunctionModel<P,D,PR,PRE>& f, const VectorFunctionModel<P,D,PR,PRE>& g) {
        return VectorFunctionModel<P,D,PR,PRE>(g._ptr->_unchecked_compose(f)); }

    friend inline VectorFunctionModel<P,D,PR,PRE> operator+(const VectorFunctionModel<P,D,PR,PRE>& f) {
        return VectorFunctionModel<P,D,PR,PRE>(f._ptr->_clone()); }
    friend inline VectorFunctionModel<P,D,PR,PRE> operator-(const VectorFunctionModel<P,D,PR,PRE>& f) {
        VectorFunctionModel<P,D,PR,PRE> r=f; for(SizeType i=0; i!=r.size(); ++i) { r[i]=-f[i]; } return r; }
    friend inline VectorFunctionModel<P,D,PR,PRE> operator+(const VectorFunctionModel<P,D,PR,PRE>& f1, const VectorFunctionModel<P,D,PR,PRE>& f2) {
        VectorFunctionModel<P,D,PR,PRE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]+f2[i]; } return r; }
    friend inline VectorFunctionModel<P,D,PR,PRE> operator-(const VectorFunctionModel<P,D,PR,PRE>& f1, const VectorFunctionModel<P,D,PR,PRE>& f2) {
        VectorFunctionModel<P,D,PR,PRE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]-f2[i]; } return r; }
    friend inline VectorFunctionModel<P,D,PR,PRE> operator+(const VectorFunctionModel<P,D,PR,PRE>& f1, const Vector<CanonicalNumericType<P,PR,PRE>>& c2) {
        VectorFunctionModel<P,D,PR,PRE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]+c2[i]; } return r; }
    friend inline VectorFunctionModel<P,D,PR,PRE> operator-(const VectorFunctionModel<P,D,PR,PRE>& f1, const Vector<CanonicalNumericType<P,PR,PRE>>& c2) {
        VectorFunctionModel<P,D,PR,PRE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]-c2[i]; } return r; }
    friend inline VectorFunctionModel<P,D,PR,PRE> operator+(const Vector<CanonicalNumericType<P,PR,PRE>>& c1, const VectorFunctionModel<P,D,PR,PRE>& f2);
    friend inline VectorFunctionModel<P,D,PR,PRE> operator-(const Vector<CanonicalNumericType<P,PR,PRE>>& c1, const VectorFunctionModel<P,D,PR,PRE>& f2);
    friend inline VectorFunctionModel<P,D,PR,PRE> operator*(const VectorFunctionModel<P,D,PR,PRE>& f1, const CanonicalNumericType<P,PR,PRE>& c2) {
        VectorFunctionModel<P,D,PR,PRE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]*c2; } return r; }
    friend inline VectorFunctionModel<P,D,PR,PRE> operator*(const CanonicalNumericType<P,PR,PRE>& c1, const VectorFunctionModel<P,D,PR,PRE>& f2) {
        VectorFunctionModel<P,D,PR,PRE> r=f2; for(SizeType i=0; i!=r.size(); ++i) { r[i]=c1*f2[i]; } return r; }

    friend inline VectorFunctionModel<P,D,PR,PRE> operator+(const VectorFunctionModel<P,D,PR,PRE>& f1, const VectorMultivariateFunction<P>& f2) {
        return f1+factory(f1).create(f2); }
    friend inline VectorFunctionModel<P,D,PR,PRE> operator-(const VectorFunctionModel<P,D,PR,PRE>& f1, const VectorMultivariateFunction<P>& f2) {
        return f1-factory(f1).create(f2); }
    friend inline VectorFunctionModel<P,D,PR,PRE> operator*(const VectorFunctionModel<P,D,PR,PRE>& f1, const ScalarMultivariateFunction<P>& f2) {
        return f1*factory(f1).create(f2); }
    friend inline VectorFunctionModel<P,D,PR,PRE> operator/(const VectorFunctionModel<P,D,PR,PRE>& f1, const ScalarMultivariateFunction<P>& f2) {
        return f1/factory(f1).create(f2); }
    friend inline VectorFunctionModel<P,D,PR,PRE> operator+(const VectorMultivariateFunction<P>& f1, const VectorFunctionModel<P,D,PR,PRE>& f2) {
        return factory(f2).create(f1)+f2; }
    friend inline VectorFunctionModel<P,D,PR,PRE> operator-(const VectorMultivariateFunction<P>& f1, const VectorFunctionModel<P,D,PR,PRE>& f2) {
        return factory(f2).create(f1)-f2; }
    friend inline VectorFunctionModel<P,D,PR,PRE> operator*(const ScalarMultivariateFunction<P>& f1, const VectorFunctionModel<P,D,PR,PRE>& f2) {
        return factory(f2).create(f1)*f2; }
    friend inline VectorFunctionModel<P,D,PR,PRE> operator/(const ScalarMultivariateFunction<P>& f1, const VectorFunctionModel<P,D,PR,PRE>& f2) {
        return factory(f2).create(f1)/f2; }


  public:
    friend NormType norm(const VectorFunctionModel<P,D,PR,PRE>& f) {
        return f._ptr->_norm(); }
    friend VectorFunctionModel<P,D,PR,PRE> embed(const DomainType& d1, const VectorFunctionModel<P,D,PR,PRE>& f, const DomainType& d2) {
        return VectorFunctionModel<P,D,PR,PRE>(f._ptr->_embed(d1,d2)); }
    friend VectorFunctionModel<P,D,PR,PRE> embed(const DomainType& d, const VectorFunctionModel<P,D,PR,PRE>& f) {
        return embed(d,f,DomainType()); }
    friend VectorFunctionModel<P,D,PR,PRE> embed(const VectorFunctionModel<P,D,PR,PRE>& f, const BoxDomainType& d) {
        return embed(DomainType(),f,d); }
    friend VectorFunctionModel<P,D,PR,PRE> embed(const VectorFunctionModel<P,D,PR,PRE>& f, const IntervalDomainType& d) {
        return embed(f,DomainType(1,d)); }
    friend VectorFunctionModel<P,D,PR,PRE> restriction(const VectorFunctionModel<P,D,PR,PRE>& f, const DomainType& d) {
        return VectorFunctionModel<P,D,PR,PRE>(f._ptr->_restriction(d)); }

    friend Vector<CanonicalNumericType<P,PR,PRE>> evaluate(const VectorFunctionModel<P,D,PR,PRE>& f, const Vector<CanonicalNumericType<P,PR,PRE>>& x) {
        return f._ptr->_evaluate(x); }
    friend Vector<Number<P>> evaluate(const VectorFunctionModel<P,D,PR,PRE>& f, const Vector<Number<P>>& x) {
        return f._ptr->_evaluate(Vector<CanonicalNumericType<P,PR,PRE>>(x,f.precision())); }

    friend Vector<CanonicalNumericType<P,PR,PRE>> unchecked_evaluate(const VectorFunctionModel<P,D,PR,PRE>& f, const Vector<CanonicalNumericType<P,PR,PRE>>& x) {
        return f._ptr->_unchecked_evaluate(x); }
    friend Vector<Number<P>> unchecked_evaluate(const VectorFunctionModel<P,D,PR,PRE>& f, const Vector<Number<P>>& x) {
        return f._ptr->_unchecked_evaluate(Vector<CanonicalNumericType<P,PR,PRE>>(x,f.precision())); }

    friend VectorFunctionModel<P,D,PR,PRE> partial_evaluate(const VectorFunctionModel<P,D,PR,PRE>& f, SizeType j, const CanonicalNumericType<P,PR,PRE>& c) {
        return VectorFunctionModel<P,D,PR,PRE>(f._ptr->_partial_evaluate(j,c)); }
    friend VectorFunctionModel<P,D,PR,PRE> partial_evaluate(const VectorFunctionModel<P,D,PR,PRE>& f, SizeType j, const Number<P>& c) {
        return partial_evaluate(f,j,CanonicalNumericType<P,PR,PRE>(c,f.precision())); }

    friend OutputStream& operator<<(OutputStream& os, const VectorFunctionModel<P,D,PR,PRE>& f) {
        return os <<  f.operator VectorMultivariateFunction<P>(); }

    friend ScalarFunctionModel<P,D,PR,PRE> unchecked_compose(const ScalarMultivariateFunction<P>& f, const VectorFunctionModel<P,D,PR,PRE>& g) {
        ScalarFunctionModelInterface<P,D,PR,PRE> const* fptr = dynamic_cast<ScalarFunctionModelInterface<P,D,PR,PRE> const*>(f.raw_pointer());
        if(fptr) { return unchecked_compose(ScalarFunctionModel<P,D,PR,PRE>(*fptr),g); } else { return compose(f,g); } }
    friend VectorFunctionModel<P,D,PR,PRE> unchecked_compose(const VectorMultivariateFunction<P>& f, const VectorFunctionModel<P,D,PR,PRE>& g) {
        VectorFunctionModelInterface<P,D,PR,PRE> const* fptr = dynamic_cast<VectorFunctionModelInterface<P,D,PR,PRE> const*>(f.raw_pointer());
        if(fptr) { return unchecked_compose(VectorFunctionModel<P,D,PR,PRE>(*fptr),g); } else { return compose(f,g); } }

    //friend VectorFunctionModel<P,D,PR,PRE> join(const ScalarFunctionModel<P,D,PR,PRE>& f1, const ScalarFunctionModel<P,D,PR,PRE>& f2);
    friend VectorFunctionModel<P,D,PR,PRE> join(const ScalarFunctionModel<P,D,PR,PRE>& f1, const VectorFunctionModel<P,D,PR,PRE>& f2) {
        return join(VectorFunctionModel<P,D,PR,PRE>(1u,f1),f2); }
    friend VectorFunctionModel<P,D,PR,PRE> join(const VectorFunctionModel<P,D,PR,PRE>& f1, const ScalarFunctionModel<P,D,PR,PRE>& f2) {
        VectorFunctionModel<P,D,PR,PRE> r=VectorFunctionModel<P,D,PR,PRE>(f1._ptr->_clone()); r._ptr->_adjoin(f2); return r; }
    friend VectorFunctionModel<P,D,PR,PRE> join(const VectorFunctionModel<P,D,PR,PRE>& f1, const VectorFunctionModel<P,D,PR,PRE>& f2) {
        return VectorFunctionModel<P,D,PR,PRE>(f1._ptr->_join(f2)); }

    friend VectorFunctionModel<P,D,PR,PRE> combine(const ScalarFunctionModel<P,D,PR,PRE>& f1, const ScalarFunctionModel<P,D,PR,PRE>& f2) {
        return VectorFunctionModel<P,D,PR,PRE>(1,f1)._ptr->_combine(VectorFunctionModel<P,D,PR,PRE>(1,f2)); };
    friend VectorFunctionModel<P,D,PR,PRE> combine(const ScalarFunctionModel<P,D,PR,PRE>& f1, const VectorFunctionModel<P,D,PR,PRE>& f2) {
        return VectorFunctionModel<P,D,PR,PRE>(1,f1)._ptr->_combine(f2); };
    friend VectorFunctionModel<P,D,PR,PRE> combine(const VectorFunctionModel<P,D,PR,PRE>& f1, const ScalarFunctionModel<P,D,PR,PRE>& f2) {
        return VectorFunctionModel<P,D,PR,PRE>(f1._ptr->_combine(VectorFunctionModel<P,D,PR,PRE>(1,f2))); };
    friend VectorFunctionModel<P,D,PR,PRE> combine(const VectorFunctionModel<P,D,PR,PRE>& f1, const VectorFunctionModel<P,D,PR,PRE>& f2) {
        return VectorFunctionModel<P,D,PR,PRE>(f1._ptr->_combine(f2)); }

    friend inline VectorFunctionModel<ValidatedTag,D,PR,PRE> refinement(const VectorFunctionModel<ValidatedTag,D,PR,PRE>& f1, const VectorFunctionModel<ValidatedTag,D,PR,PRE>& f2) {
        ARIADNE_ASSERT_MSG(f1.size()==f2.size(),"refinement(f1,f2): f1="<<f1<<", f2="<<f2<<")");
        VectorFunctionModel<ValidatedTag,D,PR,PRE> r=+f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=refinement(f1[i],f2[i]); } return r; }

    friend VectorFunctionModel<P,D,PR,PRE> antiderivative(const VectorFunctionModel<P,D,PR,PRE>& f, SizeType j) {
        VectorFunctionModel<P,D,PR,PRE> r(f);
        for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j); }
        return r;
    }

    friend VectorFunctionModel<P,D,PR,PRE> antiderivative(const VectorFunctionModel<P,D,PR,PRE>& f, SizeType j, const CanonicalNumericType<P,PR,PRE>& c) {
        VectorFunctionModel<P,D,PR,PRE> r(f);
        for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j,c); }
        return r;
    }

    friend VectorFunctionModel<P,D,PR,PRE> antiderivative(const VectorFunctionModel<P,D,PR,PRE>& f, SizeType j, const Number<P>& c) {
        return antiderivative(f,j,CanonicalNumericType<P,PR,PRE>(c,f[0].value().precision()));
    }

};

// FIXME: Implement for Multiple-Precision versions
template<class P> inline CanonicalNumeric64Type<P> unchecked_evaluate(const ScalarMultivariateFunction<P>& f, const Vector<CanonicalNumeric64Type<P>>& x) {
    ScalarMultivariateFunctionModelDPInterface<P> const* fptr = dynamic_cast<ScalarMultivariateFunctionModelDPInterface<P> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(ScalarMultivariateFunctionModelDP<P>(*fptr),x); } else { return evaluate(f,x); } }

template<class P> inline Vector<CanonicalNumeric64Type<P>> unchecked_evaluate(const VectorMultivariateFunction<P>& f, const Vector<CanonicalNumeric64Type<P>>& x) {
    VectorMultivariateFunctionModelDPInterface<P> const* fptr = dynamic_cast<VectorMultivariateFunctionModelDPInterface<P> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(VectorMultivariateFunctionModelDP<P>(*fptr),x); } else { return evaluate(f,x); } }

template<class P, class D, class PR, class PRE> inline ScalarFunctionModel<P,D,PR,PRE> unchecked_compose(const ScalarMultivariateFunction<P>& f, const VectorFunctionModel<P,D,PR,PRE>& g) {
    ScalarFunctionModelInterface<P,D,PR,PRE> const* fptr = dynamic_cast<ScalarFunctionModelInterface<P,D,PR,PRE> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(ScalarFunctionModel<P,D,PR,PRE>(*fptr),g); } else { return compose(f,g); } }
template<class P, class D, class PR, class PRE> inline VectorFunctionModel<P,D,PR,PRE> unchecked_compose(const VectorMultivariateFunction<P>& f, const VectorFunctionModel<P,D,PR,PRE>& g) {
    VectorFunctionModelInterface<P,D,PR,PRE> const* fptr = dynamic_cast<VectorFunctionModelInterface<P,D,PR,PRE> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(VectorFunctionModel<P,D,PR,PRE>(*fptr),g); } else { return compose(f,g); } }

// FIXME: Should be unneeded
template<class D, class PR, class PRE> ScalarFunctionModel<ValidatedTag,D,PR,PRE> unchecked_compose(const ScalarMultivariateFunctionModel<ValidatedTag,PR,PRE>& f, const VectorFunctionModel<ValidatedTag,D,PR,PRE>& g) {
    return VectorFunctionModel<ValidatedTag,D,PR,PRE>(g._ptr->_unchecked_compose(f)); }
template<class D, class PR, class PRE> VectorFunctionModel<ValidatedTag,D,PR,PRE> unchecked_compose(const VectorMultivariateFunctionModel<ValidatedTag,PR,PRE>& f, const VectorFunctionModel<ValidatedTag,D,PR,PRE>& g) {
    return VectorFunctionModel<ValidatedTag,D,PR,PRE>(g._ptr->_unchecked_compose(f)); }


template<class P, class D, class PR, class PRE> VectorFunctionModelElement<P,D,PR,PRE>::operator const ScalarFunctionModel<P,D,PR,PRE> () const {
    return ScalarFunctionModel<P,D,PR,PRE>(_p->get(_i)); }





// Full output
template<class T> struct Representation { const T* pointer; Representation(const T& t) : pointer(&t) { } const T& reference() const { return *pointer; } };
template<class T> inline Representation<T> representation(const T& t) { return Representation<T>(t); }

template<class T> class HasRepr {
    template<class TT, class = decltype(declval<TT>().repr(declval<OutputStream&>()))> static True test(int);
    template<class TT> static False test(...);
  public:
    static const bool value = decltype(test<T>(1))::value;
};

template<class T, EnableIf<HasRepr<T>> =dummy> inline OutputStream& operator<<(OutputStream& os, const Representation<T>& obj) { obj.reference().repr(os); return os; }
template<class T, DisableIf<HasRepr<T>> =dummy> inline OutputStream& operator<<(OutputStream& os, const Representation<T>& obj) { return os << obj.reference(); }

} // namespace Ariadne

#endif
