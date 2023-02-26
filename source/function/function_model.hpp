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

#include "function/function.decl.hpp"
#include "function/function_model_interface.hpp"

#include "numeric/operators.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/operations.hpp"
#include "function/domain.hpp"

#include "function/function_interface.hpp"
#include "function/function_mixin.hpp"
#include "function/function.hpp"

#include "function_patch.hpp"

namespace Ariadne {

template<class P, class A, class FLT, class FLTE> struct AlgebraOperations<ScalarFunctionModel<P,Real(A),FLT,FLTE>>;

// FIXME: Extend with univariate case
template<class P, class FLT, class FLTE> class FunctionModelFactory
    : public Handle<const FunctionModelFactoryInterface<P,FLT,FLTE>>
{
    using PR = PrecisionType<FLT>;
    using PRE = PrecisionType<FLTE>;
  private:
    typedef RealScalar SARG;
    typedef RealVector VARG;
    typedef IntervalDomainType SD;
    typedef BoxDomainType VD;
  public:
    typedef FunctionModelFactoryInterface<P,FLT,FLTE> Interface;

    typedef FLT CoefficientType;
    typedef FLTE ErrorValueType;
    typedef Error<FLTE> ErrorType;

    typedef P Paradigm;
    typedef SD ScalarDomainType;
    typedef VD VectorDomainType;

    using Handle<const Interface>::Handle;

    CanonicalNumericType<P,FLT,FLTE> create(Number<P> const& c) const {
        return CanonicalNumericType<P,FLT,FLTE>(this->_ptr->_create(c)); }
    ScalarMultivariateFunctionModel<P,FLT,FLTE> create(VectorDomainType const& dom, ScalarMultivariateFunction<P> const& f) const {
        return ScalarMultivariateFunctionModel<P,FLT,FLTE>(this->_ptr->_create(dom,f)); }
    VectorMultivariateFunctionModel<P,FLT,FLTE> create(VectorDomainType const& dom, VectorMultivariateFunction<P> const& f) const {
        return VectorMultivariateFunctionModel<P,FLT,FLTE>(this->_ptr->_create(dom,f)); }

    CanonicalNumericType<P,FLT,FLTE> create_number(Number<P> const& c) const {
        return CanonicalNumericType<P,FLT,FLTE>(this->_ptr->_create(c)); }
    ScalarMultivariateFunctionModel<P,FLT,FLTE> create_zero(VectorDomainType const& dom) const {
        return ScalarMultivariateFunctionModel<P,FLT,FLTE>(this->_ptr->_create_zero(dom)); }
    ScalarMultivariateFunctionModel<P,FLT,FLTE> create_constant(VectorDomainType const& dom, Number<P> const& c) const {
        return ScalarMultivariateFunctionModel<P,FLT,FLTE>(this->_ptr->_create_constant(dom,c)); }
    ScalarMultivariateFunctionModel<P,FLT,FLTE> create_constant(VectorDomainType const& dom, CanonicalNumericType<P,FLT,FLTE> const& c) const {
        return ScalarMultivariateFunctionModel<P,FLT,FLTE>(this->_ptr->_create_constant(dom,c)); }
    ScalarMultivariateFunctionModel<P,FLT,FLTE> create_coordinate(VectorDomainType const& dom, SizeType index) const {
        return ScalarMultivariateFunctionModel<P,FLT,FLTE>(this->_ptr->_create_coordinate(dom,index)); }
    VectorMultivariateFunctionModel<P,FLT,FLTE> create_zeros(SizeType rsize, VectorDomainType const& dom) const {
        return VectorMultivariateFunctionModel<P,FLT,FLTE>(this->_ptr->_create_zeros(rsize,dom)); }
    VectorMultivariateFunctionModel<P,FLT,FLTE> create_constants(VectorDomainType const& dom, Vector<Number<P>> const& c) const {
        return VectorMultivariateFunctionModel<P,FLT,FLTE>(this->_ptr->_create_constants(dom,c)); }
    VectorMultivariateFunctionModel<P,FLT,FLTE> create_projection(VectorDomainType const& dom, Range indices) const {
        return ScalarMultivariateFunctionModel<P,FLT,FLTE>(this->_ptr->_create_projection(dom,indices)); }
    VectorMultivariateFunctionModel<P,FLT,FLTE> create_identity(VectorDomainType const& dom) const {
        return VectorMultivariateFunctionModel<P,FLT,FLTE>(this->_ptr->_create_identity(dom)); }

    // FIXME: Should return a univariate model
    ScalarFunctionModel<P,VARG,FLT,FLTE> create_identity(ScalarDomainType const& dom) const {
        return ScalarFunctionModel<P,VARG,FLT,FLTE>(this->_ptr->_create_coordinate(VectorDomainType(1u,dom),0u)); }
    VectorFunctionModel<P,SARG,FLT,FLTE> create_zeros(SizeType n, ScalarDomainType const& dom) const {
        ARIADNE_NOT_IMPLEMENTED; }

    friend OutputStream& operator<<(OutputStream& os, FunctionModelFactory<P,FLT,FLTE> const& factory) { return factory._ptr->_write(os); }
};

template<class P, class FLT, class FLTE> auto inline
FunctionModelFactoryInterface<P,FLT,FLTE>::create_identity(IntervalDomainType const& dom) const -> ScalarFunctionModel<P,VARG,FLT,FLTE> {
    return ScalarFunctionModel<P,VARG,FLT,FLTE>(this->create_coordinate(BoxDomainType(1u,dom),0u)); }

template<class FCTRY, class ARG> class FunctionModelCreator {
    typedef BoundedDomainType<ARG> D;
    typedef typename FCTRY::Paradigm P;
    typedef typename FCTRY::CoefficientType FLT;
    typedef typename FCTRY::ErrorType FLTE;
  public:
    typedef FCTRY FactoryType;
    typedef D DomainType;
    typedef P Paradigm;

    explicit FunctionModelCreator(DomainType domain, FactoryType factory) : _factory(factory), _domain(domain) { }

    decltype(auto) create(Number<P> const& c) const { return this->_factory.create(c); }
    decltype(auto) create(ScalarFunction<P,ARG> const& f) { return this->_factory.create(this->_domain,f); }
    decltype(auto) create(VectorFunction<P,ARG> const& f) { return this->_factory.create(this->_domain,f); }
    decltype(auto) create_zero() { return this->_factory.create_zero(this->_domain); }
    decltype(auto) create_zeros(SizeType n) { return this->_factory.create_zeros(n,this->_domain); }
    decltype(auto) create_constant(Number<P> const& c) const { return this->_factory.create_constant(this->_domain,c); }
    decltype(auto) create_constants(Vector<Number<P>> const& c) const { return this->_factory.create_constants(this->_domain,c); }
    decltype(auto) create_identity() { return this->_factory.create_identity(this->_domain); }

    decltype(auto) create(DomainType const& dom, ScalarFunction<P,ARG> const& f) { return this->_factory.create(dom,f); }
    decltype(auto) create(DomainType const& dom, VectorFunction<P,ARG> const& f) { return this->_factory.create(dom,f); }
    decltype(auto) create_zero(DomainType const& dom) { return this->_factory.create_zero(dom); }
    decltype(auto) create_zeros(DomainType const& dom, SizeType n) { return this->_factory.create_zeros(n,dom); }
    decltype(auto) create_constant(DomainType const& dom, Number<P> const& c) const { return this->_factory.create_constant(dom,c); }
    decltype(auto) create_constants(DomainType const& dom, Vector<Number<P>> const& c) const { return this->_factory.create_constants(dom,c); }
    decltype(auto) create_identity(DomainType const& dom) { return this->_factory.create_identity(dom); }

    CanonicalNumericType<P,FLT,FLTE> const& create(CanonicalNumericType<P,FLT,FLTE> const& c) const { return c; }
  protected:
    FactoryType _factory;
    DomainType _domain;
};

// FIXME: Merge with multivariate case
template<class FCTRY> class FunctionModelCreator<FCTRY,RealScalar> {
    using P = typename FCTRY::Paradigm;
    using FLT = typename FCTRY::CoefficientType;
    using FLTE = typename FCTRY::ErrorValueType;
    using D = IntervalDomainType;
    using ARG = RealScalar;
  public:
    typedef FCTRY FactoryType;
    typedef D DomainType;
    typedef P Paradigm;

    explicit FunctionModelCreator(DomainType domain, FactoryType factory) : _factory(factory), _domain(domain) { }

    CanonicalNumericType<P,FLT,FLTE> create(Number<P> const& c) const { return this->_factory.create(c); }
    ScalarFunctionModel<P,ARG,FLT,FLTE> create(ScalarFunction<P,ARG> const& f) { ARIADNE_NOT_IMPLEMENTED; }
    VectorFunctionModel<P,ARG,FLT,FLTE> create(VectorFunction<P,ARG> const& f) { ARIADNE_NOT_IMPLEMENTED; }
    ScalarFunctionModel<P,ARG,FLT,FLTE> create_zero() { ARIADNE_NOT_IMPLEMENTED; }
    VectorFunctionModel<P,ARG,FLT,FLTE> create_zeros(SizeType n) { ARIADNE_NOT_IMPLEMENTED; }
    ScalarFunctionModel<P,ARG,FLT,FLTE> create_identity() { ARIADNE_NOT_IMPLEMENTED; }

    CanonicalNumericType<P,FLT,FLTE> const& create(CanonicalNumericType<P,FLT,FLTE> const& c) const { return c; }
  protected:
    FactoryType _factory;
    DomainType _domain;
};


//! \ingroup FunctionModule
//! \ingroup FunctionModelSubModule
//! \brief Generic class representing approximations to functions on bounded domains.
//!  \tparam P	The information paradigm tag, which can be either ValidatedTag, indicating that the approximation has a known (uniform) error bound, or ApproximateTag, indicating that no error bound is available. See the \ref InformationSubModule for more details.
//!  \tparam SIG The signature, which has the standard C++ form RES(ARG), so signature Real(RealVector) indicates a function \f$f:\R^n\to\R\f$. See \ref function_signature_section for more details.
//!  \tparam PR The precision used for numerical values within the approximation.
//!  \tparam PRE The precision used for error bounds provided for the approximation.
//! \paragraph function_model_template_parameter_note Note:
//!  It is planned to change the precision paramters \c PR and \c PRE to the actual number types used.
template<class P, class SIG, class FLT, class FLTE> class FunctionModel;

//! \ingroup FunctionModelSubModule
//! \brief Generic scalar functions on bounded domains.
template<class P, class ARG, class FLT, class FLTE> class FunctionModel<P,RealScalar(ARG),FLT,FLTE>
    : public Handle<FunctionModelInterface<P,RealScalar(ARG),FLT,FLTE>>
//    , public DispatchTranscendentalAlgebraOperations<ScalarFunctionModel<P,D,FLT,FLTE>, CanonicalNumericType<P,FLT,FLTE>>
    , public DispatchElementaryAlgebraOperations<ScalarFunctionModel<P,ARG,FLT,FLTE>, CanonicalNumericType<P,FLT,FLTE>>
    , public ProvideConcreteGenericElementaryOperations<ScalarFunctionModel<P,ARG,FLT,FLTE>,ScalarMultivariateFunction<P>>
    , public ProvideConcreteGenericElementaryOperations<ScalarFunctionModel<P,ARG,FLT,FLTE>,Number<P>>
{
    static_assert(Same<ARG,RealScalar> or Same<ARG,RealVector>,"");
    using RES=RealScalar; using SIG=RES(ARG);
    using PR=PrecisionType<FLT>; using PRE=PrecisionType<FLTE>;
    using C=typename SignatureTraits<SIG>::BoundedCodomainType;
    using D=typename SignatureTraits<SIG>::BoundedDomainType;
    static_assert(Same<D,IntervalDomainType> or Same<D,BoxDomainType>,"");
  public:
    typedef FunctionModelInterface<P,SIG,FLT,FLTE> Interface;
    typedef ScalarFunction<P,ARG> GenericType;
    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef D DomainType;
    typedef C CodomainType;
    typedef typename Interface::ValueType ValueType;
    typedef typename Interface::ErrorType ErrorType;
    typedef typename Interface::NumericType NumericType;
    typedef typename Interface::GenericNumericType GenericNumericType;
    typedef typename Interface::NormType NormType;
    typedef typename Interface::RangeType RangeType;

    template<class Y> using Argument = typename SignatureTraits<SIG>::template Argument<Y>;
    template<class Y> using Result = typename SignatureTraits<SIG>::template Result<Y>;
  public:
    explicit FunctionModel(Interface* p) : Handle<Interface>(p) { }
    FunctionModel(SharedPointer<Interface> p) : Handle<Interface>(p) { }

    FunctionModel() : FunctionModel(nullptr) { }
    FunctionModel(const SharedPointer<const Interface> p) : Handle<Interface>(p->_clone()) { }
    FunctionModel(const FunctionModel<P,SIG,FLT,FLTE>& f) : Handle<Interface>(f._ptr) { }
    FunctionModel& operator=(const FunctionModel<P,SIG,FLT,FLTE>& f) { this->_ptr=f._ptr; return *this; }
    FunctionModel(const Interface& f) : Handle<Interface>(f._clone()) { }
    FunctionModel(const Function<P,SIG>& f) : FunctionModel(dynamic_cast<Interface*>(f.raw_pointer()->_clone())) { }
    operator Function<P,SIG>() const { return Function<P,SIG>(this->_ptr->_clone()); }
    operator Interface& () { return this->reference(); }
    operator const Interface& () const { return this->reference(); }

    ScalarFunctionModel<P,ARG,FLT,FLTE>& operator=(const Number<P>& c);
    ScalarFunctionModel<P,ARG,FLT,FLTE>& operator=(const CanonicalNumericType<P,FLT,FLTE>& c);
    ScalarFunctionModel<P,ARG,FLT,FLTE>& operator=(const ScalarFunction<P,ARG>& f);
    ScalarFunctionModel<P,ARG,FLT,FLTE>& operator=(const ScalarFunctionModelInterface<P,ARG,FLT,FLTE>& f);
//    ScalarFunctionModel<P,ARG,FLT,FLTE>& operator=(const ValidatedScalarMultivariateTaylorFunctionModelDP& f);

    inline PrecisionType const precision() const { return this->value().precision(); }
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    template<class X> X operator() (const Vector<X>& x) const {
        return this->_ptr->_call(x); }
    template<class X> X evaluate(const Vector<X>& x) const {
        return this->_ptr->_call(x); }
    inline DomainType const domain() const { return this->_ptr->domain(); }
    inline CodomainType const codomain() const { return this->_ptr->codomain(); }
    inline RangeType const range() const { return this->_ptr->range(); }

    inline ValueType value() const { return this->_ptr->_value(); }
    inline ErrorType error() const { return this->_ptr->_concrete_error(); }

    inline Void clobber() { this->pointer()->clobber(); }

//    inline ScalarFunctionModel<P,ARG,FLT,FLTE> apply(UnaryElementaryOperator op) const { return ScalarFunctionModel<P,ARG,FLT,FLTE>(this->_ptr->_apply(op)); }
    inline Void restrict(const DomainType& d) { *this=restriction(*this,d); }
  public:
    friend FunctionModelCreator<FunctionModelFactory<P,FLT,FLTE>,ARG> factory(ScalarFunctionModel<P,ARG,FLT,FLTE> const& f) {
        FunctionModelFactory<P,FLT,FLTE> factory(f._ptr->_factory()); return FunctionModelCreator<FunctionModelFactory<P,FLT,FLTE>,ARG>(f.domain(),factory); }
  public:
  public:
    friend CanonicalNumericType<P,FLT,FLTE> evaluate(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, const Vector<CanonicalNumericType<P,FLT,FLTE>>& x) {
        return f._ptr->_call(x); }
    friend Number<P> evaluate(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, const Vector<Number<P>>& x) {
        return f._ptr->_call(Vector<CanonicalNumericType<P,FLT,FLTE>>(x,f.precision())); }
    friend CanonicalNumericType<P,FLT,FLTE> unchecked_evaluate(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, const Vector<CanonicalNumericType<P,FLT,FLTE>>& x) {
        return f._ptr->_unchecked_evaluate(x); }
    friend Number<P> unchecked_evaluate(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, const Vector<Number<P>>& x) {
        return f._ptr->_unchecked_evaluate(Vector<CanonicalNumericType<P,FLT,FLTE>>(x,f.precision())); }

    friend ScalarFunctionModel<P,ARG,FLT,FLTE> partial_evaluate(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, SizeType j, const CanonicalNumericType<P,FLT,FLTE>& c) {
        return ScalarFunctionModel<P,ARG,FLT,FLTE>(f._ptr->_concrete_partial_evaluate(j,c)); }
    friend ScalarFunctionModel<P,ARG,FLT,FLTE> partial_evaluate(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, SizeType j, const Number<P>& c) {
        return partial_evaluate(f,j,CanonicalNumericType<P,FLT,FLTE>(c,f.precision())); }

    friend NormType norm(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f) {
        return f._ptr->_concrete_norm(); }
    friend ScalarFunctionModel<P,ARG,FLT,FLTE> derivative(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, SizeType j) {
        return ScalarFunctionModel<P,ARG,FLT,FLTE>(f._ptr->_derivative(j)); }
    friend ScalarFunctionModel<P,ARG,FLT,FLTE> antiderivative(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, SizeType j) {
        return ScalarFunctionModel<P,ARG,FLT,FLTE>(f._ptr->_antiderivative(j)); }
    friend ScalarFunctionModel<P,ARG,FLT,FLTE> antiderivative(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, SizeType j, CanonicalNumericType<P,FLT,FLTE> c) {
        return ScalarFunctionModel<P,ARG,FLT,FLTE>(f._ptr->_antiderivative(j,c)); }
    friend ScalarFunctionModel<P,ARG,FLT,FLTE> antiderivative(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, SizeType j, const Number<P>& c) {
        return antiderivative(f,j,CanonicalNumericType<P,FLT,FLTE>(c,f.value().precision())); }

    friend ScalarFunctionModel<P,ARG,FLT,FLTE> embed(const DomainType& d1, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, const DomainType& d2) {
        return ScalarFunctionModel<P,ARG,FLT,FLTE>(f._ptr->_embed(d1,d2)); }
    friend ScalarFunctionModel<P,ARG,FLT,FLTE> embed(const DomainType& d, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f) {
        return embed(d,f,DomainType()); }
    friend ScalarFunctionModel<P,ARG,FLT,FLTE> embed(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, const BoxDomainType& d) {
        return embed(DomainType(),f,d); }
    friend ScalarFunctionModel<P,ARG,FLT,FLTE> embed(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, const IntervalDomainType& d) {
        return embed(f,DomainType(1,d)); }
    friend ScalarFunctionModel<P,ARG,FLT,FLTE> restrict(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, const DomainType& d) {
        return ScalarFunctionModel<P,ARG,FLT,FLTE>(f._ptr->_restriction(d)); }
    friend ScalarFunctionModel<P,ARG,FLT,FLTE> restriction(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, const DomainType& d) {
        return ScalarFunctionModel<P,ARG,FLT,FLTE>(f._ptr->_restriction(d)); }

    friend VectorFunctionModel<P,ARG,FLT,FLTE> join(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f2) {
        return join(VectorFunctionModel<P,ARG,FLT,FLTE>(1,f1),f2); }
    friend VectorFunctionModel<P,ARG,FLT,FLTE> combine(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f2);
  public:
    typedef ValidatedTag VP;
    friend ScalarFunctionModel<VP,ARG,FLT,FLTE> refinement(const ScalarFunctionModel<VP,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<VP,ARG,FLT,FLTE>& f2) {
        return ScalarFunctionModel<VP,ARG,FLT,FLTE>(f1._ptr->_refinement(f2)); }
    friend Boolean inconsistent(const ScalarFunctionModel<VP,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<VP,ARG,FLT,FLTE>& f2) {
        return f1._ptr->_inconsistent(f2); }
    friend Boolean refines(const ScalarFunctionModel<VP,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<VP,ARG,FLT,FLTE>& f2) {
        return f1._ptr->_refines(f2); }
  public:
    friend OutputStream& operator<<(OutputStream& os, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f) {
        return os <<  f.operator ScalarMultivariateFunction<P>(); }
};

template<class P, class ARG, class FLT, class FLTE> struct AlgebraOperations<ScalarFunctionModel<P,ARG,FLT,FLTE>> {

/*
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Nul, ScalarFunctionModel<P,ARG,FLT,FLTE> f) {
        f._ptr->_imul(CanonicalNumericType<P,FLT,FLTE>(0)); return f; }
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Neg, ScalarFunctionModel<P,ARG,FLT,FLTE> f) {
        f._ptr->_imul(CanonicalNumericType<P,FLT,FLTE>(-1)); return f; }
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Add, ScalarFunctionModel<P,ARG,FLT,FLTE> f1, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f2) {
        f1._ptr->_isma(CanonicalNumericType<P,FLT,FLTE>(+1),f2); return f1; }
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Sub, ScalarFunctionModel<P,ARG,FLT,FLTE> f1, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f2) {
        f1._ptr->_isma(CanonicalNumericType<P,FLT,FLTE>(-1),f2); return f1; }
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Mul, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f2) {
        ScalarFunctionModel<P,ARG,FLT,FLTE> r=factory(f1).create_zero(); r._ptr->_ifma(f1,f2); return r; }
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Add, ScalarFunctionModel<P,ARG,FLT,FLTE> f1, const CanonicalNumericType<P,FLT,FLTE>& c2) {
        f1._ptr->_iadd(c2); return f1; }
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Mul, ScalarFunctionModel<P,ARG,FLT,FLTE> f1, const CanonicalNumericType<P,FLT,FLTE>& c2) {
        f1._ptr->_imul(c2); return f1; }
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Add, const CanonicalNumericType<P,FLT,FLTE>& c1, ScalarFunctionModel<P,ARG,FLT,FLTE> f2) {
        f2._ptr->_iadd(c1); return f2; }
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Mul, const CanonicalNumericType<P,FLT,FLTE>& c1, ScalarFunctionModel<P,ARG,FLT,FLTE> f2) {
        f2._ptr->_imul(c1); return f2; }

    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Sub, ScalarFunctionModel<P,ARG,FLT,FLTE> f1, const CanonicalNumericType<P,FLT,FLTE>& c2) {
        return add(std::move(f1),neg(c2)); }
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Div, ScalarFunctionModel<P,ARG,FLT,FLTE> f1, const CanonicalNumericType<P,FLT,FLTE>& c2) {
        return mul(std::move(f1),rec(c2)); }
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Sub, const CanonicalNumericType<P,FLT,FLTE>& c1, ScalarFunctionModel<P,ARG,FLT,FLTE> f2) {
        return add(neg(std::move(f2)),c1); }
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Div, const CanonicalNumericType<P,FLT,FLTE>& c1, ScalarFunctionModel<P,ARG,FLT,FLTE> f2) {
        return mul(rec(std::move(f2)),c1); }

    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Rec, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f) {
        return f.apply(Rec()); }
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Div, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f2) {
        return mul(f1,rec(f2)); }
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Pow, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f1, Int n2) {
        return generic_pow(f1,n2); }

    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Abs, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f) {
        return f.apply(Abs()); }
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Max, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f2) {
        return hlf(add(add(f1,f2),abs(sub(f1,f2)))); }
    static ScalarFunctionModel<P,ARG,FLT,FLTE> apply(Min, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f2) {
        return hlf(sub(add(f1,f2),abs(sub(f1,f2)))); }
*/

    typedef ScalarFunctionModel<P,ARG,FLT,FLTE> FM; typedef CanonicalNumericType<P,FLT,FLTE> X;
    typedef ScalarFunctionModelInterface<P,ARG,FLT,FLTE> FMI;
    typedef ScalarFunctionModelAlgebraInterface<P,ARG,FLT,FLTE> FMAI;

    static FM apply(BinaryElementaryOperator op, const FM& f1, const FM& f2) {
         return FM(&dynamic_cast<FMI&>(*static_cast<FMAI const&>(f1.reference())._apply(op,f2.reference()))); }
    static FM apply(BinaryElementaryOperator op, const FM& f1, const X& c2) {
         return FM(&dynamic_cast<FMI&>(*static_cast<FMAI const&>(f1.reference())._apply(op,c2))); }
    static FM apply(BinaryElementaryOperator op, const X& c1, const FM& f2) {
         return FM(&dynamic_cast<FMI&>(*static_cast<FMAI const&>(f2.reference())._rapply(op,c1))); }
    static FM apply(UnaryElementaryOperator op, const FM& f) {
         return FM(&dynamic_cast<FMI&>(*static_cast<FMAI const&>(f.reference())._apply(op))); }
    static FM apply(GradedElementaryOperator op, const FM& f, Int n) {
         return FM(&dynamic_cast<FMI&>(*static_cast<FMAI const&>(f.reference())._apply(op,n))); }

/*
    static FM apply(BinaryElementaryOperator op, const X& c1, const FM& f2) {
        FM s1=factory(f2).create_constant(c1); return op(s1,f2); }

    static FM apply(BinaryElementaryOperator op, const FM& f1, const Number<P>& c2) {
        return op(f1,f1.create_constant(c2)); }
    static FM apply(BinaryElementaryOperator, const Number<P>& c1, const FM& f2) {
        return op(f2.create_constant(c1),f2); }

*/

};


template<class P, class ARG, class FLT, class FLTE> inline
ScalarFunctionModel<P,ARG,FLT,FLTE>& ScalarFunctionModel<P,ARG,FLT,FLTE>::operator=(const CanonicalNumericType<P,FLT,FLTE>& c) {
    (*this)*=nul(c); (*this)+=c; return *this; }
template<class P, class ARG, class FLT, class FLTE> inline
ScalarFunctionModel<P,ARG,FLT,FLTE>& ScalarFunctionModel<P,ARG,FLT,FLTE>::operator=(const Number<P>& c) {
    return (*this)=factory(*this).create(c); }
template<class P, class ARG, class FLT, class FLTE> inline
ScalarFunctionModel<P,ARG,FLT,FLTE>& ScalarFunctionModel<P,ARG,FLT,FLTE>::operator=(const ScalarFunction<P,ARG>& f) {
    return (*this)=factory(*this).create(f); }
template<class P, class ARG, class FLT, class FLTE> inline
ScalarFunctionModel<P,ARG,FLT,FLTE>& ScalarFunctionModel<P,ARG,FLT,FLTE>::operator=(const ScalarFunctionModelInterface<P,ARG,FLT,FLTE>& f) {
    return (*this)=ScalarFunctionModel<P,ARG,FLT,FLTE>(f._clone()); }



template<class V> struct Element;

template<class M> class ScaledFunctionPatch;
template<class M> class VectorScaledFunctionPatch;
template<class M> struct Element<VectorScaledFunctionPatch<M>> { typedef ScaledFunctionPatch<M> Type; };

typedef ScaledFunctionPatch<ValidatedTaylorModelDP> ValidatedScalarMultivariateTaylorFunctionModelDP;

template<class P, class ARG, class FLT, class FLTE> class VectorFunctionModelElement
    : public DispatchTranscendentalAlgebraOperations<ScalarFunctionModel<P,ARG,FLT,FLTE>, CanonicalNumericType<P,FLT,FLTE>>
    , public ProvideConcreteGenericArithmeticOperators<ScalarFunctionModel<P,ARG,FLT,FLTE>>
{
    VectorFunctionModel<P,ARG,FLT,FLTE>* _p; SizeType _i;
  public:
    typedef typename ScalarFunctionModel<P,ARG,FLT,FLTE>::GenericType GenericType;
    typedef typename ScalarFunctionModel<P,ARG,FLT,FLTE>::RangeType RangeType;
    operator const ScalarFunctionModel<P,ARG,FLT,FLTE> () const;
    VectorFunctionModelElement(VectorFunctionModel<P,ARG,FLT,FLTE>* p, SizeType i) : _p(p), _i(i) { }
    VectorFunctionModelElement<P,ARG,FLT,FLTE>& operator=(const ScalarFunctionModel<P,ARG,FLT,FLTE>& sf) {
        _p->set(_i,sf); return *this; }
    VectorFunctionModelElement<P,ARG,FLT,FLTE>& operator=(const VectorFunctionModelElement<P,ARG,FLT,FLTE>& sf) {
        return this->operator=(static_cast<ScalarFunctionModel<P,ARG,FLT,FLTE>const>(sf)); }
    Void clobber() { ScalarFunctionModel<P,ARG,FLT,FLTE> sf=_p->get(_i); sf.clobber(); _p->set(_i,sf); }
    const ScalarFunctionModel<P,ARG,FLT,FLTE> model() const { ScalarFunctionModel<P,ARG,FLT,FLTE> sf=_p->get(_i); return sf.model(); }
    const CanonicalErrorType<P,FLTE> error() const { ScalarFunctionModel<P,ARG,FLT,FLTE> sf=_p->get(_i); return sf.error(); }
    Void set_error(CanonicalErrorType<P,FLTE> e) const { ScalarFunctionModel<P,ARG,FLT,FLTE> sf=_p->get(_i); sf.set_error(e); _p->set(_i,sf); }
    Void set_error(Nat e) const { ScalarFunctionModel<P,ARG,FLT,FLTE> sf=_p->get(_i); sf.set_error(e); _p->set(_i,sf); }
    const RangeType range() const { ScalarFunctionModel<P,ARG,FLT,FLTE> sf=_p->get(_i); return sf.range(); }
    friend Boolean refines(const ScalarFunctionModel<ValidatedTag,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<ValidatedTag,ARG,FLT,FLTE>& f2);
    friend Boolean inconsistent(const ScalarFunctionModel<ValidatedTag,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<ValidatedTag,ARG,FLT,FLTE>& f2);
    friend ScalarFunctionModel<P,ARG,FLT,FLTE> antiderivative(VectorFunctionModelElement<P,ARG,FLT,FLTE> const& f, SizeType k) {
        return antiderivative(ScalarFunctionModel<P,ARG,FLT,FLTE>(f),k); }
    friend inline OutputStream& operator<<(OutputStream& os, const VectorFunctionModelElement<P,ARG,FLT,FLTE>& function) {
        return os << static_cast< const ScalarFunctionModel<P,ARG,FLT,FLTE> >(function); }
};

//! \ingroup FunctionModelSubModule
//! \brief Generic vector functions on bounded domains.
template<class P, class ARG, class FLT, class FLTE> class FunctionModel<P,RealVector(ARG),FLT,FLTE>
    : public Handle<FunctionModelInterface<P,RealVector(ARG),FLT,FLTE>>
{
    static_assert(Same<ARG,RealScalar> or Same<ARG,RealVector>,"");
    static_assert(AFloat<FLTE>);
    using RES=RealVector; using SIG=RES(ARG);
    using PR=PrecisionType<FLT>; using PRE=PrecisionType<FLTE>;
    using C=typename SignatureTraits<SIG>::BoundedCodomainType;
    using D=typename SignatureTraits<SIG>::BoundedDomainType;
  public:
    typedef FunctionModelInterface<P,SIG,FLT,FLTE> Interface;
  public:
    typedef VectorFunction<P,ARG> GenericType;
    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef D DomainType;
    typedef C CodomainType;

    typedef typename Interface::RangeType RangeType;
    typedef typename Interface::NormType NormType;
    typedef typename Interface::ValueType ValueType;
    typedef typename Interface::ErrorType ErrorType;
    typedef typename Interface::NumericType NumericType;
    typedef typename Interface::GenericNumericType GenericNumericType;

    template<class Y> using Argument = typename SignatureTraits<SIG>::template Argument<Y>;
    template<class Y> using Result = typename SignatureTraits<SIG>::template Result<Y>;
  public:
    inline explicit FunctionModel(Interface* p) : Handle<Interface>(p) { }
    inline FunctionModel(SharedPointer<Interface> p) : Handle<Interface>(p) { }

    inline FunctionModel() : FunctionModel(nullptr) { }
    inline FunctionModel(SharedPointer<const Interface> vfp)
        : FunctionModel(vfp->_clone()) { }
    inline FunctionModel(SizeType n, const ScalarFunctionModelInterface<P,ARG,FLT,FLTE>& sf) : FunctionModel() {
        FunctionModelFactory<P,FLT,FLTE> factory(sf._factory()); *this=factory.create_zeros(n,sf.domain());
        for(SizeType i=0; i!=n; ++i) { (*this)[i]=sf; } }
    inline FunctionModel(Array<ScalarFunctionModel<P,ARG,FLT,FLTE>> const& asf)
        : FunctionModel(asf.size(),asf[0]) { for(SizeType i=0; i!=asf.size(); ++i) { (*this)[i]=asf[i]; } }
    inline FunctionModel(List<ScalarFunctionModel<P,ARG,FLT,FLTE>> const& lsf)
        : FunctionModel(lsf.size(),lsf[0]) { for(SizeType i=0; i!=lsf.size(); ++i) { (*this)[i]=lsf[i]; } }
    inline FunctionModel(const Interface& f) : FunctionModel(f._clone()) { }
    inline operator Function<P,SIG> () const { return Function<P,SIG>(*this->_ptr); }

    inline SizeType result_size() const { return this->_ptr->result_size(); }
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    inline SizeType size() const { return this->_ptr->result_size(); }
    template<class XX> inline Vector<XX> operator()(const Vector<XX>& v) const { return this->_ptr->_call(v); }
    template<class XX> inline Vector<XX> evaluate(const Vector<XX>& v) const { return this->_ptr->_call(v); }
    inline ScalarFunctionModel<P,ARG,FLT,FLTE> const get(SizeType i) const { return ScalarFunctionModel<P,ARG,FLT,FLTE>(this->_ptr->_concrete_get(i)); }
    inline Void set(SizeType i, ScalarFunctionModel<P,ARG,FLT,FLTE> const& sf) { this->pointer()->_concrete_set(i,sf); }
    inline ScalarFunctionModel<P,ARG,FLT,FLTE> const operator[](SizeType i) const { return this->get(i); }
    inline VectorFunctionModelElement<P,ARG,FLT,FLTE> operator[](SizeType i) { return VectorFunctionModelElement<P,ARG,FLT,FLTE>(this,i); }
    inline VectorFunctionModel<P,ARG,FLT,FLTE> operator[](Range rng) const { VectorFunctionModel<P,ARG,FLT,FLTE> r=factory(*this).create_zeros(rng.size());
        for(SizeType i=0; i!=rng.size(); ++i) { r[i]=this->operator[](rng[i]); } return r; }
    inline PrecisionType const precision() const { return this->get(0).precision(); }
    inline DomainType const domain() const { return this->_ptr->domain(); }
    inline CodomainType const codomain() const { return this->_ptr->codomain(); }
    inline RangeType const range() const { return this->_ptr->range(); }
    inline Vector<ValueType> const values() const { return this->_ptr->_concrete_values(); }
    inline Vector<ErrorType> const errors() const { return this->_ptr->_concrete_errors(); }
    inline ErrorType const error() const { return this->_ptr->_concrete_error(); }
    inline Void clobber() { this->pointer()->clobber(); }
    inline Matrix<NumericType> const jacobian(const Vector<NumericType>& x) const;
//        Vector<Differential<NumericType>> dfx=this->_ptr->_call(Differential<NumericType>::variables(1u,x));
//        return dfx.jacobian(); }

    inline Void restrict(const DomainType& d) { *this=restriction(*this,d); }
  public:
    friend FunctionModelCreator<FunctionModelFactory<P,FLT,FLTE>,ARG> factory(VectorFunctionModel<P,ARG,FLT,FLTE> const& f) {
        FunctionModelFactory<P,FLT,FLTE> factory(f._ptr->_factory()); return FunctionModelCreator<FunctionModelFactory<P,FLT,FLTE>,ARG>(f.domain(),factory); }
  public:
    friend inline ScalarFunctionModel<P,ARG,FLT,FLTE> compose(const ScalarMultivariateFunction<P>& f, const VectorFunctionModel<P,ARG,FLT,FLTE>& g) {
        return ScalarFunctionModel<P,ARG,FLT,FLTE>(g._ptr->_compose(f)); }
    friend inline ScalarFunctionModel<P,ARG,FLT,FLTE> compose(const ScalarMultivariateFunctionPatch<P>& f, const VectorFunctionModel<P,ARG,FLT,FLTE>& g) {
        return ScalarFunctionModel<P,ARG,FLT,FLTE>(g._ptr->_compose(f)); }
    friend inline ScalarFunctionModel<P,ARG,FLT,FLTE> compose(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, const VectorFunctionModel<P,ARG,FLT,FLTE>& g) {
        return ScalarFunctionModel<P,ARG,FLT,FLTE>(g._ptr->_compose(f)); }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> compose(const VectorMultivariateFunction<P>& f, const VectorFunctionModel<P,ARG,FLT,FLTE>& g) {
        return VectorFunctionModel<P,ARG,FLT,FLTE>(g._ptr->_compose(f)); }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> compose(const VectorMultivariateFunctionPatch<P>& f, const VectorFunctionModel<P,ARG,FLT,FLTE>& g) {
        return VectorFunctionModel<P,ARG,FLT,FLTE>(g._ptr->_compose(f)); }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> compose(const VectorFunctionModel<P,ARG,FLT,FLTE>& f, const VectorFunctionModel<P,ARG,FLT,FLTE>& g) {
        return VectorFunctionModel<P,ARG,FLT,FLTE>(g._ptr->_compose(f)); }

    friend inline ScalarFunctionModel<P,ARG,FLT,FLTE> unchecked_compose(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f, const VectorFunctionModel<P,ARG,FLT,FLTE>& g) {
        return ScalarFunctionModel<P,ARG,FLT,FLTE>(g._ptr->_concrete_unchecked_compose(f)); }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> unchecked_compose(const VectorFunctionModel<P,ARG,FLT,FLTE>& f, const VectorFunctionModel<P,ARG,FLT,FLTE>& g) {
        return VectorFunctionModel<P,ARG,FLT,FLTE>(g._ptr->_concrete_unchecked_compose(f)); }

    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator+(const VectorFunctionModel<P,ARG,FLT,FLTE>& f) {
        return VectorFunctionModel<P,ARG,FLT,FLTE>(f._ptr->_clone()); }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator-(const VectorFunctionModel<P,ARG,FLT,FLTE>& f) {
        VectorFunctionModel<P,ARG,FLT,FLTE> r=f; for(SizeType i=0; i!=r.size(); ++i) { r[i]=-f[i]; } return r; }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator+(const VectorFunctionModel<P,ARG,FLT,FLTE>& f1, const VectorFunctionModel<P,ARG,FLT,FLTE>& f2) {
        VectorFunctionModel<P,ARG,FLT,FLTE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]+f2[i]; } return r; }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator-(const VectorFunctionModel<P,ARG,FLT,FLTE>& f1, const VectorFunctionModel<P,ARG,FLT,FLTE>& f2) {
        VectorFunctionModel<P,ARG,FLT,FLTE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]-f2[i]; } return r; }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator*(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f1, const VectorFunctionModel<P,ARG,FLT,FLTE>& f2) {
        VectorFunctionModel<P,ARG,FLT,FLTE> r=f2; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1*f2[i]; } return r; }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator*(const VectorFunctionModel<P,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f2) {
        VectorFunctionModel<P,ARG,FLT,FLTE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]*f2; } return r; }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator/(const VectorFunctionModel<P,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f2) {
        VectorFunctionModel<P,ARG,FLT,FLTE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]/f2; } return r; }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator+(const VectorFunctionModel<P,ARG,FLT,FLTE>& f1, const Vector<CanonicalNumericType<P,FLT,FLTE>>& c2) {
        VectorFunctionModel<P,ARG,FLT,FLTE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]+c2[i]; } return r; }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator-(const VectorFunctionModel<P,ARG,FLT,FLTE>& f1, const Vector<CanonicalNumericType<P,FLT,FLTE>>& c2) {
        VectorFunctionModel<P,ARG,FLT,FLTE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]-c2[i]; } return r; }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator+(const Vector<CanonicalNumericType<P,FLT,FLTE>>& c1, const VectorFunctionModel<P,ARG,FLT,FLTE>& f2);
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator-(const Vector<CanonicalNumericType<P,FLT,FLTE>>& c1, const VectorFunctionModel<P,ARG,FLT,FLTE>& f2);
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator*(const VectorFunctionModel<P,ARG,FLT,FLTE>& f1, const CanonicalNumericType<P,FLT,FLTE>& c2) {
        VectorFunctionModel<P,ARG,FLT,FLTE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]*c2; } return r; }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator*(const CanonicalNumericType<P,FLT,FLTE>& c1, const VectorFunctionModel<P,ARG,FLT,FLTE>& f2) {
        VectorFunctionModel<P,ARG,FLT,FLTE> r=f2; for(SizeType i=0; i!=r.size(); ++i) { r[i]=c1*f2[i]; } return r; }

    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator+(const VectorFunctionModel<P,ARG,FLT,FLTE>& f1, const VectorMultivariateFunction<P>& f2) {
        return f1+factory(f1).create(f2); }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator-(const VectorFunctionModel<P,ARG,FLT,FLTE>& f1, const VectorMultivariateFunction<P>& f2) {
        return f1-factory(f1).create(f2); }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator*(const VectorFunctionModel<P,ARG,FLT,FLTE>& f1, const ScalarMultivariateFunction<P>& f2) {
        return f1*factory(f1).create(f2); }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator/(const VectorFunctionModel<P,ARG,FLT,FLTE>& f1, const ScalarMultivariateFunction<P>& f2) {
        return f1/factory(f1).create(f2); }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator+(const VectorMultivariateFunction<P>& f1, const VectorFunctionModel<P,ARG,FLT,FLTE>& f2) {
        return factory(f2).create(f1)+f2; }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator-(const VectorMultivariateFunction<P>& f1, const VectorFunctionModel<P,ARG,FLT,FLTE>& f2) {
        return factory(f2).create(f1)-f2; }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator*(const ScalarMultivariateFunction<P>& f1, const VectorFunctionModel<P,ARG,FLT,FLTE>& f2) {
        return factory(f2).create(f1)*f2; }
    friend inline VectorFunctionModel<P,ARG,FLT,FLTE> operator/(const ScalarMultivariateFunction<P>& f1, const VectorFunctionModel<P,ARG,FLT,FLTE>& f2) {
        return factory(f2).create(f1)/f2; }


  public:
    friend NormType norm(const VectorFunctionModel<P,ARG,FLT,FLTE>& f) {
        return f._ptr->_concrete_norm(); }
    friend VectorFunctionModel<P,ARG,FLT,FLTE> embed(const DomainType& d1, const VectorFunctionModel<P,ARG,FLT,FLTE>& f, const DomainType& d2) {
        return VectorFunctionModel<P,ARG,FLT,FLTE>(f._ptr->_concrete_embed(d1,d2)); }
    friend VectorFunctionModel<P,ARG,FLT,FLTE> embed(const DomainType& d, const VectorFunctionModel<P,ARG,FLT,FLTE>& f) {
        return embed(d,f,DomainType()); }
    friend VectorFunctionModel<P,ARG,FLT,FLTE> embed(const VectorFunctionModel<P,ARG,FLT,FLTE>& f, const BoxDomainType& d) {
        return embed(DomainType(),f,d); }
    friend VectorFunctionModel<P,ARG,FLT,FLTE> embed(const VectorFunctionModel<P,ARG,FLT,FLTE>& f, const IntervalDomainType& d) {
        return embed(f,DomainType(1,d)); }
    friend VectorFunctionModel<P,ARG,FLT,FLTE> restriction(const VectorFunctionModel<P,ARG,FLT,FLTE>& f, const DomainType& d) {
        return VectorFunctionModel<P,ARG,FLT,FLTE>(f._ptr->_restriction(d)); }

    friend Vector<CanonicalNumericType<P,FLT,FLTE>> evaluate(const VectorFunctionModel<P,ARG,FLT,FLTE>& f, const Vector<CanonicalNumericType<P,FLT,FLTE>>& x) {
        return f._ptr->_call(x); }
    friend Vector<Number<P>> evaluate(const VectorFunctionModel<P,ARG,FLT,FLTE>& f, const Vector<Number<P>>& x) {
        return f._ptr->_call(Vector<CanonicalNumericType<P,FLT,FLTE>>(x,f.precision())); }

    friend Vector<CanonicalNumericType<P,FLT,FLTE>> unchecked_evaluate(const VectorFunctionModel<P,ARG,FLT,FLTE>& f, const Vector<CanonicalNumericType<P,FLT,FLTE>>& x) {
        return f._ptr->_unchecked_evaluate(x); }
    friend Vector<Number<P>> unchecked_evaluate(const VectorFunctionModel<P,ARG,FLT,FLTE>& f, const Vector<Number<P>>& x) {
        return f._ptr->_unchecked_evaluate(Vector<CanonicalNumericType<P,FLT,FLTE>>(x,f.precision())); }

    friend VectorFunctionModel<P,ARG,FLT,FLTE> partial_evaluate(const VectorFunctionModel<P,ARG,FLT,FLTE>& f, SizeType j, const CanonicalNumericType<P,FLT,FLTE>& c) {
        return VectorFunctionModel<P,ARG,FLT,FLTE>(f._ptr->_concrete_partial_evaluate(j,c)); }
    friend VectorFunctionModel<P,ARG,FLT,FLTE> partial_evaluate(const VectorFunctionModel<P,ARG,FLT,FLTE>& f, SizeType j, const Number<P>& c) {
        return partial_evaluate(f,j,CanonicalNumericType<P,FLT,FLTE>(c,f.precision())); }

    friend OutputStream& operator<<(OutputStream& os, const VectorFunctionModel<P,ARG,FLT,FLTE>& f) {
        return os <<  f.operator VectorMultivariateFunction<P>(); }

    friend ScalarFunctionModel<P,ARG,FLT,FLTE> unchecked_compose(const ScalarMultivariateFunction<P>& f, const VectorFunctionModel<P,ARG,FLT,FLTE>& g) {
        ScalarFunctionModelInterface<P,ARG,FLT,FLTE> const* fptr = dynamic_cast<ScalarFunctionModelInterface<P,ARG,FLT,FLTE> const*>(f.raw_pointer());
        if(fptr) { return unchecked_compose(ScalarFunctionModel<P,ARG,FLT,FLTE>(*fptr),g); } else { return compose(f,g); } }
    friend VectorFunctionModel<P,ARG,FLT,FLTE> unchecked_compose(const VectorMultivariateFunction<P>& f, const VectorFunctionModel<P,ARG,FLT,FLTE>& g) {
        VectorFunctionModelInterface<P,ARG,FLT,FLTE> const* fptr = dynamic_cast<VectorFunctionModelInterface<P,ARG,FLT,FLTE> const*>(f.raw_pointer());
        if(fptr) { return unchecked_compose(VectorFunctionModel<P,ARG,FLT,FLTE>(*fptr),g); } else { return compose(f,g); } }

    friend VectorFunctionModel<P,ARG,FLT,FLTE> join(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f1, const VectorFunctionModel<P,ARG,FLT,FLTE>& f2) {
        return join(VectorFunctionModel<P,ARG,FLT,FLTE>(1u,f1),f2); }
    friend VectorFunctionModel<P,ARG,FLT,FLTE> join(const VectorFunctionModel<P,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f2) {
        VectorFunctionModel<P,ARG,FLT,FLTE> r=VectorFunctionModel<P,ARG,FLT,FLTE>(f1._ptr->_clone()); r._ptr->_concrete_adjoin(f2); return r; }
    friend VectorFunctionModel<P,ARG,FLT,FLTE> join(const VectorFunctionModel<P,ARG,FLT,FLTE>& f1, const VectorFunctionModel<P,ARG,FLT,FLTE>& f2) {
        return VectorFunctionModel<P,ARG,FLT,FLTE>(f1._ptr->_concrete_join(*f2._ptr)); }

    friend VectorFunctionModel<P,ARG,FLT,FLTE> combine(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f2) {
        return VectorFunctionModel<P,ARG,FLT,FLTE>(1,f1)._ptr->_concrete_combine(VectorFunctionModel<P,ARG,FLT,FLTE>(1,f2)); };
    friend VectorFunctionModel<P,ARG,FLT,FLTE> combine(const ScalarFunctionModel<P,ARG,FLT,FLTE>& f1, const VectorFunctionModel<P,ARG,FLT,FLTE>& f2) {
        return VectorFunctionModel<P,ARG,FLT,FLTE>(1,f1)._ptr->_concrete_combine(f2); };
    friend VectorFunctionModel<P,ARG,FLT,FLTE> combine(const VectorFunctionModel<P,ARG,FLT,FLTE>& f1, const ScalarFunctionModel<P,ARG,FLT,FLTE>& f2) {
        return VectorFunctionModel<P,ARG,FLT,FLTE>(f1._ptr->_concrete_combine(VectorFunctionModel<P,ARG,FLT,FLTE>(1,f2))); };
    friend VectorFunctionModel<P,ARG,FLT,FLTE> combine(const VectorFunctionModel<P,ARG,FLT,FLTE>& f1, const VectorFunctionModel<P,ARG,FLT,FLTE>& f2) {
        return VectorFunctionModel<P,ARG,FLT,FLTE>(f1._ptr->_concrete_combine(f2)); }

    friend inline VectorFunctionModel<ValidatedTag,ARG,FLT,FLTE> refinement(const VectorFunctionModel<ValidatedTag,ARG,FLT,FLTE>& f1, const VectorFunctionModel<ValidatedTag,ARG,FLT,FLTE>& f2) {
        ARIADNE_ASSERT_MSG(f1.size()==f2.size(),"refinement(f1,f2): f1="<<f1<<", f2="<<f2<<")");
        VectorFunctionModel<ValidatedTag,ARG,FLT,FLTE> r=+f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=refinement(f1[i],f2[i]); } return r; }

    friend VectorFunctionModel<P,ARG,FLT,FLTE> antiderivative(const VectorFunctionModel<P,ARG,FLT,FLTE>& f, SizeType j) {
        VectorFunctionModel<P,ARG,FLT,FLTE> r(f);
        for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j); }
        return r;
    }

    friend VectorFunctionModel<P,ARG,FLT,FLTE> antiderivative(const VectorFunctionModel<P,ARG,FLT,FLTE>& f, SizeType j, const CanonicalNumericType<P,FLT,FLTE>& c) {
        VectorFunctionModel<P,ARG,FLT,FLTE> r(f);
        for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j,c); }
        return r;
    }

    friend VectorFunctionModel<P,ARG,FLT,FLTE> antiderivative(const VectorFunctionModel<P,ARG,FLT,FLTE>& f, SizeType j, const Number<P>& c) {
        return antiderivative(f,j,CanonicalNumericType<P,FLT,FLTE>(c,f[0].value().precision()));
    }

};

// FIXME: Implement for Multiple-Precision versions
template<class P> inline CanonicalNumericDPType<P> unchecked_evaluate(const ScalarMultivariateFunction<P>& f, const Vector<CanonicalNumericDPType<P>>& x) {
    auto const* fmptr = dynamic_cast<typename ScalarMultivariateFunctionModelDP<P>::Interface const*>(f.raw_pointer());
    if(fmptr) { return unchecked_evaluate(ScalarMultivariateFunctionModelDP<P>(*fmptr),x); }
    auto const* fptr = dynamic_cast<typename ScalarMultivariateFunctionPatch<P>::Interface const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(ScalarMultivariateFunctionPatch<P>(*fptr),x); }
    return evaluate(f,x);
}

template<class P> inline Vector<CanonicalNumericDPType<P>> unchecked_evaluate(const VectorMultivariateFunction<P>& f, const Vector<CanonicalNumericDPType<P>>& x) {
    auto const* fmptr = dynamic_cast<typename VectorMultivariateFunctionModelDP<P>::Interface const*>(f.raw_pointer());
    if(fmptr) { return unchecked_evaluate(VectorMultivariateFunctionModelDP<P>(*fmptr),x); }
    auto const* fptr = dynamic_cast<typename VectorMultivariateFunctionPatch<P>::Interface const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(VectorMultivariateFunctionPatch<P>(*fptr),x); }
    return evaluate(f,x);
}

template<class P, class ARG, class FLT, class FLTE> inline ScalarFunctionModel<P,ARG,FLT,FLTE> unchecked_compose(const ScalarMultivariateFunction<P>& f, const VectorFunctionModel<P,ARG,FLT,FLTE>& g) {
    auto const* fptr = dynamic_cast<typename ScalarFunctionModel<P,ARG,FLT,FLTE>::Interface const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(ScalarFunctionModel<P,ARG,FLT,FLTE>(*fptr),g); } else { return compose(f,g); } }
template<class P, class ARG, class FLT, class FLTE> inline VectorFunctionModel<P,ARG,FLT,FLTE> unchecked_compose(const VectorMultivariateFunction<P>& f, const VectorFunctionModel<P,ARG,FLT,FLTE>& g) {
    auto const* fptr = dynamic_cast<typename VectorFunctionModel<P,ARG,FLT,FLTE>::Interface const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(VectorFunctionModel<P,ARG,FLT,FLTE>(*fptr),g); } else { return compose(f,g); } }

// FIXME: Should be unneeded
template<class ARG, class FLT, class FLTE> ScalarFunctionModel<ValidatedTag,ARG,FLT,FLTE> unchecked_compose(const ScalarMultivariateFunctionModel<ValidatedTag,FLT,FLTE>& f, const VectorFunctionModel<ValidatedTag,ARG,FLT,FLTE>& g) {
    return VectorFunctionModel<ValidatedTag,ARG,FLT,FLTE>(g._ptr->_unchecked_compose(f)); }
template<class ARG, class FLT, class FLTE> VectorFunctionModel<ValidatedTag,ARG,FLT,FLTE> unchecked_compose(const VectorMultivariateFunctionModel<ValidatedTag,FLT,FLTE>& f, const VectorFunctionModel<ValidatedTag,ARG,FLT,FLTE>& g) {
    return VectorFunctionModel<ValidatedTag,ARG,FLT,FLTE>(g._ptr->_unchecked_compose(f)); }


template<class P, class ARG, class FLT, class FLTE> VectorFunctionModelElement<P,ARG,FLT,FLTE>::operator const ScalarFunctionModel<P,ARG,FLT,FLTE> () const {
    return ScalarFunctionModel<P,ARG,FLT,FLTE>(_p->get(_i)); }



template<class P, class... ARGS> template<class FLT, class FLTE>
FunctionPatch<P,RealScalar(ARGS...)>::FunctionPatch(FunctionModel<P,RealScalar(ARGS...),FLT,FLTE> fm)
    : _ptr(fm.pointer()) { }

template<class P, class... ARGS> template<class FLT, class FLTE>
FunctionPatch<P,RealVector(ARGS...)>::FunctionPatch(FunctionModel<P,RealVector(ARGS...),FLT,FLTE> fm)
    : _ptr(fm.pointer()) { }

} // namespace Ariadne

#endif
