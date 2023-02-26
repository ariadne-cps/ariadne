/***************************************************************************
 *            function/function_model_interface.hpp
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

/*! \file function/function_model_interface.hpp
 *  \brief Interface for functions on bounded sets.
 */

#ifndef ARIADNE_FUNCTION_MODEL_INTERFACE_HPP
#define ARIADNE_FUNCTION_MODEL_INTERFACE_HPP

#include "algebra/algebra_interface.hpp"

#include "function/function.decl.hpp"
#include "function/function_interface.hpp"
#include "function/function_patch_interface.hpp"

#include "algebra/range.hpp"
#include "numeric/operators.hpp"

namespace Ariadne {

template<class P, class F, class FE> class FunctionModelFactoryInterface;
template<class P, class ARG, class F, class FE> class FunctionModelCreatorInterface;

template<class P, class SIG, class F, class FE> class FunctionModelInterface;

template<class P, class SIG, class F, class FE> class FunctionModelAlgebraInterface;

template<class P, class ARG, class F, class FE> class FunctionModelAlgebraInterface<P,RealScalar(ARG),F,FE>
    : public virtual ElementaryAlgebraInterface<CanonicalNumericType<P,F,FE>>
{
    static_assert(Same<ARG,RealScalar> or Same<ARG,RealVector>);
    using RES = RealScalar; using SIG = RES(ARG);
    using D = BoundedDomainType<ARG>; using C = BoundedDomainType<RES>;
    using PR = PrecisionType<F>; using PRE = PrecisionType<FE>;
  public:
    typedef D DomainType;
    typedef C CodomainType;
    typedef typename FunctionModelTraits<P,F,FE>::ValueType ValueType;
  public:
    virtual ValueType const _value() const = 0;
};

template<class P, class ARG, class F, class FE> using ScalarFunctionModelInterface = FunctionModelInterface<P,RealScalar(ARG),F,FE>;
template<class P, class ARG, class F, class FE> using VectorFunctionModelInterface = FunctionModelInterface<P,RealVector(ARG),F,FE>;

template<class P, class ARG, class F, class FE> using ScalarFunctionModelAlgebraInterface = FunctionModelAlgebraInterface<P,RealScalar(ARG),F,FE>;
template<class P, class ARG, class F, class FE> using VectorFunctionModelAlgebraInterface = FunctionModelAlgebraInterface<P,RealVector(ARG),F,FE>;

template<class P, class ARG, class F, class FE> class FunctionModelAlgebraInterface<P,RealVector(ARG),F,FE>
    : public virtual WritableInterface
{
  public:
    typedef typename FunctionModelTraits<P,F,FE>::ValueType ValueType;
    typedef typename FunctionModelTraits<P,F,FE>::ErrorType ErrorType;
  public:
    virtual Void _concrete_set(SizeType, ScalarFunctionModelInterface<P,ARG,F,FE> const&) = 0;
    virtual ScalarFunctionModelInterface<P,ARG,F,FE>* _concrete_get(SizeType) const = 0;
    virtual VectorFunctionModelInterface<P,ARG,F,FE>* _concrete_join(const VectorFunctionModelInterface<P,ARG,F,FE>& f2) const = 0;
    virtual VectorFunctionModelInterface<P,ARG,F,FE>* _concrete_combine(const VectorFunctionModelInterface<P,ARG,F,FE>& f2) const = 0;
    virtual Void _concrete_adjoin(const ScalarFunctionModelInterface<P,ARG,F,FE>& f2) = 0;

    virtual Vector<ValueType> const _concrete_values() const = 0;
    virtual Vector<ErrorType> const _concrete_errors() const = 0;
};

template<class P, class SIG, class F, class FE> class FunctionModelInterface
    : public virtual FunctionInterface<P,SIG>
    , public virtual FunctionPatchInterface<P,SIG>
    , public virtual FunctionModelAlgebraInterface<P,SIG,F,FE>
{
    using RES=typename SignatureTraits<SIG>::ResultKind; using ARG=typename SignatureTraits<SIG>::ArgumentKind;
    static_assert(Same<ARG,RealScalar> or Same<ARG,RealVector>,"");
    static_assert(Same<RES,RealScalar> or Same<RES,RealVector>,"");
    using C=BoundedDomainType<RES>; using D=BoundedDomainType<ARG>;
    static_assert(Same<D,IntervalDomainType> or Same<D,BoxDomainType>,"");
    static_assert(Same<C,IntervalDomainType> or Same<C,BoxDomainType>,"");

    template<class X> using Argument = typename FunctionInterface<P,SIG>::template Argument<X>;
    template<class X> using Result = typename FunctionInterface<P,SIG>::template Result<X>;
  public:
    typedef D DomainType;
    typedef C CodomainType;
    typedef typename SignatureTraits<SIG>::template ConcreteRangeType<F> RangeType;
    typedef typename FunctionModelTraits<P,F,FE>::NormType NormType;
    typedef typename FunctionModelTraits<P,F,FE>::ValueType ValueType;
    typedef typename FunctionModelTraits<P,F,FE>::ErrorType ErrorType;
    typedef typename FunctionModelTraits<P,F,FE>::NumericType NumericType;
    typedef typename FunctionModelTraits<P,F,FE>::GenericNumericType GenericNumericType;
    typedef typename SignatureTraits<SIG>::ArgumentIndexType ArgumentIndexType;
  public:
    virtual DomainType const domain() const = 0;
    virtual CodomainType const codomain() const = 0;
    virtual RangeType const range() const = 0;
    virtual NormType const _concrete_norm() const = 0;
    virtual ErrorType const _concrete_error() const = 0;

    virtual Void clobber() = 0;

    virtual Result<CanonicalNumericType<P,F,FE>> _concrete_unchecked_evaluate(const Argument<CanonicalNumericType<P,F,FE>>& x) const = 0;
    virtual FunctionModelInterface<P,SIG,F,FE>* _concrete_partial_evaluate(SizeType j, const CanonicalNumericType<P,F,FE>& c) const = 0;
    virtual FunctionModelInterface<P,SIG,F,FE>* _embed(const DomainType& d1, const DomainType& d2) const = 0;

    virtual FunctionModelFactoryInterface<P,F,FE>* _factory() const = 0;
    virtual FunctionModelInterface<P,SIG,F,FE>* _clone() const = 0;
    virtual FunctionModelInterface<P,SIG,F,FE>* _create() const = 0;
    inline FunctionModelInterface<P,SIG,F,FE>* _copy() const { return this->_clone(); }

    virtual FunctionModelInterface<P,SIG,F,FE>* _restriction(const DomainType& d) const = 0;

    virtual FunctionModelInterface<P,SIG,F,FE>* _derivative(ArgumentIndexType j) const = 0;
    virtual FunctionModelInterface<P,SIG,F,FE>* _antiderivative(ArgumentIndexType j) const = 0;
    virtual FunctionModelInterface<P,SIG,F,FE>* _antiderivative(ArgumentIndexType j, Number<P>) const = 0;
    virtual FunctionModelInterface<P,SIG,F,FE>* _concrete_antiderivative(ArgumentIndexType j, CanonicalNumericType<P,F,FE> c) const = 0;

    using FunctionPatchInterface<P,SIG>::_compose;
    using FunctionPatchInterface<P,SIG>::_unchecked_compose;
    virtual ScalarFunctionModelInterface<P,ARG,F,FE>* _compose(const ScalarFunction<P,RES>& f) const = 0;
    virtual VectorFunctionModelInterface<P,ARG,F,FE>* _compose(const VectorFunction<P,RES>& f) const = 0;
    virtual ScalarFunctionModelInterface<P,ARG,F,FE>* _unchecked_compose(const ScalarFunction<P,RES>& f) const = 0;
    virtual VectorFunctionModelInterface<P,ARG,F,FE>* _unchecked_compose(const VectorFunction<P,RES>& f) const = 0;

    virtual Boolean _concrete_refines(const FunctionModelInterface<P,SIG,F,FE>& f) const = 0;
    virtual Boolean _concrete_inconsistent(const FunctionModelInterface<P,SIG,F,FE>& f) const = 0;
    virtual FunctionModelInterface<P,SIG,F,FE>* _concrete_refinement(const FunctionModelInterface<P,SIG,F,FE>& f) const = 0;

    virtual OutputStream& _write(OutputStream&) const = 0;

    friend OutputStream& operator<<(OutputStream& os, FunctionModelInterface<P,SIG,F,FE> const& f) {
        return f._write(os); }
};



template<class P, class F, class FE> class FunctionModelFactoryInterface
    : public FunctionPatchFactoryInterface<P>
{
    typedef RealScalar SARG;
    typedef RealVector VARG;

    typedef BoxDomainType VD;
    typedef IntervalDomainType SD;
    friend class FunctionModelFactory<P,F,FE>;
  public:
    typedef SD ScalarDomainType;
    typedef VD VectorDomainType;
    virtual ~FunctionModelFactoryInterface() = default;
    virtual FunctionModelFactoryInterface<P,F,FE>* clone() const = 0;
    inline FunctionModelFactoryInterface<P,F,FE>* _copy() const { return this->clone(); }
    virtual OutputStream& _write(OutputStream& os) const = 0;
    friend OutputStream& operator<<(OutputStream& os, FunctionModelFactoryInterface<P,F,FE> const& factory) { factory._write(os); return os; }
  private:
    using FunctionPatchFactoryInterface<P>::_create;
    virtual CanonicalNumericType<P,F,FE> _create(const Number<P>& number) const = 0;
    virtual ScalarFunctionModelInterface<P,VARG,F,FE>* _create(const VectorDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const = 0;
    virtual VectorFunctionModelInterface<P,VARG,F,FE>* _create(const VectorDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const = 0;

    virtual ScalarFunctionModelInterface<P,VARG,F,FE>* _create_zero(const VectorDomainType& domain) const = 0;
    virtual ScalarFunctionModelInterface<P,VARG,F,FE>* _create_constant(const VectorDomainType& domain, const Number<P>& value) const = 0;
    virtual ScalarFunctionModelInterface<P,VARG,F,FE>* _create_coordinate(const VectorDomainType& domain, SizeType index) const = 0;
    virtual VectorFunctionModelInterface<P,VARG,F,FE>* _create_zeros(SizeType rsize, const VectorDomainType& domain) const = 0;
    virtual VectorFunctionModelInterface<P,VARG,F,FE>* _create_constants(const VectorDomainType& domain, const Vector<Number<P>>& values) const = 0;
    virtual VectorFunctionModelInterface<P,VARG,F,FE>* _create_projection(const VectorDomainType& domain, Range indices) const = 0;
    virtual VectorFunctionModelInterface<P,VARG,F,FE>* _create_identity(const VectorDomainType& domain) const = 0;
  public:
    CanonicalNumericType<P,F,FE> create(const Number<P>& number) const {
        return CanonicalNumericType<P,F,FE>(this->_create(number)); }
    ScalarMultivariateFunctionModel<P,F,FE> create(const VectorDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const {
        return ScalarMultivariateFunctionModel<P,F,FE>(this->_create(domain,function)); }
    VectorMultivariateFunctionModel<P,F,FE> create(const VectorDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const {
        return VectorMultivariateFunctionModel<P,F,FE>(this->_create(domain,function)); }

    CanonicalNumericType<P,F,FE> create_number(const Number<P>& number) const {
        return CanonicalNumericType<P,F,FE>(this->_create(number)); }

    ScalarMultivariateFunctionModel<P,F,FE> create_zero(const VectorDomainType& domain) const {
        return ScalarMultivariateFunctionModel<P,F,FE>(_create_zero(domain)); }
    ScalarMultivariateFunctionModel<P,F,FE> create_constant(const VectorDomainType& domain, const Number<P>& value) const {
        return ScalarMultivariateFunctionModel<P,F,FE>(_create_constant(domain,value)); }
    ScalarMultivariateFunctionModel<P,F,FE> create_constant(const VectorDomainType& domain, const CanonicalNumericType<P,F,FE>& value) const {
        return ScalarMultivariateFunctionModel<P,F,FE>(_create_constant(domain,Number<P>(value))); }
    ScalarMultivariateFunctionModel<P,F,FE> create_coordinate(const VectorDomainType& domain, SizeType index) const {
        return ScalarMultivariateFunctionModel<P,F,FE>(_create_coordinate(domain,index)); }
    VectorMultivariateFunctionModel<P,F,FE> create_zeros(SizeType rsize, const VectorDomainType& domain) const {
        return VectorMultivariateFunctionModel<P,F,FE>(_create_zeros(rsize,domain)); }
    VectorMultivariateFunctionModel<P,F,FE> create_constants(const VectorDomainType& domain, const Vector<Number<P>>& values) const {
        return VectorMultivariateFunctionModel<P,F,FE>(_create_constants(domain,values)); }
    VectorMultivariateFunctionModel<P,F,FE> create_projection(const VectorDomainType& domain, Range indices) const {
        return VectorMultivariateFunctionModel<P,F,FE>(_create_projection(domain,indices)); }
    VectorMultivariateFunctionModel<P,F,FE> create_identity(const VectorDomainType& domain) const {
        return VectorMultivariateFunctionModel<P,F,FE>(_create_identity(domain)); }

    ScalarMultivariateFunctionModel<P,F,FE> create_identity(const ScalarDomainType& domain) const;
};


} // namespace Ariadne

#endif
