/***************************************************************************
 *            function_model_interface.hpp
 *
 *  Copyright 2011-17  Pieter Collins
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

/*! \file function_model_interface.hpp
 *  \brief Interface for functions on bounded sets.
 */

#ifndef ARIADNE_FUNCTION_MODEL_INTERFACE_HPP
#define ARIADNE_FUNCTION_MODEL_INTERFACE_HPP

#include "../function/function.decl.hpp"
#include "../function/function_interface.hpp"

#include "../numeric/operators.hpp"

namespace Ariadne {

template<class P, class PR, class PRE> class FunctionModelFactoryInterface;
template<class P, class D, class PR, class PRE> class FunctionModelCreatorInterface;

template<class P, class D, class C, class PR, class PRE> class FunctionModelInterface;

template<class P, class D, class PR, class PRE> class FunctionModelInterface<P,D,IntervalDomainType,PR,PRE>
    : public virtual FunctionInterface<P,D,IntervalDomainType>
{
    static_assert(IsSame<D,IntervalDomainType>::value or IsSame<D,BoxDomainType>::value,"");
    typedef IntervalDomainType C;
  public:
    typedef D DomainType;
    typedef C CodomainType;
    typedef Interval<FloatUpperBound<PR>> RangeType;
    typedef FloatError<PR> NormType;
  public:
    virtual RangeType range() const = 0;

    virtual CanonicalCoefficientType<P,PR> const& value() const = 0;
    virtual CanonicalCoefficientType<P,PR> const gradient_value(SizeType i) const = 0;
    virtual CanonicalErrorType<P,PRE> const& error() const = 0;

    virtual Void set_error(const CanonicalErrorType<P,PRE>& e) = 0;
    virtual Void clobber() = 0;

    virtual NormType const _norm() const = 0;

    virtual ScalarFunctionModelInterface<P,D,PR,PRE>* _apply(OperatorCode op) const = 0;
    virtual CanonicalNumericType<P,PR,PRE> _unchecked_evaluate(const Vector<CanonicalNumericType<P,PR,PRE>>& x) const = 0;
    virtual ScalarFunctionModelInterface<P,D,PR,PRE>* _partial_evaluate(SizeType j, const CanonicalNumericType<P,PR,PRE>& c) const = 0;

    virtual FunctionModelFactoryInterface<P,PR,PRE>* _factory() const = 0;
    virtual ScalarFunctionModelInterface<P,D,PR,PRE>* _clone() const = 0;
    virtual ScalarFunctionModelInterface<P,D,PR,PRE>* _create() const = 0;
    virtual ScalarFunctionModelInterface<P,D,PR,PRE>* _embed(const DomainType& d1, const DomainType& d2) const = 0;
    virtual ScalarFunctionModelInterface<P,D,PR,PRE>* _restriction(const DomainType& d) const = 0;

    virtual ScalarFunctionModelInterface<P,D,PR,PRE>* _derivative(ElementIndexType<D> j) const = 0;
    virtual ScalarFunctionModelInterface<P,D,PR,PRE>* _antiderivative(SizeType j) const = 0;
    virtual ScalarFunctionModelInterface<P,D,PR,PRE>* _antiderivative(SizeType j, CanonicalNumericType<P,PR,PRE> c) const = 0;

    virtual Boolean _refines(const ScalarFunctionModelInterface<P,D,PR,PRE>& f) const = 0;
    virtual Boolean _inconsistent(const ScalarFunctionModelInterface<P,D,PR,PRE>& f) const = 0;
    virtual ScalarFunctionModelInterface<P,D,PR,PRE>* _refinement(const ScalarFunctionModelInterface<P,D,PR,PRE>& f) const = 0;

    virtual Void _iadd(const CanonicalNumericType<P,PR,PRE>& c) = 0;
    virtual Void _imul(const CanonicalNumericType<P,PR,PRE>& c) = 0;
    virtual Void _isma(const CanonicalNumericType<P,PR,PRE>& c, const ScalarFunctionModelInterface<P,D,PR,PRE>& f) = 0;
    virtual Void _ifma(const ScalarFunctionModelInterface<P,D,PR,PRE>& f1, const ScalarFunctionModelInterface<P,D,PR,PRE>& f2) = 0;
};


template<class P, class D, class PR, class PRE> class FunctionModelInterface<P,D,BoxDomainType,PR,PRE>
    : public virtual FunctionInterface<P,D,BoxDomainType>
{
    static_assert(IsSame<D,IntervalDomainType>::value or IsSame<D,BoxDomainType>::value,"");
    typedef BoxDomainType C;
  public:
    typedef D DomainType;
    typedef C CodomainType;
    typedef Box<Interval<FloatUpperBound<PR>>> RangeType;
    typedef FloatError<PR> NormType;
  public:
    virtual RangeType const range() const = 0;
    virtual Vector<CanonicalErrorType<P,PRE>> const errors() const = 0;
    virtual CanonicalErrorType<P,PRE> const error() const = 0;
    virtual Void clobber() = 0;

    virtual NormType const _norm() const = 0;

    virtual FunctionModelFactoryInterface<P,PR,PRE>* _factory() const = 0;
    virtual VectorFunctionModelInterface<P,D,PR,PRE>* _clone() const = 0;
    virtual Void _set(SizeType, ScalarFunctionModelInterface<P,D,PR,PRE> const&) = 0;
    virtual ScalarFunctionModelInterface<P,D,PR,PRE>* _get(SizeType) const = 0;
    virtual VectorFunctionModelInterface<P,D,PR,PRE>* _embed(const DomainType& d1, const DomainType& d2) const = 0;
    virtual VectorFunctionModelInterface<P,D,PR,PRE>* _restriction(const DomainType& d) const = 0;
    virtual VectorFunctionModelInterface<P,D,PR,PRE>* _join(const VectorFunctionModelInterface<P,D,PR,PRE>& f2) const = 0;
    virtual VectorFunctionModelInterface<P,D,PR,PRE>* _combine(const VectorFunctionModelInterface<P,D,PR,PRE>& f2) const = 0;
    virtual Void _adjoin(const ScalarFunctionModelInterface<P,D,PR,PRE>& f2) = 0;
    virtual Vector<CanonicalNumericType<P,PR,PRE>> _unchecked_evaluate(const Vector<CanonicalNumericType<P,PR,PRE>>& x) const = 0;
    virtual ScalarFunctionModelInterface<P,D,PR,PRE>* _compose(const ScalarFunctionInterface<P>& f) const = 0;
    virtual VectorFunctionModelInterface<P,D,PR,PRE>* _compose(const VectorFunctionInterface<P>& f) const = 0;
    virtual ScalarFunctionModelInterface<P,D,PR,PRE>* _unchecked_compose(const ScalarFunctionInterface<P>& f) const = 0;
    virtual VectorFunctionModelInterface<P,D,PR,PRE>* _unchecked_compose(const VectorFunctionInterface<P>& f) const = 0;
    virtual VectorFunctionModelInterface<P,D,PR,PRE>* _partial_evaluate(SizeType j, const CanonicalNumericType<P,PR,PRE>& c) const = 0;
    virtual Void restrict(const DomainType& d) = 0;
};


template<class P, class PR, class PRE> class FunctionModelFactoryInterface
{
    typedef BoxDomainType VD;
    typedef IntervalDomainType SD;
    friend class FunctionModelFactory<P,PR,PRE>;
  public:
    typedef SD ScalarDomainType;
    typedef VD VectorDomainType;
    virtual ~FunctionModelFactoryInterface<P,PR,PRE>() = default;
    virtual FunctionModelFactoryInterface<P,PR,PRE>* clone() const = 0;
    virtual OutputStream& _write(OutputStream& os) const = 0;
    friend OutputStream& operator<<(OutputStream& os, FunctionModelFactoryInterface<P,PR,PRE> const& factory) { factory._write(os); return os; }
  private:
    virtual CanonicalNumericType<P,PR,PRE> _create(const Number<P>& number) const = 0;
    virtual ScalarFunctionModelInterface<P,VD,PR,PRE>* _create(const VectorDomainType& domain, const ScalarFunctionInterface<P,VD>& function) const = 0;
    virtual VectorFunctionModelInterface<P,VD,PR,PRE>* _create(const VectorDomainType& domain, const VectorFunctionInterface<P,VD>& function) const = 0;

    virtual ScalarFunctionModelInterface<P,VD,PR,PRE>* _create_zero(const VectorDomainType& domain) const = 0;
    virtual ScalarFunctionModelInterface<P,VD,PR,PRE>* _create_constant(const VectorDomainType& domain, const Number<P>& value) const = 0;
    virtual ScalarFunctionModelInterface<P,VD,PR,PRE>* _create_coordinate(const VectorDomainType& domain, SizeType index) const = 0;
    virtual VectorFunctionModelInterface<P,VD,PR,PRE>* _create_zeros(SizeType result_size, const VectorDomainType& domain) const = 0;
    virtual VectorFunctionModelInterface<P,VD,PR,PRE>* _create_constants(const VectorDomainType& domain, const Vector<Number<P>>& values) const = 0;
    virtual VectorFunctionModelInterface<P,VD,PR,PRE>* _create_identity(const VectorDomainType& domain) const = 0;
  public:
    CanonicalNumericType<P,PR,PRE> create(const Number<P>& number) const {
        return CanonicalNumericType<P,PR,PRE>(this->_create(number)); }
    ScalarFunctionModel<P,VD,PR,PRE> create(const VectorDomainType& domain, const ScalarFunctionInterface<P,VD>& function) const {
        return ScalarFunctionModel<P,VD,PR,PRE>(this->_create(domain,function)); }
    VectorFunctionModel<P,VD,PR,PRE> create(const VectorDomainType& domain, const VectorFunctionInterface<P,VD>& function) const {
        return VectorFunctionModel<P,VD,PR,PRE>(this->_create(domain,function)); }

    CanonicalNumericType<P,PR,PRE> create_number(const Number<P>& number) const {
        return CanonicalNumericType<P,PR,PRE>(this->_create(number)); }

    ScalarFunctionModel<P,VD,PR,PRE> create_zero(const VectorDomainType& domain) const {
        return ScalarFunctionModel<P,VD,PR,PRE>(_create_zero(domain)); }
    ScalarFunctionModel<P,VD,PR,PRE> create_constant(const VectorDomainType& domain, const Number<P>& value) const {
        return ScalarFunctionModel<P,VD,PR,PRE>(_create_constant(domain,value)); }
    ScalarFunctionModel<P,VD,PR,PRE> create_constant(const VectorDomainType& domain, const CanonicalNumericType<P,PR,PRE>& value) const {
        return ScalarFunctionModel<P,VD,PR,PRE>(_create_constant(domain,Number<P>(value))); }
    ScalarFunctionModel<P,VD,PR,PRE> create_coordinate(const VectorDomainType& domain, SizeType index) const {
        return ScalarFunctionModel<P,VD,PR,PRE>(_create_coordinate(domain,index)); }
    VectorFunctionModel<P,VD,PR,PRE> create_zeros(SizeType result_size_, const VectorDomainType& domain) const {
        return VectorFunctionModel<P,VD,PR,PRE>(_create_zeros(result_size_,domain)); }
    VectorFunctionModel<P,VD,PR,PRE> create_constants(const VectorDomainType& domain, const Vector<Number<P>>& values) const {
        return VectorFunctionModel<P,VD,PR,PRE>(_create_constants(domain,values)); }
    VectorFunctionModel<P,VD,PR,PRE> create_identity(const VectorDomainType& domain) const {
        return VectorFunctionModel<P,VD,PR,PRE>(_create_identity(domain)); }

    ScalarFunctionModel<P,VD,PR,PRE> create_identity(const ScalarDomainType& domain) const;
};


} // namespace Ariadne

#endif
