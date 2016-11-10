/***************************************************************************
 *            function_model_interface.h
 *
 *  Copyright 2011  Pieter Collins
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

/*! \file function_model_interface.h
 *  \brief Interface for functions on bounded sets.
 */

#ifndef ARIADNE_FUNCTION_MODEL_INTERFACE_H
#define ARIADNE_FUNCTION_MODEL_INTERFACE_H

#include "function/function_interface.h"

#include "numeric/operators.h"
#include "geometry/box.decl.h"

namespace Ariadne {

template<class P> struct FunctionModelTraits;

template<> struct FunctionModelTraits<ValidatedTag> {
    typedef Float64Value CoefficientType;
    typedef Float64Error ErrorType;
    typedef ValidatedNumericType NumericType;
};

template<class P> using CanonicalCoefficientType = typename FunctionModelTraits<P>::CoefficientType;
template<class P> using CanonicalErrorType = typename FunctionModelTraits<P>::ErrorType;

template<class P> class ScalarFunctionModelInterface;
template<class P> class VectorFunctionModelInterface;

template<class P> class ScalarFunctionModelInterface
    : public virtual ScalarFunctionInterface<P>
{
    typedef ExactBoxType DomainType;
  public:
    virtual UpperIntervalType range() const = 0;

    virtual CanonicalCoefficientType<P> const& value() const = 0;
    virtual CanonicalCoefficientType<P> const gradient_value(SizeType i) const = 0;
    virtual CanonicalErrorType<P> const& error() const = 0;

    virtual Void set_error(const CanonicalErrorType<P>& e) = 0;
    virtual Void clobber() = 0;

    virtual NormType const _norm() const = 0;

    virtual ScalarFunctionModelInterface<P>* _apply(OperatorCode op) const = 0;
    virtual CanonicalNumericType<P> _unchecked_evaluate(const Vector<CanonicalNumericType<P>>& x) const = 0;

    virtual ScalarFunctionModelInterface<P>* _clone() const = 0;
    virtual ScalarFunctionModelInterface<P>* _create() const = 0;
    virtual VectorFunctionModelInterface<P>* _create_vector(SizeType i) const = 0;
    virtual ScalarFunctionModelInterface<P>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const = 0;
    virtual VectorFunctionModelInterface<P>* _create_identity() const = 0;
    virtual ScalarFunctionModelInterface<P>* _restriction(const ExactBoxType& d) const = 0;

    virtual ScalarFunctionModelInterface<P>* _create_zero(DomainType const& dom) const = 0;
    virtual ScalarFunctionModelInterface<P>* _create_constant(DomainType const& dom, CanonicalNumericType<P> const& c) const = 0;
    virtual ScalarFunctionModelInterface<P>* _create_coordinate(DomainType const& dom, SizeType j) const = 0;

    virtual ScalarFunctionModelInterface<P>* _derivative(SizeType j) const = 0;
    virtual ScalarFunctionModelInterface<P>* _antiderivative(SizeType j) const = 0;
    virtual ScalarFunctionModelInterface<P>* _antiderivative(SizeType j, CanonicalNumericType<P> c) const = 0;

    virtual Boolean _refines(const ScalarFunctionModelInterface<P>& f) const = 0;
    virtual Boolean _inconsistent(const ScalarFunctionModelInterface<P>& f) const = 0;
    virtual ScalarFunctionModelInterface<P>* _refinement(const ScalarFunctionModelInterface<P>& f) const = 0;

    virtual Void _iadd(const CanonicalNumericType<P>& c) = 0;
    virtual Void _imul(const CanonicalNumericType<P>& c) = 0;
    virtual Void _isma(const CanonicalNumericType<P>& c, const ScalarFunctionModelInterface<P>& f) = 0;
    virtual Void _ifma(const ScalarFunctionModelInterface<P>& f1, const ScalarFunctionModelInterface<P>& f2) = 0;
};


template<class P> class VectorFunctionModelInterface;

template<class P> class VectorFunctionModelInterface
    : public virtual VectorFunctionInterface<P>
{
    typedef ExactBoxType DomainType;
  public:
    virtual UpperBoxType const range() const = 0;
    virtual Vector<CanonicalErrorType<P>> const errors() const = 0;
    virtual CanonicalErrorType<P> const error() const = 0;
    virtual Void clobber() = 0;

    virtual NormType const _norm() const = 0;

    virtual VectorFunctionModelInterface<P>* _clone() const = 0;
    virtual ScalarFunctionModelInterface<P>* _create_zero() const = 0;
    virtual VectorFunctionModelInterface<P>* _create_identity() const = 0;
    virtual Void _set(SizeType, ScalarFunctionModelInterface<P> const&) = 0;
    virtual ScalarFunctionModelInterface<P>* _get(SizeType) const = 0;
    virtual VectorFunctionModelInterface<P>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const = 0;
    virtual VectorFunctionModelInterface<P>* _join(const VectorFunctionModelInterface<P>& f2) const = 0;
    virtual VectorFunctionModelInterface<P>* _combine(const VectorFunctionModelInterface<P>& f2) const = 0;
    virtual Void _adjoin(const ScalarFunctionModelInterface<P>& f2) = 0;
    virtual Vector<ValidatedNumericType> _unchecked_evaluate(const Vector<CanonicalNumericType<P>>& x) const = 0;
    virtual ScalarFunctionModelInterface<P>* _compose(const ScalarFunctionInterface<P>& f) const = 0;
    virtual VectorFunctionModelInterface<P>* _compose(const VectorFunctionInterface<P>& f) const = 0;
    virtual ScalarFunctionModelInterface<P>* _unchecked_compose(const ScalarFunctionInterface<P>& f) const = 0;
    virtual VectorFunctionModelInterface<P>* _unchecked_compose(const VectorFunctionInterface<P>& f) const = 0;
    virtual VectorFunctionModelInterface<P>* _partial_evaluate(SizeType j, const CanonicalNumericType<P>& c) const = 0;
    virtual Void restrict(const ExactBoxType& d) = 0;
};


template<class P> class FunctionModelFactoryInterface
{
    typedef ExactBoxType DomainType;
  public:
    virtual FunctionModelFactoryInterface<P>* clone() const = 0;
    virtual Void write(OutputStream& os) const = 0;
    ScalarFunctionModel<P> create(const ExactBoxType& domain, const ScalarFunctionInterface<P>& function) const;
    VectorFunctionModel<P> create(const ExactBoxType& domain, const VectorFunctionInterface<P>& function) const;
    ScalarFunctionModel<P> create_zero(const ExactBoxType& domain) const;
    VectorFunctionModel<P> create_zeros(SizeType result_size, const ExactBoxType& domain) const;
    ScalarFunctionModel<P> create_constant(const ExactBoxType& domain, const Number<P>& value) const;
    ScalarFunctionModel<P> create_constant(const ExactBoxType& domain, const CanonicalNumericType<P>& value) const;
    VectorFunctionModel<P> create_constants(const ExactBoxType& domain, const Vector<Number<P>>& values) const;
    VectorFunctionModel<P> create_constants(const ExactBoxType& domain, const Vector<CanonicalNumericType<P>>& values) const;
    ScalarFunctionModel<P> create_coordinate(const ExactBoxType& domain, SizeType index) const;
    ScalarFunctionModel<P> create_identity(const ExactIntervalType& domain) const;
    VectorFunctionModel<P> create_identity(const ExactBoxType& domain) const;
    CanonicalNumericType<P> create_number(const Number<P>& number) const;
    friend OutputStream& operator<<(OutputStream& os, FunctionModelFactoryInterface<P> const& factory) {
        factory.write(os); return os; }
  private:
    virtual ScalarFunctionModelInterface<P>* _create(const ExactBoxType& domain, const ScalarFunctionInterface<P>& function) const = 0;
    virtual VectorFunctionModelInterface<P>* _create(const ExactBoxType& domain, const VectorFunctionInterface<P>& function) const = 0;


};

} // namespace Ariadne

#endif
