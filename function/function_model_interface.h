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

namespace Ariadne {

typedef Float64Value ValidatedCoefficientType;
typedef Float64Error ValidatedErrorType;

template<class P> class ScalarFunctionModelInterface;
template<class P> class VectorFunctionModelInterface;

template<> class ScalarFunctionModelInterface<ValidatedTag>
    : public virtual ScalarFunctionInterface<ValidatedTag>
{
  public:
    typedef Ariadne::ValidatedCoefficientType CoefficientType;
    typedef Ariadne::ValidatedErrorType ErrorType;
  public:
    virtual UpperIntervalType range() const = 0;

    virtual CoefficientType const& value() const = 0;
    virtual CoefficientType const gradient_value(SizeType i) const = 0;
    virtual ErrorType const& error() const = 0;

    virtual Void set_error(const ErrorType& e) = 0;
    virtual Void clobber() = 0;

    virtual NormType const _norm() const = 0;

    virtual ScalarFunctionModelInterface<ValidatedTag>* _apply(OperatorCode op) const = 0;
    virtual ValidatedNumericType _unchecked_evaluate(const Vector<ValidatedNumericType>& x) const = 0;

    virtual ScalarFunctionModelInterface<ValidatedTag>* _clone() const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _create() const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _create_vector(SizeType i) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _create_identity() const = 0;
    virtual Void restrict(const ExactBoxType& d) = 0;

    virtual ScalarFunctionModelInterface<ValidatedTag>* _create_zero(DomainType const& dom) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _create_constant(DomainType const& dom, ValidatedNumericType const& c) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _create_coordinate(DomainType const& dom, SizeType j) const = 0;

    virtual ScalarFunctionModelInterface<ValidatedTag>* _derivative(SizeType j) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _antiderivative(SizeType j) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _antiderivative(SizeType j, ValidatedNumericType c) const = 0;

    virtual Boolean _refines(const ScalarFunctionModelInterface<ValidatedTag>& f) const = 0;
    virtual Boolean _inconsistent(const ScalarFunctionModelInterface<ValidatedTag>& f) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _refinement(const ScalarFunctionModelInterface<ValidatedTag>& f) const = 0;

    virtual Void _iadd(const ValidatedNumericType& c) = 0;
    virtual Void _imul(const ValidatedNumericType& c) = 0;
    virtual Void _isma(const ValidatedNumericType& c, const ScalarFunctionModelInterface<ValidatedTag>& f) = 0;
    virtual Void _ifma(const ScalarFunctionModelInterface<ValidatedTag>& f1, const ScalarFunctionModelInterface<ValidatedTag>& f2) = 0;
};


template<class P> class VectorFunctionModelInterface;

template<> class VectorFunctionModelInterface<ValidatedTag>
    : public virtual VectorFunctionInterface<ValidatedTag>
{
  public:
    virtual UpperBoxType const range() const = 0;
    virtual Vector<ErrorType> const errors() const = 0;
    virtual ErrorType const error() const = 0;
    virtual Void clobber() = 0;

    virtual NormType const _norm() const = 0;

    virtual VectorFunctionModelInterface<ValidatedTag>* _clone() const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _create_zero() const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _create_identity() const = 0;
    virtual Void _set(SizeType, ScalarFunctionModelInterface<ValidatedTag> const&) = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _get(SizeType) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _join(const VectorFunctionModelInterface<ValidatedTag>& f2) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _combine(const VectorFunctionModelInterface<ValidatedTag>& f2) const = 0;
    virtual Void _adjoin(const ScalarFunctionModelInterface<ValidatedTag>& f2) = 0;
    virtual Vector<ValidatedNumericType> _unchecked_evaluate(const Vector<ValidatedNumericType>& x) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _compose(const ScalarFunctionInterface<ValidatedTag>& f) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _compose(const VectorFunctionInterface<ValidatedTag>& f) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _unchecked_compose(const ScalarFunctionInterface<ValidatedTag>& f) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _unchecked_compose(const VectorFunctionInterface<ValidatedTag>& f) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _partial_evaluate(SizeType j, const ValidatedNumericType& c) const = 0;
    virtual Void restrict(const ExactBoxType& d) = 0;
};


template<class X> class FunctionModelFactoryInterface;

template<> class FunctionModelFactoryInterface<ValidatedTag>
{
    typedef ExactBoxType DomainType;
  public:
    virtual FunctionModelFactoryInterface<ValidatedTag>* clone() const = 0;
    virtual Void write(OutputStream& os) const = 0;
    inline ScalarFunctionModel<ValidatedTag> create(const ExactBoxType& domain, const ScalarFunctionInterface<ValidatedTag>& function) const;
    inline VectorFunctionModel<ValidatedTag> create(const ExactBoxType& domain, const VectorFunctionInterface<ValidatedTag>& function) const;
    inline ScalarFunctionModel<ValidatedTag> create_zero(const ExactBoxType& domain) const;
    inline VectorFunctionModel<ValidatedTag> create_zeros(SizeType result_size, const ExactBoxType& domain) const;
    inline ScalarFunctionModel<ValidatedTag> create_constant(const ExactBoxType& domain, const ValidatedNumber& value) const;
    inline ScalarFunctionModel<ValidatedTag> create_constant(const ExactBoxType& domain, const ValidatedNumericType& value) const;
    inline VectorFunctionModel<ValidatedTag> create_constants(const ExactBoxType& domain, const Vector<ValidatedNumericType>& values) const;
    inline ScalarFunctionModel<ValidatedTag> create_coordinate(const ExactBoxType& domain, SizeType index) const;
    inline ScalarFunctionModel<ValidatedTag> create_identity(const ExactIntervalType& domain) const;
    inline VectorFunctionModel<ValidatedTag> create_identity(const ExactBoxType& domain) const;
    inline CanonicalNumericType<ValidatedTag> create_number(const ValidatedNumber& number) const;
    friend OutputStream& operator<<(OutputStream& os, FunctionModelFactoryInterface<ValidatedTag> const& factory) {
        factory.write(os); return os; }
  private:
    virtual ScalarFunctionModelInterface<ValidatedTag>* _create(const ExactBoxType& domain, const ScalarFunctionInterface<ValidatedTag>& function) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _create(const ExactBoxType& domain, const VectorFunctionInterface<ValidatedTag>& function) const = 0;


};

} // namespace Ariadne

#endif
