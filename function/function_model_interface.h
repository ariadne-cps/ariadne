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

#include "function/function.decl.h"
#include "function/function_interface.h"

#include "numeric/operators.h"
#include "geometry/box.decl.h"

namespace Ariadne {

template<class P, class PR, class PRE> class ScalarFunctionModelInterface
    : public virtual ScalarFunctionInterface<P>
{
  public:
    typedef ExactBoxType DomainType;
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

    virtual ScalarFunctionModelInterface<P,PR,PRE>* _apply(OperatorCode op) const = 0;
    virtual CanonicalNumericType<P,PR,PRE> _unchecked_evaluate(const Vector<CanonicalNumericType<P,PR,PRE>>& x) const = 0;
    virtual ScalarFunctionModelInterface<P,PR,PRE>* _partial_evaluate(SizeType j, const CanonicalNumericType<P,PR,PRE>& c) const = 0;

    virtual ScalarFunctionModelInterface<P,PR,PRE>* _clone() const = 0;
    virtual ScalarFunctionModelInterface<P,PR,PRE>* _create() const = 0;
    virtual VectorFunctionModelInterface<P,PR,PRE>* _create_vector(SizeType i) const = 0;
    virtual ScalarFunctionModelInterface<P,PR,PRE>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const = 0;
    virtual VectorFunctionModelInterface<P,PR,PRE>* _create_identity() const = 0;
    virtual ScalarFunctionModelInterface<P,PR,PRE>* _restriction(const ExactBoxType& d) const = 0;

    virtual ScalarFunctionModelInterface<P,PR,PRE>* _create_zero(DomainType const& dom) const = 0;
    virtual ScalarFunctionModelInterface<P,PR,PRE>* _create_constant(DomainType const& dom, Number<P> const& c) const = 0;
    virtual ScalarFunctionModelInterface<P,PR,PRE>* _create_coordinate(DomainType const& dom, SizeType j) const = 0;

    virtual ScalarFunctionModelInterface<P,PR,PRE>* _derivative(SizeType j) const = 0;
    virtual ScalarFunctionModelInterface<P,PR,PRE>* _antiderivative(SizeType j) const = 0;
    virtual ScalarFunctionModelInterface<P,PR,PRE>* _antiderivative(SizeType j, CanonicalNumericType<P,PR,PRE> c) const = 0;

    virtual Boolean _refines(const ScalarFunctionModelInterface<P,PR,PRE>& f) const = 0;
    virtual Boolean _inconsistent(const ScalarFunctionModelInterface<P,PR,PRE>& f) const = 0;
    virtual ScalarFunctionModelInterface<P,PR,PRE>* _refinement(const ScalarFunctionModelInterface<P,PR,PRE>& f) const = 0;

    virtual Void _iadd(const CanonicalNumericType<P,PR,PRE>& c) = 0;
    virtual Void _imul(const CanonicalNumericType<P,PR,PRE>& c) = 0;
    virtual Void _isma(const CanonicalNumericType<P,PR,PRE>& c, const ScalarFunctionModelInterface<P,PR,PRE>& f) = 0;
    virtual Void _ifma(const ScalarFunctionModelInterface<P,PR,PRE>& f1, const ScalarFunctionModelInterface<P,PR,PRE>& f2) = 0;
};


template<class P, class PR, class PRE> class VectorFunctionModelInterface
    : public virtual VectorFunctionInterface<P>
{
    typedef ExactBoxType DomainType;
    typedef Box<Interval<FloatUpperBound<PR>>> RangeType;
    typedef FloatError<PR> NormType;
  public:
    virtual RangeType const range() const = 0;
    virtual Vector<CanonicalErrorType<P,PRE>> const errors() const = 0;
    virtual CanonicalErrorType<P,PRE> const error() const = 0;
    virtual Void clobber() = 0;

    virtual NormType const _norm() const = 0;

    virtual VectorFunctionModelInterface<P,PR,PRE>* _clone() const = 0;
    virtual ScalarFunctionModelInterface<P,PR,PRE>* _create_zero() const = 0;
    virtual VectorFunctionModelInterface<P,PR,PRE>* _create_identity() const = 0;
    virtual Void _set(SizeType, ScalarFunctionModelInterface<P,PR,PRE> const&) = 0;
    virtual ScalarFunctionModelInterface<P,PR,PRE>* _get(SizeType) const = 0;
    virtual VectorFunctionModelInterface<P,PR,PRE>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const = 0;
    virtual VectorFunctionModelInterface<P,PR,PRE>* _restriction(const ExactBoxType& d) const = 0;
    virtual VectorFunctionModelInterface<P,PR,PRE>* _join(const VectorFunctionModelInterface<P,PR,PRE>& f2) const = 0;
    virtual VectorFunctionModelInterface<P,PR,PRE>* _combine(const VectorFunctionModelInterface<P,PR,PRE>& f2) const = 0;
    virtual Void _adjoin(const ScalarFunctionModelInterface<P,PR,PRE>& f2) = 0;
    virtual Vector<CanonicalNumericType<P,PR,PRE>> _unchecked_evaluate(const Vector<CanonicalNumericType<P,PR,PRE>>& x) const = 0;
    virtual ScalarFunctionModelInterface<P,PR,PRE>* _compose(const ScalarFunctionInterface<P>& f) const = 0;
    virtual VectorFunctionModelInterface<P,PR,PRE>* _compose(const VectorFunctionInterface<P>& f) const = 0;
    virtual ScalarFunctionModelInterface<P,PR,PRE>* _unchecked_compose(const ScalarFunctionInterface<P>& f) const = 0;
    virtual VectorFunctionModelInterface<P,PR,PRE>* _unchecked_compose(const VectorFunctionInterface<P>& f) const = 0;
    virtual VectorFunctionModelInterface<P,PR,PRE>* _partial_evaluate(SizeType j, const CanonicalNumericType<P,PR,PRE>& c) const = 0;
    virtual Void restrict(const ExactBoxType& d) = 0;
};

template<class P, class PR, class PRE> class FunctionModelFactoryInterface
{
    typedef ExactBoxType DomainType;
  public:
    virtual FunctionModelFactoryInterface<P,PR,PRE>* clone() const = 0;
    virtual Void write(OutputStream& os) const = 0;
    ScalarFunctionModel<P,PR,PRE> create(const ExactBoxType& domain, const ScalarFunctionInterface<P>& function) const;
    VectorFunctionModel<P,PR,PRE> create(const ExactBoxType& domain, const VectorFunctionInterface<P>& function) const;
    ScalarFunctionModel<P,PR,PRE> create_zero(const ExactBoxType& domain) const;
    VectorFunctionModel<P,PR,PRE> create_zeros(SizeType result_size, const ExactBoxType& domain) const;
    ScalarFunctionModel<P,PR,PRE> create_constant(const ExactBoxType& domain, const Number<P>& value) const;
    ScalarFunctionModel<P,PR,PRE> create_constant(const ExactBoxType& domain, const CanonicalNumericType<P,PR,PRE>& value) const;
    VectorFunctionModel<P,PR,PRE> create_constants(const ExactBoxType& domain, const Vector<Number<P>>& values) const;
    VectorFunctionModel<P,PR,PRE> create_constants(const ExactBoxType& domain, const Vector<CanonicalNumericType<P,PR,PRE>>& values) const;
    ScalarFunctionModel<P,PR,PRE> create_coordinate(const ExactBoxType& domain, SizeType index) const;
    ScalarFunctionModel<P,PR,PRE> create_identity(const ExactIntervalType& domain) const;
    VectorFunctionModel<P,PR,PRE> create_identity(const ExactBoxType& domain) const;
    CanonicalNumericType<P,PR,PRE> create_number(const Number<P>& number) const;
    friend OutputStream& operator<<(OutputStream& os, FunctionModelFactoryInterface<P,PR,PRE> const& factory) {
        factory.write(os); return os; }
  private:
    virtual CanonicalNumericType<P,PR,PRE> _create(const Number<P>& number) const = 0;
    virtual ScalarFunctionModelInterface<P,PR,PRE>* _create(const ExactBoxType& domain, const ScalarFunctionInterface<P>& function) const = 0;
    virtual VectorFunctionModelInterface<P,PR,PRE>* _create(const ExactBoxType& domain, const VectorFunctionInterface<P>& function) const = 0;


};

} // namespace Ariadne

#endif
