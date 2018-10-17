/***************************************************************************
 *            bounder.hpp
 *
 *  Copyright  2018  Luca Geretti
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

/*! \file bounder.hpp
 *  \brief Classes for finding the bounds of the flow for differential equations.
 */

#ifndef ARIADNE_BOUNDER_HPP
#define ARIADNE_BOUNDER_HPP

#include "../utility/typedefs.hpp"
#include "../utility/attribute.hpp"
#include "../algebra/algebra.hpp"
#include "../function/domain.hpp"
#include "../function/function_model.hpp"
#include "../output/logging.hpp"
#include "integrator_interface.hpp"

namespace Ariadne {

class BounderInterface {
  public:
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& dom, StepSizeType const& hsug) const = 0;
    virtual Void write(OutputStream& os) const = 0;

    virtual BounderInterface* clone() const = 0;

    friend inline OutputStream& operator<<(OutputStream& os, BounderInterface const& bounder) {
        bounder.write(os); return os; }

    virtual ~BounderInterface() = default;
};

class BounderBase : public BounderInterface, Loggable {
  public:
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& dom, StepSizeType const& hsug) const override;
  protected:
    virtual UpperBoxType formula(BoxDomainType const& D, BoxDomainType const& V, ValidatedVectorMultivariateFunction const& f, UpperBoxType const& B, StepSizeType const& h) const = 0;
  private:
    UpperBoxType _initial(BoxDomainType const& dom, ValidatedVectorMultivariateFunction const& f, UpperBoxType const& arg, StepSizeType const& h, PositiveFloatDPValue FORMULA_WIDENING) const;
    UpperBoxType _refinement(BoxDomainType const& dom, ValidatedVectorMultivariateFunction const& f, UpperBoxType const& B, StepSizeType const& h) const;
};

class EulerBounder final : public BounderBase {
  protected:
    virtual UpperBoxType formula(BoxDomainType const& D, BoxDomainType const& V, ValidatedVectorMultivariateFunction const& f, UpperBoxType const& B, StepSizeType const& h) const override;
  public:
    virtual Void write(OutputStream& os) const override { os << "Euler"; }
    virtual EulerBounder* clone() const override { return new EulerBounder(*this); }
};

class BounderHandle {
  private:
    SharedPointer<BounderInterface> _impl;
  public:
    BounderHandle(BounderInterface const& bounder) : _impl(bounder.clone()) { }
    BounderHandle(BounderHandle const& other) : _impl(other._impl) { }
    BounderHandle& operator=(BounderHandle const& other) { _impl = other._impl; return *this; }

    operator BounderInterface const& () const { return *_impl; }

    Void write(OutputStream& os) const { _impl->write(os); }
    BounderHandle* clone() const { return new BounderHandle(*this); }

    Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& dom, StepSizeType const& hsug) const {
        return _impl->compute(f,dom,hsug);
    }
};



} // namespace Ariadne

#endif /* ARIADNE_BOUNDER_HPP */
