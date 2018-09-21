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

namespace Ariadne {

class BounderInterface {
  public:
    virtual Pair<PositiveFloatDPValue,UpperBoxType> compute(ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPApproximation hsug) const = 0;
    virtual Void write(OutputStream& os) const = 0;

    virtual BounderInterface* clone() const = 0;

    friend inline OutputStream& operator<<(OutputStream& os, const BounderInterface& bounder) {
        bounder.write(os); return os; }

    virtual ~BounderInterface() = default;
};

class BounderBase : public BounderInterface, Loggable {
  public:
    virtual Pair<PositiveFloatDPValue,UpperBoxType> compute(ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPApproximation hsug) const;
  protected:
    virtual UpperBoxType formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const = 0;
  private:
    UpperBoxType _initial(BoxDomainType dom, ValidatedVectorFunction f, UpperBoxType arg, PositiveFloatDPValue h, PositiveFloatDPValue FORMULA_WIDENING) const;
    UpperBoxType _refinement(BoxDomainType dom, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const;
  public:
    virtual ~BounderBase() = default;

};

class EulerBounder : public BounderBase {
  protected:
    virtual UpperBoxType formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const;
  public:
    virtual Void write(OutputStream& os) const { os << "Euler"; }
    virtual EulerBounder* clone() const { return new EulerBounder(*this); }
    virtual ~EulerBounder() = default;
};

class HeunBounder : public BounderBase {
  protected:
    virtual UpperBoxType formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const;
  public:
    virtual Void write(OutputStream& os) const { os << "Heun"; }
    virtual HeunBounder* clone() const { return new HeunBounder(*this); }
    virtual ~HeunBounder() = default;
};

class RalstonBounder : public BounderBase {
  protected:
    virtual UpperBoxType formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const;
  public:
    virtual Void write(OutputStream& os) const { os << "Ralston"; }
    virtual RalstonBounder* clone() const { return new RalstonBounder(*this); }
    virtual ~RalstonBounder() = default;
};

class RungeKutta4Bounder : public BounderBase {
  protected:
    virtual UpperBoxType formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const;
  public:
    virtual Void write(OutputStream& os) const { os << "RungeKutta4"; }
    virtual RungeKutta4Bounder* clone() const { return new RungeKutta4Bounder(*this); }
    virtual ~RungeKutta4Bounder() = default;
};

class BounderHandle : public BounderInterface {
    friend class BounderFactory;
  private:
    SharedPointer<BounderInterface> _impl;
  public:
    BounderHandle(BounderInterface const& bounder) : _impl(bounder.clone()) { }
    BounderHandle(BounderHandle const& other) : _impl(other._impl) { }
    BounderHandle& operator=(BounderHandle const& other) { _impl = other._impl; return *this; }

    virtual Void write(OutputStream& os) const { _impl->write(os); }
    virtual BounderHandle* clone() const { return new BounderHandle(*this); }

    virtual Pair<PositiveFloatDPValue,UpperBoxType> compute(ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPApproximation hsug) const {
        return _impl->compute(f,dom,hsug);
    }
  public:
    virtual ~BounderHandle() = default;
};



} // namespace Ariadne

#endif /* ARIADNE_BOUNDER_HPP */
