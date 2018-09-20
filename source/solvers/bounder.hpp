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

class BounderFactory;

class BounderInterface {
  public:
    virtual Pair<PositiveFloatDPValue,UpperBoxType> flow_bounds(ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPApproximation hsug) const = 0;
};

class BounderBase : public BounderInterface {
  public:
    Pair<PositiveFloatDPValue,UpperBoxType> flow_bounds(ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPApproximation hsug) const;
  protected:
    virtual UpperBoxType formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const = 0;
  private:
    UpperBoxType _initial(ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPValue h) const;
    UpperBoxType _refinement(UpperBoxType B, ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPValue h) const;
  public:
    virtual ~BounderBase() = default;
};

class EulerBounder : public BounderBase {
  protected:
    virtual UpperBoxType formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const;
  public:
    virtual ~EulerBounder() = default;
};

class HeunBounder : public BounderBase {
  protected:
    virtual UpperBoxType formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const;
  public:
    virtual ~HeunBounder() = default;
};

class RalstonBounder : public BounderBase {
  protected:
    virtual UpperBoxType formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const;
  public:
    virtual ~RalstonBounder() = default;
};

class RungeKutta4Bounder : public BounderBase {
  protected:
    virtual UpperBoxType formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const;
  public:
    virtual ~RungeKutta4Bounder() = default;
};

class BounderHandler : public BounderInterface {
    friend class BounderFactory;
  private:
    SharedPointer<BounderInterface> _impl;
    BounderHandler(SharedPointer<BounderInterface> const& impl) : _impl(impl) { }
  public:
    BounderHandler(BounderHandler const& other) : _impl(other._impl) { }
    BounderHandler& operator=(BounderHandler const& other) { _impl = other._impl; return *this; }

    virtual Pair<PositiveFloatDPValue,UpperBoxType> flow_bounds(ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPApproximation hsug) const {
        return _impl->flow_bounds(f,dom,hsug);
    }
  public:
    virtual ~BounderHandler() = default;
};

enum class Bounder : std::uint8_t { EULER, HEUN, RALSTON, RUNGEKUTTA4 };

inline std::ostream& operator << (std::ostream& os, const Bounder& kind) {
    switch (kind) {
        case Bounder::EULER: os << "EULER"; break;
        case Bounder::HEUN: os << "HEUN"; break;
        case Bounder::RALSTON: os << "RALSTON"; break;
        case Bounder::RUNGEKUTTA4: os << "RUNGEKUTTA4"; break;
        default: ARIADNE_FAIL_MSG("Unhandled bounder method for output streaming\n");
    }
    return os;
}

class BounderFactory {
public:
    static BounderHandler create(Bounder method);
};


} // namespace Ariadne

#endif /* ARIADNE_BOUNDER_HPP */
