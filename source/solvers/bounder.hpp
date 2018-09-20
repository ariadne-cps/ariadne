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

class FlowBoundsMethodHandlerFactory;

class FlowBoundsMethodHandlerInterface {
  public:
    virtual UpperBoxType initial(ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPValue h) const = 0;
    virtual UpperBoxType refinement(UpperBoxType B, ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPValue h) const = 0;
};

class FlowBoundsMethodHandlerBase : public FlowBoundsMethodHandlerInterface {
  public:
    virtual UpperBoxType initial(ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPValue h) const;
    virtual UpperBoxType refinement(UpperBoxType B, ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPValue h) const;
  protected:
    virtual UpperBoxType formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const = 0;
  public:
    virtual ~FlowBoundsMethodHandlerBase() = default;
};

class EulerFlowBoundsHandler : public FlowBoundsMethodHandlerBase {
  protected:
    virtual UpperBoxType formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const;
  public:
    virtual ~EulerFlowBoundsHandler() = default;
};

class HeunFlowBoundsHandler : public FlowBoundsMethodHandlerBase {
  protected:
    virtual UpperBoxType formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const;
  public:
    virtual ~HeunFlowBoundsHandler() = default;
};

class RalstonFlowBoundsHandler : public FlowBoundsMethodHandlerBase {
  protected:
    virtual UpperBoxType formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const;
  public:
    virtual ~RalstonFlowBoundsHandler() = default;
};

class RungeKutta4FlowBoundsHandler : public FlowBoundsMethodHandlerBase {
  protected:
    virtual UpperBoxType formula(BoxDomainType D, BoxDomainType V, ValidatedVectorFunction f, UpperBoxType B, PositiveFloatDPValue h) const;
  public:
    virtual ~RungeKutta4FlowBoundsHandler() = default;
};

class FlowBoundsMethodHandler : public FlowBoundsMethodHandlerInterface {
    friend class FlowBoundsMethodHandlerFactory;
  private:
    SharedPointer<FlowBoundsMethodHandlerInterface> _impl;
    FlowBoundsMethodHandler(SharedPointer<FlowBoundsMethodHandlerInterface> const& impl) : _impl(impl) { }
  public:
    FlowBoundsMethodHandler(FlowBoundsMethodHandler const& other) : _impl(other._impl) { }
    FlowBoundsMethodHandler& operator=(FlowBoundsMethodHandler const& other) { _impl = other._impl; return *this; }

    virtual UpperBoxType initial(ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPValue h) const {
        return _impl->initial(f,dom,h); }
    virtual UpperBoxType refinement(UpperBoxType B, ValidatedVectorFunction f, BoxDomainType dom, PositiveFloatDPValue h) const {
        return _impl->refinement(B,f,dom,h); }
  public:
    virtual ~FlowBoundsMethodHandler() = default;
};

enum class FlowBoundsMethod : std::uint8_t { EULER, HEUN, RALSTON, RUNGEKUTTA4 };

inline std::ostream& operator << (std::ostream& os, const FlowBoundsMethod& kind) {
    switch (kind) {
        case FlowBoundsMethod::EULER: os << "EULER"; break;
        case FlowBoundsMethod::HEUN: os << "HEUN"; break;
        case FlowBoundsMethod::RALSTON: os << "RALSTON"; break;
        case FlowBoundsMethod::RUNGEKUTTA4: os << "RUNGEKUTTA4"; break;
        default: ARIADNE_FAIL_MSG("Unhandled flow bounds method for output streaming\n");
    }
    return os;
}

class FlowBoundsMethodHandlerFactory {
public:
    static FlowBoundsMethodHandler create(FlowBoundsMethod method);
};


} // namespace Ariadne

#endif /* ARIADNE_BOUNDER_HPP */
