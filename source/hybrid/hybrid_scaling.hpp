/***************************************************************************
 *            hybrid_scaling.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file hybrid_scaling.hpp
 *  \brief Scalings for variables in hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_SCALING_HPP
#define ARIADNE_HYBRID_SCALING_HPP

#include <iostream>
#include <map>

#include "../utility/container.hpp"
#include "../utility/stlio.hpp"

#include "../symbolic/variables.hpp"

namespace Ariadne {

typedef void Void;
typedef std::ostream OutputStream;


//! \ingroup HybridModule
//! \brief A class which defines the state space grid to use in location \a loc given the continuous state variables \a spc.
class HybridScalingInterface
{
  public:
    //!
    virtual ~HybridScalingInterface() = default;
    virtual HybridScalingInterface* clone() const = 0;
    virtual FloatDPValue scaling(const DiscreteLocation& loc, const RealVariable& var) const = 0;
    virtual Void write(OutputStream& os) const = 0;
};
inline OutputStream& operator<<(OutputStream& os, const HybridScalingInterface& hsc) { hsc.write(os); return os; }


//! \ingroup HybridModule
//! \brief A class which defines the state space grid to use in location \a loc given the continuous state variables \a spc.
class HybridScaling
{
    std::shared_ptr<HybridScalingInterface> _ptr;
  public:
    HybridScaling(const HybridScalingInterface& ref) : _ptr(ref.clone()) { }
    operator const HybridScalingInterface& () const { return *this->_ptr; }
    operator HybridScalingInterface& () { return *this->_ptr; }
    FloatDPValue scaling(const DiscreteLocation& loc, const RealVariable& var) const { return this->_ptr->scaling(loc,var); }
  public:
    HybridScaling(const Map<Identifier,FloatDPValue>& scalings);
    HybridScaling(const InitializerList<Pair<RealVariable,FloatDP>>& scalings);
};

class SimpleHybridScaling
    : public HybridScalingInterface
{
    Map<Identifier,FloatDPValue> _scalings;
  public:
    SimpleHybridScaling() : _scalings() { }
    SimpleHybridScaling(const Map<Identifier,FloatDPValue>& scalings) : _scalings(scalings) { }
    SimpleHybridScaling(const InitializerList<Pair<RealVariable,FloatDP>>& scalings);
    Void set_scaling(const RealVariable& var, FloatDPValue res) { ARIADNE_ASSERT(decide(res>0)); _scalings[var.name()]=res; }
    virtual SimpleHybridScaling* clone() const { return new SimpleHybridScaling(*this); }
    virtual FloatDPValue scaling(const DiscreteLocation& loc, const RealVariable& var) const {
        return (this->_scalings.has_key(var.name())) ? this->_scalings[var.name()] : FloatDPValue(1.0); }
    virtual Void write(OutputStream& os) const { os << "HybridScaling( " << this->_scalings << " )"; }
};

inline Pair<RealVariable,FloatDP> operator|(const RealVariable& var, FloatDP scal) {
    return Pair<RealVariable,FloatDP>(var,scal); }

inline HybridScaling::HybridScaling(const Map<Identifier,FloatDPValue>& scalings)
    : _ptr(new SimpleHybridScaling(scalings)) { }

inline HybridScaling::HybridScaling(const InitializerList<Pair<RealVariable,FloatDP>>& scalings)
    : _ptr(new SimpleHybridScaling(scalings)) { }

inline SimpleHybridScaling::SimpleHybridScaling(const InitializerList<Pair<RealVariable,FloatDP>>& scalings) {
    for(auto iter=scalings.begin(); iter!=scalings.end(); ++iter) {
        this->_scalings.insert(iter->first.name(),FloatDPValue(iter->second));
    }
}

}

#endif // ARIADNE_HYBRID_SCALING_HPP

