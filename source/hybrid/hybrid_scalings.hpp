/***************************************************************************
 *            hybrid/hybrid_scalings.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/*! \file hybrid/hybrid_scalings.hpp
 *  \brief Scalings for variables in hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_SCALINGS_HPP
#define ARIADNE_HYBRID_SCALINGS_HPP

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
class HybridScalingsInterface
{
  public:
    //!
    virtual ~HybridScalingsInterface() = default;
    virtual HybridScalingsInterface* clone() const = 0;
    virtual FloatDPValue scaling(const DiscreteLocation& loc, const RealVariable& var) const = 0;
    virtual Void _write(OutputStream& os) const = 0;
};
inline OutputStream& operator<<(OutputStream& os, const HybridScalingsInterface& hsc) { hsc._write(os); return os; }


//! \ingroup HybridModule
//! \brief A class which defines the state space grid to use in location \a loc given the continuous state variables \a spc.
class HybridScalings
{
    std::shared_ptr<HybridScalingsInterface> _ptr;
  public:
    HybridScalings(const HybridScalingsInterface& ref) : _ptr(ref.clone()) { }
    operator const HybridScalingsInterface& () const { return *this->_ptr; }
    operator HybridScalingsInterface& () { return *this->_ptr; }
    FloatDPValue scaling(const DiscreteLocation& loc, const RealVariable& var) const { return this->_ptr->scaling(loc,var); }
    Grid grid(const DiscreteLocation& loc, const RealSpace& spc) const;
  public:
    HybridScalings(const Map<Identifier,FloatDPValue>& scalings);
    HybridScalings(const InitializerList<Pair<RealVariable,FloatDP>>& scalings);
};

class SimpleHybridScalings
    : public HybridScalingsInterface
{
    Map<Identifier,FloatDPValue> _scalings;
  public:
    SimpleHybridScalings() : _scalings() { }
    SimpleHybridScalings(const Map<Identifier,FloatDPValue>& scalings) : _scalings(scalings) { }
    SimpleHybridScalings(const InitializerList<Pair<RealVariable,FloatDP>>& scalings);
    Void set_scaling(const RealVariable& var, FloatDPValue res) { ARIADNE_ASSERT(decide(res>0)); _scalings[var.name()]=res; }
    virtual SimpleHybridScalings* clone() const { return new SimpleHybridScalings(*this); }
    virtual FloatDPValue scaling(const DiscreteLocation& loc, const RealVariable& var) const {
        return (this->_scalings.has_key(var.name())) ? this->_scalings[var.name()] : FloatDPValue(1.0); }
    virtual Void _write(OutputStream& os) const { os << "HybridScalings( " << this->_scalings << " )"; }
};

inline Pair<RealVariable,FloatDP> operator|(const RealVariable& var, FloatDP scal) {
    return Pair<RealVariable,FloatDP>(var,scal); }

inline HybridScalings::HybridScalings(const Map<Identifier,FloatDPValue>& scalings)
    : _ptr(new SimpleHybridScalings(scalings)) { }

inline HybridScalings::HybridScalings(const InitializerList<Pair<RealVariable,FloatDP>>& scalings)
    : _ptr(new SimpleHybridScalings(scalings)) { }

inline SimpleHybridScalings::SimpleHybridScalings(const InitializerList<Pair<RealVariable,FloatDP>>& scalings) {
    for(auto iter=scalings.begin(); iter!=scalings.end(); ++iter) {
        this->_scalings.insert(iter->first.name(),FloatDPValue(iter->second));
    }
}

}

#endif // ARIADNE_HYBRID_SCALINGS_HPP

