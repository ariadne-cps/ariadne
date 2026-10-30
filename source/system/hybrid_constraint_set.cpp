/***************************************************************************
 *            system/hybrid_constraint_set.cpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

#include "function/functional.hpp"
#include "config.hpp"
#include "system/hybrid_constraint_set.hpp"
#include "utility/stlio.hpp"

namespace Ariadne {

HybridConstraintSet::HybridConstraintSet()
    : _sets()
{
}

HybridConstraintSet::HybridConstraintSet(const DiscreteLocation& loc,
                                         const List<ContinuousPredicate>& cnstr)
{
    this->adjoin(loc,RealExpressionConstraintSet(cnstr));
}

HybridConstraintSet::HybridConstraintSet(const DiscreteLocation& loc,
                                         const RealExpressionConstraintSet& cs)
{
    this->adjoin(loc,cs);
}

HybridConstraintSet::HybridConstraintSet(const Set<DiscreteLocation>& locs,
                                         const RealExpressionConstraintSet& cs)
{
    for(auto loc : locs) { this->adjoin(loc,cs); }
}

HybridConstraintSet* HybridConstraintSet::clone() const {
    return new HybridConstraintSet(*this);
}

HybridConstraintSet& HybridConstraintSet::adjoin(const DiscreteLocation& loc, RealExpressionConstraintSet const& cs) {
    this->_sets.insert(loc,cs); return *this;
}

HybridConstraintSet& HybridConstraintSet::adjoin(const DiscreteLocation& loc, List<ContinuousPredicate> const& cnstr) {
    return this->adjoin(loc,RealExpressionConstraintSet(cnstr));
}

Set<DiscreteLocation> HybridConstraintSet::locations() const {
    return this->_sets.keys();
}

Bool HybridConstraintSet::has_location(DiscreteLocation loc) const {
    return this->_sets.has_key(loc);
}

Set<RealVariable> HybridConstraintSet::variables(DiscreteLocation loc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return _sets[loc].variables();
}

RealExpressionConstraintSet const& HybridConstraintSet::continuous_set(DiscreteLocation loc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return this->_sets[loc];
}

ConstraintSet const HybridConstraintSet::euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return ConstraintSet(this->_sets[loc].euclidean_set(spc));
}

EuclideanSetInterfaceBase* HybridConstraintSet::_euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    return new ConstraintSet(this->euclidean_set(loc,spc));
}

OutputStream& HybridConstraintSet::_write(OutputStream& os) const {
    return os << "HybridConstraintSet( "<< this->_sets << " )";
}



} // namespace Ariadne
