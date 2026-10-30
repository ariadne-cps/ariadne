// Mathematical operations on bounded symbolic hybrid sets.
#include "config.hpp"
#include "system/hybrid_bounded_constraint_set.hpp"
#include "set/function_set.hpp"
#include "utility/stlio.hpp"

namespace Ariadne {

HybridBoundedConstraintSet::HybridBoundedConstraintSet()
    : _sets()
{
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                                       const RealExpressionBoundedConstraintSet& set)
{
    this->adjoin(loc,set);
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const Set<DiscreteLocation>& locs,
                                                       const RealExpressionBoundedConstraintSet& set)
{
    for (auto loc : locs) { this->adjoin(loc,set); }
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                                       const InitializerList<RealVariableInterval>& bnd)
{
    this->adjoin(loc,RealExpressionBoundedConstraintSet(bnd));
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                                       const InitializerList<RealVariableInterval>& bnd,
                                                       const InitializerList<ContinuousPredicate>& cnstr)
{
    this->adjoin(loc,RealExpressionBoundedConstraintSet(bnd,cnstr));
}

HybridBoundedConstraintSet::HybridBoundedConstraintSet(const DiscreteLocation& loc,
                                                       const RealVariablesBox& bx)
{
    this->adjoin(loc,RealExpressionBoundedConstraintSet(bx));
}

HybridBoundedConstraintSet& HybridBoundedConstraintSet::adjoin(const DiscreteLocation& loc,
                                                               const RealExpressionBoundedConstraintSet& set)
{
    _sets.insert(loc,set); return *this;
}

HybridBoundedConstraintSet* HybridBoundedConstraintSet::clone() const {
    return new HybridBoundedConstraintSet(*this);
}

Set<RealVariable> HybridBoundedConstraintSet::variables(DiscreteLocation loc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return _sets[loc].variables();
}

RealExpressionBoundedConstraintSet const& HybridBoundedConstraintSet::continuous_set(DiscreteLocation loc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return this->_sets[loc];
}

BoundedConstraintSet const HybridBoundedConstraintSet::euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    ARIADNE_ASSERT(this->_sets.has_key(loc));
    return BoundedConstraintSet(this->_sets[loc].euclidean_set(spc));
}

EuclideanSetInterfaceBase* HybridBoundedConstraintSet::_euclidean_set(DiscreteLocation loc, RealSpace spc) const {
    return new BoundedConstraintSet(this->euclidean_set(loc,spc));
}

DiscreteLocation HybridBoundedConstraintSet::location() const {
    ARIADNE_ASSERT(this->_sets.size()==1);
    return this->_sets.begin()->first;
}

Set<DiscreteLocation> HybridBoundedConstraintSet::locations() const {
    return this->_sets.keys();
}

OutputStream& HybridBoundedConstraintSet::_write(OutputStream& os) const {
    return os << "HybridBoundedConstraintSet( "<< this->_sets << " )";
}

} // namespace Ariadne
