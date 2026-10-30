// Bounded symbolic sets in hybrid state spaces.
#ifndef ARIADNE_SYSTEM_HYBRID_BOUNDED_CONSTRAINT_SET_HPP
#define ARIADNE_SYSTEM_HYBRID_BOUNDED_CONSTRAINT_SET_HPP

#include "hybrid/hybrid_set_interface.hpp"
#include "hybrid/hybrid_drawable_interface.hpp"
#include "expression-set/expression_set.hpp"
#include "utility/container.hpp"

namespace Ariadne {

//! \ingroup ExpressionSetSubModule
//! \ingroup HybridSetSubModule
//! \brief A hybrid set defined by the intersection of a box and a constraint system in each location.
class HybridBoundedConstraintSet
    : public virtual EffectiveHybridSetInterface
    , public virtual HybridDrawableInterface
{
    Map<DiscreteLocation, RealExpressionBoundedConstraintSet> _sets;
  public:
    HybridBoundedConstraintSet();
    HybridBoundedConstraintSet(const DiscreteLocation& loc,
                               const RealVariablesBox& bx);
    HybridBoundedConstraintSet(const DiscreteLocation& loc,
                               const RealExpressionBoundedConstraintSet& set);
    HybridBoundedConstraintSet(const DiscreteLocation& loc,
                               const InitializerList<RealVariableInterval>& bnd);
    //! \brief Construct a set in a single \a location with a list of \a bounds on the variables and nonlinear \a constraints.
    HybridBoundedConstraintSet(const DiscreteLocation& location,
                               const InitializerList<RealVariableInterval>& bounds,
                               const InitializerList<ContinuousPredicate>& constraints);
    HybridBoundedConstraintSet(const Set<DiscreteLocation>& locations,
                               const RealExpressionBoundedConstraintSet& constraint_set);

    HybridBoundedConstraintSet& adjoin(const DiscreteLocation& loc, const RealExpressionBoundedConstraintSet& set);

    virtual HybridBoundedConstraintSet* clone() const override;

    //! \brief If the set only has a single location, returns this, otherwise fails. DEPRECATED
    DiscreteLocation location() const;

    //! \brief The set of discrete locations in which the set is nontrivial.
    virtual Set<DiscreteLocation> locations() const override;
    //! \brief The active variables in the location \a loc.
    virtual Set<RealVariable> variables(DiscreteLocation loc) const override;
    //! \brief The subset of \f$\mathbb{R}^V\f$ obtained by restricting to location \a loc.
    RealExpressionBoundedConstraintSet const& continuous_set(DiscreteLocation loc) const;
    //! \brief The subset of \f$\mathbb{R}^n\f$ obtained by restricting to location \a loc and ordering the variables as defined by \a spc.
    BoundedConstraintSet const euclidean_set(DiscreteLocation loc, RealSpace spc) const;

    virtual LowerKleenean overlaps(const HybridExactBox& bx) const override;
    virtual LowerKleenean inside(const HybridExactBoxes& bx) const override;

    virtual LowerKleenean separated(const HybridExactBox& bx) const override;
    virtual LowerKleenean covers(const HybridExactBox& bx) const override;
    virtual HybridUpperBoxes bounding_box() const override;

    virtual OutputStream& _write(OutputStream& os) const override;
    virtual Void draw(CanvasInterface&, const Set<DiscreteLocation>&, const Variables2d&) const override;
  protected:
    virtual EuclideanSetInterfaceBase* _euclidean_set(DiscreteLocation loc, RealSpace spc) const override;
};

} // namespace Ariadne

#endif // ARIADNE_SYSTEM_HYBRID_BOUNDED_CONSTRAINT_SET_HPP
