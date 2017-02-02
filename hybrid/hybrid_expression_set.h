/***************************************************************************
 *            hybrid_set.h
 *
 *  Copyright 2008  Pieter Collins
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

/*! \file hybrid_set.h
 *  \brief Sets in hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_EXPRESSION_SET_H
#define ARIADNE_HYBRID_EXPRESSION_SET_H

#include <map>

#include <boost/iterator.hpp>
#include <boost/iterator_adaptors.hpp>

#include <memory>

#include "utility/macros.h"
#include "utility/stlio.h"
#include "utility/declarations.h"
#include "utility/container.h"
#include "geometry/function_set.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/curve.h"

#include "expression/expression_set.h"

#include "hybrid/hybrid_set.decl.h"
#include "hybrid/hybrid_set_interface.h"
#include "hybrid/hybrid_space.h"
#include "hybrid/hybrid_grid.h"
#include "geometry/point.h"
#include "geometry/box.h"

#ifdef ARIADNE_ENABLE_SERIALIZATION
#include "output/serialization.h"
#endif /* ARIADNE_ENABLE_SERIALIZATION */

#include "hybrid/hybrid_graphics_interface.h"

namespace Ariadne {


//! \ingroup ExpressionSetSubModule
//! \ingroup HybridSetSubModule
//! \brief A hybrid set defined by a box in a single location.
//! \details Does not assume a canonical order of the real variables.
template<class IVL> class HybridVariablesBox
    : public Pair<DiscreteLocation,VariablesBox<IVL>>
{
    typedef Pair<DiscreteLocation,VariablesBox<IVL>> Base;
  public:
    HybridVariablesBox(const DiscreteLocation& loc, const VariablesBox<IVL>& bx)
        : Base(loc,bx) { }
    HybridVariablesBox(const DiscreteLocation& loc, const RealSpace& spc, const Box<IVL>& bx)
        : Base(loc,VariablesBox<IVL>(spc,bx)) { }

    //! \brief The location in which the box is defined.
    DiscreteLocation location() const { return this->Base::first; }
    //! \brief The location in which the box is defined.
    VariablesBox<IVL> const& continuous_set() const { return this->Base::second; }
    //! \brief The active variables in the location \a loc.
    virtual Set<RealVariable> variables() const { return this->Base::second.variables(); }
    //! \brief The subset of \f$\mathbb{R}^n\f$ obtained by ordering the variables as defined by \a spc.
    Box<IVL> euclidean_set(const RealSpace& spc) const { return this->second.euclidean_set(spc); }
};

class RealHybridVariablesBox
    : public virtual HybridSetInterface
    , public virtual HybridDrawableInterface
    , public HybridVariablesBox<RealInterval>
{
    using HybridVariablesBox<RealInterval>::variables;
    using HybridVariablesBox<RealInterval>::euclidean_set;

    virtual RealHybridVariablesBox* clone() const;
    virtual SetInterface* _euclidean_set(DiscreteLocation, RealSpace) const;
    virtual Set<DiscreteLocation> locations() const;
    virtual Set<RealVariable> variables(DiscreteLocation) const;
    virtual RealSpace space(DiscreteLocation) const;

    virtual ValidatedSierpinskian overlaps(const HybridExactBox& bx) const override;
    virtual ValidatedSierpinskian separated(const HybridExactBox& bx) const override;
    virtual ValidatedSierpinskian covers(const HybridExactBox& bx) const override;
    virtual ValidatedSierpinskian inside(const HybridExactBox& bx) const;
//    virtual HybridUpperBox bounding_box() const;

    virtual ValidatedSierpinskian inside(const HybridExactBoxes& bx) const override;
    virtual HybridUpperBoxes bounding_box() const override;

    virtual OutputStream& write(OutputStream& os) const override;
    virtual Void draw(CanvasInterface&, const Set<DiscreteLocation>&, const Variables2d&) const override;
};


//! \ingroup ExpressionSetSubModule
//! \ingroup HybridSetSubModule
//! \brief A hybrid set defined by a constraint system in each location.
class HybridConstraintSet
    : public virtual HybridRegularSetInterface
{
    Map<DiscreteLocation, RealExpressionConstraintSet> _sets;
  public:
    HybridConstraintSet();
    //! \brief Construct a set in a single \a location with a list of \a bounds on the variables and nonlinear \a constraints.
    HybridConstraintSet(const DiscreteLocation& location,
                        const List<ContinuousPredicate>& constraints);

    virtual HybridConstraintSet* clone() const override;

    //! \brief The active variables in the location \a loc.
    virtual Set<RealVariable> variables(DiscreteLocation loc) const override;
    //! \brief The subset of \f$\mathbb{R}^V\f$ obtained by restricting to location \a loc.
    RealExpressionConstraintSet const& continuous_set(DiscreteLocation loc) const;
    //! \brief The subset of \f$\mathbb{R}^n\f$ obtained by restricting to location \a loc and ordering the variables as defined by \a spc.
    ConstraintSet const euclidean_set(DiscreteLocation loc, RealSpace spc) const;

    virtual ValidatedSierpinskian overlaps(const HybridExactBox& bx) const override;
    virtual ValidatedSierpinskian separated(const HybridExactBox& bx) const override;
    virtual ValidatedSierpinskian covers(const HybridExactBox& bx) const override;

    virtual OutputStream& write(OutputStream& os) const override;
  protected:
    virtual RegularSetInterface* _euclidean_set(DiscreteLocation loc, RealSpace spc) const override;
};

//! \ingroup ExpressionSetSubModule
//! \ingroup HybridSetSubModule
//! \brief A hybrid set defined by the intersection of a box and a constraint system in each location.
class HybridBoundedConstraintSet
    : public virtual HybridSetInterface
    , public virtual HybridDrawableInterface
{
    Map<DiscreteLocation, RealExpressionBoundedConstraintSet> _sets;
  public:
    HybridBoundedConstraintSet();
    HybridBoundedConstraintSet(const DiscreteLocation& loc,
                               const RealVariablesBox& bx);
    HybridBoundedConstraintSet(const DiscreteLocation& loc,
                               const InitializerList<RealVariableInterval>& bnd);
    //! \brief Construct a set in a single \a location with a list of \a bounds on the variables and nonlinear \a constraints.
    HybridBoundedConstraintSet(const DiscreteLocation& location,
                               const InitializerList<RealVariableInterval>& bounds,
                               const InitializerList<ContinuousPredicate>& constraints);

    DiscreteLocation location() const;

    virtual HybridBoundedConstraintSet* clone() const override;

    //! \brief The set of discrete locations in which the set is nontrivial.
    virtual Set<DiscreteLocation> locations() const override;
    //! \brief The active variables in the location \a loc.
    virtual Set<RealVariable> variables(DiscreteLocation loc) const override;
    //! \brief The subset of \f$\mathbb{R}^V\f$ obtained by restricting to location \a loc.
    RealExpressionBoundedConstraintSet const& continuous_set(DiscreteLocation loc) const;
    //! \brief The subset of \f$\mathbb{R}^n\f$ obtained by restricting to location \a loc and ordering the variables as defined by \a spc.
    BoundedConstraintSet const euclidean_set(DiscreteLocation loc, RealSpace spc) const;

    virtual ValidatedSierpinskian overlaps(const HybridExactBox& bx) const override;
    virtual ValidatedSierpinskian inside(const HybridExactBoxes& bx) const override;

    virtual ValidatedSierpinskian separated(const HybridExactBox& bx) const override;
    virtual ValidatedSierpinskian covers(const HybridExactBox& bx) const override;
    virtual HybridUpperBoxes bounding_box() const override;

    virtual OutputStream& write(OutputStream& os) const override;
    virtual Void draw(CanvasInterface&, const Set<DiscreteLocation>&, const Variables2d&) const override;
  protected:
    virtual BoundedConstraintSet* _euclidean_set(DiscreteLocation loc, RealSpace spc) const override;
};

} // namespace Ariadne

#endif // ARIADNE_HYBRID_EXPRESSION_SET_H
