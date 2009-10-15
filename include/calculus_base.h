/***************************************************************************
 *            calculus_base.h
 *
 *  Copyright  2008  Pieter Collins
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

/*! \file calculus_base.h
 *  \brief Base class for dynamical calculus routines containing implementations of routings which can be naturally composed from other routines.
 */

#ifndef ARIADNE_CALCULUS_BASE_H
#define ARIADNE_CALCULUS_BASE_H

#include "tribool.h"
#include "logging.h"
#include "function_interface.h"
#include "calculus_interface.h"

#include "numeric.h"
#include "vector.h"

/* \brief Top-level namespace. */
namespace Ariadne {

using std::pair;


template<class T> class array;

class Interval;
class VectorFunction;
template<class X> class Vector;
class Box;



/*! \brief Tools for analysing dynamical systems based on function models. */
template<class Var>
class CalculusBase
    : public CalculusInterface<Var>
    , public Loggable
{
    typedef Float R;
    typedef Float A;
    typedef Interval I;
  public:
    //!
    typedef typename CalculusInterface<Var>::BaseModelType BaseModelType;
    //!
    typedef typename CalculusInterface<Var>::FunctionModelType FunctionModelType;
    //!
    typedef typename CalculusInterface<Var>::SetModelType SetModelType;
    //!
    typedef typename CalculusInterface<Var>::TimeModelType TimeModelType;
    //!
    typedef typename CalculusInterface<Var>::PredicateModelType PredicateModelType;

    typedef FunctionModelType MapModelType;
    typedef FunctionModelType FlowModelType;
    typedef SetModelType FlowSetModelType;

    typedef Float RealType;
    typedef Interval IntervalType;
    typedef Vector<Interval> BoxType;

    typedef Float TimeType;
    typedef VectorFunction VectorFunctionType;
    typedef ScalarFunction ScalarFunctionType;
    typedef SetModelType EnclosureType;
  protected:
    tribool _tribool(const IntervalType& ivl) const {
        if(ivl.lower()>0) { return true; } else if(ivl.upper()<0) { return false; } else { return indeterminate; } }
  public:
    //@{ \name Dynamical operations

    //! \brief Computes the image of the set defined by \a set_model under the approximation of the map
    //! given by \a map_model.
    virtual SetModelType reset_step(const MapModelType& map_model,
                                    const SetModelType& set_model) const = 0;


    //! \brief Test if a set satisfied the constraint given by the guard model. Returns \a true is all
    //! points in the set satisfy the constraint, \a false if all points do not satisfy the constraint,
    //! and indeterminate otherwise.
    virtual tribool active(const PredicateModelType& guard_model,
                           const SetModelType& set_model) const = 0;

    //! \brief Computes an over-approximation to the time interval for which the \a initial_set_model
    //! touch the set specified by the \a guard model under the \a flow_model. The \a minimum and \a maximum_time
    //! gives the minimum and maximum time for which the evolution is valid.
    virtual Interval touching_time_interval(const PredicateModelType& guard_model,
                                            const FlowModelType& flow_model,
                                            const SetModelType& initial_set_model) const = 0;

    //! \brief Computes an over-approximation to the time interval for which the \a initial_set_model
    //! touch the set specified by the \a guard model under the \a flow_model. The \a minimum and \a maximum_time
    //! gives the minimum and maximum time for which the evolution is valid.
    virtual Interval touching_time_interval(const ScalarFunction& guard,
                                            const FlowModelType& flow_model,
                                            const SetModelType& initial_set_model) const
    {
        PredicateModelType guard_model=this->predicate_model(guard,flow_model.range());
        return this->touching_time_interval(guard_model,flow_model,initial_set_model);
    }

    //! \brief Computes an over-approximation to the time interval for which the \a initial_set_model
    //! touch the set specified by the \a guard model under the \a flow_model. The \a minimum and \a maximum_time
    //! gives the minimum and maximum time for which the evolution is valid. Deprecated
    virtual Interval touching_time_interval(const VectorFunction& guard,
                                            const FlowModelType& flow_model,
                                            const SetModelType& initial_set_model) const
    {
        PredicateModelType guard_model=this->map_model(guard,flow_model.range())[0];
        return this->touching_time_interval(guard_model,flow_model,initial_set_model);
    }

    //! \brief Computes an over-approximation to the touching time interval scaling the flow step to [-1,+1]
    virtual Interval scaled_touching_time_interval(const BaseModelType& guard_flow_set_model) const = 0;

    //! \brief Computes an over-approximation to the time interval for which the \a initial_set_model
    //! touch the set specified by the \a guard model under the \a flow_model. The \a minimum and \a maximum_time
    //! gives the minimum and maximum time for which the evolution is valid. Deprecated
    virtual Interval scaled_touching_time_interval(const ScalarFunction& guard,
                                                   const FlowSetModelType& flow_set_model) const
    {
        BaseModelType guard_flow_set_model=apply(guard,flow_set_model);
        return this->scaled_touching_time_interval(guard_flow_set_model);
    }

    //! \brief Computes an over-approximation to the time interval for which the \a initial_set_model
    //! touch the set specified by the \a guard model under the \a flow_model. The \a minimum and \a maximum_time
    //! gives the minimum and maximum time for which the evolution is valid. Deprecated
    virtual Interval scaled_touching_time_interval(const VectorFunction& guard,
                                                     const FlowSetModelType& flow_set_model) const
    {
        BaseModelType guard_flow_set_model=apply(guard,flow_set_model)[0];
        return this->scaled_touching_time_interval(guard_flow_set_model);
    }


    //! \brief Computes the time at which points in the \a initial_set_model cross the zero-set of the
    //! the \a guard under evolution of the \a flow_model.
    //! \brief Computes the time at which points in the \a initial_set_model cross the zero-set of the
    //! the \a guard_model under evolution of the \a flow_model.
    //! The crossing must be (differentiably) transverse.
    virtual TimeModelType crossing_time(const PredicateModelType& guard_model,
                                        const FlowModelType& flow_model,
                                        const SetModelType& initial_set_model) const = 0;

    //! \brief Computes the time at which points in the \a initial_set_model cross the zero-set of the
    //! the \a guard under evolution of the \a flow_model.
    //! The crossing must be (differentiably) transverse.
    virtual TimeModelType crossing_time(const ScalarFunction& guard,
                                        const FlowModelType& flow_model,
                                        const SetModelType& initial_set_model) const
    {
        PredicateModelType guard_model=this->predicate_model(guard,flow_model.range());
        return this->crossing_time(guard_model,flow_model,initial_set_model);
    }

    //! \brief Computes the time at which points in the \a initial_set_model cross the zero-set of the
    //! the \a guard under evolution of the \a flow_model.
    //! The crossing must be (differentiably) transverse.
    virtual TimeModelType crossing_time(const VectorFunction& guard,
                                        const FlowModelType& flow_model,
                                        const SetModelType& initial_set_model) const
    {
        PredicateModelType guard_model=this->map_model(guard,flow_model.range())[0];
        return this->crossing_time(guard_model,flow_model,initial_set_model);
    }

    //! The crossing must be (differentiably) transverse.
    virtual TimeModelType scaled_crossing_time(const BaseModelType& guard_flow_set_model) const = 0;

    //! \brief Computes the time at which points in the \a initial_set_model cross the zero-set of the
    //! the \a guard under evolution of the \a flow_model.
    //! The crossing must be (differentiably) transverse.
    virtual TimeModelType scaled_crossing_time(const ScalarFunction& guard,
                                               const FlowSetModelType& flow_set_model) const
    {
        return this->scaled_crossing_time(apply(guard,flow_set_model));
    }

    //! \brief Computes the time at which points in the \a initial_set_model cross the zero-set of the
    //! the \a guard under evolution of the \a flow_model.
    //! The crossing must be (differentiably) transverse.
    virtual TimeModelType scaled_crossing_time(const VectorFunction& guard,
                                               const FlowSetModelType& flow_set_model) const
    {
        return this->scaled_crossing_time(apply(guard,flow_set_model)[0]);
    }

    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model. The \a integration_time_model \f$\tau(e)\f$ gives the time the point
    //! starting at \f$x(e)\f$ should be flowed.
    virtual SetModelType integration_step(const FlowModelType& flow_model,
                                          const SetModelType& initial_set_model,
                                          const TimeModelType& integration_time_model) const = 0;

    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times between \a initial_time_model and \a final_time_model.
    virtual SetModelType reachability_step(const FlowModelType& flow_model,
                                           const SetModelType& initial_set_model,
                                           const TimeModelType& initial_time_model,
                                           const TimeModelType& final_time_model) const = 0;

    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times given by \a reachability_time_model.
    //! The \a reachability_time_model must have one more independent variable than the
    //! \a initial_set_model.
    //!
    //! \invariant <code>reachability_time_model.argument_size()==initial_set_model.argument_size()+1</code>
    virtual SetModelType
    reachability_step(const FlowModelType& flow_model,
                      const SetModelType& initial_set_model,
                      const TimeModelType& reachability_time_model) const = 0;

    //! \brief Computed a pair \f$(h,B)\f$ such that the flow of the vector_field \a vf starting in
    //! domain \a d remains in \a B for times up to \a h. The maximum allowable \a h and maximum
    //! allowable diameter of \a B are given.
    virtual pair<TimeType,BoxType>
    flow_bounds(const VectorFunctionType& vf,
                const BoxType& d,
                const RealType& maximum_step_size,
                const RealType& maximum_bound_diameter) const = 0;

    //! \brief Computed a pair \f$(h,B)\f$ such that the flow of the vector_field \a vf starting in
    //! domain \a d remains in \a B for times up to \a h. The maximum allowable \a h is given.
    virtual pair<TimeType,BoxType>
    flow_bounds(const VectorFunctionType& vf,
                const BoxType& d,
                const RealType& maximum_step_size) const
    {
        RealType maximum_bound_diameter=std::numeric_limits<double>::max();
        return this->flow_bounds(vf,d,maximum_step_size,maximum_bound_diameter);
    }

    //@}

    //@{ \name Constructing models for functions
    //! \brief A model for the map \a f over the domain \a d.
    virtual MapModelType map_model(const VectorFunctionType& f, const BoxType& d) const = 0;

    //! \brief A model for the flow determined by the vector field \a vf over the initial domain \a d,
    //! valid for times up to \a h, assuming that the state remains in the bounding box \a b.
    virtual FlowModelType flow_model(const VectorFunctionType& vf, const BoxType& d,
                                     const TimeType& h, const BoxType& b) const = 0;

    //! \brief A model for the real-valued function \a g over the domain \a d. \deprecated
    virtual PredicateModelType predicate_model(const VectorFunctionType& g, const BoxType& d) const = 0;

    //! \brief A model for the real-valued function \a g over the domain \a d.
    virtual PredicateModelType predicate_model(const ScalarFunctionType& g, const BoxType& d) const = 0;

    //! \brief A model for the constant time \a t over the box \a d.
    virtual TimeModelType time_model(const Float& t, const BoxType& d) const = 0;
    //@}

    //@{ \name Set-based operations
    //! \brief Compute a model for the given box \a bx.
    virtual SetModelType set_model(const BoxType& bx) const = 0;
    //! \brief Compute an enclosure for the set model \a s.
    virtual EnclosureType enclosure(const SetModelType& s) const = 0;
    //! \brief Tests if the set described by the model \a s is disjoint from the box \a box.
    virtual tribool disjoint(const SetModelType& s, const BoxType& bx) const = 0;
    //! \brief A box containing the set \a s.
    virtual BoxType bounding_box(const SetModelType& s) const = 0;
    //! \brief A list of sets obtained by subdividing the set \a s into at least two smaller pieces.
    virtual array<SetModelType> subdivide(const SetModelType& s) const = 0;
    //! \brief An over-approximation to the set \a s with a simplified description.
    virtual SetModelType simplify(const SetModelType& s) const = 0;
    //@}

  public:
    //! \brief Test if a box satisfies the constraint given by the guard. Returns \a true is all points
    //! in the box satisfy the constraint, \a false if all points do not satisfy the constraint, and
    //! indeterminate otherwise.
    tribool active(const ScalarFunctionType& guard,  const BoxType& box) const {
        return this->_tribool(guard.evaluate(box)); }
    tribool active(const VectorFunctionType& guard,  const BoxType& box) const {
        Vector<Interval> range=guard.evaluate(box);
        return this->_tribool(range[0]); }

    //! \brief Test if a set satisfied the constraint given by the guard. Returns \a true is all points
    //! in the set satisfy the constraint, \a false if all points do not satisfy the constraint, and
    //! indeterminate otherwise.
    tribool active(const ScalarFunctionType& guard,  const SetModelType& set_model) const {
        return this->active(this->predicate_model(guard,set_model.bounding_box()),set_model); }
    tribool active(const VectorFunctionType& guard,  const SetModelType& set_model) const {
        TimeModelType guard_set_model = apply(guard,set_model)[0];
        Interval guard_range=guard_set_model.range();
        tribool guard_active=guard_range.lower()>0 ? tribool(true) : guard_range.upper()<0 ? tribool(false) : indeterminate;
        return guard_active;
    }

    //! \brief Computes the image of the set defined by \a set_model under the \a map.
    SetModelType reset_step(const VectorFunctionType& map,
                            const SetModelType& set_model) const {
        return this->reset_step(this->map_model(map,set_model.range()),set_model); }

    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model. The \a integration_time gives the time all points should be flowed.
    SetModelType integration_step(const FlowModelType& flow_model,
                                  const SetModelType& initial_set_model,
                                  const TimeType& integration_time) const {
        return this->integration_step(flow_model,initial_set_model,this->time_model(integration_time,initial_set_model.domain())); }

    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times between \a initial_time and \a final_time.
    SetModelType reachability_step(const FlowModelType& flow_model,
                                   const SetModelType& initial_set_model,
                                   const TimeType& initial_time,
                                   const TimeType& final_time) const {
        return this->reachability_step(flow_model,initial_set_model,this->time_model(initial_time,initial_set_model.domain()),this->time_model(final_time,initial_set_model.domain())); };

    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times between \a initial_time and \a final_time_model.
    SetModelType reachability_step(const FlowModelType& flow_model,
                                   const SetModelType& initial_set_model,
                                   const TimeType& initial_time,
                                   const TimeModelType& final_time_model) const {
        return this->reachability_step(flow_model,initial_set_model,this->time_model(initial_time,initial_set_model.domain()),final_time_model); }

    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times between \a initial_time_model and \a final_time.
    SetModelType reachability_step(const FlowModelType& flow_model,
                                   const SetModelType& initial_set_model,
                                   const TimeModelType& initial_time_model,
                                   const TimeType& final_time) const {
        return this->reachability_step(flow_model,initial_set_model,initial_time_model,this->time_model(final_time,initial_set_model.domain())); }

    //! \brief Gives the extended time model for the reachability step between the
    //! \a initial_time and the \a final_time_model. The new time is given by
    //! \f$\tau'(e,s) = (1-s)\tau_0+s\tau_1(e)\f$.
    TimeModelType reachability_time(const TimeType& initial_time,
                                    const TimeModelType& final_time_model) const {
        return static_cast<const CalculusInterface<Var>*>(this)->reachability_time(this->time_model(initial_time,final_time_model.domain()),final_time_model); }

    //! \brief Gives the extended time model for the reachability step between the
    //! \a initial_time_model and the \a final_time. The new time is given by
    //! \f$\tau'(e,s) = (1-s)\tau_0(e)+s\tau_1\f$.
    TimeModelType reachability_time(const TimeModelType& initial_time_model,
                                    const TimeType& final_time) const {
        return static_cast<const CalculusInterface<Var>*>(this)->reachability_time(initial_time_model,this->time_model(final_time,initial_time_model.domain())); };

};



}




#endif /* ARIADNE_DYNAMICAL_CALCULUS_H */
