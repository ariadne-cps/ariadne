/***************************************************************************
 *            taylor_calculus.h
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
 
/*! \file taylor_calculus.h
 *  \brief Methods of taylor calculus based on the TaylorModel class.
 */


#ifndef ARIADNE_TAYLOR_CALCULUS_H
#define ARIADNE_TAYLOR_CALCULUS_H

#include "tribool.h"
#include "logging.h"
#include "calculus_interface.h"
#include "calculus_base.h"

/* \brief Top-level namespace. */
namespace Ariadne {
 
using std::pair;


template<class T> class array;

class Interval;
class ExpressionInterface;
class FunctionInterface;
template<class X> class Vector;
class Box;
class TaylorModel;
class TaylorSet;
class TaylorFunction;

/*! \brief Tools for analysing dynamical systems based on function models. */
class TaylorCalculus
    : public CalculusBase<TaylorModel>
{
    typedef Float R;
    typedef Float A;
    typedef Interval I;
  private:
    ushort _spacial_order;
    ushort _temporal_order;
    double _sweep_threshold;
    ushort _order;
    ushort _smoothness;
    double _maximum_step_size;

    shared_ptr<TaylorModel::Accuracy> _spacial_accuracy_ptr;
  public:
    //!
    //!
    typedef TaylorModel::Accuracy AccuracyType;
    typedef TaylorModel VariableType;
    typedef TaylorModel TimeModelType;
    typedef TaylorSet SetModelType;
    typedef TaylorFunction MapModelType;
    typedef TaylorFunction FlowModelType;
    typedef TaylorExpression PredicateModelType;
    typedef Float TimeType;
    typedef Float RealType;
    typedef Interval IntervalType;
    typedef Vector<Interval> BoxType;
    typedef FunctionInterface FunctionType;
    typedef ExpressionInterface ExpressionType;

    typedef SetModelType EnclosureType;
  public:
    using CalculusBase<TaylorModel>::verbosity;
    using CalculusBase<TaylorModel>::active;
    using CalculusBase<TaylorModel>::reset_step;
    using CalculusBase<TaylorModel>::integration_step;
    using CalculusBase<TaylorModel>::reachability_step;

  public:
    //! \brief Default constructor.
    TaylorCalculus();

    //! \brief Test if a set satisfied the constraint given by the guard model. Returns \a true is all 
    //! points in the set satisfy the constraint, \a false if all points do not satisfy the constraint, 
    //! and indeterminate otherwise.
    tribool active(const PredicateModelType& guard_model, 
                   const SetModelType& _set_model) const;

    //! \brief Computes an over-approximation to the time interval for which the \a initial_set_model 
    //! touch the set specified by the \a guard model under the \a flow_model. The \a minimum and \a maximum_time 
    //! gives the minimum and maximum time for which the evolution is valid.
    IntervalType
    touching_time_interval(const PredicateModelType& guard_model, 
                           const FlowModelType& flow_model, 
                           const SetModelType& initial_set_model) const;

    using CalculusBase<TaylorModel>::crossing_time;

    //! \brief Computes the time at which points in the \a initial_set_model cross the zero-set of the
    //! the \a guard_model under evolution of the \a flow_model, for times between the \a minimum_time and \a maximum_time.
    //! The crossing must be (differentiably) transverse.
    TimeModelType crossing_time(const PredicateModelType& guard_model,
                                const FlowModelType& flow_model, 
                                const SetModelType& initial_set_model) const;

  
    //! \brief Computes the image of the set defined by \a set_model under the map \a map. 
    SetModelType reset_step(const FunctionType& map, 
                            const SetModelType& set_model) const;
  
    //! \brief Computes the image of the set defined by \a set_model under the approximation of the map 
    //! given by \a map_model.
    SetModelType reset_step(const MapModelType& map_model, 
                            const SetModelType& set_model) const;
  
    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model. The \a integration_time_model \f$\tau(e)\f$ gives the time the point 
    //! starting at \f$x(e)\f$ should be flowed.
    SetModelType integration_step(const FlowModelType& flow_model, 
                                  const SetModelType& initial_set_model, 
                                  const TimeModelType& integration_time_model) const;
  
    //! \brief Gives the extended time model for the reachability step between the
    //! \a initial_time_model and the \a final_time_model. The new time is given by
    //! \f$\tau'(e,s) = (1-s)\tau_0(e)+s\tau_1(e)\f$.
    TimeModelType reachability_time(const TimeModelType& initial_time_model, 
                                    const TimeModelType& final_time_model) const;

    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times given by \a reachability_time_model. 
    //! The \a reachability_time_model must have one more independent variable than the 
    //! \a initial_set_model.
    //! 
    //! \invariant <code>reachability_time_model.argument_size()==initial_set_model.argument_size()+1</code>
    virtual SetModelType 
    reachability_step(const FlowModelType& flow_model, 
                      const SetModelType& initial_set_model, 
                      const TimeModelType& reachability_time_model) const;
  
    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times between \a initial_time_model and \a final_time_model.
    SetModelType reachability_step(const FlowModelType& flow_model, 
                                   const SetModelType& initial_set_model, 
                                   const TimeType& initial_timel, 
                                   const TimeType& final_time) const;
  
    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times between \a initial_time_model and \a final_time_model.
    SetModelType reachability_step(const FlowModelType& flow_model, 
                                   const SetModelType& initial_set_model, 
                                   const TimeModelType& initial_time_model, 
                                   const TimeModelType& final_time_model) const;
  
    //! \brief Computed a pair \f$(h,B)\f$ such that the flow of the vector_field \a vf starting in
    //! domain \a d remains in \a B for times up to \a h. The maximum allowable \a h and maximum
    //! allowable diameter of \a B are given.
    pair<TimeType,BoxType> 
    flow_bounds(const FunctionType& vf, 
                const BoxType& d, 
                const RealType& maximum_step_size, 
                const RealType& maximum_bound_diameter) const;



    //! \brief A model for the map \a f over the domain \a d.
    MapModelType map_model(const FunctionType& f, const BoxType& d) const;

    //! \brief A model for the flow determined by the vector field \a vf over the initial domain \a d,
    //! valid for times up to \a h, assuming that the state remains in the bounding box \a b.
    FlowModelType flow_model(const FunctionType& vf, const BoxType& d, 
                             const TimeType& h, const BoxType& b) const;

    //! \brief A model for the real-valued function \a g over the domain \a d.
    PredicateModelType predicate_model(const ExpressionType& g, const BoxType& d) const;

    //! \brief A model for the real-valued function \a g over the domain \a d. \deprecated
    PredicateModelType predicate_model(const FunctionType& g, const BoxType& d) const;

    //! \brief A model for the constant time \a t over the box \a d.
    TimeModelType time_model(const Float& t, const BoxType& d) const;


    //@{ \name Set-based operations
    //! \brief Compute a model for the given box \a bx.
    SetModelType set_model(const BoxType& bx) const;
    //! \brief Compute a model for the given enclosure \a e.
    SetModelType set_model(const EnclosureType& e) const;
    //! \brief Compute an enclosure for the set model \a s.
    virtual EnclosureType enclosure(const SetModelType& s) const;
    //! \brief Tests if the set described by the model \a s is disjoint from the box \a box.
    tribool disjoint(const SetModelType& s, const BoxType& bx) const;
    //! \brief A box containing the set \a s.
    BoxType bounding_box(const SetModelType& s) const;
    //! \brief A list of sets obtained by subdividing the set \a s into at least two smaller pieces.
    array<SetModelType> subdivide(const SetModelType& s) const;
    //! \brief An over-approximation to the set \a s with a simplified description.
    SetModelType simplify(const SetModelType& s) const;
    //@}

};


std::pair<Float, Vector<Interval> > flow_bounds(FunctionInterface const&,Vector<Interval> const&,Float const&);

}


#endif /* ARIADNE_TAYLOR_CALCULUS_H */
