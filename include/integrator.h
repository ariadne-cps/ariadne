/***************************************************************************
 *            integrator.h
 *
 *  Copyright  2006-10  Pieter Collins
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

/*! \file integrator.h
 *  \brief Solver classes for differential equations.
 */

#ifndef ARIADNE_INTEGRATOR_H
#define ARIADNE_INTEGRATOR_H

#include <exception>
#include <stdexcept>
#include <string>

#include "integrator_interface.h"
#include "function_interface.h"

#include "attribute.h"
#include "logging.h"
#include "pointer.h"
#include "affine.h"

namespace Ariadne {

class Real;
class Interval;
template<class X> class Vector;
template<class X> class Differential;
template<class X> class Procedure;
template<class X> class Polynomial;
typedef Vector< Procedure<Interval> > IntervalVectorProcedure;
typedef Vector< Differential<Interval> > IntervalDifferentialVector;
template<class X> class FunctionModelFactoryInterface;
typedef FunctionModelFactoryInterface<Interval> IntervalFunctionModelFactoryInterface;
typedef shared_ptr<const IntervalFunctionModelFactoryInterface> FunctionFactoryPointer;

class Sweeper;

struct LipschitzConstant : Attribute<double> { LipschitzConstant(double v) : Attribute<double>(v) { } };
struct StepMaximumError : Attribute<double> { StepMaximumError(double v) : Attribute<double>(v) { } };
struct StepSweepThreshold : Attribute<double> { StepSweepThreshold(double v) : Attribute<double>(v) { } };
struct SpacialOrder : Attribute<uint> { SpacialOrder(uint v) : Attribute<uint>(v) { } };
struct TemporalOrder : Attribute<uint> { TemporalOrder(uint v) : Attribute<uint>(v) { } };
struct MinimumSpacialOrder : Attribute<uint> { MinimumSpacialOrder(uint v) : Attribute<uint>(v) { } };
struct MinimumTemporalOrder : Attribute<uint> { MinimumTemporalOrder(uint v) : Attribute<uint>(v) { } };
struct MaximumSpacialOrder : Attribute<uint> { MaximumSpacialOrder(uint v) : Attribute<uint>(v) { } };
struct MaximumTemporalOrder : Attribute<uint> { MaximumTemporalOrder(uint v) : Attribute<uint>(v) { } };

static const Generator<LipschitzConstant> lipschitz_constant = Generator<LipschitzConstant>();
static const Generator<StepMaximumError> step_maximum_error = Generator<StepMaximumError>();
static const Generator<StepSweepThreshold> step_sweep_threshold = Generator<StepSweepThreshold>();
static const Generator<SpacialOrder> spacial_order = Generator<SpacialOrder>();
static const Generator<TemporalOrder> temporal_order = Generator<TemporalOrder>();
static const Generator<MinimumTemporalOrder> minimum_temporal_order = Generator<MinimumTemporalOrder>();
static const Generator<MinimumSpacialOrder> minimum_spacial_order = Generator<MinimumSpacialOrder>();
static const Generator<MaximumSpacialOrder> maximum_spacial_order = Generator<MaximumSpacialOrder>();
static const Generator<MaximumTemporalOrder> maximum_temporal_order = Generator<MaximumTemporalOrder>();

class IntegratorBase
    : public IntegratorInterface
    , public Loggable
{
  protected:
    //! \brief Construct from an error bound for a single step, a constant describing the maximum Lh allowed, and a sweep threshold for the global evolution.
    IntegratorBase(MaximumError e, SweepThreshold swp, LipschitzConstant l);
    IntegratorBase(MaximumError e, LipschitzConstant l);
  public:
    //! \brief A threshold for the error estimate of the approximation.
    virtual void set_maximum_error(double e) { assert(e>0.0); this->_maximum_error=e; }
    virtual double maximum_error() const  { return this->_maximum_error; }
    //! \brief The fraction L(f)*h used for a time step.
    //! The convergence of the Picard iteration is approximately Lf*h.
    void set_lipschitz_tolerance(double lt) { _lipschitz_tolerance = lt; }
    double lipschitz_tolerance() const { return this->_lipschitz_tolerance; }

    //! \brief The class which constructs functions for representing the flow.
    const IntervalFunctionModelFactoryInterface& function_factory() const;
    //! \brief Set the class which constructs functions for representing the flow.
    void set_function_factory(const IntervalFunctionModelFactoryInterface& factory);


    virtual Pair<Float,IntervalVector>
    flow_bounds(const IntervalVectorFunction& vector_field,
                const IntervalVector& state_domain,
                const Float& suggested_time_step) const;

    virtual IntervalVectorFunctionModel
    flow_step(const IntervalVectorFunction& vector_field,
              const IntervalVector& state_domain,
              Float& suggested_time_step) const;

    virtual IntervalVectorFunctionModel
    flow(const IntervalVectorFunction& vector_field,
         const IntervalVector& state_domain,
         const Real& time) const;

    virtual IntervalVectorFunctionModel
    flow(const IntervalVectorFunction& vector_field,
         const IntervalVector& state_domain,
         const Interval& time_domain) const;

    virtual IntervalVectorFunctionModel
    flow_step(const IntervalVectorFunction& vector_field,
              const IntervalVector& state_domain,
              const Float& suggested_time_step,
              const IntervalVector& bounding_box) const = 0;

  public:
    double _maximum_error;
    double _lipschitz_tolerance;
    FunctionFactoryPointer _function_factory_ptr;
};

//! \brief An integrator which uses a validated Picard iteration on Taylor models.
class TaylorPicardIntegrator
    : public IntegratorBase
{
    double _step_maximum_error;
    double _step_sweep_threshold;
    uint _maximum_temporal_order;
  public:
    //! \brief Default constructor.
    TaylorPicardIntegrator(MaximumError err)
        : IntegratorBase(err,SweepThreshold(err/1024),LipschitzConstant(0.5))
        , _step_maximum_error(err/128), _step_sweep_threshold(err/(1024*128)), _maximum_temporal_order(16) { }

    //! \brief Constructor.
    TaylorPicardIntegrator(MaximumError err, SweepThreshold swp, LipschitzConstant lip,
                           StepMaximumError lerr, StepSweepThreshold lswp, MaximumTemporalOrder maxto)
        : IntegratorBase(err,swp,lip), _step_maximum_error(lerr), _step_sweep_threshold(lswp), _maximum_temporal_order(maxto) { }

    //! \brief The order of the method in time.
    uint maximum_temporal_order() const { return this->_maximum_temporal_order; }
    void set_maximum_temporal_order(uint m) { this->_maximum_temporal_order=m; }
    //! \brief  Set the sweep threshold of the Taylor model for a single step.
    double step_sweep_threshold() const { return this->_step_sweep_threshold; }
    void set_step_sweep_threshold(double lt) { _step_sweep_threshold = lt; }
    //! \brief  Set the maximum error of a single step.
    double step_maximum_error() const { return this->_step_maximum_error; }
    void set_step_maximum_error(double e) { _step_maximum_error = e; }

    virtual TaylorPicardIntegrator* clone() const { return new TaylorPicardIntegrator(*this); }
    virtual void write(std::ostream& os) const;

    virtual IntervalVectorFunctionModel
    flow_step(const IntervalVectorFunction& vector_field,
              const IntervalVector& state_domain,
              const Float& time_step,
              const IntervalVector& bounding_box) const;

    using IntegratorBase::flow_step;
};

//! \brief An integrator which computes the Taylor series of the flow function with remainder term.
class TaylorSeriesIntegrator
    : public IntegratorBase
{
    double _step_maximum_error;
    double _step_sweep_threshold;
    uint _minimum_spacial_order;
    uint _minimum_temporal_order;
    uint _maximum_spacial_order;
    uint _maximum_temporal_order;
  public:
    //! \brief Constructor.
    TaylorSeriesIntegrator(MaximumError err)
        : IntegratorBase(err,SweepThreshold(err/1024),LipschitzConstant(0.5)), _step_maximum_error(err/128), _step_sweep_threshold(err/(128*1024))
        , _minimum_spacial_order(1), _minimum_temporal_order(4), _maximum_spacial_order(4), _maximum_temporal_order(12) { }

    //! \brief Constructor.
    TaylorSeriesIntegrator(MaximumError err, SweepThreshold swp)
        : IntegratorBase(err,swp,LipschitzConstant(0.5)), _step_maximum_error(err/128), _step_sweep_threshold(swp/1024)
        , _minimum_spacial_order(1), _minimum_temporal_order(4), _maximum_spacial_order(4), _maximum_temporal_order(12) { }

    //! \brief Constructor.
    TaylorSeriesIntegrator(MaximumError err, SweepThreshold gswp, LipschitzConstant lip,
                           StepMaximumError lerr, StepSweepThreshold lswp,
                           MinimumSpacialOrder minso, MinimumTemporalOrder minto,
                           MaximumSpacialOrder maxso, MaximumTemporalOrder maxto)
        : IntegratorBase(err,gswp,lip), _step_maximum_error(lerr), _step_sweep_threshold(lswp)
        , _minimum_spacial_order(minso), _minimum_temporal_order(minto), _maximum_spacial_order(maxso), _maximum_temporal_order(maxto) { }

    //! \brief The order of the method in space.
    uint minimum_spacial_order() const { return this->_minimum_spacial_order; }
    void set_minimum_spacial_order(uint n) { this->_minimum_spacial_order=n; }
    //! \brief The order of the method in space.
    uint minimum_temporal_order() const { return this->_minimum_temporal_order; }
    void set_minimum_temporal_order(uint m) { this->_minimum_temporal_order=m; }
    //! \brief The maximum order of the method in time.
    uint maximum_spacial_order() const { return this->_maximum_spacial_order; }
    void set_maximum_spacial_order(uint n) { this->_maximum_spacial_order=n; }
    //! \brief The maximum order of the method in time.
    uint maximum_temporal_order() const { return this->_maximum_temporal_order; }
    void set_maximum_temporal_order(uint m) { this->_maximum_temporal_order=m; }
    //! \brief  Set the sweep threshold of the Taylor model representing a single step.
    double step_sweep_threshold() const { return this->_step_sweep_threshold; }
    void set_step_sweep_threshold(double lswp) { _step_sweep_threshold = lswp; }
    //! \brief  Set the sweep threshold of the Taylor model.
    double step_maximum_error() const { return this->_step_maximum_error; }
    void set_step_maximum_error(double e) { _step_maximum_error = e; }

    virtual TaylorSeriesIntegrator* clone() const { return new TaylorSeriesIntegrator(*this); }
    virtual void write(std::ostream& os) const;

    virtual Pair<Float,IntervalVector>
    flow_bounds(const IntervalVectorFunction& vector_field,
                const IntervalVector& state_domain,
                const Float& suggested_time_step) const;

    virtual IntervalVectorFunctionModel
    flow_step(const IntervalVectorFunction& vector_field,
              const IntervalVector& state_domain,
              const Float& time_step,
              const IntervalVector& bounding_box) const;

    using IntegratorBase::flow_step;
};


//! \brief An integrator computes a approximation to the flow which is affine in space.
//! \internal This code is written to allow higher-spacial order approximations.
class AffineIntegrator
    : public IntegratorBase
{
    uint _spacial_order;
    uint _temporal_order;
  public:
    AffineIntegrator(MaximumError maximum_error, TemporalOrder temporal_order)
        : IntegratorBase(maximum_error,lipschitz_constant=0.5), _spacial_order(1u), _temporal_order(temporal_order) { }
    AffineIntegrator(MaximumError maximum_error, SpacialOrder spacial_order, TemporalOrder temporal_order)
        : IntegratorBase(maximum_error,lipschitz_constant=0.5), _spacial_order(spacial_order), _temporal_order(temporal_order) { }

    //! \brief The order of the method in space.
    uint spacial_order() const { return this->_spacial_order; }
    //! \brief The order of the method in time.
    uint temporal_order() const { return this->_temporal_order; }
    virtual AffineIntegrator* clone() const { return new AffineIntegrator(*this); }
    virtual void write(std::ostream& os) const;

    virtual IntervalVectorFunctionModel
    flow_step(const IntervalVectorFunction& vector_field,
              const IntervalVector& state_domain,
              const Float& time_step,
              const IntervalVector& bounding_box) const;

    using IntegratorBase::flow_step;

    //! \brief Compute the derivative of the flow of f at time zero within \a dom.
    IntervalDifferentialVector
    flow_derivative(const IntervalVectorFunction& f,
                    const IntervalVector& dom) const;
};


} // namespace Ariadne

#endif /* ARIADNE_INTEGRATOR_H */
