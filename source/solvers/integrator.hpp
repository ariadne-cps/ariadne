/***************************************************************************
 *            integrator.hpp
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

/*! \file integrator.hpp
 *  \brief Solver classes for differential equations.
 */

#ifndef ARIADNE_INTEGRATOR_HPP
#define ARIADNE_INTEGRATOR_HPP

#include <exception>
#include <stdexcept>
#include <string>

#include "solvers/integrator_interface.hpp"
#include "function/function_interface.hpp"

#include "utility/declarations.hpp"
#include "utility/attribute.hpp"
#include "utility/logging.hpp"
#include "utility/pointer.hpp"
#include "function/affine.hpp"

namespace Ariadne {

class Real;
template<class X> class Vector;
template<class X> class Differential;
template<class X> class Procedure;
template<class X> class Polynomial;
typedef Differential<ValidatedNumericType> ValidatedDifferential;
typedef Vector< Procedure<ValidatedNumericType> > ValidatedVectorProcedure;
typedef FunctionModelFactoryInterface<ValidatedTag> ValidatedFunctionModelDPFactoryInterface;
typedef std::shared_ptr<const ValidatedFunctionModelDPFactoryInterface> FunctionFactoryPointer;

struct LipschitzConstant : Attribute<double> { LipschitzConstant(double v) : Attribute<double>(v) { } };
struct StepMaximumError : Attribute<double> { StepMaximumError(double v) : Attribute<double>(v) { } };
struct StepSweepThreshold : Attribute<double> { StepSweepThreshold(double v) : Attribute<double>(v) { } };
struct SpacialOrder : Attribute<Nat> { SpacialOrder(Nat v) : Attribute<Nat>(v) { } };
struct TemporalOrder : Attribute<Nat> { TemporalOrder(Nat v) : Attribute<Nat>(v) { } };
struct MinimumSpacialOrder : Attribute<Nat> { MinimumSpacialOrder(Nat v) : Attribute<Nat>(v) { } };
struct MinimumTemporalOrder : Attribute<Nat> { MinimumTemporalOrder(Nat v) : Attribute<Nat>(v) { } };
struct MaximumSpacialOrder : Attribute<Nat> { MaximumSpacialOrder(Nat v) : Attribute<Nat>(v) { } };
struct MaximumTemporalOrder : Attribute<Nat> { MaximumTemporalOrder(Nat v) : Attribute<Nat>(v) { } };

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
    virtual Void set_maximum_error(double e) { assert(e>0.0); this->_maximum_error=e; }
    virtual double maximum_error() const  { return this->_maximum_error; }
    //! \brief The fraction L(f)*h used for a time step.
    //! The convergence of the Picard iteration is approximately Lf*h.
    Void set_lipschitz_tolerance(double lt) { _lipschitz_tolerance = lt; }
    double lipschitz_tolerance() const { return this->_lipschitz_tolerance; }
    //! \brief  Set maximum size used for a single step.
    double maximum_step_size() const { return this->_maximum_step_size; }
    Void set_maximum_step_size(double hmax) { this->_maximum_step_size = hmax; }

    //! \brief The class which constructs functions for representing the flow.
    const ValidatedFunctionModelDPFactoryInterface& function_factory() const;
    //! \brief Set the class which constructs functions for representing the flow.
    Void set_function_factory(const ValidatedFunctionModelDPFactoryInterface& factory);


    virtual Pair<FloatDPValue,UpperBoxType>
    flow_bounds(const ValidatedVectorFunction& vector_field,
                const ExactBoxType& state_domain,
                const RawFloatDP& maximum_time_step) const;

    virtual ValidatedVectorFunctionModelDP
    flow_step(const ValidatedVectorFunction& vector_field,
              const ExactBoxType& state_domain,
              RawFloatDP& suggested_time_step) const;

    virtual ValidatedVectorFunctionModelDP
    flow_to(const ValidatedVectorFunction& vector_field,
         const ExactBoxType& state_domain,
         const Real& time) const;

    //! \brief Solve \f$\der{\phi}(x,t)=f(\phi(x,t))\f$ for \f$t\in[0,T_{\max}]\f$.
    virtual List<ValidatedVectorFunctionModelDP>
    flow(const ValidatedVectorFunction& vector_field,
         const ExactBoxType& state_domain,
         const Real& minimum_time,
         const Real& maximum_time) const;

    //! \brief Solve \f$\der{\phi}(x,t)=f(\phi(x,t))\f$ for \f$t\in[0,T_{\max}]\f$.
    virtual List<ValidatedVectorFunctionModelDP>
    flow(const ValidatedVectorFunction& vector_field,
         const ExactBoxType& state_domain,
         const Real& maximum_time) const;

    virtual ValidatedVectorFunctionModelDP
    flow_step(const ValidatedVectorFunction& vector_field,
              const ExactBoxType& state_domain,
              const FloatDPValue& suggested_time_step,
              const UpperBoxType& bounding_box) const = 0;

  public:
    double _maximum_error;
    double _lipschitz_tolerance;
    double _maximum_step_size;
    FunctionFactoryPointer _function_factory_ptr;
};

//! \brief An integrator which uses a validated Picard iteration on Taylor models.
class TaylorPicardIntegrator
    : public IntegratorBase
{
    double _step_maximum_error;
    double _step_sweep_threshold;
    Nat _maximum_temporal_order;
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
    Nat maximum_temporal_order() const { return this->_maximum_temporal_order; }
    Void set_maximum_temporal_order(Nat m) { this->_maximum_temporal_order=m; }
    //! \brief  Set the sweep threshold of the Taylor model for a single step.
    double step_sweep_threshold() const { return this->_step_sweep_threshold; }
    Void set_step_sweep_threshold(double lt) { _step_sweep_threshold = lt; }
    //! \brief  Set the maximum error of a single step.
    double step_maximum_error() const { return this->_step_maximum_error; }
    Void set_step_maximum_error(double e) { _step_maximum_error = e; }

    virtual TaylorPicardIntegrator* clone() const { return new TaylorPicardIntegrator(*this); }
    virtual Void write(OutputStream& os) const;

    virtual ValidatedVectorFunctionModelDP
    flow_step(const ValidatedVectorFunction& vector_field,
              const ExactBoxType& state_domain,
              const FloatDPValue& time_step,
              const UpperBoxType& bounding_box) const;

    using IntegratorBase::flow_step;
};

//! \brief An integrator which computes the Taylor series of the flow function with remainder term.
class TaylorSeriesIntegrator
    : public IntegratorBase
{
    double _step_maximum_error;
    double _step_sweep_threshold;
    Nat _minimum_spacial_order;
    Nat _minimum_temporal_order;
    Nat _maximum_spacial_order;
    Nat _maximum_temporal_order;
  public:
    //! \brief Constructor.
    TaylorSeriesIntegrator(MaximumError err);

    //! \brief Constructor.
    TaylorSeriesIntegrator(MaximumError err, SweepThreshold swp, LipschitzConstant lip=0.5);

    //! \brief Constructor.
    TaylorSeriesIntegrator(MaximumError err, SweepThreshold gswp, LipschitzConstant lip,
                           StepMaximumError lerr, StepSweepThreshold lswp, MaximumTemporalOrder maxto);

    //! \brief Constructor.
    TaylorSeriesIntegrator(MaximumError err, SweepThreshold gswp, LipschitzConstant lip,
                           StepMaximumError lerr, StepSweepThreshold lswp,
                           MinimumSpacialOrder minso, MinimumTemporalOrder minto,
                           MaximumSpacialOrder maxso, MaximumTemporalOrder maxto);

    //! \brief The order of the method in space.
    Nat minimum_spacial_order() const { return this->_minimum_spacial_order; }
    Void set_minimum_spacial_order(Nat n) { this->_minimum_spacial_order=n; }
    //! \brief The order of the method in space.
    Nat minimum_temporal_order() const { return this->_minimum_temporal_order; }
    Void set_minimum_temporal_order(Nat m) { this->_minimum_temporal_order=m; }
    //! \brief The maximum order of the method in time.
    Nat maximum_spacial_order() const { return this->_maximum_spacial_order; }
    Void set_maximum_spacial_order(Nat n) { this->_maximum_spacial_order=n; }
    //! \brief The maximum order of the method in time.
    Nat maximum_temporal_order() const { return this->_maximum_temporal_order; }
    Void set_maximum_temporal_order(Nat m) { this->_maximum_temporal_order=m; }
    //! \brief  Set the sweep threshold of the Taylor model representing a single step.
    double step_sweep_threshold() const { return this->_step_sweep_threshold; }
    Void set_step_sweep_threshold(double lswp) { _step_sweep_threshold = lswp; }
    //! \brief  Set the sweep threshold of the Taylor model.
    double step_maximum_error() const { return this->_step_maximum_error; }
    Void set_step_maximum_error(double e) { _step_maximum_error = e; }

    virtual TaylorSeriesIntegrator* clone() const { return new TaylorSeriesIntegrator(*this); }
    virtual Void write(OutputStream& os) const;

    virtual Pair<FloatDPValue,UpperBoxType>
    flow_bounds(const ValidatedVectorFunction& vector_field,
                const ExactBoxType& state_domain,
                const RawFloatDP& suggested_time_step) const;

    virtual ValidatedVectorFunctionModelDP
    flow_step(const ValidatedVectorFunction& vector_field,
              const ExactBoxType& state_domain,
              const FloatDPValue& time_step,
              const UpperBoxType& bounding_box) const;

    using IntegratorBase::flow_step;
};


//! \brief An integrator computes a approximation to the flow which is affine in space.
//! \internal This code is written to allow higher-spacial order approximations.
class AffineIntegrator
    : public IntegratorBase
{
    Nat _spacial_order;
    Nat _temporal_order;
  public:
    AffineIntegrator(MaximumError maximum_error, TemporalOrder temporal_order)
        : IntegratorBase(maximum_error,lipschitz_constant=0.5), _spacial_order(1u), _temporal_order(temporal_order) { }
    AffineIntegrator(MaximumError maximum_error, SpacialOrder spacial_order, TemporalOrder temporal_order)
        : IntegratorBase(maximum_error,lipschitz_constant=0.5), _spacial_order(spacial_order), _temporal_order(temporal_order) { }

    //! \brief The order of the method in space.
    Nat spacial_order() const { return this->_spacial_order; }
    //! \brief The order of the method in time.
    Nat temporal_order() const { return this->_temporal_order; }
    virtual AffineIntegrator* clone() const { return new AffineIntegrator(*this); }
    virtual Void write(OutputStream& os) const;

    virtual ValidatedVectorFunctionModelDP
    flow_step(const ValidatedVectorFunction& vector_field,
              const ExactBoxType& state_domain,
              const FloatDPValue& time_step,
              const UpperBoxType& bounding_box) const;

    using IntegratorBase::flow_step;

    //! \brief Compute the derivative of the flow of f at time zero within \a dom.
    Vector<ValidatedDifferential>
    flow_derivative(const ValidatedVectorFunction& f,
                    const Vector<ValidatedNumericType>& dom) const;
};


} // namespace Ariadne

#endif /* ARIADNE_INTEGRATOR_HPP */
