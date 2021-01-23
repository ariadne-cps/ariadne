/***************************************************************************
 *            solvers/integrator.hpp
 *
 *  Copyright  2006-20  Pieter Collins
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

/*! \file solvers/integrator.hpp
 *  \brief Solver classes for differential equations.
 */

#ifndef ARIADNE_INTEGRATOR_HPP
#define ARIADNE_INTEGRATOR_HPP

#include <exception>
#include <stdexcept>
#include <string>

#include "../solvers/integrator_interface.hpp"
#include "../solvers/bounder.hpp"
#include "../function/function_interface.hpp"

#include "../utility/declarations.hpp"
#include "../utility/attribute.hpp"
#include "../numeric/dyadic.hpp"
#include "../output/logging.hpp"
#include "../utility/pointer.hpp"
#include "../function/affine.hpp"
#include "../algebra/sweeper.hpp"

#include "../configuration/searchable_configuration.hpp"
#include "../configuration/configurable.hpp"
#include "../configuration/configurable.tpl.hpp"
#include "../configuration/configuration_property.hpp"
#include "../configuration/configuration_property.tpl.hpp"

namespace Ariadne {

class Real;
template<class X> class Vector;
template<class X> class Differential;
template<class X> class Procedure;
template<class I, class X> class Polynomial;
template<class X> using MultivariatePolynomial = Polynomial<MultiIndex,X>;
typedef FunctionModelFactoryInterface<ValidatedTag,DoublePrecision> ValidatedFunctionModelDPFactoryInterface;
typedef SharedPointer<const ValidatedFunctionModelDPFactoryInterface> ValidatedFunctionModelDPFactoryPointer;
typedef SharedPointer<const ValidatedFunctionModelDPFactoryInterface> FunctionFactoryPointer;
typedef SharedPointer<const BounderInterface> BounderPointer;

struct LipschitzConstant : Attribute<ApproximateDouble> { using Attribute<ApproximateDouble>::Attribute; };
struct StepMaximumError : Attribute<ApproximateDouble> { using Attribute<ApproximateDouble>::Attribute; };
struct StepSweepThreshold : Attribute<ApproximateDouble> { using Attribute<ApproximateDouble>::Attribute; };
struct Order : Attribute<DegreeType> { Order(DegreeType v) : Attribute<DegreeType>(v) { } };
struct SpacialOrder : Attribute<DegreeType> { SpacialOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct TemporalOrder : Attribute<DegreeType> { TemporalOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct MinimumSpacialOrder : Attribute<DegreeType> { MinimumSpacialOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct MinimumTemporalOrder : Attribute<DegreeType> { MinimumTemporalOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct MaximumSpacialOrder : Attribute<DegreeType> { MaximumSpacialOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct MaximumTemporalOrder : Attribute<DegreeType> { MaximumTemporalOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct StartingStepSizeNumRefinements : Attribute<DegreeType> { StartingStepSizeNumRefinements(DegreeType v) : Attribute<DegreeType>(v) { } };

static const Generator<LipschitzConstant> lipschitz_constant = Generator<LipschitzConstant>();
static const Generator<StepMaximumError> step_maximum_error = Generator<StepMaximumError>();
static const Generator<StepSweepThreshold> step_sweep_threshold = Generator<StepSweepThreshold>();
static const Generator<Order> order = Generator<Order>();
static const Generator<SpacialOrder> spacial_order = Generator<SpacialOrder>();
static const Generator<TemporalOrder> temporal_order = Generator<TemporalOrder>();
static const Generator<MinimumTemporalOrder> minimum_temporal_order = Generator<MinimumTemporalOrder>();
static const Generator<MinimumSpacialOrder> minimum_spacial_order = Generator<MinimumSpacialOrder>();
static const Generator<MaximumSpacialOrder> maximum_spacial_order = Generator<MaximumSpacialOrder>();
static const Generator<MaximumTemporalOrder> maximum_temporal_order = Generator<MaximumTemporalOrder>();

//! \brief Class used for storing the result of a flow step.
class FlowStepModelType : public ValidatedVectorMultivariateFunctionModelDP {
  public:
    using ValidatedVectorMultivariateFunctionModelDP::ValidatedVectorMultivariateFunctionModelDP;
    FlowStepModelType(ValidatedVectorMultivariateFunctionModelDP const& fsmt) : ValidatedVectorMultivariateFunctionModelDP(fsmt) { }
    friend OutputStream& operator<<(OutputStream& os, FlowStepModelType const& flwstpm);
};

//! \brief Class used for storing the result of a flow tube.
class FlowModelType : public List<ValidatedVectorMultivariateFunctionModelDP> {
  public:
    using List<ValidatedVectorMultivariateFunctionModelDP>::List;
    friend OutputStream& operator<<(OutputStream& os, FlowModelType const& flwm);
};

class IntegratorBase : public IntegratorInterface, public Configurable<IntegratorBase>
{
  protected:
    //! \brief Construct from an error bound for a single step, a constant describing the maximum Lh allowed, and a sweep threshold for the global evolution.
    IntegratorBase(MaximumError e, Sweeper<FloatDP> sweeper, LipschitzConstant l, StartingStepSizeNumRefinements nr);
    IntegratorBase(MaximumError e, LipschitzConstant l, StartingStepSizeNumRefinements nr);
  public:
    //! \brief A threshold for the error estimate of the approximation.
    virtual Void set_maximum_error(ApproximateDouble e) { assert(cast_exact(e)>0.0_x); this->_maximum_error=cast_exact(e); }
    virtual ExactDouble maximum_error() const  { return this->_maximum_error; }
    //! \brief The fraction L(f)*h used for a time step.
    //! The convergence of the Picard iteration is approximately Lf*h.
    Void set_lipschitz_tolerance(ApproximateDouble lt) { _lipschitz_tolerance = cast_exact(lt); }
    ExactDouble lipschitz_tolerance() const { return this->_lipschitz_tolerance; }

    //! \brief The number of times the starting step size obtained from the Lipschitz step is divided by two
    Void set_starting_step_size_num_refinements(StartingStepSizeNumRefinements nr) { _starting_step_size_num_refinements = nr; }
    DegreeType starting_step_size_num_refinements() const { return _starting_step_size_num_refinements; }

    //! \brief The class which constructs functions for representing the flow.
    const ValidatedFunctionModelDPFactoryInterface& function_factory() const;
    //! \brief Set the class which constructs functions for representing the flow.
    Void set_function_factory(const ValidatedFunctionModelDPFactoryInterface& factory);

    //! \brief The class that computes bounds.
    const BounderInterface& bounder() const;
    //! \brief Set the class that computes bounds.
    Void set_bounder(const BounderInterface& bounder);

    virtual StepSizeType
    starting_time_step_size(const ValidatedVectorMultivariateFunction& vector_field,
                            const BoxDomainType& state_domain,
                            const StepSizeType& evaluation_time_step,
                            const StepSizeType& maximum_time_step) const;

    virtual StepSizeType
    starting_time_step_size(const ValidatedVectorMultivariateFunction& vector_field,
                            const BoxDomainType& state_domain,
                            const BoxDomainType& parameter_domain,
                            const StepSizeType& evaluation_time_step,
                            const StepSizeType& maximum_time_step) const;


    virtual Pair<StepSizeType,UpperBoxType>
    flow_bounds(const ValidatedVectorMultivariateFunction& vector_field,
                const ExactBoxType& state_domain,
                const StepSizeType& starting_time_step) const;


    virtual Pair<StepSizeType,UpperBoxType>
    flow_bounds(const ValidatedVectorMultivariateFunction& vector_field,
                const ExactBoxType& state_domain,
                const ExactBoxType& parameter_domain,
                const StepSizeType& starting_time_step) const;


    virtual Pair<StepSizeType,UpperBoxType>
    flow_bounds(const ValidatedVectorMultivariateFunction& differential_equation,
                const ExactBoxType& state_domain,
                const StepSizeType& starting_time,
                const ExactBoxType& parameter_domain,
                const StepSizeType& starting_time_step) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& evaluation_time_step,
              StepSizeType& suggested_time_step) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const = 0;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const = 0;

  public:
    ExactDouble _maximum_error;
    ExactDouble _lipschitz_tolerance;
    DegreeType _starting_step_size_num_refinements;
    FunctionFactoryPointer _function_factory_ptr;
    BounderPointer _bounder_ptr;
};

template<> class Configuration<IntegratorBase> : public SearchableConfiguration {
  public:
    typedef ExactDouble RealType;
    typedef ApproximateDouble ApproximateRealType;
    typedef RangeConfigurationProperty<DegreeType> DegreeTypeProperty;
    typedef RangeConfigurationProperty<RealType> RealTypeProperty;
    typedef ListConfigurationProperty<ValidatedFunctionModelDPFactoryInterface> FunctionFactoryProperty;
    typedef ListConfigurationProperty<BounderInterface> BounderProperty;

    Configuration();

    //! \brief The number of times the starting step size obtained from the Lipschitz step is divided by two
    DegreeType const& starting_step_size_num_refinements() const { return dynamic_cast<DegreeTypeProperty const&>(*properties().get("starting_step_size_num_refinements")).get(); }
    void set_starting_step_size_num_refinements(DegreeType const& value) { dynamic_cast<DegreeTypeProperty&>(*properties().get("starting_step_size_num_refinements")).set(value); }
    void set_starting_step_size_num_refinements(DegreeType const& lower, DegreeType const& upper) { dynamic_cast<DegreeTypeProperty&>(*properties().get("starting_step_size_num_refinements")).set(lower,upper); }

    //! \brief The fraction L(f)*h used for a time step.
    //! \details The convergence of the Picard iteration is approximately Lf*h.
    RealType const& lipschitz_tolerance() const { return dynamic_cast<RealTypeProperty const&>(*properties().get("lipschitz_tolerance")).get(); }
    void set_lipschitz_tolerance(ApproximateRealType const& value) { dynamic_cast<RealTypeProperty&>(*properties().get("lipschitz_tolerance")).set(cast_exact(value)); }
    void set_lipschitz_tolerance(ApproximateRealType const& lower, ApproximateRealType const& upper) { dynamic_cast<RealTypeProperty&>(*properties().get("lipschitz_tolerance")).set(cast_exact(lower),cast_exact(upper)); }

    //! \brief A threshold for the error estimate of the approximation.
    RealType const& maximum_error() const { return dynamic_cast<RealTypeProperty const&>(*properties().get("maximum_error")).get(); }
    void set_maximum_error(ApproximateRealType const& value) { dynamic_cast<RealTypeProperty&>(*properties().get("maximum_error")).set(cast_exact(value)); }
    void set_maximum_error(ApproximateRealType const& lower, ApproximateRealType const& upper) { dynamic_cast<RealTypeProperty&>(*properties().get("maximum_error")).set(cast_exact(lower),cast_exact(upper)); }

    //! \brief The function factory to be used.
    ValidatedFunctionModelDPFactoryInterface const& function_factory() const { return dynamic_cast<FunctionFactoryProperty const&>(*properties().get("function_factory")).get(); }
    void set_function_factory(ValidatedFunctionModelDPFactoryInterface const& function_factory) { dynamic_cast<FunctionFactoryProperty&>(*properties().get("function_factory")).set(function_factory); }
    void set_function_factory(SharedPointer<ValidatedFunctionModelDPFactoryInterface> const& function_factory) { dynamic_cast<FunctionFactoryProperty&>(*properties().get("function_factory")).set(function_factory); }
    void set_function_factory(List<SharedPointer<ValidatedFunctionModelDPFactoryInterface>> const& function_factories) { dynamic_cast<FunctionFactoryProperty&>(*properties().get("function_factory")).set(function_factories); }

    //! \brief The bounder to be used.
    BounderInterface const& bounder() const { return dynamic_cast<BounderProperty const&>(*properties().get("bounder")).get(); }
    void set_bounder(BounderInterface const& bounder) { dynamic_cast<BounderProperty&>(*properties().get("bounder")).set(bounder); }
    void set_bounder(SharedPointer<BounderInterface> const& bounder) { dynamic_cast<BounderProperty&>(*properties().get("bounder")).set(bounder); }
};

//! \brief An integrator which uses a validated Picard iteration on Taylor models.
class TaylorPicardIntegrator
    : public IntegratorBase
{
    ExactDouble _step_maximum_error;
    Sweeper<FloatDP> _sweeper;
    DegreeType _minimum_temporal_order;
    DegreeType _maximum_temporal_order;
  public:
    //! \brief Default constructor.
    TaylorPicardIntegrator(MaximumError err);

    //! \brief Constructor.
    TaylorPicardIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip);

    //! \brief Constructor.
    TaylorPicardIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip, StartingStepSizeNumRefinements nr);

    //! \brief Constructor.
    TaylorPicardIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip,
                                                   StepMaximumError lerr, MinimumTemporalOrder minto, MaximumTemporalOrder maxto);

    //! \brief Constructor.
    TaylorPicardIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip, StartingStepSizeNumRefinements nr,
                           StepMaximumError lerr, MinimumTemporalOrder minto, MaximumTemporalOrder maxto);

    //! \brief The order of the method in time.
    DegreeType minimum_temporal_order() const { return this->_minimum_temporal_order; }
    Void set_minimum_temporal_order(Nat m) { this->_minimum_temporal_order=m; }
    DegreeType maximum_temporal_order() const { return this->_maximum_temporal_order; }
    Void set_maximum_temporal_order(Nat m) { this->_maximum_temporal_order=m; }
    //! \brief  Set the sweep threshold of the Taylor model.
    Sweeper<FloatDP> const& sweeper() const { return this->_sweeper; }
    Void set_sweeper(Sweeper<FloatDP> const& sweeper) { _sweeper = sweeper; }
    //! \brief  Set the maximum error of a single step.
    ExactDouble step_maximum_error() const { return this->_step_maximum_error; }
    Void set_step_maximum_error(ApproximateDouble e) { _step_maximum_error = cast_exact(e); }

    virtual TaylorPicardIntegrator* clone() const { return new TaylorPicardIntegrator(*this); }
    virtual Void _write(OutputStream& os) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const;

    using IntegratorBase::flow_step;

  private:
    FlowStepModelType
    _flow_step(const ValidatedVectorMultivariateFunction& vector_field_or_differential_equation,
               const ExactBoxType& state_domain,
               const ExactIntervalType& time_domain,
               const ExactBoxType& parameter_domain,
               const UpperBoxType& bounding_box) const;
};



//! \brief An integrator which computes the Taylor series of the flow function with remainder term.
class TaylorSeriesIntegrator
    : public IntegratorBase
{
    Sweeper<FloatDP> _sweeper;
    DegreeType _order;
  public:
    //! \brief Constructor.
    TaylorSeriesIntegrator(MaximumError err, Order order);

    //! \brief Constructor.
    TaylorSeriesIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip,
                           Order order);

    //! \brief The order of the method in space and time.
    DegreeType order() const { return this->_order; }
    Void set_order(DegreeType n) { this->_order=n; }
    //! \brief  Set the sweep threshold of the Taylor model.
    Sweeper<FloatDP> const& sweeper() const { return this->_sweeper; }
    Void set_sweeper(Sweeper<FloatDP> const& sweeper) { _sweeper = sweeper; }

    virtual TaylorSeriesIntegrator* clone() const { return new TaylorSeriesIntegrator(*this); }
    virtual Void _write(OutputStream& os) const;

    virtual Pair<StepSizeType,UpperBoxType>
    flow_bounds(const ValidatedVectorMultivariateFunction& vector_field,
                const ExactBoxType& state_domain,
                const StepSizeType& suggested_time_step) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const;

    using IntegratorBase::flow_step;

  private:

};


//! \brief An integrator which computes the Taylor series of the flow function with remainder term.
class GradedTaylorSeriesIntegrator
    : public IntegratorBase
{
    ExactDouble _step_maximum_error;
    Sweeper<FloatDP> _sweeper;
    DegreeType _minimum_spacial_order;
    DegreeType _minimum_temporal_order;
    DegreeType _maximum_spacial_order;
    DegreeType _maximum_temporal_order;
  public:
    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(MaximumError err);

    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper);

    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip);

    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip,
                           StepMaximumError lerr, MaximumTemporalOrder maxto);

    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(MaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzConstant lip,
                           StepMaximumError lerr,
                           MinimumSpacialOrder minso, MinimumTemporalOrder minto,
                           MaximumSpacialOrder maxso, MaximumTemporalOrder maxto);

    //! \brief The order of the method in space.
    DegreeType minimum_spacial_order() const { return this->_minimum_spacial_order; }
    Void set_minimum_spacial_order(DegreeType n) { this->_minimum_spacial_order=n; }
    //! \brief The order of the method in space.
    DegreeType minimum_temporal_order() const { return this->_minimum_temporal_order; }
    Void set_minimum_temporal_order(DegreeType m) { this->_minimum_temporal_order=m; }
    //! \brief The maximum order of the method in time.
    DegreeType maximum_spacial_order() const { return this->_maximum_spacial_order; }
    Void set_maximum_spacial_order(DegreeType n) { this->_maximum_spacial_order=n; }
    //! \brief The maximum order of the method in time.
    DegreeType maximum_temporal_order() const { return this->_maximum_temporal_order; }
    Void set_maximum_temporal_order(DegreeType m) { this->_maximum_temporal_order=m; }
    //! \brief  Set the sweep threshold of the Taylor model.
    Sweeper<FloatDP> const& sweeper() const { return this->_sweeper; }
    Void set_sweeper(Sweeper<FloatDP> const& sweeper) { _sweeper = sweeper; }
    //! \brief  Set the sweep threshold of the Taylor model.
    ExactDouble step_maximum_error() const { return this->_step_maximum_error; }
    Void set_step_maximum_error(ApproximateDouble e) { _step_maximum_error = cast_exact(e); }

    virtual GradedTaylorSeriesIntegrator* clone() const { return new GradedTaylorSeriesIntegrator(*this); }
    virtual Void _write(OutputStream& os) const;

    virtual Pair<StepSizeType,UpperBoxType>
    flow_bounds(const ValidatedVectorMultivariateFunction& vector_field,
                const ExactBoxType& state_domain,
                const StepSizeType& suggested_time_step) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const;

    using IntegratorBase::flow_step;

  private:

};


//! \brief An integrator computes a approximation to the flow which is affine in space.
//! \internal This code is written to allow higher-spacial order approximations.
class AffineIntegrator
    : public IntegratorBase
{
    DegreeType _spacial_order;
    DegreeType _temporal_order;
  public:
    AffineIntegrator(MaximumError maximum_error, TemporalOrder temporal_order);
    AffineIntegrator(MaximumError maximum_error, SpacialOrder spacial_order, TemporalOrder temporal_order);

    //! \brief The order of the method in space.
    DegreeType spacial_order() const { return this->_spacial_order; }
    //! \brief The order of the method in time.
    DegreeType temporal_order() const { return this->_temporal_order; }
    virtual AffineIntegrator* clone() const { return new AffineIntegrator(*this); }
    virtual Void _write(OutputStream& os) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const;

    using IntegratorBase::flow_step;

    //! \brief Compute the derivative of the flow of f at time zero within \a dom.
    Vector<Differential<FloatDPBounds>>
    flow_derivative(const ValidatedVectorMultivariateFunction& f,
                    const Vector<FloatDPBounds>& dom) const;
};

template<class P> class Sweeper;

ValidatedVectorMultivariateFunctionModelDP
series_flow_step(const ValidatedVectorMultivariateFunction& f,
                 const ExactBoxType& domx,
                 const Interval<StepSizeType>& domt,
                 const ExactBoxType& doma,
                 const UpperBoxType& bndbx,
                 DegreeType deg,
                 Sweeper<FloatDP> swp);


} // namespace Ariadne

#endif /* ARIADNE_INTEGRATOR_HPP */
