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

#include "solvers/integrator_interface.hpp"
#include "solvers/bounder.hpp"
#include "function/function_interface.hpp"
#include "function/taylor_function.hpp"

#include "utility/declarations.hpp"
#include "utility/attribute.hpp"
#include "numeric/dyadic.hpp"
#include "conclog/logging.hpp"
#include "utility/pointer.hpp"
#include "function/affine.hpp"
#include "algebra/sweeper.hpp"

#include "function/function_patch.hpp"

#include "pronest/searchable_configuration.hpp"
#include "pronest/configurable.hpp"
#include "pronest/configuration_interface.hpp"
#include "pronest/configurable.tpl.hpp"
#include "pronest/configuration_property.hpp"
#include "pronest/configuration_property.tpl.hpp"

using namespace ConcLog;

namespace Ariadne {

using ProNest::Configurable;
using ProNest::Configuration;

class Real;
template<class X> class Vector;
template<class X> class Differential;
template<class X> class Procedure;
template<class I, class X> class Polynomial;
template<class X> using MultivariatePolynomial = Polynomial<MultiIndex,X>;

template<class P> class FunctionPatchFactoryInterface;
typedef FunctionPatchFactoryInterface<ValidatedTag> ValidatedFunctionPatchFactoryInterface;
typedef SharedPointer<const ValidatedFunctionPatchFactoryInterface> ValidatedFunctionPatchFactoryPointer;
typedef SharedPointer<const ValidatedFunctionPatchFactoryInterface> FunctionFactoryPointer;
typedef SharedPointer<const BounderInterface> BounderPointer;

struct StepMaximumError : Attribute<ApproximateDouble> {
    StepMaximumError(ApproximateDouble x) : Attribute<ApproximateDouble>(x) { }
    StepMaximumError(double x) : Attribute<ApproximateDouble>(x) { }
};

struct StepSweepThreshold : Attribute<ApproximateDouble> { using Attribute<ApproximateDouble>::Attribute; };
struct Order : Attribute<DegreeType> { Order(DegreeType v) : Attribute<DegreeType>(v) { } };
struct SpacialOrder : Attribute<DegreeType> { SpacialOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct TemporalOrder : Attribute<DegreeType> { TemporalOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct MinimumSpacialOrder : Attribute<DegreeType> { MinimumSpacialOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct MinimumTemporalOrder : Attribute<DegreeType> { MinimumTemporalOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct MaximumSpacialOrder : Attribute<DegreeType> { MaximumSpacialOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct MaximumTemporalOrder : Attribute<DegreeType> { MaximumTemporalOrder(DegreeType v) : Attribute<DegreeType>(v) { } };

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
class FlowStepModelType : public ValidatedVectorMultivariateFunctionPatch {
  public:
    using ValidatedVectorMultivariateFunctionPatch::ValidatedVectorMultivariateFunctionPatch;
    FlowStepModelType(ValidatedVectorMultivariateFunctionPatch const& fsmt) : ValidatedVectorMultivariateFunctionPatch(fsmt) { }
    friend OutputStream& operator<<(OutputStream& os, FlowStepModelType const& flwstpm);
};

//! \brief Class used for storing the result of a flow tube.
class FlowModelType : public List<ValidatedVectorMultivariateFunctionPatch> {
  public:
    using List<ValidatedVectorMultivariateFunctionPatch>::List;
    friend OutputStream& operator<<(OutputStream& os, FlowModelType const& flwm);
};

class IntegratorBase
    : public IntegratorInterface, public Configurable<IntegratorBase>
{
  protected:
    //! \brief Construct with a configuration
    IntegratorBase(Configuration<IntegratorBase> const& config);
  public:

    //! \brief The class which constructs functions for representing the flow.
    const ValidatedFunctionPatchFactory& function_factory() const;
    //! \brief Set the class which constructs functions for representing the flow.
    Void set_function_factory(const ValidatedFunctionPatchFactory& factory);

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
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
};

class BoundedIntegratorBase : public IntegratorBase {
  protected:
    BoundedIntegratorBase(Configuration<BoundedIntegratorBase> const& config);
  public:

    Configuration<BoundedIntegratorBase> const& configuration() const;

    virtual Pair<StepSizeType,UpperBoxType>
    flow_bounds(const ValidatedVectorMultivariateFunction& vector_field,
                const ExactBoxType& state_domain,
                const StepSizeType& maximum_time_step) const;

    virtual Pair<StepSizeType,UpperBoxType>
    flow_bounds(const ValidatedVectorMultivariateFunction& vector_field,
                const ExactBoxType& state_domain,
                const ExactBoxType& parameter_domain,
                const StepSizeType& maximum_time_step) const;

    virtual Pair<StepSizeType,UpperBoxType>
    flow_bounds(const ValidatedVectorMultivariateFunction& differential_equation,
                const ExactBoxType& state_domain,
                const StepSizeType& starting_time,
                const ExactBoxType& parameter_domain,
                const StepSizeType& maximum_time_step) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              StepSizeType& suggested_time_step) const override;

    using IntegratorBase::flow_step;
};

//! \ingroup DifferentialEquationSubModule
//! \brief An integrator which uses a validated Picard iteration on Taylor models.
class TaylorPicardIntegrator : public BoundedIntegratorBase
{
  public:

    //! \brief Construct from a configuration
    TaylorPicardIntegrator(Configuration<TaylorPicardIntegrator> const& config);

    Configuration<TaylorPicardIntegrator> const& configuration() const;

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

    using BoundedIntegratorBase::flow_step;

  private:
    FlowStepModelType
    _flow_step(const ValidatedVectorMultivariateFunction& vector_field_or_differential_equation,
               const ExactBoxType& state_domain,
               const ExactIntervalType& time_domain,
               const ExactBoxType& parameter_domain,
               const UpperBoxType& bounding_box) const;
};

/*
//! \ingroup DifferentialEquationSubModule
class GradedTaylorPicardIntegrator
        : public IntegratorBase
{
    ExactDouble _step_maximum_error;
    GradedSweeper<FloatDP> _sweeper;
    ExactDouble _error_refinement_minimum_improvement_percentage;
    DegreeType _order;
public:
    //! \brief Default constructor.
    GradedTaylorPicardIntegrator(StepMaximumError err, Order order);

    //! \brief The order of the method.
    DegreeType order() const { return this->_order; }
    Void set_order(Nat m) { this->_order=m; this->_sweeper = GradedSweeper<FloatDP>(DoublePrecision(),m); }
    //! \brief  Set the maximum error of a single step.
    ExactDouble step_maximum_error() const { return this->_step_maximum_error; }
    Void set_step_maximum_error(ApproximateDouble e) { _step_maximum_error = cast_exact(e); }
    ExactDouble error_refinement_minimum_improvement_percentage() const { return this->_error_refinement_minimum_improvement_percentage; }
    Void set_error_refinement_minimum_improvement_percentage(ApproximateDouble e) { _error_refinement_minimum_improvement_percentage = cast_exact(e); }

    virtual GradedTaylorPicardIntegrator* clone() const { return new GradedTaylorPicardIntegrator(*this); }
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

private:
    FlowStepModelType
    _flow_step(const ValidatedVectorMultivariateFunction& vector_field_or_differential_equation,
               const ExactBoxType& state_domain,
               const ExactIntervalType& time_domain,
               const ExactBoxType& parameter_domain,
               const UpperBoxType& bounding_box) const;
};



//! \ingroup DifferentialEquationSubModule
//! \brief An integrator which computes the Taylor series of the flow function with remainder term.
class TaylorSeriesIntegrator
    : public BoundedIntegratorBase
{
    Sweeper<FloatDP> _sweeper;
    DegreeType _order;
  public:
    //! \brief Constructor.
    TaylorSeriesIntegrator(StepMaximumError err, Order order);
    //! \brief Constructor.
    TaylorSeriesIntegrator(Sweeper<FloatDP> const& sweeper, LipschitzTolerance lip, Order order);

    //! \brief The order of the method in space and time.
    DegreeType order() const { return this->_order; }
    Void set_order(DegreeType n) { this->_order=n; }
    //! \brief  Set the sweep threshold of the Taylor model.
    Sweeper<FloatDP> const& sweeper() const { return this->_sweeper; }
    Void set_sweeper(Sweeper<FloatDP> const& sweeper) { _sweeper = sweeper; }

    virtual TaylorSeriesIntegrator* clone() const { return new TaylorSeriesIntegrator(*this); }
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

    using BoundedIntegratorBase::flow_step;

  private:
};

*/

//! \brief An integrator which computes the Taylor series of the flow function with remainder term.
class GradedTaylorSeriesIntegrator
    : public BoundedIntegratorBase
{
  public:
    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(Configuration<GradedTaylorSeriesIntegrator> const& config);

    Configuration<GradedTaylorSeriesIntegrator> const& configuration() const;

    virtual GradedTaylorSeriesIntegrator* clone() const { return new GradedTaylorSeriesIntegrator(*this); }
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

    using BoundedIntegratorBase::flow_step;
};

//! \brief An integrator computes a approximation to the flow which is affine in space.
//! \internal This code is written to allow higher-spacial order approximations.
class AffineIntegrator
    : public BoundedIntegratorBase
{
  public:
    AffineIntegrator(Configuration<AffineIntegrator> const& config);

    Configuration<AffineIntegrator> const& configuration() const;

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

    using BoundedIntegratorBase::flow_step;

    //! \brief Compute the derivative of the flow of f at time zero within \a dom.
    Vector<Differential<FloatDPBounds>>
    flow_derivative(const ValidatedVectorMultivariateFunction& f,
                    const Vector<FloatDPBounds>& dom) const;
};

template<class P> class Sweeper;

ValidatedVectorMultivariateFunctionPatch
series_flow_step(const ValidatedVectorMultivariateFunction& f,
                 const ExactBoxType& domx,
                 const Interval<StepSizeType>& domt,
                 const ExactBoxType& doma,
                 const UpperBoxType& bndbx,
                 DegreeType deg,
                 Sweeper<FloatDP> swp);


} // namespace Ariadne

namespace ProNest {

using Ariadne::IntegratorBase;
using Ariadne::BoundedIntegratorBase;
using Ariadne::TaylorPicardIntegrator;
using Ariadne::GradedTaylorSeriesIntegrator;
using Ariadne::AffineIntegrator;
using Ariadne::BounderInterface;
using Ariadne::DegreeType;
using Ariadne::EulerBounder;
using Ariadne::FloatDP;
using Ariadne::Sweeper;
using Ariadne::ThresholdSweeper;
using Ariadne::TaylorFunctionFactory;
using Ariadne::DoublePrecision;
using Ariadne::ValidatedFunctionPatchFactoryInterface;
using Ariadne::cast_exact;

template<> struct LinearSearchSpaceConverter<unsigned short> : ConfigurationSearchSpaceConverterInterface<unsigned short> {
    int to_int(unsigned short const& value) const override { return static_cast<int>(value); }
    unsigned short from_int(int i) const override { return static_cast<unsigned short>(i); }
    ConfigurationSearchSpaceConverterInterface* clone() const override { return new LinearSearchSpaceConverter(*this); }
};

template<> struct Configuration<IntegratorBase> : public SearchableConfiguration {
    typedef Configuration<IntegratorBase> C;
    typedef InterfaceListConfigurationProperty<ValidatedFunctionPatchFactoryInterface> FunctionFactoryProperty;

    Configuration() {
        add_property("function_factory",FunctionFactoryProperty(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(DoublePrecision(),1e-7))));
    }

    //! \brief The function factory to be used.
    ValidatedFunctionPatchFactoryInterface const& function_factory() const { return at<FunctionFactoryProperty>("function_factory").get(); }
};

template<> struct Configuration<BoundedIntegratorBase> : public Configuration<IntegratorBase> {
public:
    typedef double RealType;
    typedef RangeConfigurationProperty<DegreeType> DegreeTypeProperty;
    typedef RangeConfigurationProperty<RealType> RealTypeProperty;
    typedef InterfaceListConfigurationProperty<BounderInterface> BounderProperty;

    Configuration() {
        add_property("lipschitz_tolerance",RealTypeProperty(0.5,Log2SearchSpaceConverter<RealType>()));
        add_property("bounder",BounderProperty(EulerBounder()));
    }

    //! \brief The fraction L(f)*h used for a time step.
    //! \details The convergence of the Picard iteration is approximately Lf*h.
    RealType const& lipschitz_tolerance() const { return at<RealTypeProperty>("lipschitz_tolerance").get(); }

    //! \brief The bounder to be used.
    BounderInterface const& bounder() const { return at<BounderProperty>("bounder").get(); }
};

template<> struct Configuration<TaylorPicardIntegrator> : public Configuration<BoundedIntegratorBase> {
    typedef Configuration<TaylorPicardIntegrator> C;
    typedef double RealType;
    typedef RangeConfigurationProperty<DegreeType> DegreeTypeProperty;
    typedef RangeConfigurationProperty<RealType> RealTypeProperty;
    typedef InterfaceListConfigurationProperty<ValidatedFunctionPatchFactoryInterface> FunctionFactoryProperty;
    typedef InterfaceListConfigurationProperty<BounderInterface> BounderProperty;
    typedef HandleListConfigurationProperty<Sweeper<FloatDP>> SweeperProperty;

    Configuration() {
        add_property("sweeper",SweeperProperty(ThresholdSweeper<FloatDP>(DoublePrecision(),1e-7)));
        add_property("step_maximum_error",RealTypeProperty(1e-6,Log10SearchSpaceConverter<RealType>()));
        add_property("minimum_temporal_order",DegreeTypeProperty(0u));
        add_property("maximum_temporal_order",DegreeTypeProperty(12u));
    }

    //! \brief The sweeper to be used when creating a flow function.
    Sweeper<FloatDP> const& sweeper() const { return at<SweeperProperty>("sweeper").get(); }
    C& set_sweeper(Sweeper<FloatDP> const& sweeper) { at<SweeperProperty>("sweeper").set(sweeper); return *this; }
    C& set_sweeper(List<Sweeper<FloatDP>> const& sweepers) { at<SweeperProperty>("sweeper").set(sweepers); return *this; }

    //! \brief The maximum error produced on a single step of integration.
    RealType const& step_maximum_error() const { return at<RealTypeProperty>("step_maximum_error").get(); }
    C& set_step_maximum_error(double const& value) { at<RealTypeProperty>("step_maximum_error").set(value); return *this; }
    C& set_step_maximum_error(double const& lower, double const& upper) { at<RealTypeProperty>("step_maximum_error").set(lower,upper); return *this; }

    //! \brief The minimum temporal order to be used for the flow function.
    DegreeType const& minimum_temporal_order() const { return at<DegreeTypeProperty>("minimum_temporal_order").get(); }
    C& set_minimum_temporal_order(DegreeType const& value) { at<DegreeTypeProperty>("minimum_temporal_order").set(value); return *this; }
    C& set_minimum_temporal_order(DegreeType const& lower, DegreeType const& upper) { at<DegreeTypeProperty>("minimum_temporal_order").set(lower,upper); return *this; }

    //! \brief The maximum temporal order to be used for the flow function.
    DegreeType const& maximum_temporal_order() const { return at<DegreeTypeProperty>("maximum_temporal_order").get(); }
    C& set_maximum_temporal_order(DegreeType const& value) { at<DegreeTypeProperty>("maximum_temporal_order").set(value); return *this; }
    C& set_maximum_temporal_order(DegreeType const& lower, DegreeType const& upper) { at<DegreeTypeProperty>("maximum_temporal_order").set(lower,upper); return *this; }

    //! Base properties

    ValidatedFunctionPatchFactoryInterface const& function_factory() const { return at<FunctionFactoryProperty>("function_factory").get(); }
    C& set_function_factory(ValidatedFunctionPatchFactoryInterface const& factory) { at<FunctionFactoryProperty>("function_factory").set(factory); return *this; }
    C& set_function_factory(SharedPointer<ValidatedFunctionPatchFactoryInterface> const& factory) { at<FunctionFactoryProperty>("function_factory").set(factory); return *this; }

    //! \brief The fraction L(f)*h used for a time step.
    //! \details The convergence of the Picard iteration is approximately Lf*h.
    RealType const& lipschitz_tolerance() const { return at<RealTypeProperty>("lipschitz_tolerance").get(); }
    C& set_lipschitz_tolerance(double const& value) { at<RealTypeProperty>("lipschitz_tolerance").set(value); return *this; }
    C& set_lipschitz_tolerance(double const& lower, double const& upper) { at<RealTypeProperty>("lipschitz_tolerance").set(lower,upper); return *this; }

    //! \brief The bounder to be used.
    BounderInterface const& bounder() const { return at<BounderProperty>("bounder").get(); }
    C& set_bounder(BounderInterface const& bounder) { at<BounderProperty>("bounder").set(bounder); return *this; }
    C& set_bounder(SharedPointer<BounderInterface> const& bounder) { at<BounderProperty>("bounder").set(bounder); return *this; }
};

template<> struct Configuration<GradedTaylorSeriesIntegrator> : public Configuration<BoundedIntegratorBase> {

    typedef Configuration<GradedTaylorSeriesIntegrator> C;
    typedef double RealType;
    typedef RangeConfigurationProperty<DegreeType> DegreeTypeProperty;
    typedef RangeConfigurationProperty<RealType> RealTypeProperty;
    typedef InterfaceListConfigurationProperty<ValidatedFunctionPatchFactoryInterface> FunctionFactoryProperty;
    typedef InterfaceListConfigurationProperty<BounderInterface> BounderProperty;
    typedef HandleListConfigurationProperty<Sweeper<FloatDP>> SweeperProperty;

    Configuration() {
        add_property("sweeper",SweeperProperty(ThresholdSweeper<FloatDP>(DoublePrecision(),1e-7)));
        add_property("step_maximum_error",RealTypeProperty(1e-6,Log10SearchSpaceConverter<RealType>()));
        add_property("minimum_spacial_order",DegreeTypeProperty(1u));
        add_property("maximum_spacial_order",DegreeTypeProperty(4u));
        add_property("minimum_temporal_order",DegreeTypeProperty(4u));
        add_property("maximum_temporal_order",DegreeTypeProperty(12u));
    }

    //! \brief The sweeper to be used when creating a flow function.
    Sweeper<FloatDP> const& sweeper() const { return at<SweeperProperty>("sweeper").get(); }
    C& set_sweeper(Sweeper<FloatDP> const& sweeper) { at<SweeperProperty>("sweeper").set(sweeper); return *this; }
    C& set_sweeper(List<Sweeper<FloatDP>> const& sweepers) { at<SweeperProperty>("sweeper").set(sweepers); return *this; }

    //! \brief The maximum error produced on a single step of integration.
    RealType const& step_maximum_error() const { return at<RealTypeProperty>("step_maximum_error").get(); }
    C& set_step_maximum_error(double const& value) { at<RealTypeProperty>("step_maximum_error").set(value); return *this; }
    C& set_step_maximum_error(double const& lower, double const& upper) { at<RealTypeProperty>("step_maximum_error").set(lower,upper); return *this; }

    //! \brief The minimum spacial order to be used for the flow function.
    DegreeType const& minimum_spacial_order() const { return at<DegreeTypeProperty>("minimum_spacial_order").get(); }
    C& set_minimum_spacial_order(DegreeType const& value) { at<DegreeTypeProperty>("minimum_spacial_order").set(value); return *this; }
    C& set_minimum_spacial_order(DegreeType const& lower, DegreeType const& upper) { at<DegreeTypeProperty>("minimum_spacial_order").set(lower,upper); return *this; }

    //! \brief The maximum spacial order to be used for the flow function.
    DegreeType const& maximum_spacial_order() const { return at<DegreeTypeProperty>("maximum_spacial_order").get(); }
    C& set_maximum_spacial_order(DegreeType const& value) { at<DegreeTypeProperty>("maximum_spacial_order").set(value); return *this; }
    C& set_maximum_spacial_order(DegreeType const& lower, DegreeType const& upper) { at<DegreeTypeProperty>("maximum_spacial_order").set(lower,upper); return *this; }

    //! \brief The maximum temporal order to be used for the flow function.
    DegreeType const& minimum_temporal_order() const { return at<DegreeTypeProperty>("minimum_temporal_order").get(); }
    C& set_minimum_temporal_order(DegreeType const& value) { at<DegreeTypeProperty>("minimum_temporal_order").set(value); return *this; }
    C& set_minimum_temporal_order(DegreeType const& lower, DegreeType const& upper) { at<DegreeTypeProperty>("minimum_temporal_order").set(lower,upper); return *this; }

    //! \brief The maximum temporal order to be used for the flow function.
    DegreeType const& maximum_temporal_order() const { return at<DegreeTypeProperty>("maximum_temporal_order").get(); }
    C& set_maximum_temporal_order(DegreeType const& value) { at<DegreeTypeProperty>("maximum_temporal_order").set(value); return *this; }
    C& set_maximum_temporal_order(DegreeType const& lower, DegreeType const& upper) { at<DegreeTypeProperty>("maximum_temporal_order").set(lower,upper); return *this; }

    //! Base properties

    ValidatedFunctionPatchFactoryInterface const& function_factory() const { return at<FunctionFactoryProperty>("function_factory").get(); }
    C& set_function_factory(ValidatedFunctionPatchFactoryInterface const& factory) { at<FunctionFactoryProperty>("function_factory").set(factory); return *this; }
    C& set_function_factory(SharedPointer<ValidatedFunctionPatchFactoryInterface> const& factory) { at<FunctionFactoryProperty>("function_factory").set(factory); return *this; }

    //! \brief The fraction L(f)*h used for a time step.
    //! \details The convergence of the Picard iteration is approximately Lf*h.
    RealType const& lipschitz_tolerance() const { return at<RealTypeProperty>("lipschitz_tolerance").get(); }
    C& set_lipschitz_tolerance(double const& value) { at<RealTypeProperty>("lipschitz_tolerance").set(value); return *this; }
    C& set_lipschitz_tolerance(double const& lower, double const& upper) { at<RealTypeProperty>("lipschitz_tolerance").set(lower,upper); return *this; }

    //! \brief The bounder to be used.
    BounderInterface const& bounder() const { return at<BounderProperty>("bounder").get(); }
    C& set_bounder(BounderInterface const& bounder) { at<BounderProperty>("bounder").set(bounder); return *this; }
    C& set_bounder(SharedPointer<BounderInterface> const& bounder) { at<BounderProperty>("bounder").set(bounder); return *this; }
};

template<> struct Configuration<AffineIntegrator> : public Configuration<BoundedIntegratorBase> {
    typedef Configuration<AffineIntegrator> C;
    typedef RangeConfigurationProperty<DegreeType> DegreeTypeProperty;

    Configuration() {
        add_property("spacial_order",DegreeTypeProperty(2u));
        add_property("temporal_order",DegreeTypeProperty(1u));
    }

    //! \brief The temporal order to be used for the flow function.
    DegreeType const& spacial_order() const { return at<DegreeTypeProperty>("spacial_order").get(); }
    C& set_spacial_order(DegreeType const& value) { at<DegreeTypeProperty>("spacial_order").set(value); return *this; }
    C& set_spacial_order(DegreeType const& lower, DegreeType const& upper) { at<DegreeTypeProperty>("spacial_order").set(lower,upper); return *this; }

    //! \brief The temporal order to be used for the flow function.
    DegreeType const& temporal_order() const { return at<DegreeTypeProperty>("temporal_order").get(); }
    C& set_temporal_order(DegreeType const& value) { at<DegreeTypeProperty>("temporal_order").set(value); return *this; }
    C& set_temporal_order(DegreeType const& lower, DegreeType const& upper) { at<DegreeTypeProperty>("temporal_order").set(lower,upper); return *this; }

    //! Base properties

    ValidatedFunctionPatchFactoryInterface const& function_factory() const { return at<FunctionFactoryProperty>("function_factory").get(); }
    C& set_function_factory(ValidatedFunctionPatchFactoryInterface const& factory) { at<FunctionFactoryProperty>("function_factory").set(factory); return *this; }
    C& set_function_factory(SharedPointer<ValidatedFunctionPatchFactoryInterface> const& factory) { at<FunctionFactoryProperty>("function_factory").set(factory); return *this; }

    //! \brief The fraction L(f)*h used for a time step.
    //! \details The convergence of the Picard iteration is approximately Lf*h.
    RealType const& lipschitz_tolerance() const { return at<RealTypeProperty>("lipschitz_tolerance").get(); }
    C& set_lipschitz_tolerance(double const& value) { at<RealTypeProperty>("lipschitz_tolerance").set(value); return *this; }
    C& set_lipschitz_tolerance(double const& lower, double const& upper) { at<RealTypeProperty>("lipschitz_tolerance").set(lower,upper); return *this; }

    //! \brief The bounder to be used.
    BounderInterface const& bounder() const { return at<BounderProperty>("bounder").get(); }
    C& set_bounder(BounderInterface const& bounder) { at<BounderProperty>("bounder").set(bounder); return *this; }
    C& set_bounder(SharedPointer<BounderInterface> const& bounder) { at<BounderProperty>("bounder").set(bounder); return *this; }
};

}

#endif /* ARIADNE_INTEGRATOR_HPP */
