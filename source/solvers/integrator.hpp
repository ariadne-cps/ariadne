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
#include "pronest/configuration_property.hpp"

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

    IntegratorInterface* clone() const override;
    Void _write(OutputStream& os) const override;

    FlowStepModelType flow_step(const ValidatedVectorMultivariateFunction& vector_field, const ExactBoxType& state_domain,
                                const StepSizeType& time_step, const UpperBoxType& bounding_box) const override;

    FlowStepModelType flow_step(const ValidatedVectorMultivariateFunction& differential_equation, const ExactBoxType& state_domain,
                                const Interval<StepSizeType>& time_domain, const ExactBoxType& parameter_domain, const UpperBoxType& bounding_box) const override;

    using BoundedIntegratorBase::flow_step;

  private:
    FlowStepModelType
    _flow_step(const ValidatedVectorMultivariateFunction& vector_field_or_differential_equation,
               const ExactBoxType& state_domain,
               const ExactIntervalType& time_domain,
               const ExactBoxType& parameter_domain,
               const UpperBoxType& bounding_box) const;
};

//! \ingroup DifferentialEquationSubModule
class GradedTaylorPicardIntegrator
        : public IntegratorBase
{
  public:
    //! \brief Default constructor.
    GradedTaylorPicardIntegrator(Configuration<GradedTaylorPicardIntegrator> const& config);

    Configuration<GradedTaylorPicardIntegrator> const& configuration() const;

    IntegratorInterface* clone() const override;
    Void _write(OutputStream& os) const override;

    FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const override;

    FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const override;

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
  public:
    TaylorSeriesIntegrator(Configuration<TaylorSeriesIntegrator> const& config);

    Configuration<TaylorSeriesIntegrator> const& configuration() const;

    IntegratorInterface* clone() const override;
    Void _write(OutputStream& os) const override;

    FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const override;

    FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const override;

    using BoundedIntegratorBase::flow_step;
};

//! \brief An integrator which computes the Taylor series of the flow function with remainder term.
class GradedTaylorSeriesIntegrator
    : public BoundedIntegratorBase
{
  public:
    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(Configuration<GradedTaylorSeriesIntegrator> const& config);

    Configuration<GradedTaylorSeriesIntegrator> const& configuration() const;

    IntegratorInterface* clone() const override;

    Void _write(OutputStream& os) const override;

    FlowStepModelType flow_step(const ValidatedVectorMultivariateFunction& vector_field, const ExactBoxType& state_domain,
                                const StepSizeType& time_step, const UpperBoxType& bounding_box) const override;

    FlowStepModelType flow_step(const ValidatedVectorMultivariateFunction& differential_equation, const ExactBoxType& state_domain,
                                const Interval<StepSizeType>& time_domain, const ExactBoxType& parameter_domain, const UpperBoxType& bounding_box) const override;

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

    IntegratorInterface* clone() const override;
    Void _write(OutputStream& os) const override;

    FlowStepModelType flow_step(const ValidatedVectorMultivariateFunction& vector_field, const ExactBoxType& state_domain,
                                const StepSizeType& time_step, const UpperBoxType& bounding_box) const override;

    FlowStepModelType flow_step(const ValidatedVectorMultivariateFunction& differential_equation, const ExactBoxType& state_domain,
                                const Interval<StepSizeType>& time_domain, const ExactBoxType& parameter_domain, const UpperBoxType& bounding_box) const override;

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
using Ariadne::GradedTaylorPicardIntegrator;
using Ariadne::TaylorSeriesIntegrator;
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
    typedef HandleListConfigurationProperty<Sweeper<FloatDP>> SweeperProperty;

    Configuration() {
        add_property("sweeper",SweeperProperty(ThresholdSweeper<FloatDP>(dp,Configuration<ThresholdSweeper<FloatDP>>().set_sweep_threshold(1e-8))));
    }

    Sweeper<FloatDP> const& sweeper() const { return at<SweeperProperty>("sweeper").get(); }
};

template<> struct Configuration<BoundedIntegratorBase> : public Configuration<IntegratorBase> {
public:
    typedef double RealType;
    typedef RangeConfigurationProperty<DegreeType> DegreeTypeProperty;
    typedef RangeConfigurationProperty<RealType> RealTypeProperty;
    typedef InterfaceListConfigurationProperty<BounderInterface> BounderProperty;

    Configuration() {
        add_property("bounder",BounderProperty(EulerBounder(Configuration<EulerBounder>())));
    }

    //! \brief The bounder to be used.
    BounderInterface const& bounder() const { return at<BounderProperty>("bounder").get(); }
};

template<> struct Configuration<TaylorPicardIntegrator> : public Configuration<BoundedIntegratorBase> {
    typedef Configuration<TaylorPicardIntegrator> C;
    typedef double RealType;
    typedef RangeConfigurationProperty<DegreeType> DegreeTypeProperty;
    typedef RangeConfigurationProperty<RealType> RealTypeProperty;
    typedef InterfaceListConfigurationProperty<BounderInterface> BounderProperty;
    typedef HandleListConfigurationProperty<Sweeper<FloatDP>> SweeperProperty;

    Configuration() {
        add_property("step_maximum_error",RealTypeProperty(1e-6,Log10SearchSpaceConverter<RealType>()));
        add_property("minimum_temporal_order",DegreeTypeProperty(0u));
        add_property("maximum_temporal_order",DegreeTypeProperty(12u));
    }

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

    //! \brief The sweeper to be used when creating a flow function.
    Sweeper<FloatDP> const& sweeper() const { return at<SweeperProperty>("sweeper").get(); }
    C& set_sweeper(Sweeper<FloatDP> const& sweeper) { at<SweeperProperty>("sweeper").set(sweeper); return *this; }
    C& set_sweeper(List<Sweeper<FloatDP>> const& sweepers) { at<SweeperProperty>("sweeper").set(sweepers); return *this; }

    //! \brief The bounder to be used.
    BounderInterface const& bounder() const { return at<BounderProperty>("bounder").get(); }
    C& set_bounder(BounderInterface const& bounder) { at<BounderProperty>("bounder").set(bounder); return *this; }
    C& set_bounder(SharedPointer<BounderInterface> const& bounder) { at<BounderProperty>("bounder").set(bounder); return *this; }
    C& set_bounder(List<SharedPointer<BounderInterface>> const& bounders) { at<BounderProperty>("bounder").set(bounders); return *this; }
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
        add_property("step_maximum_error",RealTypeProperty(1e-6,Log10SearchSpaceConverter<RealType>()));
        add_property("minimum_spacial_order",DegreeTypeProperty(1u));
        add_property("maximum_spacial_order",DegreeTypeProperty(4u));
        add_property("minimum_temporal_order",DegreeTypeProperty(4u));
        add_property("maximum_temporal_order",DegreeTypeProperty(12u));
    }

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

    //! \brief The sweeper to be used when creating a flow function.
    Sweeper<FloatDP> const& sweeper() const { return at<SweeperProperty>("sweeper").get(); }
    C& set_sweeper(Sweeper<FloatDP> const& sweeper) { at<SweeperProperty>("sweeper").set(sweeper); return *this; }
    C& set_sweeper(List<Sweeper<FloatDP>> const& sweepers) { at<SweeperProperty>("sweeper").set(sweepers); return *this; }

    //! \brief The bounder to be used.
    BounderInterface const& bounder() const { return at<BounderProperty>("bounder").get(); }
    C& set_bounder(BounderInterface const& bounder) { at<BounderProperty>("bounder").set(bounder); return *this; }
    C& set_bounder(SharedPointer<BounderInterface> const& bounder) { at<BounderProperty>("bounder").set(bounder); return *this; }
};

template<> struct Configuration<GradedTaylorPicardIntegrator> : public Configuration<IntegratorBase> {

    typedef Configuration<GradedTaylorPicardIntegrator> C;
    typedef double RealType;
    typedef RangeConfigurationProperty<DegreeType> DegreeTypeProperty;
    typedef RangeConfigurationProperty<RealType> RealTypeProperty;
    typedef HandleListConfigurationProperty<Sweeper<FloatDP>> SweeperProperty;

    Configuration() {
        add_property("step_maximum_error",RealTypeProperty(1e-6,Log10SearchSpaceConverter<RealType>()));
        add_property("temporal_order",DegreeTypeProperty(4u));
        add_property("error_refinement_minimum_improvement_percentage",DegreeTypeProperty(2u));
    }

    //! \brief The maximum error produced on a single step of integration.
    RealType const& step_maximum_error() const { return at<RealTypeProperty>("step_maximum_error").get(); }
    C& set_step_maximum_error(double const& value) { at<RealTypeProperty>("step_maximum_error").set(value); return *this; }
    C& set_step_maximum_error(double const& lower, double const& upper) { at<RealTypeProperty>("step_maximum_error").set(lower,upper); return *this; }

    //! \brief The temporal order to be used for the flow function.
    DegreeType const& temporal_order() const { return at<DegreeTypeProperty>("temporal_order").get(); }
    C& set_temporal_order(DegreeType const& value) { at<DegreeTypeProperty>("temporal_order").set(value); return *this; }
    C& set_temporal_order(DegreeType const& lower, DegreeType const& upper) { at<DegreeTypeProperty>("temporal_order").set(lower,upper); return *this; }

    //! \brief The minimum percentage of improvement required to continue refining the error, between 1 and 100.
    DegreeType const& error_refinement_minimum_improvement_percentage() const { return at<DegreeTypeProperty>("error_refinement_minimum_improvement_percentage").get(); }
    C& set_error_refinement_minimum_improvement_percentage(DegreeType const& value) { at<DegreeTypeProperty>("error_refinement_minimum_improvement_percentage").set(value); return *this; }
    C& set_error_refinement_minimum_improvement_percentage(DegreeType const& lower, DegreeType const& upper) { at<DegreeTypeProperty>("error_refinement_minimum_improvement_percentage").set(lower,upper); return *this; }

    //! Base properties

    //! \brief The sweeper to be used when creating a flow function.
    Sweeper<FloatDP> const& sweeper() const { return at<SweeperProperty>("sweeper").get(); }
    C& set_sweeper(Sweeper<FloatDP> const& sweeper) { at<SweeperProperty>("sweeper").set(sweeper); return *this; }
    C& set_sweeper(List<Sweeper<FloatDP>> const& sweepers) { at<SweeperProperty>("sweeper").set(sweepers); return *this; }
};

template<> struct Configuration<TaylorSeriesIntegrator> : public Configuration<BoundedIntegratorBase> {

    typedef Configuration<TaylorSeriesIntegrator> C;
    typedef double RealType;
    typedef RangeConfigurationProperty<DegreeType> DegreeTypeProperty;
    typedef RangeConfigurationProperty<RealType> RealTypeProperty;
    typedef InterfaceListConfigurationProperty<BounderInterface> BounderProperty;
    typedef HandleListConfigurationProperty<Sweeper<FloatDP>> SweeperProperty;

    Configuration() {
        add_property("order",DegreeTypeProperty(4u));
    }

    //! \brief The order to be used for the flow function.
    DegreeType const& order() const { return at<DegreeTypeProperty>("order").get(); }
    C& set_order(DegreeType const& value) { at<DegreeTypeProperty>("order").set(value); return *this; }
    C& set_order(DegreeType const& lower, DegreeType const& upper) { at<DegreeTypeProperty>("order").set(lower,upper); return *this; }

    //! Inherited properties

    //! \brief The sweeper to be used when creating a flow function.
    Sweeper<FloatDP> const& sweeper() const { return at<SweeperProperty>("sweeper").get(); }
    C& set_sweeper(Sweeper<FloatDP> const& sweeper) { at<SweeperProperty>("sweeper").set(sweeper); return *this; }
    C& set_sweeper(List<Sweeper<FloatDP>> const& sweepers) { at<SweeperProperty>("sweeper").set(sweepers); return *this; }

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

    //! Inherited properties

    //! \brief The sweeper to be used when creating a flow function.
    Sweeper<FloatDP> const& sweeper() const { return at<SweeperProperty>("sweeper").get(); }
    C& set_sweeper(Sweeper<FloatDP> const& sweeper) { at<SweeperProperty>("sweeper").set(sweeper); return *this; }
    C& set_sweeper(List<Sweeper<FloatDP>> const& sweepers) { at<SweeperProperty>("sweeper").set(sweepers); return *this; }

    //! \brief The bounder to be used.
    BounderInterface const& bounder() const { return at<BounderProperty>("bounder").get(); }
    C& set_bounder(BounderInterface const& bounder) { at<BounderProperty>("bounder").set(bounder); return *this; }
    C& set_bounder(SharedPointer<BounderInterface> const& bounder) { at<BounderProperty>("bounder").set(bounder); return *this; }
};

}

#endif /* ARIADNE_INTEGRATOR_HPP */
