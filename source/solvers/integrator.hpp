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

#include "utility/declarations.hpp"
#include "utility/attribute.hpp"
#include "numeric/dyadic.hpp"
#include "conclog/logging.hpp"
#include "utility/pointer.hpp"
#include "function/affine.hpp"
#include "algebra/sweeper.hpp"
#include "algebra/matrix.hpp"

#include "function/function_patch.hpp"
#include "function/taylor_function.hpp"

using namespace ConcLog;

namespace Ariadne {

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
    : public IntegratorInterface
{
  protected:
    //! \brief Construct with a sweeper for the function factory
    IntegratorBase(Sweeper<FloatDP> sweeper);
  public:

    //! \brief The class which constructs functions for representing the flow.
    const ValidatedFunctionPatchFactory& function_factory() const;
    //! \brief Set the class which constructs functions for representing the flow.
    Void set_function_factory(const ValidatedFunctionPatchFactory& factory);

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const Suggestion<StepSizeType>& suggested_time_step) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step) const;

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

  private:
    ValidatedFunctionPatchFactory _function_factory;
};

class BoundedIntegratorBase : public IntegratorBase {
  protected:
    BoundedIntegratorBase(Sweeper<FloatDP> sweeper, LipschitzTolerance lipschitz_tolerance);
  public:

    //! \brief The class that computes bounds.
    const BounderInterface& bounder() const;
    //! \brief Set the class that computes bounds.
    Void set_bounder(const BounderInterface& bounder);

    virtual Pair<StepSizeType,UpperBoxType>
    flow_bounds(const ValidatedVectorMultivariateFunction& vector_field,
                const ExactBoxType& state_domain,
                const Suggestion<StepSizeType>& maximum_time_step) const;


    virtual Pair<StepSizeType,UpperBoxType>
    flow_bounds(const ValidatedVectorMultivariateFunction& vector_field,
                const ExactBoxType& state_domain,
                const ExactBoxType& parameter_domain,
                const Suggestion<StepSizeType>& maximum_time_step) const;


    virtual Pair<StepSizeType,UpperBoxType>
    flow_bounds(const ValidatedVectorMultivariateFunction& differential_equation,
                const ExactBoxType& state_domain,
                const StepSizeType& starting_time,
                const ExactBoxType& parameter_domain,
                const Suggestion<StepSizeType>& maximum_time_step) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step) const override;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const Suggestion<StepSizeType>& suggested_time_step) const override;

    using IntegratorBase::flow_step;

  private:
    BounderPointer _bounder_ptr;
};

//! \ingroup DifferentialEquationSubModule
//! \brief An integrator which uses a validated Picard iteration on Taylor models.
class TaylorPicardIntegrator
    : public BoundedIntegratorBase
{
    ExactDouble _step_maximum_error;
    Sweeper<FloatDP> _sweeper;
    DegreeType _minimum_temporal_order;
    DegreeType _maximum_temporal_order;
  public:
    //! \brief Default constructor.
    TaylorPicardIntegrator(StepMaximumError err);

    //! \brief Constructor.
    TaylorPicardIntegrator(StepMaximumError lerr, Sweeper<FloatDP> const& sweeper, LipschitzTolerance lip,
                           MinimumTemporalOrder minto, MaximumTemporalOrder maxto);

    //! \brief The order of the method in time.
    DegreeType minimum_temporal_order() const { return this->_minimum_temporal_order; }
    Void set_minimum_temporal_order(Nat m) { this->_minimum_temporal_order=m; }
    DegreeType maximum_temporal_order() const { return this->_maximum_temporal_order; }
    Void set_maximum_temporal_order(Nat m) { this->_maximum_temporal_order=m; }
    //! \brief  Set the sweep threshold of the Taylor model.
    Sweeper<FloatDP> const& sweeper() const { return this->_sweeper; }
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

//! \ingroup DifferentialEquationSubModule
class GradedTaylorPicardIntegrator
        : public IntegratorBase
{
    ExactDouble _step_maximum_error;
    GradedThresholdSweeper<FloatDP> _sweeper;
    ExactDouble _error_refinement_minimum_improvement_percentage;
    Nat _maximum_error_refinement_iterations;
    DegreeType _order;
    ApproximateDouble _sweep_threshold;
    Bool _diagnostics;
  public:
    //! \brief Construct with degree truncation only.
    GradedTaylorPicardIntegrator(StepMaximumError err, Order order);
    //! \brief Construct with degree and coefficient-magnitude truncation.
    GradedTaylorPicardIntegrator(StepMaximumError err, Order order, StepSweepThreshold sweep_threshold);

    //! \brief The order of the method.
    DegreeType order() const { return this->_order; }
    Void set_order(Nat m) {
        this->_order=m;
        this->_sweeper = GradedThresholdSweeper<FloatDP>(DoublePrecision(),m,FloatDP(cast_exact(this->_sweep_threshold),DoublePrecision()));
        this->set_function_factory(ValidatedFunctionPatchFactory(make_taylor_function_patch_factory(this->_sweeper)));
    }
    //! \brief  Set the maximum error of a single step.
    ExactDouble step_maximum_error() const { return this->_step_maximum_error; }
    Void set_step_maximum_error(ApproximateDouble e) { _step_maximum_error = cast_exact(e); }
    ExactDouble error_refinement_minimum_improvement_percentage() const { return this->_error_refinement_minimum_improvement_percentage; }
    Void set_error_refinement_minimum_improvement_percentage(ApproximateDouble e) { _error_refinement_minimum_improvement_percentage = cast_exact(e); }
    //! \brief Limit error-refinement iterations; zero preserves the historical unlimited behaviour.
    Void set_maximum_error_refinement_iterations(Nat n) { _maximum_error_refinement_iterations=n; }
    Nat maximum_error_refinement_iterations() const { return _maximum_error_refinement_iterations; }
    Void set_diagnostics(Bool diagnostics) { _diagnostics=diagnostics; }

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

    using IntegratorBase::flow_step;

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
    Void set_sweeper(Sweeper<FloatDP> const& sweeper) {
        _sweeper = sweeper;
        this->set_function_factory(ValidatedFunctionPatchFactory(make_taylor_function_patch_factory(_sweeper)));
    }

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

//! \ingroup DifferentialEquationSubModule
//! \brief An integrator which computes the Taylor series of the flow function with remainder term.
class TaylorSeriesBounderIntegrator
    : public TaylorSeriesIntegrator
{
  public:
    //! \brief Constructor.
    TaylorSeriesBounderIntegrator(StepMaximumError err, Order order);

    //! \brief Constructor.
    TaylorSeriesBounderIntegrator(StepMaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzTolerance lip,
                                   Order order);

    virtual TaylorSeriesBounderIntegrator* clone() const { return new TaylorSeriesBounderIntegrator(*this); }
    virtual Void _write(OutputStream& os) const;

    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const Suggestion<StepSizeType>& suggested_time_step) const;

    using TaylorSeriesIntegrator::flow_step;

  private:
    ExactDouble _step_maximum_error;
};




//! \brief An integrator which computes the Taylor series of the flow function with remainder term.
class GradedTaylorSeriesIntegrator
    : public BoundedIntegratorBase
{
    ExactDouble _step_maximum_error;
    Sweeper<FloatDP> _sweeper;
    DegreeType _minimum_spacial_order;
    DegreeType _minimum_temporal_order;
    DegreeType _maximum_spacial_order;
    DegreeType _maximum_temporal_order;
  public:
    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(StepMaximumError err);

    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(StepMaximumError err, Sweeper<FloatDP> const& sweeper);

    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(StepMaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzTolerance lip);

    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(StepMaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzTolerance lip, MaximumTemporalOrder maxto);

    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(StepMaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzTolerance lip,
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
    Void set_sweeper(Sweeper<FloatDP> const& sweeper) {
        _sweeper = sweeper;
        this->set_function_factory(ValidatedFunctionPatchFactory(make_taylor_function_patch_factory(_sweeper)));
    }
    //! \brief  Set the sweep threshold of the Taylor model.
    ExactDouble step_maximum_error() const { return this->_step_maximum_error; }
    Void set_step_maximum_error(ApproximateDouble e) { _step_maximum_error = cast_exact(e); }

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

//! \brief A factorisation of a parameterised state as x(s)=c+A*y(s).
//!
//! The normalised mapping y(s) retains the Taylor-model dependence on the
//! original enclosure parameters.  Keeping this object separate from the local
//! flow is the first step towards propagating Flow*-style local initial sets
//! without storing mutable state in the integrator itself.
class PreconditionedTaylorSeriesState {
  private:
    Vector<FloatDP> _centre;
    Matrix<FloatDP> _linear_map;
    ExactBoxType _local_domain;
    ValidatedVectorMultivariateFunctionPatch _normalised_mapping;
  public:
    PreconditionedTaylorSeriesState(
        Vector<FloatDP> centre,
        Matrix<FloatDP> linear_map,
        ExactBoxType local_domain,
        ValidatedVectorMultivariateFunctionPatch normalised_mapping)
        : _centre(std::move(centre)),
          _linear_map(std::move(linear_map)),
          _local_domain(std::move(local_domain)),
          _normalised_mapping(std::move(normalised_mapping)) { }

    Vector<FloatDP> const& centre() const { return _centre; }
    Matrix<FloatDP> const& linear_map() const { return _linear_map; }
    ExactBoxType const& local_domain() const { return _local_domain; }
    ValidatedVectorMultivariateFunctionPatch const& normalised_mapping() const { return _normalised_mapping; }
    ExactBoxType const parameter_domain() const { return _normalised_mapping.domain(); }
};

//! \brief Result of one preconditioned graded Taylor-series step.
class PreconditionedTaylorSeriesStep {
  private:
    StepSizeType _time_step;
    ValidatedVectorMultivariateFunctionPatch _flowpipe_mapping;
    ValidatedVectorMultivariateFunctionPatch _evolved_mapping;
    PreconditionedTaylorSeriesState _final_state;
  public:
    PreconditionedTaylorSeriesStep(
        StepSizeType time_step,
        ValidatedVectorMultivariateFunctionPatch flowpipe_mapping,
        ValidatedVectorMultivariateFunctionPatch evolved_mapping,
        PreconditionedTaylorSeriesState final_state)
        : _time_step(time_step),
          _flowpipe_mapping(std::move(flowpipe_mapping)),
          _evolved_mapping(std::move(evolved_mapping)),
          _final_state(std::move(final_state)) { }

    StepSizeType const& time_step() const { return _time_step; }
    ValidatedVectorMultivariateFunctionPatch const& flowpipe_mapping() const { return _flowpipe_mapping; }
    ValidatedVectorMultivariateFunctionPatch const& evolved_mapping() const { return _evolved_mapping; }
    PreconditionedTaylorSeriesState const& final_state() const { return _final_state; }
};

enum class TaylorSeriesPreconditioning { IDENTITY, QR };

//! \brief A graded Taylor-series integrator with explicit preconditioning support.
class PreconditionedGradedTaylorSeriesIntegrator
    : public GradedTaylorSeriesIntegrator
{
  private:
    TaylorSeriesPreconditioning _preconditioning=TaylorSeriesPreconditioning::QR;
    Bool _diagnostics=false;
    Bool _carried_expansion_diagnostics=false;
  public:
    //! \brief Construct with QR-compatible fixed default orders.
    PreconditionedGradedTaylorSeriesIntegrator(StepMaximumError err);

    using GradedTaylorSeriesIntegrator::GradedTaylorSeriesIntegrator;
    using GradedTaylorSeriesIntegrator::flow_step;

    TaylorSeriesPreconditioning preconditioning() const { return _preconditioning; }
    Void set_preconditioning(TaylorSeriesPreconditioning value) { _preconditioning=value; }
    Bool diagnostics() const { return _diagnostics; }
    Void set_diagnostics(Bool value) { _diagnostics=value; }
    Bool carried_expansion_diagnostics() const { return _carried_expansion_diagnostics; }
    Void set_carried_expansion_diagnostics(Bool value) { _carried_expansion_diagnostics=value; }

    virtual PreconditionedGradedTaylorSeriesIntegrator* clone() const override {
        return new PreconditionedGradedTaylorSeriesIntegrator(*this);
    }
    virtual Void _write(OutputStream& os) const override;

    //! \brief Factor a parameterised state as x(s)=c+A*y(s), retaining the
    //! complete Taylor-model dependence in y(s).
    PreconditionedTaylorSeriesState
    precondition(const ValidatedVectorMultivariateFunctionPatch& state) const;

    //! \brief Propagate a preconditioned local initial set for one validated
    //! time step, retaining its parameter dependence.
    PreconditionedTaylorSeriesStep
    step(const ValidatedVectorMultivariateFunction& vector_field,
         const PreconditionedTaylorSeriesState& state,
         const Suggestion<StepSizeType>& suggested_time_step) const;

    //! \brief Compute a graded Taylor flow after diagonal affine
    //! preconditioning of the state domain to the unit box.
    //!
    //! The returned model is expressed again on the original physical state
    //! domain, so it remains compatible with IntegratorInterface and existing
    //! evolvers.
    virtual FlowStepModelType
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const override;
};



//! \brief An integrator computes a approximation to the flow which is affine in space.
//! \internal This code is written to allow higher-spacial order approximations.
class AffineIntegrator
    : public BoundedIntegratorBase
{
    DegreeType _spacial_order;
    DegreeType _temporal_order;
  public:
    AffineIntegrator(SpacialOrder spacial_order, TemporalOrder temporal_order);

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

#endif /* ARIADNE_INTEGRATOR_HPP */
