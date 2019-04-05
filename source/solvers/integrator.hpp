/***************************************************************************
 *            integrator.hpp
 *
 *  Copyright  2006-10  Pieter Collins
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

/*! \file integrator.hpp
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

namespace Ariadne {

class Real;
template<class X> class Vector;
template<class X> class Differential;
template<class X> class Procedure;
template<class I, class X> class Polynomial;
template<class X> using MultivariatePolynomial = Polynomial<MultiIndex,X>;
typedef Differential<ValidatedNumericType> ValidatedDifferential;
typedef Vector< Procedure<ValidatedNumericType> > ValidatedVectorProcedure;
typedef FunctionModelFactoryInterface<ValidatedTag> ValidatedFunctionModelDPFactoryInterface;
typedef SharedPointer<const ValidatedFunctionModelDPFactoryInterface> FunctionFactoryPointer;
typedef SharedPointer<const BounderInterface> BounderPointer;

struct LipschitzConstant : Attribute<double> { LipschitzConstant(double v) : Attribute<double>(v) { } };
struct StepMaximumError : Attribute<double> { StepMaximumError(double v) : Attribute<double>(v) { } };
struct StepSweepThreshold : Attribute<double> { StepSweepThreshold(double v) : Attribute<double>(v) { } };
struct Order : Attribute<DegreeType> { Order(DegreeType v) : Attribute<DegreeType>(v) { } };
struct SpacialOrder : Attribute<DegreeType> { SpacialOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct TemporalOrder : Attribute<DegreeType> { TemporalOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct MinimumSpacialOrder : Attribute<DegreeType> { MinimumSpacialOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct MinimumTemporalOrder : Attribute<DegreeType> { MinimumTemporalOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct MaximumSpacialOrder : Attribute<DegreeType> { MaximumSpacialOrder(DegreeType v) : Attribute<DegreeType>(v) { } };
struct MaximumTemporalOrder : Attribute<DegreeType> { MaximumTemporalOrder(DegreeType v) : Attribute<DegreeType>(v) { } };

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
    StepSizeType maximum_step_size() const { return this->_maximum_step_size; }
    Void set_maximum_step_size(StepSizeType hmax) { this->_maximum_step_size = hmax; }

    //! \brief The class which constructs functions for representing the flow.
    const ValidatedFunctionModelDPFactoryInterface& function_factory() const;
    //! \brief Set the class which constructs functions for representing the flow.
    Void set_function_factory(const ValidatedFunctionModelDPFactoryInterface& factory);

    //! \brief The class that computes bounds.
    const BounderInterface& bounder() const;
    //! \brief Set the class that computes bounds.
    Void set_bounder(const BounderInterface& bounder);

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

    virtual ValidatedVectorMultivariateFunctionModelDP
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              StepSizeType& suggested_time_step) const;

    virtual ValidatedVectorMultivariateFunctionModelDP
    flow_to(const ValidatedVectorMultivariateFunction& vector_field,
         const ExactBoxType& state_domain,
         const Real& time) const;

    //! \brief Solve \f$\der{\phi}(x,t)=f(\phi(x,t))\f$ for \f$t\in[0,T_{\max}]\f$.
    virtual List<ValidatedVectorMultivariateFunctionModelDP>
    flow(const ValidatedVectorMultivariateFunction& vector_field,
         const ExactBoxType& state_domain,
         const Real& minimum_time,
         const Real& maximum_time) const;

    //! \brief Solve \f$\der{\phi}(x,t)=f(\phi(x,t))\f$ for \f$t\in[0,T_{\max}]\f$.
    virtual List<ValidatedVectorMultivariateFunctionModelDP>
    flow(const ValidatedVectorMultivariateFunction& vector_field,
         const ExactBoxType& state_domain,
         const Real& maximum_time) const;


    virtual ValidatedVectorMultivariateFunctionModelDP
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const = 0;

    virtual ValidatedVectorMultivariateFunctionModelDP
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const = 0;

  public:
    double _maximum_error;
    double _lipschitz_tolerance;
    StepSizeType _maximum_step_size;
    FunctionFactoryPointer _function_factory_ptr;
    BounderPointer _bounder_ptr;
};

//! \brief An integrator which uses a validated Picard iteration on Taylor models.
class TaylorPicardIntegrator
    : public IntegratorBase
{
    double _step_maximum_error;
    double _step_sweep_threshold;
    DegreeType _maximum_temporal_order;
  public:
    //! \brief Default constructor.
    TaylorPicardIntegrator(MaximumError err)
        : IntegratorBase(err,SweepThreshold(err/1024),LipschitzConstant(0.5))
        , _step_maximum_error(err/128), _step_sweep_threshold(err/(1024*1024)), _maximum_temporal_order(12) { }

    //! \brief Constructor.
    TaylorPicardIntegrator(MaximumError err, SweepThreshold swp, LipschitzConstant lip,
                           StepMaximumError lerr, StepSweepThreshold lswp, MaximumTemporalOrder maxto)
        : IntegratorBase(err,swp,lip), _step_maximum_error(lerr), _step_sweep_threshold(lswp), _maximum_temporal_order(maxto) { }

    //! \brief The order of the method in time.
    DegreeType maximum_temporal_order() const { return this->_maximum_temporal_order; }
    Void set_maximum_temporal_order(DegreeType m) { this->_maximum_temporal_order=m; }
    //! \brief  Set the sweep threshold of the Taylor model for a single step.
    double step_sweep_threshold() const { return this->_step_sweep_threshold; }
    Void set_step_sweep_threshold(double lt) { _step_sweep_threshold = lt; }
    //! \brief  Set the maximum error of a single step.
    double step_maximum_error() const { return this->_step_maximum_error; }
    Void set_step_maximum_error(double e) { _step_maximum_error = e; }

    virtual TaylorPicardIntegrator* clone() const { return new TaylorPicardIntegrator(*this); }
    virtual Void write(OutputStream& os) const;

    virtual ValidatedVectorMultivariateFunctionModelDP
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const;

    virtual ValidatedVectorMultivariateFunctionModelDP
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const;

    using IntegratorBase::flow_step;

  private:
    ValidatedVectorMultivariateFunctionModelDP
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
    double _step_sweep_threshold;
    DegreeType _order;
  public:
    //! \brief Constructor.
    TaylorSeriesIntegrator(MaximumError err, Order order);

    //! \brief Constructor.
    TaylorSeriesIntegrator(MaximumError err, SweepThreshold gswp, LipschitzConstant lip,
                           StepSweepThreshold lswp, Order order);

    //! \brief The order of the method in space and time.
    DegreeType order() const { return this->_order; }
    Void set_order(DegreeType n) { this->_order=n; }
    //! \brief  Set the sweep threshold of the Taylor model representing a single step.
    double step_sweep_threshold() const { return this->_step_sweep_threshold; }
    Void set_step_sweep_threshold(double lswp) { _step_sweep_threshold = lswp; }

    virtual TaylorSeriesIntegrator* clone() const { return new TaylorSeriesIntegrator(*this); }
    virtual Void write(OutputStream& os) const;

    virtual Pair<StepSizeType,UpperBoxType>
    flow_bounds(const ValidatedVectorMultivariateFunction& vector_field,
                const ExactBoxType& state_domain,
                const StepSizeType& suggested_time_step) const;

    virtual ValidatedVectorMultivariateFunctionModelDP
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const;

    virtual ValidatedVectorMultivariateFunctionModelDP
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const;

    using IntegratorBase::flow_step;

  private:

};

//! \brief An integrator which computes the Taylor polynomial of the flow function and estimates a high-order remainder term.
class TaylorPolynomialIntegrator
    : public IntegratorBase
{
    DegreeType _order;
  public:
    //! \brief Constructor.
    TaylorPolynomialIntegrator(MaximumError e, LipschitzConstant l, Order order);

    virtual ValidatedVectorMultivariateFunctionModelDP
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              StepSizeType& time_step) const override;

    virtual ValidatedVectorMultivariateFunctionModelDP
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const override;

    virtual ValidatedVectorMultivariateFunctionModelDP
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const override { ARIADNE_NOT_IMPLEMENTED; }

    virtual TaylorPolynomialIntegrator* clone() const override { return new TaylorPolynomialIntegrator(*this); }

    virtual Void write(OutputStream& os) const override { os << "TaylorPolynomialIntegrator(order="<<this->_order<<")"; }
};




//! \brief An integrator which computes the Taylor series of the flow function with remainder term.
class GradedTaylorSeriesIntegrator
    : public IntegratorBase
{
    double _step_maximum_error;
    double _step_sweep_threshold;
    DegreeType _minimum_spacial_order;
    DegreeType _minimum_temporal_order;
    DegreeType _maximum_spacial_order;
    DegreeType _maximum_temporal_order;
  public:
    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(MaximumError err);

    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(MaximumError err, SweepThreshold swp, LipschitzConstant lip=0.5);

    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(MaximumError err, SweepThreshold gswp, LipschitzConstant lip,
                           StepMaximumError lerr, StepSweepThreshold lswp, MaximumTemporalOrder maxto);

    //! \brief Constructor.
    GradedTaylorSeriesIntegrator(MaximumError err, SweepThreshold gswp, LipschitzConstant lip,
                           StepMaximumError lerr, StepSweepThreshold lswp,
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
    //! \brief  Set the sweep threshold of the Taylor model representing a single step.
    double step_sweep_threshold() const { return this->_step_sweep_threshold; }
    Void set_step_sweep_threshold(double lswp) { _step_sweep_threshold = lswp; }
    //! \brief  Set the sweep threshold of the Taylor model.
    double step_maximum_error() const { return this->_step_maximum_error; }
    Void set_step_maximum_error(double e) { _step_maximum_error = e; }

    virtual GradedTaylorSeriesIntegrator* clone() const { return new GradedTaylorSeriesIntegrator(*this); }
    virtual Void write(OutputStream& os) const;

    virtual Pair<StepSizeType,UpperBoxType>
    flow_bounds(const ValidatedVectorMultivariateFunction& vector_field,
                const ExactBoxType& state_domain,
                const StepSizeType& suggested_time_step) const;

    virtual ValidatedVectorMultivariateFunctionModelDP
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const;

    virtual ValidatedVectorMultivariateFunctionModelDP
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
    virtual Void write(OutputStream& os) const;

    virtual ValidatedVectorMultivariateFunctionModelDP
    flow_step(const ValidatedVectorMultivariateFunction& vector_field,
              const ExactBoxType& state_domain,
              const StepSizeType& time_step,
              const UpperBoxType& bounding_box) const;

    virtual ValidatedVectorMultivariateFunctionModelDP
    flow_step(const ValidatedVectorMultivariateFunction& differential_equation,
              const ExactBoxType& state_domain,
              const Interval<StepSizeType>& time_domain,
              const ExactBoxType& parameter_domain,
              const UpperBoxType& bounding_box) const;

    using IntegratorBase::flow_step;

    //! \brief Compute the derivative of the flow of f at time zero within \a dom.
    Vector<ValidatedDifferential>
    flow_derivative(const ValidatedVectorMultivariateFunction& f,
                    const Vector<ValidatedNumericType>& dom) const;
};


} // namespace Ariadne

#endif /* ARIADNE_INTEGRATOR_HPP */
