/***************************************************************************
 *            solvers/bounder.hpp
 *
 *  Copyright  2018-20  Luca Geretti
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

/*! \file solvers/bounder.hpp
 *  \brief Classes for finding the bounds of the flow for differential equations.
 */

#ifndef ARIADNE_BOUNDER_HPP
#define ARIADNE_BOUNDER_HPP

#include "utility/typedefs.hpp"
#include "utility/attribute.hpp"
#include "algebra/algebra.hpp"
#include "function/domain.hpp"
#include "function/function_model.hpp"
#include "conclog/logging.hpp"

#include "pronest/searchable_configuration.hpp"
#include "pronest/configurable.hpp"
#include "pronest/configuration_interface.hpp"
#include "pronest/configuration_property.hpp"

#include "integrator_interface.hpp"

using namespace ConcLog;

namespace Ariadne {

using ProNest::Configuration;
using ProNest::Configurable;

class BoundingNotFoundException : public std::runtime_error {
  public:
    BoundingNotFoundException(const String& str) : std::runtime_error(str) { }
};

struct LipschitzTolerance : Attribute<ExactDouble> { using Attribute<ExactDouble>::Attribute; };
static const Generator<LipschitzTolerance> lipschitz_tolerance = Generator<LipschitzTolerance>();
static const LipschitzTolerance DEFAULT_LIPSCHITZ_TOLERANCE(0.5_x);
struct MinimumStepSize : Attribute<ExactDouble> { using Attribute<ExactDouble>::Attribute; };
static const Generator<MinimumStepSize> minimum_step_size = Generator<MinimumStepSize>();
static const MinimumStepSize DEFAULT_MINIMUM_STEP_SIZE(0.00000095367431640625_x);

//! \ingroup DifferentialEquationSubModule
//! \brief Interface for classes calculating the bounds of a flow.
class BounderInterface {
  public:
    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x)\f$ starting in \f$dom\f$  for time step \f$h\leq h_{sug}\f$.
    //! <br>
    //! Arguments: \f$f\f$ is the \a vector_field, \f$dom\f$ is the \a state_domain and \f$h_{sug}\f$ is the suggested time step.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& hsug) const = 0;

    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x,t)\f$ starting in \f$dom\f$  for time step \f$h\leq h_{sug}\f$.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& differential_equation, BoxDomainType const& state_domain, StepSizeType const& initial_time, StepSizeType const& suggested_time_step) const = 0;

    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x,a)\f$ starting in \f$D\f$ with \f$a\in A\f$ for time step \f$h\leq h_{sug}\f$.
    //! <br>
    //! Arguments: \f$f\f$ is the \a vector_field, \f$D\f$ is the \a state_domain and \f$h_{sug}\f$ is the suggested time step.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& vector_field, BoxDomainType const& state_domain, BoxDomainType const& parameter_domain, StepSizeType const& suggested_time_step) const = 0;

    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x,t,a)\f$ starting in \f$D\f$ with \f$a\in A\f$ for time step \f$h\leq h_{sug}\f$.
    //! Arguments: \f$f\f$ is the \a vector_field, \f$dom\f$ is the \a state_domain and \f$h_{sug}\f$ is the suggested time step.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& differential_equation, BoxDomainType const& state_domain, StepSizeType const& initial_time, BoxDomainType const& parameter_domain, StepSizeType const& suggested_time_step) const = 0;

  public:
    virtual Void _write(OutputStream& os) const = 0;
    virtual BounderInterface* clone() const = 0;

    friend inline OutputStream& operator<<(OutputStream& os, BounderInterface const& bounder) {
        bounder._write(os); return os; }

    virtual ~BounderInterface() = default;
};

class BounderBase : public BounderInterface, public Configurable<BounderBase> {
  public:
    BounderBase(Configuration<BounderBase> const& config);

    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& hsug) const override;
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, StepSizeType const& hsug) const override;
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, BoxDomainType const& A, StepSizeType const& hsug) const override = 0;
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, StepSizeType const& hsug) const override = 0;
};

//! \ingroup DifferentialEquationSubModule
//! \brief Compute bounds on the flow of a differential equation using a set-based Euler method.
class EulerBounder final : public BounderBase {
  public:
    EulerBounder(Configuration<EulerBounder> const& config);
    Configuration<EulerBounder> const& configuration() const;

    using BounderBase::compute;

    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, BoxDomainType const& A, StepSizeType const& hsug) const override;
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, StepSizeType const& hsug) const override;
    virtual Void _write(OutputStream& os) const override;
    virtual BounderInterface* clone() const override;
  private:
    // Compute the bounds on the reach set of dx/dt = f(x,t,a) starting at time t, for a suggested step size of h.
    Pair<StepSizeType,UpperBoxType> _compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, StepSizeType const& hsug) const;
    // Compute a bounding box for the flow of f starting in D at T.lower() to time T.upper() with parameters A, assuming flow remains in B,
    // but widening D by INITIAL_WIDENING and the flow vectors by VECTOR_WIDENING.
    UpperBoxType _formula(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, IntervalDomainType const& T, BoxDomainType const& A, UpperBoxType const& B, PositiveFloatDP INITIAL_WIDENING, PositiveFloatDP VECTOR_WIDENING) const;
    UpperBoxType _refinement(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, IntervalDomainType const& T, BoxDomainType const& A, UpperBoxType const& B) const;
};

} // namespace Ariadne

namespace ProNest {

using Ariadne::BounderBase;
using Ariadne::EulerBounder;
using Ariadne::DegreeType;

template<> struct Configuration<BounderBase> : public SearchableConfiguration {
public:
    typedef Configuration<BounderBase> C;
    typedef double RealType;
    typedef RangeConfigurationProperty<RealType> RealTypeProperty;

    Configuration() {
        add_property("lipschitz_tolerance",RealTypeProperty(0.5,Log2SearchSpaceConverter<RealType>()));
        add_property("minimum_step_size",RealTypeProperty(1e-8,Log2SearchSpaceConverter<RealType>()));
    }

    //! \brief The fraction L(f)*h used for a time step.
    //! \details The convergence of the Picard iteration is approximately Lf*h.
    RealType const& lipschitz_tolerance() const { return at<RealTypeProperty>("lipschitz_tolerance").get(); }

    //! \brief The minimum allowable step size for finding a bound.
    //! Increasing this value prevents getting too low (and consequently slow) steps.
    RealType const& minimum_step_size() const { return at<RealTypeProperty>("minimum_step_size").get(); }
};

template<> struct Configuration<EulerBounder> : public Configuration<BounderBase> {
public:
    typedef Configuration<EulerBounder> C;
    typedef double RealType;
    typedef RangeConfigurationProperty<RealType> RealTypeProperty;

    Configuration() {
    }

    //! Inherited properties

    //! \brief The fraction L(f)*h used for a time step.
    //! \details The convergence of the Picard iteration is approximately Lf*h.
    RealType const& lipschitz_tolerance() const { return at<RealTypeProperty>("lipschitz_tolerance").get(); }
    C& set_lipschitz_tolerance(double const& value) { at<RealTypeProperty>("lipschitz_tolerance").set(value); return *this; }
    C& set_lipschitz_tolerance(double const& lower, double const& upper) { at<RealTypeProperty>("lipschitz_tolerance").set(lower,upper); return *this; }

    //! \brief The minimum allowable step size for finding a bound.
    //! Increasing this value prevents getting too low (and consequently slow) steps.
    RealType const& minimum_step_size() const { return at<RealTypeProperty>("minimum_step_size").get(); }
    C& set_minimum_step_size(RealType const& value) { at<RealTypeProperty>("minimum_step_size").set(value); return *this; }
    C& set_minimum_step_size(RealType const& lower, RealType const& upper) { at<RealTypeProperty>("minimum_step_size").set(lower,upper); return *this; }
};

}

#endif /* ARIADNE_BOUNDER_HPP */
