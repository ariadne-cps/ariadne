/***************************************************************************
 *            differential_inclusion.hpp
 *
 *  Copyright  2008-18  Luca Geretti, Pieter Collins, Sanja Zivanovic
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

/*! \file differential_inclusion.hpp
 *  \brief Methods for computing solutions of differential inclusions.
 */

#ifndef ARIADNE_DIFFERENTIAL_INCLUSION_HPP
#define ARIADNE_DIFFERENTIAL_INCLUSION_HPP

#include "../utility/typedefs.hpp"
#include "../utility/attribute.hpp"
#include "../algebra/sweeper.hpp"
#include "../algebra/algebra.hpp"
#include "../function/domain.hpp"
#include "../function/function_model.hpp"
#include "../function/formula.hpp"
#include "../symbolic/expression_set.hpp"
#include "../output/logging.hpp"
#include "../solvers/integrator_interface.hpp"
#include "inclusion_vector_field.hpp"

namespace Ariadne {

class Real;

struct StepSize : public Attribute<FloatDP> { };
struct NumberOfStepsBetweenSimplifications : public Attribute<Nat> { };
struct NumberOfVariablesToKeep : public Attribute<Nat> { };

Generator<StepSize> step_size;
Generator<NumberOfStepsBetweenSimplifications> number_of_steps_between_simplifications;
Generator<NumberOfVariablesToKeep> number_of_variables_to_keep;

using ValidatedScalarMultivariateFunctionModelType = ValidatedScalarMultivariateFunctionModelDP;
using ValidatedVectorMultivariateFunctionModelType = ValidatedVectorMultivariateFunctionModelDP;

using ThresholdSweeperDP = ThresholdSweeper<FloatDP>;
using GradedSweeperDP = GradedSweeper<FloatDP>;
using GradedThresholdSweeperDP = GradedThresholdSweeper<FloatDP>;
using SweeperDP = Sweeper<FloatDP>;

using ApproximateTimeStepType = PositiveFloatDPApproximation;
using ExactTimeStepType = PositiveFloatDPValue;

BoxDomainType initial_ranges_to_box(RealVariablesBox const& var_ranges);

inline Vector<FloatDPValue> const& cast_exact(Vector<FloatDPError> const& v) {
    return reinterpret_cast<Vector<FloatDPValue>const&>(v); }

FloatDP volume(Vector<ApproximateIntervalType> const& box);

inline Bool refines(Vector<UpperIntervalType> const& v1, UpperBoxType const& bx2) {
    return refines(v1,static_cast<Vector<UpperIntervalType>const&>(bx2)); }

Box<Interval<FloatDPValue>> over_approximation(Box<Interval<Real>> const&);

ValidatedVectorMultivariateFunctionModelDP add_errors(ValidatedVectorMultivariateFunctionModelDP phi, Vector<ErrorType> const& e);

ValidatedVectorMultivariateFunction build_Fw(ValidatedVectorMultivariateFunction const& F, Vector<ValidatedScalarMultivariateFunction> const& w);

template<class F1, class F2, class F3, class... FS> decltype(auto) combine(F1 const& f1, F2 const& f2, F3 const& f3, FS const& ... fs) {
    return combine(combine(f1,f2),f3,fs...); }
template<class F1, class F2, class F3, class... FS> decltype(auto) join(F1 const& f1, F2 const& f2, F3 const& f3, FS const& ... fs) {
    return join(join(f1,f2),f3,fs...); }

//! \brief The function dexp(x)=(exp(x)-1)/x.
//! Note that the function is positive and monotone increasing.
template<class F> PositiveUpperBound<F> dexp(UpperBound<F> const& x) {
    if(x.raw()>=0) { return cast_positive(exp(x)-1)/cast_positive(cast_exact(x)); }
    else { return cast_positive(cast_exact(1-exp(x)))/cast_positive(-x); }
}

template<class F> PositiveLowerBound<F> dexp(LowerBound<F> const& x) {
    if(x.raw()>=0) { return cast_positive(exp(x)-1)/cast_positive(cast_exact(x)); }
    else { return cast_positive(cast_exact(1-exp(x)))/cast_positive(-x); }
}

template<class F> PositiveBounds<F> dexp(Bounds<F> const& x) {
    return PositiveBounds<F>(dexp(x.lower()),dexp(x.upper()));
}

struct C1Norms {
    FloatDPError K;
    Vector<FloatDPError> Kj;
    FloatDPError pK;
    Vector<FloatDPError> pKj;
    FloatDPError L;
    Vector<FloatDPError> Lj;
    FloatDPError pL;
    Vector<FloatDPError> pLj;
    FloatDPError H;
    Vector<FloatDPError> Hj;
    FloatDPError pH;
    Vector<FloatDPError> pHj;
    FloatDPError expLambda;
    FloatDPError expL;

    C1Norms(FloatDPError const&,Vector<FloatDPError> const&,FloatDPError const&,Vector<FloatDPError> const&,FloatDPError const&,Vector<FloatDPError> const&,FloatDPError const&,Vector<FloatDPError> const&,FloatDPError const&,Vector<FloatDPError> const&,FloatDPError const&,Vector<FloatDPError> const&,FloatDPError const&,FloatDPError const&);
    Tuple<FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,FloatDPError> values() const;

    SizeType dimension() const { return _dimension; }
private:
    SizeType _dimension;
};

inline std::ostream& operator << (std::ostream& os, const C1Norms& n) {
    os << "K=" << n.K << ", Kj=" << n.Kj << ", K'=" << n.pK << ", Kj'=" << n.Kj <<
          ", L=" << n.L << ", Lj=" << n.Lj << ", L'=" << n.pL << ", Lj'=" << n.Lj <<
          ", H=" << n.H << ", Hj=" << n.Hj << ", H'=" << n.pH << ", Hj'=" << n.Hj <<
          ", expLambda=" << n.expLambda << ", expL=" << n.expL;
    return os;
}

C1Norms compute_norms(InclusionVectorField const&, PositiveFloatDPValue const&, UpperBoxType const&);

enum class InputApproximationKind : std::uint8_t { ZERO, CONSTANT, AFFINE, SINUSOIDAL, PIECEWISE };

inline std::ostream& operator << (std::ostream& os, const InputApproximationKind& kind) {
    switch (kind) {
        case InputApproximationKind::ZERO: os << "ZERO"; break;
        case InputApproximationKind::CONSTANT: os << "CONSTANT"; break;
        case InputApproximationKind::AFFINE: os << "AFFINE"; break;
        case InputApproximationKind::SINUSOIDAL: os << "SINUSOIDAL"; break;
        case InputApproximationKind::PIECEWISE: os << "PIECEWISE"; break;
        default: ARIADNE_FAIL_MSG("Unhandled input approximation for output streaming\n");
    }
    return os;
}

template<InputApproximationKind A> struct InputApproximationKindTrait {
    static InputApproximationKind kind() { return A; }
};

class ZeroApproximation : public InputApproximationKindTrait<InputApproximationKind::ZERO> { };
class ConstantApproximation : public InputApproximationKindTrait<InputApproximationKind::CONSTANT> { };
class AffineApproximation : public InputApproximationKindTrait<InputApproximationKind::AFFINE> { };
class SinusoidalApproximation : public InputApproximationKindTrait<InputApproximationKind::SINUSOIDAL> { };
class PiecewiseApproximation : public InputApproximationKindTrait<InputApproximationKind::PIECEWISE> { };

template<class A> ErrorType r_value();
template<> ErrorType r_value<ZeroApproximation>() { return ErrorType(0u); }
template<> ErrorType r_value<ConstantApproximation>() { return ErrorType(1u); }
template<> ErrorType r_value<AffineApproximation>() { return ErrorType(5.0/3u); }
template<> ErrorType r_value<SinusoidalApproximation>() { return ErrorType(5.0/4u); }
template<> ErrorType r_value<PiecewiseApproximation>() { return ErrorType(1.3645_upper); }

template<class A> constexpr Nat num_params_per_input();
template<> constexpr Nat num_params_per_input<ZeroApproximation>() { return 0u; }
template<> constexpr Nat num_params_per_input<ConstantApproximation>() { return 1u; }
template<> constexpr Nat num_params_per_input<AffineApproximation>() { return 2u; }
template<> constexpr Nat num_params_per_input<SinusoidalApproximation>() { return 2u; }
template<> constexpr Nat num_params_per_input<PiecewiseApproximation>() { return 2u; }

enum class InputsRelationKind : std::uint8_t { AFFINE, SINGULAR, ADDITIVE};

template<InputsRelationKind R> struct InputsRelationKindTrait {
    static InputsRelationKind kind() { return R; }
};

class AffineInputs : public InputsRelationKindTrait<InputsRelationKind::AFFINE> { };
class SingularInput : public InputsRelationKindTrait<InputsRelationKind::SINGULAR> { };
class AdditiveInputs : public InputsRelationKindTrait<InputsRelationKind::ADDITIVE> { };

inline std::ostream& operator<<(std::ostream& os, const InputsRelationKind& kind) {
    switch (kind) {
        case InputsRelationKind::AFFINE: os << "AFFINE"; break;
        case InputsRelationKind::SINGULAR: os << "SINGULAR"; break;
        case InputsRelationKind::ADDITIVE: os << "ADDITIVE"; break;
        default: ARIADNE_FAIL_MSG("Unhandled InputsRoles for output streaming\n");
    }
    return os;
}

ErrorType zeroparam_worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h);
ErrorType zeroparam_component_error(C1Norms const& n, PositiveFloatDPValue const& h, SizeType j);
ErrorType oneparam_worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h);
ErrorType oneparam_component_error(C1Norms const& n, PositiveFloatDPValue const& h, SizeType j);
template<class R> ErrorType twoparam_worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r);
template<class R> ErrorType twoparam_component_error(C1Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r, SizeType j);

template<class A, class R, EnableIf<IsSame<A,ZeroApproximation>> = dummy> ErrorType worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h) { return zeroparam_worstcase_error(n,h); }
template<class A, class R, EnableIf<IsSame<A,ConstantApproximation>> = dummy> ErrorType worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h) { return oneparam_worstcase_error(n,h); }
template<class A, class R, EnableIf<Not<Or<IsSame<A,ZeroApproximation>,IsSame<A,ConstantApproximation>>>> = dummy> ErrorType worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h) { return twoparam_worstcase_error<R>(n,h,r_value<A>()); }
template<class A, class R, EnableIf<IsSame<A,ZeroApproximation>> = dummy> ErrorType component_error(C1Norms const& n, PositiveFloatDPValue const& h, SizeType j) { return zeroparam_component_error(n,h,j); }
template<class A, class R, EnableIf<IsSame<A,ConstantApproximation>> = dummy> ErrorType component_error(C1Norms const& n, PositiveFloatDPValue const& h, SizeType j) { return oneparam_component_error(n,h,j); }
template<class A, class R, EnableIf<Not<Or<IsSame<A,ZeroApproximation>,IsSame<A,ConstantApproximation>>>> = dummy> ErrorType component_error(C1Norms const& n, PositiveFloatDPValue const& h, SizeType j) { return twoparam_component_error<R>(n,h,r_value<A>(),j); }


class InputApproximator;
class InputApproximatorInterface;

class InputApproximatorFactory {
public:
    InputApproximator create(InclusionVectorField const& di, InputApproximationKind kind, SweeperDP sweeper) const;
};

template<class A> class ApproximationErrorProcessorInterface {
public:
    virtual Vector<ErrorType> process(PositiveFloatDPValue const& h, UpperBoxType const& B) const = 0;
};

template<class A, class R>
class ApproximationErrorProcessor : public ApproximationErrorProcessorInterface<A>, public Loggable {
  public:
    ApproximationErrorProcessor(InclusionVectorField const& ivf) : _ivf(ivf), _enable_componentwise_error(false) { }
    virtual Vector<ErrorType> process(PositiveFloatDPValue const& h, UpperBoxType const& B) const override;
  private:
    InclusionVectorField const& _ivf;
  protected:
    Boolean _enable_componentwise_error; // TODO: remove such option as soon as the DI paper is completed
  private:
    Vector<ErrorType> process(C1Norms const& n, PositiveFloatDPValue const& h) const;
  public:
    virtual ~ApproximationErrorProcessor() = default;
};


template<class A>
class ApproximationErrorProcessorFactory {
    typedef ApproximationErrorProcessorInterface<A> Processor;
public:
    SharedPointer<Processor> create(InclusionVectorField const& ivf) const {
        if (ivf.is_input_additive()) return SharedPointer<Processor>(new ApproximationErrorProcessor<A,AdditiveInputs>(ivf));
        else if (ivf.number_of_inputs() == 1) return SharedPointer<Processor>(new ApproximationErrorProcessor<A,SingularInput>(ivf));
        else return SharedPointer<Processor>(new ApproximationErrorProcessor<A,AffineInputs>(ivf));
    }
};

class InputApproximatorInterface {
  public:
    virtual InputApproximationKind kind() const = 0;
    virtual Vector<ErrorType> compute_errors(PositiveFloatDPValue h, UpperBoxType const& B) const = 0;
    virtual BoxDomainType build_flow_domain(BoxDomainType D, BoxDomainType V, PositiveFloatDPValue h) const = 0;
    virtual Vector<ValidatedScalarMultivariateFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const = 0;
};


class InputApproximator : public InputApproximatorInterface {
    friend class InputApproximatorFactory;
  private:
    SharedPointer<InputApproximatorInterface> _impl;
    InputApproximator(SharedPointer<InputApproximatorInterface> const& impl) : _impl(impl) { }
  public:
    InputApproximator(InputApproximator const& other) : _impl(other._impl) { }
    InputApproximator& operator=(InputApproximator const& other) { _impl = other._impl; return *this; }
  public:
    virtual InputApproximationKind kind() const override { return _impl->kind(); }
    virtual Vector<ErrorType> compute_errors(PositiveFloatDPValue h, UpperBoxType const& B) const override { return _impl->compute_errors(h,B); }
    virtual BoxDomainType build_flow_domain(BoxDomainType D, BoxDomainType V, PositiveFloatDPValue h) const override { return _impl->build_flow_domain(D,V,h); }
    virtual Vector<ValidatedScalarMultivariateFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const override { return _impl->build_w_functions(DVh,n,m); }
    virtual ~InputApproximator() = default;
};


template<class A>
class InputApproximatorBase : public InputApproximatorInterface {
    friend class InputApproximatorFactory;
  protected:
    InclusionVectorField const& _ivf;
    SweeperDP _sweeper;
    SharedPointer<ApproximationErrorProcessorInterface<A>> _processor;
    InputApproximatorBase(InclusionVectorField const& ivf, SweeperDP const& sweeper) :
        _ivf(ivf), _sweeper(sweeper), _processor(ApproximationErrorProcessorFactory<A>().create(ivf)), _kind(A::kind()), _num_params_per_input(num_params_per_input<A>()) { }
  private:
    const InputApproximationKind _kind;
    const Nat _num_params_per_input;
  public:
    virtual InputApproximationKind kind() const override { return _kind; }
    virtual Vector<ErrorType> compute_errors(PositiveFloatDPValue h, UpperBoxType const& B) const override { return _processor->process(h,B); }
    virtual BoxDomainType build_flow_domain(BoxDomainType D, BoxDomainType V, PositiveFloatDPValue h) const override;
    virtual Vector<ValidatedScalarMultivariateFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const override;
    virtual ~InputApproximatorBase() = default;
};


class Reconditioner {
  public:
    virtual Void simplify(ValidatedVectorMultivariateFunctionModelType& phi) const = 0;
    virtual ValidatedVectorMultivariateFunctionModelType expand_errors(ValidatedVectorMultivariateFunctionModelType Phi) const = 0;
};


class LohnerReconditioner : public Reconditioner, public Loggable {
    SweeperDP _sweeper;
    Nat _number_of_variables_to_keep;
public:
    LohnerReconditioner(SweeperDP sweeper, Nat number_of_variables_to_keep);
    void set_sweeper(SweeperDP sweeper) { _sweeper = sweeper; }
    void set_number_of_variables_to_keep(Nat num_variables_to_keep) { _number_of_variables_to_keep = num_variables_to_keep; }
    virtual ValidatedVectorMultivariateFunctionModelType expand_errors(ValidatedVectorMultivariateFunctionModelType f) const override;
    virtual Void simplify(ValidatedVectorMultivariateFunctionModelType& f) const override;
    virtual ~LohnerReconditioner() = default;
};


class InclusionIntegratorInterface {
  public:
    virtual List<ValidatedVectorMultivariateFunctionModelType> flow(InclusionVectorField const& ivf, BoxDomainType const& initial, Real T) = 0;
    virtual Pair<StepSizeType,UpperBoxType> flow_bounds(ValidatedVectorMultivariateFunction f, BoxDomainType dom, StepSizeType hsug) const = 0;
    virtual ValidatedVectorMultivariateFunctionModelType reach(InclusionVectorField const& ivf, BoxDomainType D, ValidatedVectorMultivariateFunctionModelType evolve_function, UpperBoxType B, PositiveFloatDPValue t, PositiveFloatDPValue h) const = 0;
};

class InclusionIntegrator : public virtual InclusionIntegratorInterface, public Loggable {
  protected:
    List<InputApproximationKind> _approximations;
    SharedPointer<InputApproximator> _approximator;
    SharedPointer<Reconditioner> _reconditioner;
    SweeperDP _sweeper;
    FloatDP _step_size;
    Nat _number_of_steps_between_simplifications;
    Nat _number_of_variables_to_keep;
  public:
    InclusionIntegrator(List<InputApproximationKind> approximations, SweeperDP sweeper, StepSize step_size);
    template<class... AS> InclusionIntegrator(List<InputApproximationKind> approximations, SweeperDP sweeper, StepSize step_size, AS... attributes);
  public:
    InclusionIntegrator& set(NumberOfStepsBetweenSimplifications n) { _number_of_steps_between_simplifications=n; return *this; }
    InclusionIntegrator& set(NumberOfVariablesToKeep n) { _number_of_variables_to_keep=n; return *this; }
    template<class A, class... AS> InclusionIntegrator& set(A a, AS... as) { this->set(a); this->set(as...); return *this; }

    virtual List<ValidatedVectorMultivariateFunctionModelType> flow(InclusionVectorField const& ivf, BoxDomainType const& initial, Real T) override;

    virtual Pair<StepSizeType,UpperBoxType> flow_bounds(ValidatedVectorMultivariateFunction f, BoxDomainType dom, StepSizeType hsug) const override;
    virtual ValidatedVectorMultivariateFunctionModelType reach(InclusionVectorField const& ivf, BoxDomainType D, ValidatedVectorMultivariateFunctionModelType evolve_function, UpperBoxType B, PositiveFloatDPValue t, PositiveFloatDPValue h) const override;
  private:
    ValidatedVectorMultivariateFunctionModelType compute_flow_function(ValidatedVectorMultivariateFunction const& dyn, BoxDomainType const& domain, UpperBoxType const& B) const;
    ValidatedVectorMultivariateFunctionModelDP build_reach_function(ValidatedVectorMultivariateFunctionModelDP evolve_function, ValidatedVectorMultivariateFunctionModelDP Phi, PositiveFloatDPValue t, PositiveFloatDPValue new_t) const;
    ValidatedVectorMultivariateFunctionModelDP evaluate_evolve_function(ValidatedVectorMultivariateFunctionModelDP reach_function, PositiveFloatDPValue t) const;
    ValidatedVectorMultivariateFunctionModelDP build_secondhalf_piecewise_reach_function(ValidatedVectorMultivariateFunctionModelDP evolve_function, ValidatedVectorMultivariateFunctionModelDP Phi, PositiveFloatDPValue t, PositiveFloatDPValue new_t) const;
    Vector<ValidatedScalarMultivariateFunction> build_secondhalf_piecewise_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const;
  private:
    Bool must_recondition(Nat step) const;
};

template<class... AS> InclusionIntegrator::InclusionIntegrator(List<InputApproximationKind> approximations, SweeperDP sweeper, StepSize step_size_, AS... attributes)
                : InclusionIntegrator::InclusionIntegrator(approximations, sweeper,step_size_) {
    this->set(attributes...);
    _reconditioner.reset(new LohnerReconditioner(_sweeper,_number_of_variables_to_keep));
}

BoxDomainType build_flow_domain(BoxDomainType D, BoxDomainType V, PositiveFloatDPValue h, Nat num_params);


} // namespace Ariadne;

#endif // ARIADNE_DIFFERENTIAL_INCLUSION_HPP
