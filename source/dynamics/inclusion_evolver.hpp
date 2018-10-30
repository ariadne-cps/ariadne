/***************************************************************************
 *            inclusion_evolver.hpp
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

/*! \file inclusion_evolver.hpp
 *  \brief Evolver for differential inclusion dynamics.
 */

#ifndef ARIADNE_INCLUSION_EVOLVER_HPP
#define ARIADNE_INCLUSION_EVOLVER_HPP

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

struct StepSize : public Attribute<StepSizeType> { };
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

using TimeStepType = Dyadic;

BoxDomainType initial_ranges_to_box(RealVariablesBox const& var_ranges);

inline Vector<FloatDPValue> const& cast_exact(Vector<FloatDPError> const& v) {
    return reinterpret_cast<Vector<FloatDPValue>const&>(v); }

FloatDP volume(Vector<ApproximateIntervalType> const& box);

inline Bool refines(Vector<UpperIntervalType> const& v1, UpperBoxType const& bx2) {
    return refines(v1,static_cast<Vector<UpperIntervalType>const&>(bx2)); }

Box<Interval<FloatDPValue>> over_approximation(Box<Interval<Real>> const&);

Void add_errors(ValidatedVectorMultivariateFunctionModelDP& phi, Vector<ErrorType> const& e);

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


class InclusionApproximatorHandle;
class InclusionApproximatorInterface;

class InclusionApproximatorFactory {
public:
    InclusionApproximatorHandle create(InclusionVectorField const& di, InputApproximationKind const& kind, SweeperDP const& sweeper) const;
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
    Boolean _enable_componentwise_error; // TODO: remove such option as soon as a paper presenting component-wise error is published
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

class InclusionApproximatorInterface {
  public:
    virtual InputApproximationKind kind() const = 0;
    virtual ValidatedVectorMultivariateFunctionModelType reach(InclusionVectorField const& ivf, BoxDomainType const& D, ValidatedVectorMultivariateFunctionModelType const& evolve_function, UpperBoxType const& B, TimeStepType const& t, StepSizeType const& h) const = 0;
    virtual ValidatedVectorMultivariateFunctionModelDP evolve(ValidatedVectorMultivariateFunctionModelDP const& reach_function, TimeStepType const& t) const = 0;
    virtual Pair<StepSizeType,UpperBoxType> flow_bounds(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& domx, BoxDomainType const& doma, StepSizeType const& hsug) const = 0;
};


class InclusionApproximatorHandle {
    friend class InclusionApproximatorFactory;
  private:
    SharedPointer<InclusionApproximatorInterface> _impl;
    InclusionApproximatorHandle(SharedPointer<InclusionApproximatorInterface> const& impl) : _impl(impl) { }
  public:
    InclusionApproximatorHandle(InclusionApproximatorHandle const& other) : _impl(other._impl) { }
    InclusionApproximatorHandle& operator=(InclusionApproximatorHandle const& other) { _impl = other._impl; return *this; }

    operator InclusionApproximatorInterface const& () const { return *_impl; }
    InclusionApproximatorHandle* clone() const { return new InclusionApproximatorHandle(*this); }

    InputApproximationKind kind() const { return _impl->kind(); }
    Pair<StepSizeType,UpperBoxType> flow_bounds(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& domx, BoxDomainType const& doma, StepSizeType const& hsug) const { return _impl->flow_bounds(f,domx,doma,hsug); }
    ValidatedVectorMultivariateFunctionModelType reach(InclusionVectorField const& ivf, BoxDomainType const& D, ValidatedVectorMultivariateFunctionModelType const& evolve_function, UpperBoxType const& B, TimeStepType const& t, StepSizeType const& h) const { return _impl->reach(ivf,D,evolve_function,B,t,h); }
    ValidatedVectorMultivariateFunctionModelDP evolve(ValidatedVectorMultivariateFunctionModelDP const& reach_function, TimeStepType const& t) const { return _impl->evolve(reach_function,t); }
  public:
    virtual ~InclusionApproximatorHandle() = default;
};


template<class A>
class InputApproximatorBase : public InclusionApproximatorInterface, Loggable {
    friend class InclusionApproximatorFactory;
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

    virtual ValidatedVectorMultivariateFunctionModelType reach(InclusionVectorField const& ivf, BoxDomainType const& D, ValidatedVectorMultivariateFunctionModelType const& evolve_function, UpperBoxType const& B, TimeStepType const& t, StepSizeType const& h) const override;
    virtual ValidatedVectorMultivariateFunctionModelDP evolve(ValidatedVectorMultivariateFunctionModelDP const& reach_function, TimeStepType const& t) const override;
    virtual Pair<StepSizeType,UpperBoxType> flow_bounds(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& domx, BoxDomainType const& doma, StepSizeType const& hsug) const override;

    virtual ~InputApproximatorBase() = default;
  private:
    Vector<ErrorType> compute_errors(StepSizeType const& h, UpperBoxType const& B) const { return _processor->process(PositiveFloatDPValue(h,DoublePrecision()),B); };
    BoxDomainType build_parameter_domain(BoxDomainType const& V) const;
    ValidatedVectorMultivariateFunctionModelType compute_flow_function(ValidatedVectorMultivariateFunction const& dyn, BoxDomainType const& domx, Interval<TimeStepType> const& domt, BoxDomainType const& doma, UpperBoxType const& B) const;
    ValidatedVectorMultivariateFunctionModelDP build_reach_function(ValidatedVectorMultivariateFunctionModelDP const& evolve_function, ValidatedVectorMultivariateFunctionModelDP const& Phi, TimeStepType const& t, TimeStepType const& new_t) const;
    ValidatedVectorMultivariateFunctionModelDP build_secondhalf_piecewise_reach_function(ValidatedVectorMultivariateFunctionModelDP const& evolve_function, ValidatedVectorMultivariateFunctionModelDP const& Phi, TimeStepType const& t, TimeStepType const& new_t) const;
    Vector<ValidatedScalarMultivariateFunction> build_w_functions(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) const;
    Vector<ValidatedScalarMultivariateFunction> build_secondhalf_piecewise_w_functions(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) const;
};


class Reconditioner {
  public:
    virtual Void simplify(ValidatedVectorMultivariateFunctionModelType& phi) const = 0;
    virtual ValidatedVectorMultivariateFunctionModelType expand_errors(ValidatedVectorMultivariateFunctionModelType const& Phi) const = 0;
};


class LohnerReconditioner : public Reconditioner, public Loggable {
    SweeperDP _sweeper;
    Nat _number_of_variables_to_keep;
public:
    LohnerReconditioner(SweeperDP sweeper, Nat number_of_variables_to_keep);
    void set_sweeper(SweeperDP sweeper) { _sweeper = sweeper; }
    void set_number_of_variables_to_keep(Nat num_variables_to_keep) { _number_of_variables_to_keep = num_variables_to_keep; }
    virtual ValidatedVectorMultivariateFunctionModelType expand_errors(ValidatedVectorMultivariateFunctionModelType const& f) const override;
    virtual Void simplify(ValidatedVectorMultivariateFunctionModelType& f) const override;
    virtual ~LohnerReconditioner() = default;
};

class InclusionEvolver : public Loggable {
  protected:
    List<InputApproximationKind> _approximations;
    SharedPointer<Reconditioner> _reconditioner;
    SweeperDP _sweeper;
    StepSizeType _step_size;
    Nat _number_of_steps_between_simplifications;
    Nat _number_of_variables_to_keep;
  public:
    InclusionEvolver(List<InputApproximationKind> approximations, SweeperDP sweeper, StepSizeType step_size);
    template<class... AS> InclusionEvolver(List<InputApproximationKind> approximations, SweeperDP sweeper, StepSizeType step_size_, AS... attributes)
                        : InclusionEvolver(approximations, sweeper,step_size_) {
            this->set(attributes...); _reconditioner.reset(new LohnerReconditioner(_sweeper,_number_of_variables_to_keep)); }
  public:
    InclusionEvolver& set(NumberOfStepsBetweenSimplifications n) { _number_of_steps_between_simplifications=n; return *this; }
    InclusionEvolver& set(NumberOfVariablesToKeep n) { _number_of_variables_to_keep=n; return *this; }
    template<class A, class... AS> InclusionEvolver& set(A a, AS... as) { this->set(a); this->set(as...); return *this; }

    List<ValidatedVectorMultivariateFunctionModelType> flow(InclusionVectorField const& ivf, BoxDomainType const& initial, Real const& T);

  private:
    Bool must_recondition(Nat step) const;
};

} // namespace Ariadne;

#endif // ARIADNE_INCLUSION_EVOLVER_HPP