/***************************************************************************
 *            inclusion_evolver.hpp
 *
 *  Copyright  2008-20  Luca Geretti, Pieter Collins, Sanja Zivanovic
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
#include "../numeric/numeric.hpp"
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

template<class C, class T> Bool instance_of(T* const obj) {
    return (dynamic_cast<const C*>(obj) != nullptr);
}

BoxDomainType initial_ranges_to_box(RealVariablesBox const& var_ranges);

inline Vector<FloatDPValue> const& cast_exact(Vector<FloatDPError> const& v) {
    return reinterpret_cast<Vector<FloatDPValue>const&>(v); }

inline Bool refines(Vector<UpperIntervalType> const& v1, UpperBoxType const& bx2) {
    return refines(v1,cast_vector(bx2)); }

Box<Interval<FloatDPValue>> over_approximation(Box<Interval<Real>> const&);

Void add_errors(ValidatedVectorMultivariateFunctionModelDP& phi, Vector<ErrorType> const& e);

ValidatedVectorMultivariateFunction build_Fw(ValidatedVectorMultivariateFunction const& F, Vector<ValidatedScalarMultivariateFunction> const& w);

template<class F1, class F2, class F3, class... FS> decltype(auto) combine(F1 const& f1, F2 const& f2, F3 const& f3, FS const& ... fs) {
    return combine(combine(f1,f2),f3,fs...); }
template<class F1, class F2, class F3, class... FS> decltype(auto) join(F1 const& f1, F2 const& f2, F3 const& f3, FS const& ... fs) {
    return join(join(f1,f2),f3,fs...); }
template<class F1, class F2, class F3, class... FS> decltype(auto) product(F1 const& f1, F2 const& f2, F3 const& f3, FS const& ... fs) {
    return product(product(f1,f2),f3,fs...); }

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

struct InputApproximationInterface {
    virtual Void write(OutputStream& os) const = 0;
    virtual InputApproximationInterface* clone() const = 0;
    virtual ~InputApproximationInterface() = default;
    virtual Nat index() const = 0;
};

inline std::ostream& operator << (std::ostream& os, const InputApproximationInterface& a) { a.write(os); return os; }

struct ZeroApproximation : public InputApproximationInterface {
    virtual Void write(OutputStream& os) const override { os << "ZERO"; }
    virtual Nat index() const override { return 0; }
    virtual ZeroApproximation* clone() const override { return new ZeroApproximation(*this); }
};
struct ConstantApproximation : public InputApproximationInterface {
    virtual Void write(OutputStream& os) const override { os << "CONSTANT"; }
    virtual Nat index() const override { return 1; }
    virtual ConstantApproximation* clone() const override { return new ConstantApproximation(*this); }
};
struct AffineApproximation : public InputApproximationInterface {
    virtual Void write(OutputStream& os) const override { os << "AFFINE"; }
    virtual Nat index() const override { return 2; }
    virtual AffineApproximation* clone() const override { return new AffineApproximation(*this); }
};
struct SinusoidalApproximation : public InputApproximationInterface {
    virtual Void write(OutputStream& os) const override { os << "SINUSOIDAL"; }
    virtual Nat index() const override { return 3; }
    virtual SinusoidalApproximation* clone() const override { return new SinusoidalApproximation(*this); }
};
struct PiecewiseApproximation : public InputApproximationInterface {
    virtual Void write(OutputStream& os) const override { os << "PIECEWISE"; }
    virtual Nat index() const override { return 4; }
    virtual PiecewiseApproximation* clone() const override { return new PiecewiseApproximation(*this); }
};

class InputApproximation {
  private:
    SharedPointer<InputApproximationInterface> _impl;
  public:
    InputApproximation(InputApproximationInterface const& impl) : _impl(impl.clone()) { }
    InputApproximation(InputApproximation const& other) : _impl(other._impl) { }

    InputApproximation& operator=(InputApproximation const& other) { _impl = other._impl; return *this; }

    template<class A> Bool handles(A const& a) const { return instance_of<A>(&*_impl); }

    operator InputApproximationInterface const& () const { return *_impl; }
    InputApproximation* clone() const { return new InputApproximation(*this); }

    Void write(OutputStream& os) const { os << *_impl; }

    virtual ~InputApproximation() = default;
};

template<class A> ErrorType r_value();
template<> ErrorType r_value<AffineApproximation>() { return ErrorType(5.0/3u); }
template<> ErrorType r_value<SinusoidalApproximation>() { return ErrorType(5.0/4u); }
template<> ErrorType r_value<PiecewiseApproximation>() { return ErrorType(1.3645_upper); }

template<class A> constexpr Nat const_num_params_per_input();
template<> constexpr Nat const_num_params_per_input<ZeroApproximation>() { return 0u; }
template<> constexpr Nat const_num_params_per_input<ConstantApproximation>() { return 1u; }
template<> constexpr Nat const_num_params_per_input<AffineApproximation>() { return 2u; }
template<> constexpr Nat const_num_params_per_input<SinusoidalApproximation>() { return 2u; }
template<> constexpr Nat const_num_params_per_input<PiecewiseApproximation>() { return 2u; }

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

FloatDPError zeroparam_worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h);
FloatDPError zeroparam_component_error(C1Norms const& n, PositiveFloatDPValue const& h, SizeType j);
FloatDPError oneparam_worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h);
FloatDPError oneparam_component_error(C1Norms const& n, PositiveFloatDPValue const& h, SizeType j);
template<class R> FloatDPError twoparam_worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h, FloatDPError const& r);
template<class R> FloatDPError twoparam_component_error(C1Norms const& n, PositiveFloatDPValue const& h, FloatDPError const& r, SizeType j);

template<class A, class R, EnableIf<IsSame<A,ZeroApproximation>> = dummy> FloatDPError worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h) { return zeroparam_worstcase_error(n,h); }
template<class A, class R, EnableIf<IsSame<A,ConstantApproximation>> = dummy> FloatDPError worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h) { return oneparam_worstcase_error(n,h); }
template<class A, class R, EnableIf<Not<Or<IsSame<A,ZeroApproximation>,IsSame<A,ConstantApproximation>>>> = dummy> FloatDPError worstcase_error(C1Norms const& n, PositiveFloatDPValue const& h) { return twoparam_worstcase_error<R>(n,h,r_value<A>()); }
template<class A, class R, EnableIf<IsSame<A,ZeroApproximation>> = dummy> FloatDPError component_error(C1Norms const& n, PositiveFloatDPValue const& h, SizeType j) { return zeroparam_component_error(n,h,j); }
template<class A, class R, EnableIf<IsSame<A,ConstantApproximation>> = dummy> FloatDPError component_error(C1Norms const& n, PositiveFloatDPValue const& h, SizeType j) { return oneparam_component_error(n,h,j); }
template<class A, class R, EnableIf<Not<Or<IsSame<A,ZeroApproximation>,IsSame<A,ConstantApproximation>>>> = dummy> FloatDPError component_error(C1Norms const& n, PositiveFloatDPValue const& h, SizeType j) { return twoparam_component_error<R>(n,h,r_value<A>(),j); }


class InclusionApproximatorHandle;
class InclusionApproximatorInterface;

class InclusionApproximatorFactory {
public:
    InclusionApproximatorHandle create(InclusionVectorField const& di, InputApproximation const& approximation, SweeperDP const& sweeper) const;
};

template<class A> class ApproximationErrorProcessorInterface {
public:
    virtual Vector<FloatDPError> process(PositiveFloatDPValue const& h, UpperBoxType const& B) const = 0;
};

template<class A, class R>
class ApproximationErrorProcessor : public ApproximationErrorProcessorInterface<A>, public Loggable {
  public:
    ApproximationErrorProcessor(InclusionVectorField const& ivf) : _ivf(ivf), _enable_componentwise_error(false) { }
    virtual Vector<FloatDPError> process(PositiveFloatDPValue const& h, UpperBoxType const& B) const override;
  private:
    InclusionVectorField const& _ivf;
  protected:
    Boolean _enable_componentwise_error; // TODO: remove such option as soon as a paper presenting component-wise error is published
  private:
    Vector<FloatDPError> process(C1Norms const& n, PositiveFloatDPValue const& h) const;
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
    virtual Void write(OutputStream& os) const = 0;
    virtual Bool operator==(const InclusionApproximatorInterface& rhs) const = 0;
    virtual Bool operator<(const InclusionApproximatorInterface& rhs) const = 0;
    virtual Nat index() const = 0;
    virtual Nat num_params_per_input() const = 0;
    virtual ValidatedVectorMultivariateFunctionModelType reach(InclusionVectorField const& ivf, BoxDomainType const& D, ValidatedVectorMultivariateFunctionModelType const& evolve_function, UpperBoxType const& B, TimeStepType const& t, StepSizeType const& h) const = 0;
    virtual ValidatedVectorMultivariateFunctionModelDP evolve(ValidatedVectorMultivariateFunctionModelDP const& reach_function, TimeStepType const& t) const = 0;
    virtual Pair<StepSizeType,UpperBoxType> flow_bounds(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& domx, BoxDomainType const& doma, StepSizeType const& hsug) const = 0;

    friend std::ostream& operator<<(std::ostream& os, const InclusionApproximatorInterface& approximator) { approximator.write(os); return os; }
};


class InclusionApproximatorHandle {
    friend class InclusionApproximatorFactory;
  private:
    SharedPointer<InclusionApproximatorInterface> _impl;
    InclusionApproximatorHandle(SharedPointer<InclusionApproximatorInterface> const& other) : _impl(other) { }
  public:
    InclusionApproximatorHandle(InclusionApproximatorHandle const& other) : _impl(other._impl) { }
    InclusionApproximatorHandle& operator=(InclusionApproximatorHandle const& other) { _impl = other._impl; return *this; }

    Bool operator==(const InclusionApproximatorInterface& rhs) const { return *_impl == rhs; }
    Bool operator<(const InclusionApproximatorHandle& rhs) const { return *_impl < *rhs._impl; }

    template<class A> Bool handles(A const& a) const { return instance_of<A>(&*_impl); }

    virtual Nat num_params_per_input() const { return _impl->num_params_per_input(); }

    operator InclusionApproximatorInterface const& () const { return *_impl; }
    InclusionApproximatorHandle* clone() const { return new InclusionApproximatorHandle(*this); }

    friend std::ostream& operator<<(std::ostream& os, const InclusionApproximatorHandle& approximator) { os << *approximator._impl; return os; }

    Pair<StepSizeType,UpperBoxType> flow_bounds(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& domx, BoxDomainType const& doma, StepSizeType const& hsug) const { return _impl->flow_bounds(f,domx,doma,hsug); }
    ValidatedVectorMultivariateFunctionModelType reach(InclusionVectorField const& ivf, BoxDomainType const& D, ValidatedVectorMultivariateFunctionModelType const& evolve_function, UpperBoxType const& B, TimeStepType const& t, StepSizeType const& h) const { return _impl->reach(ivf,D,evolve_function,B,t,h); }
    ValidatedVectorMultivariateFunctionModelDP evolve(ValidatedVectorMultivariateFunctionModelDP const& reach_function, TimeStepType const& t) const { return _impl->evolve(reach_function,t); }
  public:
    virtual ~InclusionApproximatorHandle() = default;
};


template<class A>
class InclusionApproximatorBase : public InclusionApproximatorInterface, Loggable {
    friend class InclusionApproximatorFactory;
  protected:
    InclusionVectorField const& _ivf;
    SweeperDP _sweeper;
    SharedPointer<ApproximationErrorProcessorInterface<A>> _processor;
    InclusionApproximatorBase(InclusionVectorField const& ivf, SweeperDP const& sweeper) :
        _ivf(ivf), _sweeper(sweeper), _processor(ApproximationErrorProcessorFactory<A>().create(ivf)), _num_params_per_input(const_num_params_per_input<A>()) { }
  private:
    const Nat _num_params_per_input;
  public:
    virtual Void write(OutputStream& os) const override { os << A(); }
    virtual ValidatedVectorMultivariateFunctionModelType reach(InclusionVectorField const& ivf, BoxDomainType const& D, ValidatedVectorMultivariateFunctionModelType const& evolve_function, UpperBoxType const& B, TimeStepType const& t, StepSizeType const& h) const override;
    virtual ValidatedVectorMultivariateFunctionModelDP evolve(ValidatedVectorMultivariateFunctionModelDP const& reach_function, TimeStepType const& t) const override;
    virtual Pair<StepSizeType,UpperBoxType> flow_bounds(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& domx, BoxDomainType const& doma, StepSizeType const& hsug) const override;

    virtual Bool operator==(const InclusionApproximatorInterface& rhs) const override;
    virtual Bool operator<(const InclusionApproximatorInterface& rhs) const override;

    virtual Nat index() const override { return A().index(); }

    virtual Nat num_params_per_input() const { return _num_params_per_input; }

    virtual ~InclusionApproximatorBase() = default;
  private:
    Vector<ErrorType> compute_errors(StepSizeType const& h, UpperBoxType const& B) const { return _processor->process(PositiveFloatDPValue(h,DoublePrecision()),B); };
    BoxDomainType build_parameter_domain(BoxDomainType const& V) const;
    ValidatedVectorMultivariateFunctionModelType compute_flow_function(ValidatedVectorMultivariateFunction const& dyn, BoxDomainType const& domx, Interval<TimeStepType> const& domt, BoxDomainType const& doma, UpperBoxType const& B) const;
    ValidatedVectorMultivariateFunctionModelDP build_reach_function(ValidatedVectorMultivariateFunctionModelDP const& evolve_function, ValidatedVectorMultivariateFunctionModelDP const& Phi, TimeStepType const& t, TimeStepType const& new_t) const;
    ValidatedVectorMultivariateFunctionModelDP build_secondhalf_piecewise_reach_function(ValidatedVectorMultivariateFunctionModelDP const& evolve_function, ValidatedVectorMultivariateFunctionModelDP const& Phi, TimeStepType const& t, TimeStepType const& new_t) const;
    Vector<ValidatedScalarMultivariateFunction> build_w_functions(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) const;
    Vector<ValidatedScalarMultivariateFunction> build_secondhalf_piecewise_w_functions(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) const;
};


class ReconditionerInterface {
  public:
    virtual Void simplify(ValidatedVectorMultivariateFunctionModelType& phi) const = 0;
    virtual ValidatedVectorMultivariateFunctionModelType expand_errors(ValidatedVectorMultivariateFunctionModelType const& Phi) const = 0;
    virtual Void set_number_of_variables_to_keep(Nat num_variables_to_keep) = 0;
    virtual Nat number_of_steps_between_simplifications() const = 0;
    virtual Bool must_recondition(Nat step) const = 0;
    virtual ReconditionerInterface* clone() const = 0;
    virtual ~ReconditionerInterface() = default;
};


class LohnerReconditioner : public ReconditionerInterface, public Loggable {
    SweeperDP _sweeper;
    Nat _number_of_steps_between_simplifications;
    Nat _number_of_variables_to_keep;
public:
    LohnerReconditioner(SweeperDP sweeper, Nat number_of_steps_between_simplifications_, Nat number_of_variables_to_keep_)
        : _sweeper(sweeper), _number_of_steps_between_simplifications(number_of_steps_between_simplifications_), _number_of_variables_to_keep(number_of_variables_to_keep_) { }
    virtual LohnerReconditioner* clone() const override { return new LohnerReconditioner(*this); }
    Void set_sweeper(SweeperDP sweeper) { _sweeper = sweeper; }
    virtual Void set_number_of_variables_to_keep(Nat num_variables_to_keep) override { _number_of_variables_to_keep = num_variables_to_keep; }
    Void set_number_of_steps_between_simplifications(Nat number_of_steps_between_simplifications) { _number_of_steps_between_simplifications = number_of_steps_between_simplifications; }
    Nat number_of_steps_between_simplifications() const { return _number_of_steps_between_simplifications; }
    virtual ValidatedVectorMultivariateFunctionModelType expand_errors(ValidatedVectorMultivariateFunctionModelType const& f) const override;
    virtual Void simplify(ValidatedVectorMultivariateFunctionModelType& f) const override;
    virtual Bool must_recondition(Nat step) const override {
        return (step%_number_of_steps_between_simplifications == _number_of_steps_between_simplifications-1);
    }
};

class Reconditioner {
private:
    SharedPointer<ReconditionerInterface> _impl;
public:
    Reconditioner(ReconditionerInterface const& other) : _impl(other.clone()) { }
    Reconditioner(Reconditioner const& other) : _impl(other._impl) { }

    Reconditioner& operator=(Reconditioner const& other) { _impl = other._impl; return *this; }

    Void simplify(ValidatedVectorMultivariateFunctionModelType& phi) const { _impl->simplify(phi); }
    ValidatedVectorMultivariateFunctionModelType expand_errors(ValidatedVectorMultivariateFunctionModelType const& Phi) const { return _impl->expand_errors(Phi); }
    Void set_number_of_variables_to_keep(Nat num_variables_to_keep) { _impl->set_number_of_variables_to_keep(num_variables_to_keep); }
    Nat number_of_steps_between_simplifications() const { return _impl->number_of_steps_between_simplifications(); }
    Bool must_recondition(Nat step) const { return _impl->must_recondition(step); }

    virtual ~Reconditioner() = default;
};


class InclusionEvolver : public Loggable {
  protected:
    List<InputApproximation> _approximations;
    SweeperDP _sweeper;
    StepSizeType _step_size;
    Reconditioner _reconditioner;

  public:
    InclusionEvolver(List<InputApproximation> const& approximations, SweeperDP const& sweeper, StepSizeType const& step_size, Reconditioner const& reconditioner);
  public:
    List<ValidatedVectorMultivariateFunctionModelType> flow(InclusionVectorField const& ivf, BoxDomainType const& initial, Real const& T);
};

} // namespace Ariadne;

#endif // ARIADNE_INCLUSION_EVOLVER_HPP
