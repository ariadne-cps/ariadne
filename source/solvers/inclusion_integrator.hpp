/***************************************************************************
 *            inclusion_integrator.hpp
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

#ifndef ARIADNE_INCLUSION_INTEGRATOR_HPP
#define ARIADNE_INCLUSION_INTEGRATOR_HPP

#include "../utility/typedefs.hpp"
#include "../utility/attribute.hpp"
#include "../algebra/sweeper.hpp"
#include "../algebra/algebra.hpp"
#include "../function/domain.hpp"
#include "../function/function_model.hpp"
#include "../function/formula.hpp"
#include "../function/symbolic_function.hpp"
#include "../symbolic/expression_set.hpp"
#include "../output/logging.hpp"
#include "../solvers/integrator_interface.hpp"

namespace Ariadne {

class Real;

using ValidatedScalarMultivariateFunctionModelType = ValidatedScalarMultivariateFunctionModelDP;
using ValidatedVectorMultivariateFunctionModelType = ValidatedVectorMultivariateFunctionModelDP;

using TimeStepType = Dyadic;

template<class C, class T> Bool instance_of(T* const obj) {
    return (dynamic_cast<const C*>(obj) != nullptr);
}

inline TimeStepType lower_bound(TimeStepType const& t) {
    return TimeStepType(FloatDPLowerBound(t,DoublePrecision()).raw()); }
inline TimeStepType lower_bound(Real const& t) {
    return TimeStepType(FloatDPLowerBound(t,DoublePrecision()).raw()); }

inline Vector<FloatDPValue> const& cast_exact(Vector<FloatDPError> const& v) {
    return reinterpret_cast<Vector<FloatDPValue>const&>(v); }

Void add_errors(ValidatedVectorMultivariateFunctionModelDP& phi, Vector<ErrorType> const& e);

EffectiveVectorMultivariateFunction substitute_v_with_w(EffectiveVectorMultivariateFunction const& F, Vector<EffectiveScalarMultivariateFunction> const& w);

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

C1Norms compute_norms(EffectiveVectorMultivariateFunction const&, Vector<EffectiveVectorMultivariateFunction> const&, BoxDomainType const&, PositiveFloatDPValue const&, UpperBoxType const&);

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

template<class A> Vector<EffectiveScalarMultivariateFunction> build_w_functions(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m);

class InclusionIntegratorHandle;
class InclusionIntegratorInterface;

class InclusionIntegratorFactory {
    SharedPointer<IntegratorInterface> _integrator;
public:
    InclusionIntegratorFactory(IntegratorInterface const& integrator) : _integrator(integrator.clone()) { }
    InclusionIntegratorHandle create(EffectiveVectorMultivariateFunction const& f, BoxDomainType const& inputs, InputApproximation const& approximation) const;
    InclusionIntegratorFactory* clone() const { return new InclusionIntegratorFactory(*this); }
};

template<class A> class ApproximationErrorProcessorInterface {
public:
    virtual Vector<ErrorType> process(PositiveFloatDPValue const& h, UpperBoxType const& B) const = 0;
};

template<class A, class R>
class ApproximationErrorProcessor : public ApproximationErrorProcessorInterface<A>, public Loggable {
  public:
    ApproximationErrorProcessor(EffectiveVectorMultivariateFunction const& f, BoxDomainType const& inputs) : _f(f), _inputs(inputs), _enable_componentwise_error(true) { }
    virtual Vector<ErrorType> process(PositiveFloatDPValue const& h, UpperBoxType const& B) const override;
  private:
    EffectiveVectorMultivariateFunction const& _f;
    BoxDomainType const& _inputs;
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
    SharedPointer<Processor> create(EffectiveVectorMultivariateFunction const& f, BoxDomainType const& inputs) const {
        Nat n = f.result_size();
        Nat m = inputs.size();
        ARIADNE_ASSERT_MSG(f.argument_size()-n == m, "ApproximationErrorProcessorFactory was given an incompatible f argument space in respect to the inputs box");
        Set<Nat> input_idx;
        for (Nat i : range(n,n+m)) { input_idx.insert(i); }
        if (is_additive_in(f,input_idx)) return SharedPointer<Processor>(new ApproximationErrorProcessor<A,AdditiveInputs>(f,inputs));
        else if (m == 1) return SharedPointer<Processor>(new ApproximationErrorProcessor<A,SingularInput>(f,inputs));
        else return SharedPointer<Processor>(new ApproximationErrorProcessor<A,AffineInputs>(f,inputs));
    }
};

class InclusionIntegratorInterface {
  public:
    virtual Void write(OutputStream& os) const = 0;
    virtual Bool operator==(const InclusionIntegratorInterface& rhs) const = 0;
    virtual Bool operator<(const InclusionIntegratorInterface& rhs) const = 0;
    virtual Nat index() const = 0;
    virtual Nat num_params_per_input() const = 0;
    virtual List<ValidatedVectorMultivariateFunctionModelType> reach(BoxDomainType const& D, ValidatedVectorMultivariateFunctionModelType const& evolve_function, UpperBoxType const& B, TimeStepType const& t, StepSizeType const& h) const = 0;
    virtual ValidatedVectorMultivariateFunctionModelDP evolve(ValidatedVectorMultivariateFunctionModelDP const& reach_function, TimeStepType const& t) const = 0;
    virtual Pair<StepSizeType,UpperBoxType> flow_bounds(BoxDomainType const& domx, BoxDomainType const& doma, StepSizeType const& hsug) const = 0;

    friend std::ostream& operator<<(std::ostream& os, const InclusionIntegratorInterface& approximator) { approximator.write(os); return os; }
};


template<class A>
class InclusionIntegrator : public InclusionIntegratorInterface, Loggable {
    friend class InclusionIntegratorFactory;
  protected:
    EffectiveVectorMultivariateFunction const& _f;
    BoxDomainType const& _inputs;
    SharedPointer<IntegratorInterface> _integrator;
    SharedPointer<ApproximationErrorProcessorInterface<A>> _processor;
    InclusionIntegrator(EffectiveVectorMultivariateFunction const& f, BoxDomainType const& inputs, SharedPointer<IntegratorInterface> const& integrator) :
        _f(f), _inputs(inputs), _integrator(integrator), _processor(ApproximationErrorProcessorFactory<A>().create(f,inputs)), _num_params_per_input(const_num_params_per_input<A>()) { }
  private:
    const Nat _num_params_per_input;
  public:
    virtual Void write(OutputStream& os) const override { os << A(); }
    virtual List<ValidatedVectorMultivariateFunctionModelType> reach(BoxDomainType const& D, ValidatedVectorMultivariateFunctionModelType const& evolve_function, UpperBoxType const& B, TimeStepType const& t, StepSizeType const& h) const override;
    virtual ValidatedVectorMultivariateFunctionModelDP evolve(ValidatedVectorMultivariateFunctionModelDP const& reach_function, TimeStepType const& t) const override;
    virtual Pair<StepSizeType,UpperBoxType> flow_bounds(BoxDomainType const& domx, BoxDomainType const& doma, StepSizeType const& hsug) const override;

    virtual Bool operator==(const InclusionIntegratorInterface& rhs) const override;
    virtual Bool operator<(const InclusionIntegratorInterface& rhs) const override;

    virtual Nat index() const override { return A().index(); }

    virtual Nat num_params_per_input() const override { return _num_params_per_input; }

    virtual ~InclusionIntegrator() = default;
  private:
    Vector<ErrorType> compute_errors(StepSizeType const& h, UpperBoxType const& B) const { return _processor->process(PositiveFloatDPValue(h,DoublePrecision()),B); };
    BoxDomainType build_parameter_domain(BoxDomainType const& V) const;
    ValidatedVectorMultivariateFunctionModelDP build_reach_function(ValidatedVectorMultivariateFunctionModelDP const& evolve_function, ValidatedVectorMultivariateFunctionModelDP const& Phi, TimeStepType const& t, TimeStepType const& new_t) const;
    ValidatedVectorMultivariateFunctionModelDP build_secondhalf_piecewise_reach_function(ValidatedVectorMultivariateFunctionModelDP const& evolve_function, ValidatedVectorMultivariateFunctionModelDP const& Phi, TimeStepType const& t, TimeStepType const& new_t) const;
    Vector<EffectiveScalarMultivariateFunction> build_secondhalf_piecewise_w_functions(Interval<TimeStepType> const& domt, BoxDomainType const& doma, SizeType n, SizeType m) const;
};


class InclusionIntegratorHandle {
    friend class InclusionIntegratorFactory;
  private:
    SharedPointer<InclusionIntegratorInterface> _impl;
    InclusionIntegratorHandle(SharedPointer<InclusionIntegratorInterface> const& other) : _impl(other) { }
  public:
    InclusionIntegratorHandle(InclusionIntegratorHandle const& other) : _impl(other._impl) { }
    InclusionIntegratorHandle& operator=(InclusionIntegratorHandle const& other) { _impl = other._impl; return *this; }

    Bool operator==(const InclusionIntegratorInterface& rhs) const { return *_impl == rhs; }
    Bool operator<(const InclusionIntegratorHandle& rhs) const { return *_impl < *rhs._impl; }

    template<class A> Bool handles(A const& a) const { return instance_of<A>(&*_impl); }

    virtual Nat num_params_per_input() const { return _impl->num_params_per_input(); }

    operator InclusionIntegratorInterface const& () const { return *_impl; }
    InclusionIntegratorHandle* clone() const { return new InclusionIntegratorHandle(*this); }

    friend std::ostream& operator<<(std::ostream& os, const InclusionIntegratorHandle& approximator) { os << *approximator._impl; return os; }

    Pair<StepSizeType,UpperBoxType> flow_bounds(BoxDomainType const& domx, BoxDomainType const& doma, StepSizeType const& hsug) const { return _impl->flow_bounds(domx,doma,hsug); }
    List<ValidatedVectorMultivariateFunctionModelType> reach(BoxDomainType const& D, ValidatedVectorMultivariateFunctionModelType const& evolve_function, UpperBoxType const& B, TimeStepType const& t, StepSizeType const& h) const { return _impl->reach(D,evolve_function,B,t,h); }
    ValidatedVectorMultivariateFunctionModelDP evolve(ValidatedVectorMultivariateFunctionModelDP const& reach_function, TimeStepType const& t) const { return _impl->evolve(reach_function,t); }
  public:
    virtual ~InclusionIntegratorHandle() = default;
};


} // namespace Ariadne;

#endif // ARIADNE_INCLUSION_INTEGRATOR_HPP
