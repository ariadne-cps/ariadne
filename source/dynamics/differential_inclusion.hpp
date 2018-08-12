/***************************************************************************
 *            differential_inclusion.hpp
 *
 *  Copyright  2008-17  Pieter Collins, Sanja Zivanovic
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file differential_inclusion.hpp
 *  \brief Methods for computing solutions of differential inclusions.
 */

#ifndef ARIADNE_DIFFERENTIAL_INCLUSION_HPP
#define ARIADNE_DIFFERENTIAL_INCLUSION_HPP

#include "utility/typedefs.hpp"
#include "utility/attribute.hpp"
#include "algebra/sweeper.hpp"
#include "algebra/algebra.hpp"
#include "function/domain.hpp"
#include "function/function_model.hpp"
#include "function/formula.hpp"
#include "expression/expression_set.hpp"

namespace Ariadne {

class Real;

struct StepSize : public Attribute<FloatDP> { };
struct NumberOfStepsBetweenSimplifications : public Attribute<Nat> { };
struct NumberOfVariablesToKeep : public Attribute<Nat> { };

Generator<StepSize> step_size;
Generator<NumberOfStepsBetweenSimplifications> number_of_steps_between_simplifications;
Generator<NumberOfVariablesToKeep> number_of_variables_to_keep;

using ValidatedScalarFunctionModelType = ValidatedScalarFunctionModelDP;
using ValidatedVectorFunctionModelType = ValidatedVectorFunctionModelDP;

using ThresholdSweeperDP = ThresholdSweeper<FloatDP>;
using GradedSweeperDP = GradedSweeper<FloatDP>;
using GradedThresholdSweeperDP = GradedThresholdSweeper<FloatDP>;
using SweeperDP = Sweeper<FloatDP>;

using ApproximateTimeStepType = PositiveFloatDPApproximation;
using ExactTimeStepType = PositiveFloatDPValue;

Pair<RealAssignment,RealInterval> centered_variable_transformation(RealVariable const& v, RealInterval const& bounds);
Pair<RealAssignments,RealVariablesBox> centered_variables_transformation(RealVariablesBox const& inputs);
Tuple<ValidatedVectorFunction,Vector<ValidatedVectorFunction>,BoxDomainType> expression_to_function(DottedRealAssignments const& dynamics, const RealVariablesBox& inputs);
BoxDomainType bounds_to_domain(RealVariablesBox const& var_box);

Vector<FloatDPValue> const& cast_exact(Vector<FloatDPError> const& v) {
    return reinterpret_cast<Vector<FloatDPValue>const&>(v); }

FloatDP volume(Vector<ApproximateIntervalType> const& box);

Bool refines(Vector<UpperIntervalType> const& v1, UpperBoxType const& bx2) {
    return refines(v1,static_cast<Vector<UpperIntervalType>const&>(bx2)); }

Box<Interval<FloatDPValue>> over_approximation(Box<Interval<Real>> const&);

ValidatedVectorFunctionModelDP add_errors(ValidatedVectorFunctionModelDP phi, Vector<ErrorType> const& e);

Boolean inputs_are_additive(Vector<ValidatedVectorFunction> const &g);

ValidatedVectorFunction construct_f_plus_gu(ValidatedVectorFunction const &f, Vector<ValidatedVectorFunction> const &g);

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

struct DifferentialInclusion {
public:
    friend class DifferentialInclusionIVP;
    const ValidatedVectorFunction f;
    const Vector<ValidatedVectorFunction> g;
    const BoxDomainType V;
private:
    DifferentialInclusion(Tuple<ValidatedVectorFunction,Vector<ValidatedVectorFunction>,BoxDomainType> const& components) : f(std::get<0>(components)), g(std::get<1>(components)), V(std::get<2>(components)) { }
};

struct DifferentialInclusionIVP {
public:
    const DifferentialInclusion di;
    const BoxDomainType X0;
private:
    DifferentialInclusionIVP(DifferentialInclusion const& di, BoxDomainType const& X0) : di(di), X0(X0) { }
public:
    DifferentialInclusionIVP(DottedRealAssignments const& dynamics, const RealVariablesBox& inputs, const RealVariablesBox& initial)
        : di(DifferentialInclusion(expression_to_function(dynamics,inputs))), X0(bounds_to_domain(initial)) { }
};

std::ostream& operator << (std::ostream& os, const DifferentialInclusionIVP& ivp) {
    os << "IVP: \nf: " << ivp.di.f << "\ng: " << ivp.di.g << "\nV: " << ivp.di.V << "\nX0: " << ivp.X0 << "\n";
    return os;
}

struct Norms {
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

    Norms(FloatDPError const&,Vector<FloatDPError> const&,FloatDPError const&,Vector<FloatDPError> const&,FloatDPError const&,Vector<FloatDPError> const&,FloatDPError const&,Vector<FloatDPError> const&,FloatDPError const&,Vector<FloatDPError> const&,FloatDPError const&,Vector<FloatDPError> const&,FloatDPError const&,FloatDPError const&);
    Tuple<FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,Vector<FloatDPError>,FloatDPError,FloatDPError> values() const;

    SizeType dimension() const { return _dimension; }
private:
    SizeType _dimension;
};

Norms compute_norms(DifferentialInclusion const&, PositiveFloatDPValue const&, UpperBoxType const&);

enum class InputApproximation { ZERO, CONSTANT, AFFINE, SINUSOIDAL, PIECEWISE };

std::ostream& operator << (std::ostream& os, const InputApproximation& kind) {
    switch (kind) {
        case InputApproximation::ZERO: os << "ZERO"; break;
        case InputApproximation::CONSTANT: os << "CONSTANT"; break;
        case InputApproximation::AFFINE: os << "AFFINE"; break;
        case InputApproximation::SINUSOIDAL: os << "SINUSOIDAL"; break;
        case InputApproximation::PIECEWISE: os << "PIECEWISE"; break;
    }
    return os;
}

class ParametricInputApproximation {
protected:
    ParametricInputApproximation(InputApproximation kind) : _kind(kind) { }
public:
    InputApproximation kind() { return _kind; }
private:
    InputApproximation _kind;
};

class ZeroApproximation : public ParametricInputApproximation {
public: ZeroApproximation() : ParametricInputApproximation(InputApproximation::ZERO) { }
};
class ConstantApproximation : public ParametricInputApproximation {
public: ConstantApproximation() : ParametricInputApproximation(InputApproximation::CONSTANT) { }
};
class AffineApproximation : public ParametricInputApproximation {
public: AffineApproximation() : ParametricInputApproximation(InputApproximation::AFFINE) { }
};
class SinusoidalApproximation : public ParametricInputApproximation {
public: SinusoidalApproximation() : ParametricInputApproximation(InputApproximation::SINUSOIDAL) { }
};
class PiecewiseApproximation : public ParametricInputApproximation {
public: PiecewiseApproximation() : ParametricInputApproximation(InputApproximation::PIECEWISE) { }
};

template<class A> ErrorType r_value();
template<> ErrorType r_value<ZeroApproximation>() { return ErrorType(0u); }
template<> ErrorType r_value<ConstantApproximation>() { return ErrorType(1u); }
template<> ErrorType r_value<AffineApproximation>() { return ErrorType(5.0/3u); }
template<> ErrorType r_value<SinusoidalApproximation>() { return ErrorType(5.0/4u); }
template<> ErrorType r_value<PiecewiseApproximation>() { return ErrorType(1.3645_upper); }

enum class InputsRoles { AFFINE, SINGULAR, ADDITIVE};

class InputsRole {
protected:
    InputsRole(InputsRoles kind) : _kind(kind) { }
public:
    InputsRoles kind() { return _kind; }
private:
    InputsRoles _kind;
};

class AffineInputs : public InputsRole {
public: AffineInputs() : InputsRole(InputsRoles::AFFINE) { }
};
class SingularInput : public InputsRole {
public: SingularInput() : InputsRole(InputsRoles::SINGULAR) { }
};
class AdditiveInputs : public InputsRole {
public: AdditiveInputs() : InputsRole(InputsRoles::ADDITIVE) { }
};

std::ostream& operator << (std::ostream& os, const InputsRoles& kind) {
    switch (kind) {
        case InputsRoles::AFFINE: os << "AFFINE"; break;
        case InputsRoles::SINGULAR: os << "SINGULAR"; break;
        case InputsRoles::ADDITIVE: os << "ADDITIVE"; break;
    }
    return os;
}

template<class R> ErrorType error(Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r);
template<class R> ErrorType component_error(Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r, SizeType j);
template<class A, class R> ErrorType error(Norms const& n, PositiveFloatDPValue const& h) {
    return error<R>(n,h,r_value<A>()); }
template<class A, class R> ErrorType component_error(Norms const& n, PositiveFloatDPValue const& h, SizeType j) {
    return component_error<R>(n,h,r_value<A>(),j); }

template<> ErrorType error<ZeroApproximation,AffineInputs>(Norms const& n, PositiveFloatDPValue const& h) {
    return min(n.pK*n.expLambda*h,(n.K*2u+n.pK)*h); }
template<> ErrorType error<ZeroApproximation,AdditiveInputs>(Norms const& n, PositiveFloatDPValue const& h) {
    return error<ZeroApproximation,AffineInputs>(n,h); }
template<> ErrorType error<ZeroApproximation,SingularInput>(Norms const& n, PositiveFloatDPValue const& h) {
    return error<ZeroApproximation,AffineInputs>(n,h); }
template<> ErrorType component_error<ZeroApproximation,AffineInputs>(Norms const& n, PositiveFloatDPValue const& h, SizeType j) {
    return min(n.pK*n.expL*h,(n.Kj[j]*2u+n.pKj[j])*h); }
template<> ErrorType component_error<ZeroApproximation,AdditiveInputs>(Norms const& n, PositiveFloatDPValue const& h, SizeType j) {
    return component_error<ZeroApproximation,AffineInputs>(n,h,j); }
template<> ErrorType component_error<ZeroApproximation,SingularInput>(Norms const& n, PositiveFloatDPValue const& h, SizeType j) {
    return component_error<ZeroApproximation,AffineInputs>(n,h,j); }

template<> ErrorType error<ConstantApproximation,AffineInputs>(Norms const& n, PositiveFloatDPValue const& h) {
    return pow(h,2u)*((n.K+n.pK)*n.pL/3u + n.pK*2u*(n.L+n.pL)*n.expLambda); }
template<> ErrorType error<ConstantApproximation,SingularInput>(Norms const& n, PositiveFloatDPValue const& h) {
    return error<ConstantApproximation,AffineInputs>(n,h); }
template<> ErrorType error<ConstantApproximation,AdditiveInputs>(Norms const& n, PositiveFloatDPValue const& h) {
    return error<ConstantApproximation,AffineInputs>(n,h); }
template<> ErrorType component_error<ConstantApproximation,AffineInputs>(Norms const& n, PositiveFloatDPValue const& h, SizeType j) {
    return n.pLj[j]*(n.K+n.pK)*pow(h,2u)/3u + ((n.Lj[j]+n.pLj[j])*2u*n.pK)*cast_positive(cast_exact((n.L*n.expL*h+1u-n.expL)/pow(n.L,2u))); }
template<> ErrorType component_error<ConstantApproximation,SingularInput>(Norms const& n, PositiveFloatDPValue const& h, SizeType j) {
    return component_error<ConstantApproximation,AffineInputs>(n,h,j); }
template<> ErrorType component_error<ConstantApproximation,AdditiveInputs>(Norms const& n, PositiveFloatDPValue const& h, SizeType j) {
    return component_error<ConstantApproximation,AffineInputs>(n,h,j); }

template<> ErrorType error<AffineInputs>(Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r) {
    return ((r*r+1u)*n.pL*n.pK + (r+1u)*h*n.pK*((n.pH*2u*r + n.H)*(n.K+r*n.pK)+n.L*n.L+(n.L*3u*r+n.pL*r*r*2u)*n.pL)*n.expLambda + (r+1u)/6u*h*(n.K+n.pK)*((n.H*n.pK+n.L*n.pL)*3u+(n.pH*n.K+n.L*n.pL)*4u))/cast_positive(+1u-h*n.L/2u-h*n.pL*r)*pow(h,2u)/4u; }
template<> ErrorType error<SingularInput>(Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r) {
    return ((r+1u)*n.pK*((n.pH*2u*r+n.H)*(n.K+r*n.pK)+pow(n.L,2)+(n.L*3u*r+pow(r,2)*2u*n.pL)*n.pL)*n.expLambda + (n.K+n.pK)/6u*((r+1u)*((n.H*n.pK+n.L*n.pL)*3u +(n.pH*n.K+n.L*n.pL)*4u) + (n.pH*n.pK+n.pL*n.pL)*8u*(r*r+1u)))*pow(h,3u)/4u/cast_positive(+1u-h*n.L/2u-h*n.pL*r); }
template<> ErrorType error<AdditiveInputs>(Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r) {
    return (n.H*(n.K+n.pK)/2u + (n.L*n.L+n.H*(n.K+r*n.pK))*n.expLambda)/cast_positive(+1u-h*n.L/2u)*(r+1u)*n.pK*pow(h,3u)/4u; }
template<> ErrorType component_error<AffineInputs>(Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r, SizeType j) {
    return ErrorType(0u); }
template<> ErrorType component_error<SingularInput>(Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r, SizeType j) {
    return ErrorType(0u); }
template<> ErrorType component_error<AdditiveInputs>(Norms const& n, PositiveFloatDPValue const& h, ErrorType const& r, SizeType j) {
    return ErrorType(0u); }

class InputApproximator;
class InputApproximatorInterface;

class InputApproximatorFactory {
public:
    InputApproximator create(DifferentialInclusion const& di, InputApproximation kind, SweeperDP sweeper) const;
};

class ApproximationErrorProcessorInterface {
public:
    virtual Vector<ErrorType> process(PositiveFloatDPValue const& h, UpperBoxType const& B) const = 0;
};

template<class A, class R>
class ApproximationErrorProcessor : public ApproximationErrorProcessorInterface {
  public:
    ApproximationErrorProcessor(DifferentialInclusion const& di) : _di(di), _enable_componentwise_error(false) { }
    virtual Vector<ErrorType> process(PositiveFloatDPValue const& h, UpperBoxType const& B) const override;
  protected:
    Boolean _enable_componentwise_error; // TODO: remove such option as soon as the DI paper is completed
  private:
    Vector<ErrorType> process(Norms const& n, PositiveFloatDPValue const& h) const;
    DifferentialInclusion const& _di;
};

template<class A, class R> Vector<ErrorType> ApproximationErrorProcessor<A,R>::process(Norms const& n, PositiveFloatDPValue const& h) const {
    Vector<ErrorType> result(n.dimension(),error<A,R>(n,h));
    if (_enable_componentwise_error)
        for (auto j: range(n.dimension()))
            result[j] = min(result[j],component_error<A,R>(n,h,j));
    return result;
}

template<class A, class R> Vector<ErrorType> ApproximationErrorProcessor<A,R>::process(PositiveFloatDPValue const& h, UpperBoxType const& B) const {
    Norms norms = compute_norms(_di,h,B);
    if (inputs_are_additive(_di.g))
        norms.pK=mag(norm(_di.V));
    return process(norms,h);
}

template<class A>
class ApproximationErrorProcessorFactory {
    typedef ApproximationErrorProcessorInterface Processor;
public:
    SharedPointer<Processor> create(DifferentialInclusion const& di) const {
        if (inputs_are_additive(di.g)) return SharedPointer<Processor>(new ApproximationErrorProcessor<A,AdditiveInputs>(di));
        else if (di.g.size() == 1) return SharedPointer<Processor>(new ApproximationErrorProcessor<A,SingularInput>(di));
        else return SharedPointer<Processor>(new ApproximationErrorProcessor<A,AffineInputs>(di));
    }
};

class InputApproximatorInterface {
  public:
    virtual InputApproximation kind() const = 0;
    virtual Vector<ErrorType> compute_errors(PositiveFloatDPValue h, UpperBoxType const& B) const = 0;
    virtual BoxDomainType build_flow_domain(BoxDomainType D, BoxDomainType V, PositiveFloatDPValue h) const = 0;
    virtual Vector<ValidatedScalarFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const = 0;
};

class InputApproximatorBase : public InputApproximatorInterface {
  protected:
    DifferentialInclusion const& _di;
    SweeperDP _sweeper;
    SharedPointer<ApproximationErrorProcessorInterface> _processor;
    InputApproximatorBase(DifferentialInclusion const& di, SharedPointer<ApproximationErrorProcessorInterface> const& processor,
                       SweeperDP sweeper, InputApproximation kind, Nat num_params_per_input) :
        _di(di), _processor(processor), _sweeper(sweeper), _kind(kind), _num_params_per_input(num_params_per_input) { }
  private:
    InputApproximation _kind;
    Nat _num_params_per_input;
  public:
    virtual InputApproximation kind() const override { return _kind; }
    virtual Vector<ErrorType> compute_errors(PositiveFloatDPValue h, UpperBoxType const& B) const override { return _processor->process(h,B); }
    virtual BoxDomainType build_flow_domain(BoxDomainType D, BoxDomainType V, PositiveFloatDPValue h) const override;
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
    virtual InputApproximation kind() const override { return _impl->kind(); }
    virtual Vector<ErrorType> compute_errors(PositiveFloatDPValue h, UpperBoxType const& B) const override { return _impl->compute_errors(h,B); }
    virtual BoxDomainType build_flow_domain(BoxDomainType D, BoxDomainType V, PositiveFloatDPValue h) const override { return _impl->build_flow_domain(D,V,h); }
    virtual Vector<ValidatedScalarFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const override { return _impl->build_w_functions(DVh,n,m); }
    Vector<ValidatedScalarFunction> build_secondhalf_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const;
};

class ZeroInputApproximator : public InputApproximatorBase {
    friend class InputApproximatorFactory;
protected:
    ZeroInputApproximator(DifferentialInclusion const& di, SweeperDP sweeper)
            : InputApproximatorBase(di,ApproximationErrorProcessorFactory<ZeroApproximation>().create(di),sweeper,InputApproximation::ZERO,0u) { }
    virtual Vector<ValidatedScalarFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const override;
};

class ConstantInputApproximator : public InputApproximatorBase {
    friend class InputApproximatorFactory;
protected:
    ConstantInputApproximator(DifferentialInclusion const& di, SweeperDP sweeper)
            : InputApproximatorBase(di,ApproximationErrorProcessorFactory<ConstantApproximation>().create(di),sweeper,InputApproximation::CONSTANT,1u) { }
    virtual Vector<ValidatedScalarFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const override;
};

class AffineInputApproximator : public InputApproximatorBase {
    friend class InputApproximatorFactory;
protected:
    AffineInputApproximator(DifferentialInclusion const& di, SweeperDP sweeper)
            : InputApproximatorBase(di,ApproximationErrorProcessorFactory<AffineApproximation>().create(di),sweeper,InputApproximation::AFFINE,2u) { }
    virtual Vector<ValidatedScalarFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const override;
};

class SinusoidalInputApproximator : public InputApproximatorBase {
    friend class InputApproximatorFactory;
protected:
    SinusoidalInputApproximator(DifferentialInclusion const& di, SweeperDP sweeper)
            : InputApproximatorBase(di,ApproximationErrorProcessorFactory<SinusoidalApproximation>().create(di),sweeper,InputApproximation::SINUSOIDAL,2u) { }
    virtual Vector<ValidatedScalarFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const override;
};

class PiecewiseInputApproximator : public InputApproximatorBase {
    friend class InputApproximatorFactory;
protected:
    PiecewiseInputApproximator(DifferentialInclusion const& di, SweeperDP sweeper)
            : InputApproximatorBase(di,ApproximationErrorProcessorFactory<PiecewiseApproximation>().create(di),sweeper,InputApproximation::PIECEWISE,2u) { }
    virtual Vector<ValidatedScalarFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const override;
    Vector<ValidatedScalarFunction> build_firsthalf_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const;
public:
    Vector<ValidatedScalarFunction> build_secondhalf_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const;
};

class Reconditioner {
  public:
    virtual Void simplify(ValidatedVectorFunctionModelType& phi) const = 0;
    virtual ValidatedVectorFunctionModelType expand_errors(ValidatedVectorFunctionModelType Phi) const = 0;
};


class LohnerReconditioner : public Reconditioner, public Loggable {
    SweeperDP _sweeper;
    Nat _number_of_variables_to_keep;
public:
    LohnerReconditioner(SweeperDP sweeper, Nat number_of_variables_to_keep);
    void set_sweeper(SweeperDP sweeper) { _sweeper = sweeper; }
    void set_number_of_variables_to_keep(Nat number_of_variables_to_keep) { _number_of_variables_to_keep = number_of_variables_to_keep; }
    virtual ValidatedVectorFunctionModelType expand_errors(ValidatedVectorFunctionModelType f) const override;
    virtual Void simplify(ValidatedVectorFunctionModelType& f) const override;
};


class InclusionIntegratorInterface {
  public:
    virtual List<ValidatedVectorFunctionModelType> flow(DifferentialInclusionIVP const& di_ivp, Real T) = 0;
    virtual Pair<ExactTimeStepType,UpperBoxType> flow_bounds(ValidatedVectorFunction f, BoxDomainType V, BoxDomainType D, ApproximateTimeStepType hsug) const = 0;
    virtual ValidatedVectorFunctionModelType reach(DifferentialInclusion const& di, BoxDomainType D, ValidatedVectorFunctionModelType evolve_function, UpperBoxType B, PositiveFloatDPValue t, PositiveFloatDPValue h) const = 0;
};

class InclusionIntegrator : public virtual InclusionIntegratorInterface, public Loggable {
  protected:
    List<InputApproximation> _approximations;
    SharedPointer<InputApproximator> _approximator;
    SharedPointer<Reconditioner> _reconditioner;
    SweeperDP _sweeper;
    FloatDP _step_size;
    Nat _number_of_steps_between_simplifications;
    Nat _number_of_variables_to_keep;
  public:
    InclusionIntegrator(List<InputApproximation> approximations, SweeperDP sweeper, StepSize step_size);
    template<class... AS> InclusionIntegrator(List<InputApproximation> approximations, SweeperDP sweeper, StepSize step_size, AS... attributes)
        : InclusionIntegrator(approximations, sweeper,step_size) {
        this->set(attributes...);
        _reconditioner.reset(new LohnerReconditioner(_sweeper,_number_of_variables_to_keep));
    }
  public:
    InclusionIntegrator& set(NumberOfStepsBetweenSimplifications n) { _number_of_steps_between_simplifications=n; return *this; }
    InclusionIntegrator& set(NumberOfVariablesToKeep n) { _number_of_variables_to_keep=n; return *this; }
    template<class A, class... AS> InclusionIntegrator& set(A a, AS... as) { this->set(a); this->set(as...); return *this; }

    virtual List<ValidatedVectorFunctionModelType> flow(DifferentialInclusionIVP const& di_ivp, Real T) override;

    virtual Pair<ExactTimeStepType,UpperBoxType> flow_bounds(ValidatedVectorFunction f, BoxDomainType V, BoxDomainType D, ApproximateTimeStepType hsug) const override;
    virtual ValidatedVectorFunctionModelType reach(DifferentialInclusion const& di, BoxDomainType D, ValidatedVectorFunctionModelType evolve_function, UpperBoxType B, PositiveFloatDPValue t, PositiveFloatDPValue h) const override;
  private:
    ValidatedVectorFunctionModelType compute_flow_function(ValidatedVectorFunction const& dyn, BoxDomainType const& domain, UpperBoxType const& B) const;
    ValidatedVectorFunctionModelDP build_reach_function(ValidatedVectorFunctionModelDP evolve_function, ValidatedVectorFunctionModelDP Phi, PositiveFloatDPValue t, PositiveFloatDPValue new_t) const;
    ValidatedVectorFunctionModelDP evaluate_evolve_function(ValidatedVectorFunctionModelDP reach_function, PositiveFloatDPValue t) const;
    ValidatedVectorFunctionModelDP build_secondhalf_piecewise_reach_function(ValidatedVectorFunctionModelDP evolve_function, ValidatedVectorFunctionModelDP Phi, SizeType m, PositiveFloatDPValue t, PositiveFloatDPValue new_t) const;
};


BoxDomainType build_flow_domain(BoxDomainType D, BoxDomainType V, PositiveFloatDPValue h, Nat num_params);


} // namespace Ariadne;

#endif // ARIADNE_DIFFERENTIAL_INCLUSION_HPP
