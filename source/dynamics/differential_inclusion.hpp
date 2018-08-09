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

enum class InputApproximationKind { ZERO, CONSTANT, AFFINE, SINUSOIDAL, PIECEWISE};

std::ostream& operator << (std::ostream& os, const InputApproximationKind& kind) {
    switch (kind) {
        case InputApproximationKind::ZERO: os << "ZERO"; break;
        case InputApproximationKind::CONSTANT: os << "CONSTANT"; break;
        case InputApproximationKind::AFFINE: os << "AFFINE"; break;
        case InputApproximationKind::SINUSOIDAL: os << "SINUSOIDAL"; break;
        case InputApproximationKind::PIECEWISE: os << "PIECEWISE"; break;
    }
    return os;
}

BoxDomainType bounds_to_domain(RealVariablesBox const& var_box);
Pair<RealAssignment,RealInterval> centered_variable_transformation(RealVariable const& v, RealInterval const& bounds);
Pair<RealAssignments,RealVariablesBox> centered_variables_transformation(RealVariablesBox const& inputs);

Vector<FloatDPValue> const& cast_exact(Vector<FloatDPError> const& v) {
    return reinterpret_cast<Vector<FloatDPValue>const&>(v); }

FloatDP volume(Vector<ApproximateIntervalType> const& box);

Bool refines(Vector<UpperIntervalType> const& v1, UpperBoxType const& bx2) {
    return refines(v1,static_cast<Vector<UpperIntervalType>const&>(bx2)); }

Box<Interval<FloatDPValue>> over_approximation(Box<Interval<Real>> const&);


Boolean inputs_are_additive(Vector<ValidatedVectorFunction> const &g, UpperBoxType const &B);

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

Norms compute_norms(ValidatedVectorFunction const&, Vector<ValidatedVectorFunction> const&, BoxDomainType const&, PositiveFloatDPValue const&, UpperBoxType const&);

class ApproximationErrorProcessor {
  public:
    ApproximationErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, InputApproximationKind kind) : _f(f), _g(g), _V(V), _kind(kind), _enable_componentwise_error(false) { }
    Vector<ErrorType> process(PositiveFloatDPValue const& h, UpperBoxType const& B) const;
  protected:
    virtual Vector<ErrorType> compute_errors(Norms const&, PositiveFloatDPValue const&) const = 0;
    InputApproximationKind _kind;
    Boolean _enable_componentwise_error; // TODO: remove such option as soon as the DI paper is completed
  private:
    ValidatedVectorFunction const& _f;
    Vector<ValidatedVectorFunction> const& _g;
    BoxDomainType const& _V;
};

class ZeroErrorProcessor : public ApproximationErrorProcessor {
public:
    ZeroErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V) : ApproximationErrorProcessor(f,g,V,InputApproximationKind::ZERO) { }
    virtual Vector<ErrorType> compute_errors(Norms const&, PositiveFloatDPValue const&) const override;
};

class ConstantErrorProcessor : public ApproximationErrorProcessor {
public:
    ConstantErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V) : ApproximationErrorProcessor(f,g,V,InputApproximationKind::CONSTANT) { }
    virtual Vector<ErrorType> compute_errors(Norms const&, PositiveFloatDPValue const&) const override;
};

class AffineErrorProcessor : public ApproximationErrorProcessor {
public:
    AffineErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V) : ApproximationErrorProcessor(f,g,V,InputApproximationKind::AFFINE) { }
    virtual Vector<ErrorType> compute_errors(Norms const&, PositiveFloatDPValue const&) const override;
};

class AdditiveAffineErrorProcessor : public ApproximationErrorProcessor {
  public:
    AdditiveAffineErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V) : ApproximationErrorProcessor(f,g,V,InputApproximationKind::AFFINE) { }
    virtual Vector<ErrorType> compute_errors(Norms const&, PositiveFloatDPValue const&) const override;
};

class SingleInputAffineErrorProcessor : public ApproximationErrorProcessor {
public:
    SingleInputAffineErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V) : ApproximationErrorProcessor(f,g,V,InputApproximationKind::AFFINE) { }
    virtual Vector<ErrorType> compute_errors(Norms const&, PositiveFloatDPValue const&) const override;
};

class SinusoidalErrorProcessor : public ApproximationErrorProcessor {
public:
    SinusoidalErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V) : ApproximationErrorProcessor(f,g,V,InputApproximationKind::SINUSOIDAL) { }
    virtual Vector<ErrorType> compute_errors(Norms const&, PositiveFloatDPValue const&) const override;
};

class AdditiveSinusoidalErrorProcessor : public ApproximationErrorProcessor {
public:
    AdditiveSinusoidalErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V) : ApproximationErrorProcessor(f,g,V,InputApproximationKind::SINUSOIDAL) { }
    virtual Vector<ErrorType> compute_errors(Norms const&, PositiveFloatDPValue const&) const override;
};

class SingleInputSinusoidalErrorProcessor : public ApproximationErrorProcessor {
public:
    SingleInputSinusoidalErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V) : ApproximationErrorProcessor(f,g,V,InputApproximationKind::SINUSOIDAL) { }
    virtual Vector<ErrorType> compute_errors(Norms const&, PositiveFloatDPValue const&) const override;
};

class PiecewiseErrorProcessor : public ApproximationErrorProcessor {
public:
    PiecewiseErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V) : ApproximationErrorProcessor(f,g,V,InputApproximationKind::PIECEWISE) { }
    virtual Vector<ErrorType> compute_errors(Norms const&, PositiveFloatDPValue const&) const override;
};

class AdditivePiecewiseErrorProcessor : public ApproximationErrorProcessor {
public:
    AdditivePiecewiseErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V) : ApproximationErrorProcessor(f,g,V,InputApproximationKind::PIECEWISE) { }
    virtual Vector<ErrorType> compute_errors(Norms const&, PositiveFloatDPValue const&) const override;
};

class SingleInputPiecewiseErrorProcessor : public ApproximationErrorProcessor {
public:
    SingleInputPiecewiseErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V) : ApproximationErrorProcessor(f,g,V,InputApproximationKind::PIECEWISE) { }
    virtual Vector<ErrorType> compute_errors(Norms const&, PositiveFloatDPValue const&) const override;
};

class InputApproximator {
  protected:
    ValidatedVectorFunction const& _f;
    Vector<ValidatedVectorFunction> const& _g;
    BoxDomainType const& _V;
    SweeperDP _sweeper;
    SharedPointer<ApproximationErrorProcessor> _additive_processor;
    SharedPointer<ApproximationErrorProcessor> _single_input_processor;
    SharedPointer<ApproximationErrorProcessor> _generic_processor;
    InputApproximator(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, const BoxDomainType& V,
                       SweeperDP sweeper, InputApproximationKind kind, Nat num_params_per_input) :
        _f(f), _g(g), _V(V), _sweeper(sweeper), _kind(kind), _num_params_per_input(num_params_per_input) { }
  private:
    InputApproximationKind _kind;
    Nat _num_params_per_input;
  public:
    InputApproximationKind getKind() const { return _kind; }
    Vector<ErrorType> compute_errors(PositiveFloatDPValue h, UpperBoxType const& B) const;
    BoxDomainType build_flow_domain(BoxDomainType D, BoxDomainType V, PositiveFloatDPValue h) const;
    virtual Vector<ValidatedScalarFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const = 0;
};

class ZeroInputApproximator : public InputApproximator {
public:
    ZeroInputApproximator(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, const BoxDomainType& V, SweeperDP sweeper)
            : InputApproximator(f,g,V,sweeper,InputApproximationKind::ZERO,0u) {
        _additive_processor.reset(new ZeroErrorProcessor(_f,_g,_V));
        _single_input_processor.reset(new ZeroErrorProcessor(_f,_g,_V));
        _generic_processor.reset(new ZeroErrorProcessor(_f,_g,_V));
    }
protected:
    virtual Vector<ValidatedScalarFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const override;
};

class ConstantInputApproximator : public InputApproximator {
public:
    ConstantInputApproximator(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, const BoxDomainType& V, SweeperDP sweeper)
            : InputApproximator(f,g,V,sweeper,InputApproximationKind::CONSTANT,1u) {
        _additive_processor.reset(new ConstantErrorProcessor(_f,_g,_V));
        _single_input_processor.reset(new ConstantErrorProcessor(_f,_g,_V));
        _generic_processor.reset(new ConstantErrorProcessor(_f,_g,_V));
    }
protected:
    virtual Vector<ValidatedScalarFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const override;
};

class AffineInputApproximator : public InputApproximator {
public:
    AffineInputApproximator(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, const BoxDomainType& V, SweeperDP sweeper)
            : InputApproximator(f,g,V,sweeper,InputApproximationKind::AFFINE,2u) {
        _additive_processor.reset(new AdditiveAffineErrorProcessor(_f,_g,_V));
        _single_input_processor.reset(new SingleInputAffineErrorProcessor(_f,_g,_V));
        _generic_processor.reset(new AffineErrorProcessor(_f,_g,_V));
    }
protected:
    virtual Vector<ValidatedScalarFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const override;
};

class SinusoidalInputApproximator : public InputApproximator {
public:
    SinusoidalInputApproximator(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, const BoxDomainType& V, SweeperDP sweeper)
            : InputApproximator(f,g,V,sweeper,InputApproximationKind::SINUSOIDAL,2u) {
        _additive_processor.reset(new AdditiveSinusoidalErrorProcessor(_f,_g,_V));
        _single_input_processor.reset(new SingleInputSinusoidalErrorProcessor(_f,_g,_V));
        _generic_processor.reset(new SinusoidalErrorProcessor(_f,_g,_V));
    }
protected:
    virtual Vector<ValidatedScalarFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const override;
};

class PiecewiseInputApproximator : public InputApproximator {
public:
    PiecewiseInputApproximator(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, const BoxDomainType& V, SweeperDP sweeper)
            : InputApproximator(f,g,V,sweeper,InputApproximationKind::PIECEWISE,2u) {
        _additive_processor.reset(new AdditivePiecewiseErrorProcessor(_f,_g,_V));
        _single_input_processor.reset(new SingleInputPiecewiseErrorProcessor(_f,_g,_V));
        _generic_processor.reset(new PiecewiseErrorProcessor(_f,_g,_V));
    }
public:
    virtual Vector<ValidatedScalarFunction> build_w_functions(BoxDomainType DVh, SizeType n, SizeType m) const override;
    Vector<ValidatedScalarFunction> build_firsthalf_approximating_function(BoxDomainType DVh, SizeType n, SizeType m) const;
    Vector<ValidatedScalarFunction> build_secondhalf_approximating_function(BoxDomainType DVh, SizeType n, SizeType m) const;
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
    virtual List<ValidatedVectorFunctionModelType> flow(const List<DottedRealAssignment>& dynamics, const RealVariablesBox& inputs, const RealVariablesBox& initial, ValidatedVectorFunction f, Vector<ValidatedVectorFunction> g, BoxDomainType V, BoxDomainType X0, Real T) = 0;

    virtual Pair<ExactTimeStepType,UpperBoxType> flow_bounds(ValidatedVectorFunction f, BoxDomainType V, BoxDomainType D, ApproximateTimeStepType hsug) const = 0;
    virtual ValidatedVectorFunctionModelType compute_flow_function(const List<DottedRealAssignment>& dynamics, const RealVariablesBox& inputs, const RealVariablesBox& initial, ValidatedVectorFunction f,
                                                                   Vector<ValidatedVectorFunction> g, BoxDomainType V,
                                                                   BoxDomainType D, ExactTimeStepType h, UpperBoxType B) const = 0;
};

class InclusionIntegrator : public virtual InclusionIntegratorInterface, public Loggable {
  protected:
    List<SharedPointer<InputApproximator>> _approximations;
    SharedPointer<InputApproximator> _approximation;
    SharedPointer<Reconditioner> _reconditioner;
    SweeperDP _sweeper;
    FloatDP _step_size;
    Nat _number_of_steps_between_simplifications;
    Nat _number_of_variables_to_keep;
  public:
    InclusionIntegrator(List<SharedPointer<InputApproximator>> approximations, SweeperDP sweeper, StepSize step_size);
    template<class... AS> InclusionIntegrator(List<SharedPointer<InputApproximator>> approximations, SweeperDP sweeper, StepSize step_size, AS... attributes)
        : InclusionIntegrator(approximations, sweeper,step_size) {
        this->set(attributes...);
        _reconditioner.reset(new LohnerReconditioner(_sweeper,_number_of_variables_to_keep));
    }
  public:
    InclusionIntegrator& set(NumberOfStepsBetweenSimplifications n) { _number_of_steps_between_simplifications=n; return *this; }
    InclusionIntegrator& set(NumberOfVariablesToKeep n) { _number_of_variables_to_keep=n; return *this; }
    template<class A, class... AS> InclusionIntegrator& set(A a, AS... as) { this->set(a); this->set(as...); return *this; }

    virtual List<ValidatedVectorFunctionModelType> flow(const List<DottedRealAssignment>& dynamics, const RealVariablesBox& inputs, const RealVariablesBox& initial, ValidatedVectorFunction f, Vector<ValidatedVectorFunction> g, BoxDomainType V, BoxDomainType X0, Real T) override;

    virtual Pair<ExactTimeStepType,UpperBoxType> flow_bounds(ValidatedVectorFunction f, BoxDomainType V, BoxDomainType D, ApproximateTimeStepType hsug) const override;

    virtual ValidatedVectorFunctionModelType compute_flow_function(const List<DottedRealAssignment>& dynamics, const RealVariablesBox& inputs, const RealVariablesBox& initial, ValidatedVectorFunction f,
                                                                   Vector<ValidatedVectorFunction> g, BoxDomainType V,
                                                                   BoxDomainType D, ExactTimeStepType h, UpperBoxType B) const override;
  private:
    ValidatedVectorFunctionModelDP build_reach_function(ValidatedVectorFunctionModelDP evolve_function, ValidatedVectorFunctionModelDP Phi, PositiveFloatDPValue t, PositiveFloatDPValue new_t) const;
    ValidatedVectorFunctionModelDP build_secondhalf_piecewise_reach_function(ValidatedVectorFunctionModelDP evolve_function, ValidatedVectorFunctionModelDP Phi, SizeType m, PositiveFloatDPValue t, PositiveFloatDPValue new_t) const;
};


BoxDomainType build_flow_domain(BoxDomainType D, BoxDomainType V, PositiveFloatDPValue h, Nat num_params);


} // namespace Ariadne;

#endif // ARIADNE_DIFFERENTIAL_INCLUSION_HPP
