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
using SweeperDP = Sweeper<FloatDP>;

using ApproximateTimeStepType = PositiveFloatDPApproximation;
using ExactTimeStepType = PositiveFloatDPValue;

Vector<FloatDPValue> const& cast_exact(Vector<FloatDPError> const& v) {
    return reinterpret_cast<Vector<FloatDPValue>const&>(v); }

Bool refines(Vector<UpperIntervalType> const& v1, UpperBoxType const& bx2) {
    return refines(v1,static_cast<Vector<UpperIntervalType>const&>(bx2)); }

Box<Interval<FloatDPValue>> over_approximation(Box<Interval<Real>> const&);

template<class F1, class F2, class F3, class... FS> decltype(auto) combine(F1 const& f1, F2 const& f2, F3 const& f3, FS const& ... fs) {
    return combine(combine(f1,f2),f3,fs...); }
template<class F1, class F2, class F3, class... FS> decltype(auto) join(F1 const& f1, F2 const& f2, F3 const& f3, FS const& ... fs) {
    return join(join(f1,f2),f3,fs...); }

BoxDomainType error_domain(SizeType n, FloatDPError e) {
    FloatDPValue const& v=cast_exact(e);
    return BoxDomainType(n,IntervalDomainType(-v,+v));
}

//! \brief The function dexp(x)=(exp(x)-1)/x.
//! Note that the function is positive and monotone increasing.
template<class F> PositiveUpperBound<F> dexp(UpperBound<F> const& x) {
    if(x.raw()>=0) { return cast_positive(exp(x)-1)/cast_positive(cast_exact(x)); }
    else { return cast_positive(cast_exact(1-exp(x)))/cast_positive(-x); }
}


class Reconditioner {
  public:
    virtual Void simplify(ValidatedVectorFunctionModelType& phi) const = 0;
    virtual ValidatedVectorFunctionModelType expand_errors(ValidatedVectorFunctionModelType Phi) const = 0;
};

Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPUpperBound> compute_norms(ValidatedVectorFunction const&, Vector<ValidatedVectorFunction> const&, BoxDomainType const&, UpperBoxType const& B);
Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPError,FloatDPUpperBound> compute_norms_additive(ValidatedVectorFunction const&, Vector<ValidatedVectorFunction> const&, BoxDomainType const&, UpperBoxType const& B);

class InclusionErrorProcessor {
  public:
    InclusionErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B);
    ErrorType process() const;
  protected:
    virtual ErrorType compute_error(FloatDPError const&,FloatDPError const&,FloatDPError const&,FloatDPError const&,FloatDPError const&,FloatDPError const&,FloatDPUpperBound const&,PositiveFloatDPValue const&) const = 0;
  private:
    ValidatedVectorFunction const& _f;
    Vector<ValidatedVectorFunction> const& _g;
    BoxDomainType const& _V;
    PositiveFloatDPValue const& _h;
    UpperBoxType const& _B;
};

class AffineErrorProcessor : public InclusionErrorProcessor {
  public:
    AffineErrorProcessor(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType const& V, PositiveFloatDPValue const& h, UpperBoxType const& B);
  public:
    virtual ErrorType compute_error(FloatDPError const&,FloatDPError const&,FloatDPError const&,FloatDPError const&,FloatDPError const&,FloatDPError const&,FloatDPUpperBound const&,PositiveFloatDPValue const&) const;
};


class LohnerReconditioner : public Reconditioner, public Loggable {
    SweeperDP _sweeper;
    Nat _number_of_variables_to_keep;
public:
    LohnerReconditioner(SweeperDP sweeper, Nat number_of_variables_to_keep);
    virtual ValidatedVectorFunctionModelType expand_errors(ValidatedVectorFunctionModelType Phi) const override;
    virtual Void simplify(ValidatedVectorFunctionModelType& phi) const override;
};

class InclusionIntegratorInterface {
  public:
    virtual List<ValidatedVectorFunctionModelType> flow(ValidatedVectorFunction f, Vector<ValidatedVectorFunction> g, BoxDomainType V, BoxDomainType X0, Real T) const = 0;

    virtual Pair<ExactTimeStepType,UpperBoxType> flow_bounds(ValidatedVectorFunction f, Vector<ValidatedVectorFunction> g, UpperBoxType V, ExactBoxType D, ApproximateTimeStepType hsug) const = 0;
    virtual ValidatedVectorFunctionModelType compute_flow_function(ValidatedVectorFunction f,
                                                                   Vector<ValidatedVectorFunction> g, BoxDomainType V,
                                                                   BoxDomainType D, ExactTimeStepType h, UpperBoxType B) const = 0;
};

class InclusionIntegratorBase : public virtual InclusionIntegratorInterface, public Loggable {
  protected:
    SharedPointer<Reconditioner> _reconditioner;
    SweeperDP _sweeper;
    FloatDP _step_size;
    Nat _number_of_steps_between_simplifications;
    Nat _number_of_variables_to_keep;
    InclusionIntegratorBase(SweeperDP sweeper, StepSize step_size);
    template<class... AS> InclusionIntegratorBase(SweeperDP sweeper, StepSize step_size, AS... attributes)
        : InclusionIntegratorBase(sweeper,step_size) {
        this->set(attributes...);
        _reconditioner.reset(new LohnerReconditioner(_sweeper,_number_of_variables_to_keep));
    }
  public:
    InclusionIntegratorBase& set(NumberOfStepsBetweenSimplifications n) { _number_of_steps_between_simplifications=n; return *this; }
    InclusionIntegratorBase& set(NumberOfVariablesToKeep n) { _number_of_variables_to_keep=n; return *this; }
    template<class A, class... AS> InclusionIntegratorBase& set(A a, AS... as) { this->set(a); this->set(as...); return *this; }

    ValidatedVectorFunctionModelType expand_errors(ValidatedVectorFunctionModelType Phi) const;
    ValidatedVectorFunctionModelType simplify(ValidatedVectorFunctionModelType Phi) const;
    virtual List<ValidatedVectorFunctionModelType> flow(ValidatedVectorFunction f, Vector<ValidatedVectorFunction> g, BoxDomainType V, BoxDomainType X0, Real T) const override;

    virtual Pair<ExactTimeStepType,UpperBoxType> flow_bounds(ValidatedVectorFunction f, Vector<ValidatedVectorFunction> g, UpperBoxType V, ExactBoxType D, ApproximateTimeStepType hsug) const override;

    virtual ValidatedVectorFunctionModelType compute_flow_function(ValidatedVectorFunction f,
                                                                   Vector<ValidatedVectorFunction> g, BoxDomainType V,
                                                                   BoxDomainType D, ExactTimeStepType h, UpperBoxType B) const override;
  protected:
    virtual ErrorType compute_error(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType V, PositiveFloatDPValue h, UpperBoxType const& B) const = 0;
    virtual BoxDomainType compute_flow_domain(BoxDomainType D, PositiveFloatDPValue h, BoxDomainType V, ErrorType e) const = 0;
    virtual ValidatedVectorFunctionModelType compute_approximating_function(BoxDomainType DHPE, SizeType n) const = 0;
  private:
    ValidatedVectorFunctionModelDP compute_reach_function(ValidatedVectorFunctionModelDP evolve_function, ValidatedVectorFunctionModelDP Phi, PositiveFloatDPValue t, PositiveFloatDPValue new_t) const;
};

class InclusionIntegratorAffineW : public InclusionIntegratorBase {
  public:
    template<class... AS> InclusionIntegratorAffineW(SweeperDP sweeper, StepSize step_size, AS... attributes)
        : InclusionIntegratorBase(sweeper,step_size,attributes...) { }
  protected:
    virtual ErrorType compute_error(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType V, PositiveFloatDPValue h, UpperBoxType const& B) const override;
    virtual BoxDomainType compute_flow_domain(BoxDomainType D, PositiveFloatDPValue h, BoxDomainType V, ErrorType e) const override;
    virtual ValidatedVectorFunctionModelType compute_approximating_function(BoxDomainType DHPE, SizeType n) const override;
  private:
    Tuple<FloatDPError,FloatDPError,FloatDPError,FloatDPUpperBound> compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, UpperBoxType const& B) const;
};



class InclusionIntegratorConstantW : public InclusionIntegratorBase {
  public:
    template<class... AS> InclusionIntegratorConstantW(SweeperDP sweeper, StepSize step_size, AS... attributes)
        : InclusionIntegratorBase(sweeper,step_size,attributes...) {  }
  protected:
    virtual ErrorType compute_error(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, BoxDomainType V, PositiveFloatDPValue h, UpperBoxType const& B) const override;
    virtual BoxDomainType compute_flow_domain(BoxDomainType D, PositiveFloatDPValue h, BoxDomainType V, ErrorType e) const override;
    virtual ValidatedVectorFunctionModelType compute_approximating_function(BoxDomainType DHPE, SizeType n) const override;
  private:
    Tuple<FloatDPError,FloatDPError,FloatDPUpperBound> compute_norms(ValidatedVectorFunction const& f, Vector<ValidatedVectorFunction> const& g, UpperBoxType const& B) const;
};

} // namespace Ariadne;

#endif // ARIADNE_DIFFERENTIAL_INCLUSION_HPP
