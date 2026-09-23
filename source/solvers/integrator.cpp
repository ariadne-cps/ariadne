/***************************************************************************
 *            solvers/integrator.cpp
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

#include "function/functional.hpp"
#include "config.hpp"

#include <iomanip>

#include "solvers/integrator.hpp"
#include "solvers/bounder.hpp"

#include "conclog/logging.hpp"
#include "utility/container.hpp"
#include "utility/tuple.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/differential.hpp"
#include "algebra/sweeper.hpp"
#include "algebra/algebra.hpp"
#include "function/function.hpp"
#include "function/function_patch.hpp"
#include "function/function_model.hpp"
#include "function/formula.hpp"
#include "function/scaling.hpp"
#include "function/taylor_function.hpp"
#include "function/taylor_model.hpp"

#include "function/polynomial.hpp"
#include "geometry/interval.hpp"

#include "algebra/expansion.inl.hpp"

using namespace ConcLog;

namespace Ariadne {

typedef ValidatedVectorMultivariateTaylorFunctionModelDP FlowStepTaylorModelType;


// TODO: Move this functionality to Numeric
inline ApproximateNumber convert_to_approximate_number(ValidatedUpperNumber y) { return Approximation<FloatDP>(y.get(dp)); }
inline Bool refines(ValidatedUpperNumber y1, ValidatedUpperNumber y2) { return refines(y1.get(dp),y2.get(dp)); }
inline ApproximateNumber operator-(ValidatedUpperNumber y1, ValidatedUpperNumber y2) {
    return convert_to_approximate_number(y1)-convert_to_approximate_number(y2); }
inline ApproximateNumber operator/(ApproximateNumber y1, ValidatedUpperNumber y2) {
    return y1/convert_to_approximate_number(y2); }
inline ExactNumber cast_exact(ApproximateNumber y) { return cast_exact(y.get(dp)); }

OutputStream& operator<<(OutputStream& os, FlowStepModelType const& fsm) {
    return os << static_cast<ValidatedVectorMultivariateFunctionPatch const&>(fsm);
}

OutputStream& operator<<(OutputStream& os, FlowModelType const& fm) {
    os << "[ "; for(SizeType i=0; i!=fm.size(); ++i) { if (i!=0u) { os << ",\n"; } os << "  " << fm[i]; } os << " ]"; return os ;
}

static const FloatDP zero=FloatDP(0,dp);

inline UpperBoxType operator+(Vector<ExactIntervalType> bx, Vector<UpperIntervalType> const& ex) {
    return Vector<UpperIntervalType>(bx) + ex;
}

inline UpperBoxType operator+(Vector<UpperIntervalType> bx, Vector<FloatDPBounds> const& v) {
    return bx + Vector<UpperIntervalType>(v);
}

inline UpperBoxType operator+(Vector<ExactIntervalType> bx, Vector<FloatDPBounds> const& v) {
    return Vector<UpperIntervalType>(bx) + Vector<UpperIntervalType>(v);
}

inline ExactDouble cast_exact_double(Attribute<ApproximateDouble> a) { return cast_exact(static_cast<ApproximateDouble>(a)); }

IncompleteFlowException::IncompleteFlowException(const StringType& what, FlowStepModelType const& model)
    : std::runtime_error(what), _computed_model(new FlowStepModelType(model)) {
}

IntegratorBase::IntegratorBase(Sweeper<FloatDP> s)
    : _function_factory(make_taylor_function_patch_factory(s)) { }

Void
IntegratorBase::set_function_factory(const ValidatedFunctionPatchFactory& factory)
{
    this->_function_factory=factory;
}

const ValidatedFunctionPatchFactory&
IntegratorBase::function_factory() const
{
    return this->_function_factory;
}

BoundedIntegratorBase::BoundedIntegratorBase(Sweeper<FloatDP> sweeper, LipschitzTolerance lipschitz) :
        IntegratorBase(sweeper), _bounder_ptr(new EulerBounder(lipschitz))
{
    ARIADNE_PRECONDITION(lipschitz>0.0_x)
}

const BounderInterface&
BoundedIntegratorBase::bounder() const
{
    return *this->_bounder_ptr;
}

Void
BoundedIntegratorBase::set_bounder(const BounderInterface& bounder)
{
    this->_bounder_ptr=BounderPointer(bounder.clone());
}

Pair<StepSizeType,UpperBoxType>
BoundedIntegratorBase::flow_bounds(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& D, const Suggestion<StepSizeType>& hsug) const {
    return this->_bounder_ptr->compute(vf,D,hsug);
}

Pair<StepSizeType,UpperBoxType>
BoundedIntegratorBase::flow_bounds(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& D, const ExactBoxType& A, const Suggestion<StepSizeType>& hsug) const {
    return this->_bounder_ptr->compute(vf,D,A,hsug);
}

Pair<StepSizeType,UpperBoxType>
BoundedIntegratorBase::flow_bounds(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& D, StepSizeType const& t, const ExactBoxType& A, const Suggestion<StepSizeType>& hsug) const {
    return this->_bounder_ptr->compute(vf,D,t,A,hsug);
}


FlowStepModelType
IntegratorBase::flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx) const
{
    StepSizeType hsug = 1.0_dy;
    return this->flow_step(vf,dx,suggest(hsug));
}

FlowStepModelType
IntegratorBase::flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx, const Suggestion<StepSizeType>& hsug) const
{
    StepSizeType h=static_cast<const StepSizeType&>(hsug);
    StepSizeType hprev=h*1.5_dy;
    while(true) {
        try {
            return flow_step(vf,dx,h,dx);
        } catch(const FlowTimeStepException& e) {
            StepSizeType hnew=hlf(hprev);
            hprev=h;
            h=StepSizeType(hnew.get_d());
            CONCLOG_PRINTLN_AT(1,"Reduced h to "<<h);
        }
    }
}

FlowStepModelType
IntegratorBase::flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx, const StepSizeType& h) const
{
    // The Picard implementation works on the centred auxiliary time domain
    // [-h,+h] and only restricts the returned model to [0,h] at the end.
    // Hence its a-priori state bound must enclose both the forward and backward
    // displacement from dx.  A one-sided +h*f(dx) bound is invalid here.
    auto const radius = 1.5_dy * (dx-dx.centre());
    auto const displacement = (1.5_dy * h) * cast_singleton(image(dx,vf));
    UpperBoxType const forward_box = dx + radius + displacement;
    UpperBoxType const backward_box = dx + radius + (-displacement);
    UpperBoxType bx = hull(forward_box,backward_box);
    bx = hull(UpperBoxType(dx),bx);
    StepSizeType hred=h;
    FlowStepModelType phi = this->flow_step(vf,dx,hred,bx);
    Nat reduction_count=0u;
    while (not definitely(subset(phi.range(),cast_exact_box(bx)))) {
        std::cerr << "[FixedStepBoundCheck]"
                  << " reduction=" << reduction_count
                  << " h=" << hred
                  << " bx=" << bx
                  << " flow_range=" << phi.range()
                  << std::endl;
        if(reduction_count>=32u) {
            std::stringstream msg;
            msg << "IntegratorBase::flow_step(vf,dx,h): flow range did not fit the "
                << "fixed bounding box after " << reduction_count
                << " step halvings; bx=" << bx << ", last_range=" << phi.range();
            throw IncompleteFlowException(msg.str(),phi);
        }
        hred=hlf(hred);
        ++reduction_count;
        phi = this->flow_step(vf,dx,hred,bx);
    }
    if(reduction_count!=0u) {
        std::cerr << "[FixedStepBoundCheck]"
                  << " accepted=true reductions=" << reduction_count
                  << " h=" << hred
                  << " bx=" << bx
                  << " flow_range=" << phi.range()
                  << std::endl;
    }
    if (hred==h) {
        return phi;
    } else {
        std::stringstream msg;
        msg << "BoundedIntegratorBase::flow_step(vf,dx,h): vf=" << vf << ", dx=" << dx << ", h=" << h << "; Could not bound flow of vf over state domain dx for time h";
        throw IncompleteFlowException(msg.str(), phi);
    }
}


FlowStepModelType
BoundedIntegratorBase::flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx, const Suggestion<StepSizeType>& hsug) const
{
    StepSizeType h;
    UpperBoxType bx;
    make_lpair(h,bx)=this->flow_bounds(vf,dx,hsug);
    StepSizeType hprev=h*1.5_dy;
    while(true) {
        try {
            return this->flow_step(vf,dx,h,bx);
        } catch(const FlowTimeStepException& e) {
            StepSizeType hnew=hlf(hprev);
            hprev=h;
            h=StepSizeType(hnew.get_d());
            CONCLOG_PRINTLN_AT(1,"Reduced h to "<<h);
        }
    }
}

FlowStepModelType
BoundedIntegratorBase::flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx, const StepSizeType& h) const
{
    Suggestion<StepSizeType> hsug = suggest(h);
    StepSizeType hbnd;
    UpperBoxType bx;
    make_lpair(hbnd,bx)=this->flow_bounds(vf,dx,hsug);
    if (hbnd < h) {
        std::stringstream msg;
        msg << "BoundedIntegratorBase::flow_step(vf,dx,h): vf=" << vf << ", dx=" << dx << ", h=" << h << "; Could not bound flow of vf over state domain dx for time h";
        FlowStepModelType phi = this->flow_step(vf,dx,hbnd,bx);
        throw IncompleteFlowException(msg.str(), phi);
    }
    FlowStepModelType phi = this->flow_step(vf,dx,h,bx);
    return phi;
}


inline ExactDouble operator*(ExactDouble, TwoExp);
inline ExactDouble operator/(ExactDouble, TwoExp);

TaylorPicardIntegrator::TaylorPicardIntegrator(StepMaximumError err)
    : TaylorPicardIntegrator(err,ThresholdSweeper<FloatDP>(DP(),err.value()/8),DEFAULT_LIPSCHITZ_TOLERANCE,MinimumTemporalOrder(0),MaximumTemporalOrder(12)) { }

TaylorPicardIntegrator::TaylorPicardIntegrator(StepMaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzTolerance lip,
                                               MinimumTemporalOrder minto, MaximumTemporalOrder maxto)
    : BoundedIntegratorBase(sweeper, lip), _step_maximum_error(cast_exact(err.value())), _sweeper(sweeper)
    , _minimum_temporal_order(minto), _maximum_temporal_order(maxto) { }

FlowStepModelType
TaylorPicardIntegrator::flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx, const StepSizeType& h, const UpperBoxType& bx) const
{
    ARIADNE_PRECONDITION(vf.result_size()==dx.dimension());
    ARIADNE_PRECONDITION(vf.argument_size()==dx.dimension());
    ARIADNE_PRECONDITION(bx.dimension()==dx.dimension());
    return this->_flow_step(vf,dx,IntervalDomainType(0,h),BoxDomainType(0u),bx);
}

FlowStepModelType
TaylorPicardIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& D, const Interval<StepSizeType>& T, const ExactBoxType& A, const UpperBoxType& B) const
{
    ARIADNE_PRECONDITION(f.result_size()==D.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==D.dimension()+T.dimension()+A.dimension());
    ARIADNE_PRECONDITION(B.dimension()==D.dimension());
    return this->_flow_step(f,D,IntervalDomainType(T),A,B);
}

FlowStepModelType
TaylorPicardIntegrator::_flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& D, const ExactIntervalType& T, const ExactBoxType& A, const UpperBoxType& B) const
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("f="<<f);
    CONCLOG_PRINTLN("D="<<D<<" T="<<T<<", A="<<A<<", B="<<B);

    const bool is_autonomous = (f.argument_size()==D.dimension()+A.dimension());

    const SizeType nx=D.size();
    const SizeType na=A.size();

    Range tarng = is_autonomous ? Range(nx+1u,nx+1u+na) : Range(nx,nx+1u+na);

    StepSizeType t=static_cast<StepSizeType>(T.lower_bound());
    StepSizeType h=static_cast<StepSizeType>(T.upper_bound())-t;

    // Time interval centred on initial time, which will make the antiderivative more efficient
    ExactIntervalType wT(t-h,t+h);
    ARIADNE_ASSERT(t==med(wT));

    ExactBoxType dom=join(D,T,A);
    ExactBoxType wdom=join(D,wT,A);
    UpperBoxType const& bx=B;
    CONCLOG_PRINTLN_AT(2,"dom="<<dom<<", wdom="<<wdom);

    FlowStepModelType phi0=this->function_factory().create_projection(wdom,range(0,nx));
    CONCLOG_PRINTLN_AT(1,"phi0="<<phi0);
    FlowStepModelType phi=this->function_factory().create_constants(wdom,cast_singleton(bx));
    FlowStepModelType ta=this->function_factory().create_projection(wdom,tarng);

    CONCLOG_PRINTLN_AT(1,"phi="<<phi);
    FlowStepModelType fphi=compose(f,join(phi0,ta));
    for(DegreeType k=0; k!=this->_maximum_temporal_order; ++k) {
        Bool below_maximum_error=definitely(phi.error()<this->step_maximum_error());
        try {
            fphi=compose(f,join(std::move(phi),ta));
            CONCLOG_PRINTLN_AT(2,"fphi="<<fphi);
        } catch(...) {
            ARIADNE_THROW(FlowTimeStepException,"TaylorPicardIntegrator::flow_step","Could not evaluate f="<<f<<" over model join(phi,t,a) with phi="<<phi);
        }
        // NOTE: In principle safer to use antiderivative(fphi,nx,t) here,
        // but since t is the midpoint of wdom, the (standard) antiderivative works
        // TODO: Change based antiderivative to be efficient when t is midpoint of domain
        phi=antiderivative(fphi,nx)+phi0;
        CONCLOG_PRINTLN_AT(2,"phi="<<phi);
        if(below_maximum_error && k>=this->_minimum_temporal_order) { break; }
    }
    if (possibly(phi.error()>this->step_maximum_error())) {
        ARIADNE_THROW(FlowTimeStepException,"TaylorPicardIntegrator::flow_step","Integration of "<<f<<" starting in "<<D<<" over time interval "<<T<<" of length "<<h<<" has error "<<phi.error()<<" after "<<this->_maximum_temporal_order<<" iterations, which exceeds step maximum error "<<this->step_maximum_error());
    }

    return restriction(phi,dom);
}

Void TaylorPicardIntegrator::_write(OutputStream& os) const {
    os << "TaylorPicardIntegrator"
       << ", function_factory = " << this->function_factory()
       << ", step_maximum_error = " << this->step_maximum_error()
       << ", sweeper = " << this->sweeper()
       << ", minimum_temporal_order = " << this->minimum_temporal_order()
       << ", maximum_temporal_order = " << this->maximum_temporal_order()
       << " )";
}

GradedTaylorPicardIntegrator::GradedTaylorPicardIntegrator(StepMaximumError err, Order order)
        : GradedTaylorPicardIntegrator(err,order,step_sweep_threshold=0.0) { }

GradedTaylorPicardIntegrator::GradedTaylorPicardIntegrator(StepMaximumError err, Order order, StepSweepThreshold threshold)
        : IntegratorBase(GradedThresholdSweeper<FloatDP>(
              DoublePrecision(), order, FloatDP(cast_exact(threshold.value()),DoublePrecision()))),
          _step_maximum_error(cast_exact(err.value())),
          _sweeper(GradedThresholdSweeper<FloatDP>(
              DoublePrecision(), order, FloatDP(cast_exact(threshold.value()),DoublePrecision()))),
          _error_refinement_minimum_improvement_percentage(cast_exact(0.02)),
          _maximum_error_refinement_iterations(0u),
          _order(order), _sweep_threshold(threshold.value()), _diagnostics(false) { }

FlowStepModelType
GradedTaylorPicardIntegrator::flow_step(const ValidatedVectorMultivariateFunction& vf, const ExactBoxType& dx, const StepSizeType& h, const UpperBoxType& bx) const
{
    ARIADNE_PRECONDITION(vf.result_size()==dx.dimension());
    ARIADNE_PRECONDITION(vf.argument_size()==dx.dimension());
    ARIADNE_PRECONDITION(bx.dimension()==dx.dimension());
    return this->_flow_step(vf,dx,IntervalDomainType(0,h),BoxDomainType(0u),bx);
}

FlowStepModelType
GradedTaylorPicardIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& D, const Interval<StepSizeType>& T, const ExactBoxType& A, const UpperBoxType& B) const
{
    ARIADNE_PRECONDITION(f.result_size()==D.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==D.dimension()+T.dimension()+A.dimension());
    ARIADNE_PRECONDITION(B.dimension()==D.dimension());
    return this->_flow_step(f,D,IntervalDomainType(T),A,B);
}

FlowStepModelType
GradedTaylorPicardIntegrator::_flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& D, const ExactIntervalType& T, const ExactBoxType& A, const UpperBoxType& B) const
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN("f="<<f);
    CONCLOG_PRINTLN("D="<<D<<" T="<<T<<", A="<<A<<", B="<<B);

    const bool is_autonomous = (f.argument_size()==D.dimension()+A.dimension());

    const SizeType nx=D.size();
    const SizeType na=A.size();

    Range tarng = is_autonomous ? Range(nx+1u,nx+1u+na) : Range(nx,nx+1u+na);

    StepSizeType t=static_cast<StepSizeType>(T.lower_bound());
    StepSizeType h=static_cast<StepSizeType>(T.upper_bound())-t;

    // Time interval centred on initial time, which will make the antiderivative more efficient
    ExactIntervalType wT(t-h,t+h);
    ARIADNE_ASSERT(t==med(wT));

    ExactBoxType dom=join(D,T,A);
    ExactBoxType wdom=join(D,wT,A);
    UpperBoxType const& bx=B;
    CONCLOG_PRINTLN_AT(2,"dom="<<dom<<", wdom="<<wdom);

    FlowStepModelType phi0=this->function_factory().create_projection(wdom,range(0,nx));
    CONCLOG_PRINTLN_AT(1,"phi0="<<phi0);
    FlowStepModelType phi=this->function_factory().create_constants(wdom,cast_singleton(bx));
    FlowStepModelType ta=this->function_factory().create_projection(wdom,tarng);

    CONCLOG_PRINTLN_AT(1,"phi="<<phi);

    auto diagnostic_nnz = [](FlowStepModelType const& model) {
        auto const& taylor_model =
            dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(model.reference());
        SizeType nnz=0u;
        for(SizeType i=0; i!=taylor_model.size(); ++i) {
            nnz+=taylor_model[i].number_of_nonzeros();
        }
        return nnz;
    };
    auto diagnostic_print = [&](const char* phase, DegreeType iteration, FlowStepModelType const& model) {
        if(this->_diagnostics) {
            std::cerr << "[GradedStepBreakdown]"
                      << " phase=" << phase
                      << " iteration=" << iteration
                      << " nnz=" << diagnostic_nnz(model)
                      << " errors=" << model.errors()
                      << " error=" << model.error()
                      << std::endl;
        }
    };

    FlowStepModelType fphi=compose(f,join(phi0,ta));
    diagnostic_print("initial_fphi",0u,fphi);
    for (DegreeType k=0; k!=this->_order; ++k) {
        try {
            fphi=compose(f,join(std::move(phi),ta));
            CONCLOG_PRINTLN_AT(2,"fphi="<<fphi);
        } catch (...) {
            ARIADNE_THROW(FlowTimeStepException,"GradedTaylorPicardIntegrator::flow_step","Could not evaluate f="<<f<<" over model join(phi,t,a) with phi="<<phi);
        }
        phi=antiderivative(fphi,nx)+phi0;
        CONCLOG_PRINTLN_AT(2,"phi="<<phi);
        diagnostic_print("picard",k+1u,phi);
    }
    auto errors = phi.errors();
    CONCLOG_PRINTLN_AT(2,"initial errors to validate=" << errors);
    fphi=compose(f,join(std::move(phi),ta));
    diagnostic_print("validation_fphi",1u,fphi);
    phi=antiderivative(fphi,nx)+phi0;
    diagnostic_print("validation_phi",1u,phi);
    auto new_errors = phi.errors();
    for (SizeType i=0; i<errors.size(); ++i) {
        if (not refines(new_errors[i],errors[i])) {
            ARIADNE_THROW(FlowTimeStepException,"GradedTaylorPicardIntegrator::flow_step","Integration of "<<f<<" starting in "<<D<<" over time interval "<<T<<" of length "<<h<<" has errors "<<new_errors<<", which are not smaller than previous errors " << errors);
        }
    }
    errors = new_errors;
    CONCLOG_PRINTLN_AT(2,"validated errors=" << errors);
    DegreeType diagnostic_refinement_iteration=0u;
    while (true) {
        ++diagnostic_refinement_iteration;
        fphi=compose(f,join(std::move(phi),ta));
        diagnostic_print("refinement_fphi",diagnostic_refinement_iteration,fphi);
        phi=antiderivative(fphi,nx)+phi0;
        diagnostic_print("refinement_phi",diagnostic_refinement_iteration,phi);
        new_errors = phi.errors();
        if(this->_maximum_error_refinement_iterations!=0u
           && diagnostic_refinement_iteration>=this->_maximum_error_refinement_iterations) {
            errors=new_errors;
            break;
        }
        Bool has_improved = false;
        for (SizeType i=0; i<errors.size(); ++i) {
            if (possibly(errors[i] > 0)) {
                auto error_improvement = cast_exact((errors[i]-new_errors[i])/errors[i]);
                if (error_improvement >= this->_error_refinement_minimum_improvement_percentage) {
                    has_improved = true;
                    break;
                }
            }
        }
        if (not has_improved) break;
        errors = new_errors;
        CONCLOG_PRINTLN_VAR_AT(2,errors);
    }

    diagnostic_print("final",diagnostic_refinement_iteration,phi);
    if (possibly(phi.error()>this->step_maximum_error())) {
        ARIADNE_THROW(FlowTimeStepException,"GradedTaylorPicardIntegrator::flow_step","Integration of "<<f<<" starting in "<<D<<" over time interval "<<T<<" of length "<<h<<" has error "<<phi.error()<<", which exceeds step maximum error "<<this->step_maximum_error());
    }

    return restriction(phi,dom);
}

Void GradedTaylorPicardIntegrator::_write(OutputStream& os) const {
    os << "GradedTaylorPicardIntegrator"
       << "(step_maximum_error = " << this->step_maximum_error()
       << ", function_factory = " << this->function_factory()
       << ", order = " << this->order()
       << " )";
}


} // namespace Ariadne


#include "algebra/graded.hpp"
#include "function/procedure.hpp"

namespace Ariadne {

typedef Procedure<ValidatedNumber> ValidatedProcedure;
typedef Differential<FloatDPBounds> ValidatedDifferential;
typedef Graded<ValidatedDifferential> GradedValidatedDifferential;
typedef FloatDPBounds ValidatedNumericType;
Bool operator<(const MultiIndex& a1, const MultiIndex& a2);

static const TwoExp SWEEP_THRESHOLD_RATIO=TwoExp(-10);

TaylorSeriesIntegrator::TaylorSeriesIntegrator(Sweeper<FloatDP> const& sweeper, LipschitzTolerance lip, Order ord)
    : BoundedIntegratorBase(sweeper, lip), _sweeper(sweeper), _order(ord)
{ }

TaylorSeriesIntegrator::TaylorSeriesIntegrator(StepMaximumError err, Order ord)
    : TaylorSeriesIntegrator(ThresholdSweeper<FloatDP>(DP(),err.value()*SWEEP_THRESHOLD_RATIO),DEFAULT_LIPSCHITZ_TOLERANCE,ord)
{ }


GradedTaylorSeriesIntegrator::GradedTaylorSeriesIntegrator(StepMaximumError err)
    : GradedTaylorSeriesIntegrator(err,ThresholdSweeper<FloatDP>(DP(),err.value()*SWEEP_THRESHOLD_RATIO),DEFAULT_LIPSCHITZ_TOLERANCE)
{ }

TaylorSeriesBounderIntegrator::TaylorSeriesBounderIntegrator(StepMaximumError err, Sweeper<FloatDP> const& swp, LipschitzTolerance lip, Order ord)
    : TaylorSeriesIntegrator(swp,lip,ord), _step_maximum_error(cast_exact(err.value()))
{ }

TaylorSeriesBounderIntegrator::TaylorSeriesBounderIntegrator(StepMaximumError err, Order ord)
    : TaylorSeriesIntegrator(err,ord), _step_maximum_error(cast_exact(err.value()))
{ }

GradedTaylorSeriesIntegrator::GradedTaylorSeriesIntegrator(StepMaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzTolerance lip)
    : GradedTaylorSeriesIntegrator(err,sweeper,lip,MaximumTemporalOrder(12))
{ }

GradedTaylorSeriesIntegrator::GradedTaylorSeriesIntegrator(StepMaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzTolerance lip,
                                                           MaximumTemporalOrder maxto)
    : GradedTaylorSeriesIntegrator(err,sweeper,lip,MinimumSpacialOrder(1),MinimumTemporalOrder(4),MaximumSpacialOrder(4),maxto)
{
}

GradedTaylorSeriesIntegrator::GradedTaylorSeriesIntegrator(StepMaximumError err, Sweeper<FloatDP> const& sweeper, LipschitzTolerance lip,
                                                           MinimumSpacialOrder minso, MinimumTemporalOrder minto,
                                                           MaximumSpacialOrder maxso, MaximumTemporalOrder maxto)
    : BoundedIntegratorBase(sweeper, lip), _step_maximum_error(cast_exact(err.value())), _sweeper(sweeper)
    , _minimum_spacial_order(minso), _minimum_temporal_order(minto), _maximum_spacial_order(maxso), _maximum_temporal_order(maxto)
{ }


namespace {

ExactIntervalType forwards_backwards_time_domain(ExactIntervalType domt) {
    Dyadic t0(domt.lower_bound());
    Dyadic tf(domt.upper_bound());
    return ExactIntervalType(t0-(tf-t0),tf);
}

template<class F> GradedValidatedDifferential flow(const F& f, const ExactIntervalType& c, DegreeType so, DegreeType to) {
    ValidatedDifferential x=make_differential_variable(1u,so,cast_singleton(c),0u);
    GradedValidatedDifferential y=make_graded(x);
    GradedValidatedDifferential t=create_graded(x);

    for(DegreeType n=0; n!=to; ++n) {
        t=f(y);
        y=antidifferential(t);
    }

    return y;
}

template<class X> Void append_join(Expansion<MultiIndex,X>& e, const MultiIndex& a1, const DegreeType a2, const X& c) {
    MultiIndex a(a1.size()+1);
    for(SizeType i=0; i!=a1.size(); ++i) { a[i]=a1[i]; }
    a[a1.size()]=a2;
    e.append(a,c);
}



Void graded_flow_init(const Vector<ValidatedProcedure>& f,
               Vector<GradedValidatedDifferential>& fy, List<GradedValidatedDifferential>& tmp,
               Vector<GradedValidatedDifferential>& yta,
               const Vector<ValidatedNumericType>& x, const ValidatedNumericType& t0, const Vector<ValidatedNumericType>& a,
               DegreeType so, DegreeType to)
{
    CONCLOG_SCOPE_CREATE;
    const SizeType xs=x.size();
    const SizeType as=a.size();
    const SizeType ress=f.result_size();
    const SizeType tmps=f.temporaries_size();
    const SizeType args=f.argument_size();

    ARIADNE_ASSERT(args==xs+1u+as or args==xs+as);
    const bool is_autonomous = (args==xs+as);

    CharacteristicsType<ValidatedDifferential> prs(xs+as,so,x.element_characteristics());
    GradedValidatedDifferential null(prs);
    GradedValidatedDifferential zero(0u,prs);
    fy=Vector< GradedValidatedDifferential >(ress,null);
    tmp=List< GradedValidatedDifferential >(tmps,null);
    yta=Vector< GradedValidatedDifferential >(args,zero);
    for(SizeType i=0; i!=xs; ++i) {
        yta[i]=GradedValidatedDifferential(ValidatedDifferential::variable(xs+as,so,x[i],i));
    }
    if (is_autonomous) {
        for(SizeType i=0; i!=as; ++i) {
            yta[xs+i]=GradedValidatedDifferential(ValidatedDifferential::variable(xs+as,so,a[i],xs+i));
        }
    } else {
        yta[xs]=GradedValidatedDifferential(ValidatedDifferential::constant(xs+as,so,t0));
        for(SizeType i=0; i!=as; ++i) {
            yta[xs+1u+i]=GradedValidatedDifferential(ValidatedDifferential::variable(xs+as,so,a[i],xs+i));
        }
    }
    CONCLOG_PRINTLN_AT(1,"fy="<<fy<<", tmp="<<tmp<<", yta="<<yta);
}


Void graded_flow_iterate(const Vector<ValidatedProcedure>& p,
                         Vector<GradedValidatedDifferential>& fy, List<GradedValidatedDifferential>& tmp, Vector<GradedValidatedDifferential>& yta,
                         SizeType diagnostic_call=std::numeric_limits<SizeType>::max(),
                         const char* diagnostic_branch="",
                         DegreeType diagnostic_iteration=0u)
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN_AT(1,"degree="<<yta[0].degree());
    const bool is_autonomous = (p.argument_size()==yta[0][0].argument_size());
    const SizeType n=p.result_size();

    ValidatedDifferential z=nul(yta[0][0]);
    Ariadne::compute_procedure(p,fy,tmp,yta);

    // Temporary diagnostic for the same-state second-step comparison.
    // Calls 1 and 2 are respectively the IDENTITY and QR probes.  Expose the
    // magnitude after each Procedure instruction so we can identify which
    // operation first amplifies the bounding graded differential.
    if((diagnostic_call==1u || diagnostic_call==2u)
       && diagnostic_iteration>=1u) {
        auto differential_coefficient_mag =
            [](ValidatedDifferential const& d) {
                auto r=mag(d.value());
                for(auto const& term : d.expansion()) {
                    r=max(r,mag(term.coefficient()));
                }
                return r;
            };
        auto graded_mag =
            [&](GradedValidatedDifferential const& g) {
                auto r=mag(z.value());
                for(SizeType k=0u; k!=g.size(); ++k) {
                    r=max(r,differential_coefficient_mag(g[k]));
                }
                return r;
            };
        for(SizeType j=0u; j!=tmp.size(); ++j) {
            std::cerr << "[GradedProcedureDiagnostic]"
                      << " call=" << diagnostic_call
                      << " branch=" << diagnostic_branch
                      << " iteration=" << diagnostic_iteration
                      << " instruction=" << j
                      << " op=" << p._instructions[j]
                      << " coeff_mag=" << graded_mag(tmp[j])
                      << std::endl;
        }
    }

    for(SizeType i=0; i!=n; ++i) {
        yta[i]=antidifferential(fy[i]);
    }

    if (is_autonomous) {
        for(SizeType i=n; i!=yta.size(); ++i) { yta[i].append(z); }
    } else {
        GradedValidatedDifferential& t=yta[n];
        if(t.degree()==0 && not is_autonomous) { t.append(z+1); } else { t.append(z); }
        for(SizeType i=n+1u; i!=yta.size(); ++i) { yta[i].append(z); }
    }
}






Vector<GradedValidatedDifferential>
graded_flow_differential(Vector<GradedValidatedDifferential> const& dphic, Vector<GradedValidatedDifferential> const& dphib,
                         DegreeType so, DegreeType to)
{
    const SizeType rs=dphic.size();
    const SizeType as=dphic[0][0].argument_size();
    auto z=dphic[0][0].zero_coefficient();

    Vector<GradedValidatedDifferential> gdphi(rs,GradedValidatedDifferential(List<ValidatedDifferential>(to+1u,ValidatedDifferential(as,so,z))));
    for(SizeType i=0; i!=rs; ++i) {
        for(DegreeType j=0; j!=to; ++j) {
            for(ValidatedDifferential::ConstIterator iter=dphic[i][j].begin(); iter!=dphic[i][j].end(); ++iter) {
                if(iter->index().degree()<so) { gdphi[i][j].expansion().append(iter->index(),iter->coefficient()); }
            }
            for(ValidatedDifferential::ConstIterator iter=dphib[i][j].begin(); iter!=dphib[i][j].end(); ++iter) {
                if(iter->index().degree()==so) { gdphi[i][j].expansion().append(iter->index(),iter->coefficient()); }
            }
        }
        DegreeType j=to;
        for(ValidatedDifferential::ConstIterator iter=dphib[i][j].begin(); iter!=dphib[i][j].end(); ++iter) {
            gdphi[i][j].expansion().append(iter->index(),iter->coefficient());
        }
    }
    CONCLOG_PRINTLN_AT(1,"gdphi="<<gdphi);

    return gdphi;
}

Vector<ValidatedDifferential>
differential(Vector<GradedValidatedDifferential> const& gdphi, SizeType gind,
             DegreeType so, DegreeType to)
{
    SizeType rs=gdphi.size();
    SizeType as=gdphi[0][0].argument_size();
    auto z=gdphi[0][0].zero_coefficient();

    Vector<ValidatedDifferential> dphi(rs,as+1u,so+to,z);
    MultiIndex a(as+1u);
    for(SizeType i=0; i!=rs; ++i) {
        Expansion<MultiIndex,FloatDPBounds>& component=dphi[i].expansion();
        for(DegreeType j=0; j<=to; ++j) {
            a[gind]=j;
            const Expansion<MultiIndex,FloatDPBounds>& expansion=gdphi[i][j].expansion();
            for(auto term : expansion) {
                for(SizeType k=0; k!=gind; ++k) { a[k]=term.index()[k]; }
                for(SizeType k=gind; k!=as; ++k) { a[k+1u]=term.index()[k]; }
                component.append(a,term.coefficient());
            }
        }
    }
    CONCLOG_PRINTLN_AT(1,"dphi="<<dphi);
    return dphi;
}

Vector<ValidatedDifferential>
flow_differential(Vector<GradedValidatedDifferential> const& dphic, Vector<GradedValidatedDifferential> const& dphib,
                  DegreeType so, DegreeType to)
{
    Vector<GradedValidatedDifferential> gdphi=graded_flow_differential(dphic,dphib,so,to);
    return differential(gdphi, dphic.size(),so,to);
}


FlowStepTaylorModelType make_taylor_function_model(const Vector<Differential<FloatBounds<DP>>>& df, const ExactBoxType& dom, Sweeper<FloatDP> swp) {
    ARIADNE_ASSERT(df.argument_size()==dom.dimension());
    const SizeType rs=df.size();
    const SizeType as=dom.dimension();
    const DegreeType deg = df.degree();
    FlowStepTaylorModelType tf(rs,dom,swp);

    Vector<Differential<FloatBounds<DP>>> ds=scale(Differential<FloatBounds<DP>>::variables(deg,Vector<FloatBounds<DP>>(as,dp)),dom);
    CONCLOG_PRINTLN_AT(1,"ds="<<ds<<"\rs");
    Vector<Differential<FloatBounds<DP>>> dfs = compose(df,ds);

    for(SizeType i=0; i!=rs; ++i) {
        ValidatedTaylorModelDP& model=tf.model(i);
        Expansion<MultiIndex,FloatDP>& expansion=model.expansion();
        FloatDPError& error=model.error();
        error=0u;
        expansion.reserve(dfs[i].expansion().number_of_nonzeros());

        typename Differential<FloatDPBounds>::ConstIterator iter=dfs[i].begin();
        while(iter!=dfs[i].end()) {
            MultiIndex const a=iter->index();
            FloatDPBounds coef=iter->coefficient();
            FloatDP x=coef.value();
            error+=coef.error();
            expansion.append(a,x);
            ++iter;
        }
        model.cleanup();
    }
    return tf;
}

FlowStepTaylorModelType flow_function(const Vector<Differential<FloatBounds<DP>>>& dphi, const ExactBoxType& domx, const ExactIntervalType& domt, const ExactBoxType& doma, Sweeper<FloatDP> swp) {
    StepSizeType t=static_cast<StepSizeType>(domt.lower_bound());
    StepSizeType h=static_cast<StepSizeType>(domt.upper_bound())-t;
    ExactIntervalType wdt(t-h,t+h);

    return restriction(make_taylor_function_model(dphi,join(domx,wdt,doma),swp),join(domx,domt,doma));
}

// Experimental evaluator for an affine-preconditioned vector field which keeps
// the original physical Procedure intact.  Instead of first constructing the
// dense local function g(y)=A^{-1}f(c+Ay), transform the graded arguments to
// physical coordinates, evaluate the original sparse Procedure, and transform
// only its result back to local coordinates.  This is a diagnostic path used
// to separate expression densification from the graded-differential
// representation itself.
Void graded_flow_iterate_affine_procedure(
        const Vector<ValidatedProcedure>& physical_p,
        const Vector<FloatDP>& centre,
        const Matrix<FloatDP>& A,
        const Matrix<FloatDPBounds>& inverse_A,
        Vector<GradedValidatedDifferential>& local_fy,
        Vector<GradedValidatedDifferential>& physical_fy,
        List<GradedValidatedDifferential>& physical_tmp,
        Vector<GradedValidatedDifferential>& local_yta,
        const Vector<ValidatedNumericType>* direct_physical_values,
        const char* diagnostic_branch,
        DegreeType diagnostic_iteration)
{
    const SizeType n=physical_p.result_size();
    ARIADNE_ASSERT(physical_p.argument_size()==n);
    ARIADNE_ASSERT(local_yta.size()==n);
    ARIADNE_ASSERT(centre.size()==n);
    ARIADNE_ASSERT(A.row_size()==n && A.column_size()==n);
    ARIADNE_ASSERT(inverse_A.row_size()==n && inverse_A.column_size()==n);

    Vector<GradedValidatedDifferential> physical_yta(local_yta);

    // Build c+A*y coefficient by coefficient.  Doing this directly on the
    // Differential coefficients preserves the affine dependence instead of
    // compiling it into a denser Procedure in the rotated variables.
    for(SizeType i=0u; i!=n; ++i) {
        for(SizeType k=0u; k!=local_yta[i].size(); ++k) {
            ValidatedDifferential d=nul(local_yta[0u][k]);
            if(k==0u && direct_physical_values!=nullptr) {
                // Keep the already validated physical flow-box value instead
                // of interval-mapping the local axis-aligned box back through
                // A.  Retain the local-coordinate gradient A so only the
                // zero-order box wrapping is removed by this diagnostic.
                Covector<FloatDPBounds> gradient(n,FloatDPBounds(0,dp));
                for(SizeType j=0u; j!=n; ++j) {
                    gradient[j]=FloatDPBounds(A[i][j]);
                }
                d=ValidatedDifferential::affine(
                    n,d.degree(),(*direct_physical_values)[i],gradient);
            } else {
                if(k==0u) {
                    d+=FloatDPBounds(centre[i]);
                }
                for(SizeType j=0u; j!=n; ++j) {
                    d+=local_yta[j][k]*FloatDPBounds(A[i][j]);
                }
            }
            physical_yta[i][k]=d;
        }
    }

    Ariadne::compute_procedure(
        physical_p,physical_fy,physical_tmp,physical_yta);

    // Transform only the current temporal coefficient of f back with A^{-1}.
    // local_fy stores the previous temporal coefficients and receives one new
    // coefficient on each call, exactly like compute_procedure does in the
    // ordinary graded iterator.
    for(SizeType i=0u; i!=n; ++i) {
        ValidatedDifferential gi=nul(physical_fy[0u].back());
        for(SizeType j=0u; j!=n; ++j) {
            gi+=physical_fy[j].back()*inverse_A[i][j];
        }
        local_fy[i].append(gi);
    }

    if(diagnostic_iteration>=1u) {
        auto differential_coefficient_mag =
            [](ValidatedDifferential const& d) {
                auto r=mag(d.value());
                for(auto const& term : d.expansion()) {
                    r=max(r,mag(term.coefficient()));
                }
                return r;
            };
        auto graded_mag =
            [&](GradedValidatedDifferential const& g) {
                auto r=differential_coefficient_mag(g[0u]);
                for(SizeType k=1u; k!=g.size(); ++k) {
                    r=max(r,differential_coefficient_mag(g[k]));
                }
                return r;
            };
        auto vector_mag =
            [&](Vector<GradedValidatedDifferential> const& w) {
                auto r=graded_mag(w[0u]);
                for(SizeType i=1u; i!=w.size(); ++i) {
                    r=max(r,graded_mag(w[i]));
                }
                return r;
            };
        std::cerr << "[AffineProcedureIterationDiagnostic]"
                  << " branch=" << diagnostic_branch
                  << " iteration=" << diagnostic_iteration
                  << " physical_f_coeff_mag=" << vector_mag(physical_fy)
                  << " local_f_coeff_mag=" << vector_mag(local_fy)
                  << std::endl;
    }

    for(SizeType i=0u; i!=n; ++i) {
        local_yta[i]=antidifferential(local_fy[i]);
    }
}


FlowStepTaylorModelType
graded_series_centre_polynomial_step(
        const Vector<ValidatedProcedure>& p,
        const ExactBoxType& domx,
        const ExactIntervalType& domt,
        Sweeper<FloatDP> const& sweeper,
        DegreeType so,
        DegreeType to)
{
    const SizeType n=domx.dimension();
    ExactBoxType doma;
    Vector<ValidatedNumericType> dx=cast_singleton(domx);
    Vector<ValidatedNumericType> mdx=midpoint(dx);
    Vector<ValidatedNumericType> da;
    StepSizeType t=static_cast<StepSizeType>(domt.lower_bound());
    StepSizeType h=static_cast<StepSizeType>(domt.upper_bound())-t;
    ExactIntervalType widt(t-h,t+h);
    Scalar<ValidatedNumericType> mdt=midpoint(cast_singleton(widt));

    ValidatedDifferential dzero(n,so,dx.element_characteristics());
    GradedValidatedDifferential null(0u,dzero);
    Vector<GradedValidatedDifferential> dphic(0u,null),fdphic(0u,null);
    List<GradedValidatedDifferential> tmpdphic;

    Ariadne::graded_flow_init(
        p,fdphic,tmpdphic,dphic,mdx,mdt,da,so,to);
    for(DegreeType i=0u; i!=to; ++i) {
        graded_flow_iterate(p,fdphic,tmpdphic,dphic);
    }

    // Use the centre branch for every retained coefficient, including the
    // highest temporal/spatial terms.  This is intentionally not a validated
    // flow enclosure; it is the polynomial candidate whose residual we want
    // to measure before attaching a separate validated remainder.
    Vector<ValidatedDifferential> dphi=
        flow_differential(dphic,dphic,so,to);
    return flow_function(dphi,domx,domt,doma,sweeper);
}


FlowStepTaylorModelType
graded_series_flow_step_affine_procedure(
        const Vector<ValidatedProcedure>& physical_p,
        const Vector<FloatDP>& centre,
        const Matrix<FloatDP>& A,
        const Matrix<FloatDPBounds>& inverse_A,
        const ExactBoxType& domy,
        const ExactIntervalType& domt,
        const UpperBoxType& bndy,
        const UpperBoxType& physical_bounding_box,
        Sweeper<FloatDP> const& sweeper,
        DegreeType so,
        DegreeType to)
{
    const SizeType n=domy.dimension();
    ARIADNE_PRECONDITION(physical_p.result_size()==n);
    ARIADNE_PRECONDITION(physical_p.argument_size()==n);

    Vector<ValidatedNumericType> dy=cast_singleton(domy);
    Scalar<ValidatedNumericType> dt=cast_singleton(domt);
    Vector<ValidatedNumericType> by=cast_singleton(bndy);
    Vector<ValidatedNumericType> physical_bounding_values=
        cast_singleton(physical_bounding_box);
    Vector<ValidatedNumericType> mdy=midpoint(dy);
    ExactBoxType doma;
    Vector<ValidatedNumericType> da;
    StepSizeType t=static_cast<StepSizeType>(domt.lower_bound());
    StepSizeType h=static_cast<StepSizeType>(domt.upper_bound())-t;
    ExactIntervalType widt(t-h,t+h);
    Scalar<ValidatedNumericType> mdt=midpoint(cast_singleton(widt));

    ValidatedDifferential dzero(n,so,dy.element_characteristics());
    GradedValidatedDifferential null(0u,dzero);
    Vector<GradedValidatedDifferential> dphic(0u,null),fdphic(0u,null);
    Vector<GradedValidatedDifferential> dphib(0u,null),fdphib(0u,null);
    List<GradedValidatedDifferential> unused_tmp_c,unused_tmp_b;

    Ariadne::graded_flow_init(
        physical_p,fdphic,unused_tmp_c,dphic,mdy,mdt,da,so,to);
    Ariadne::graded_flow_init(
        physical_p,fdphib,unused_tmp_b,dphib,by,dt,da,so,to);

    GradedValidatedDifferential physical_null(dphic[0u].characteristics());
    Vector<GradedValidatedDifferential> physical_f_c(n,physical_null);
    Vector<GradedValidatedDifferential> physical_f_b(n,physical_null);
    List<GradedValidatedDifferential> physical_tmp_c(
        physical_p.temporaries_size(),physical_null);
    List<GradedValidatedDifferential> physical_tmp_b(
        physical_p.temporaries_size(),physical_null);

    for(DegreeType i=0u; i!=to; ++i) {
        graded_flow_iterate_affine_procedure(
            physical_p,centre,A,inverse_A,
            fdphic,physical_f_c,physical_tmp_c,dphic,
            nullptr,"centre",i+1u);
        graded_flow_iterate_affine_procedure(
            physical_p,centre,A,inverse_A,
            fdphib,physical_f_b,physical_tmp_b,dphib,
            &physical_bounding_values,"bounding",i+1u);
    }

    Vector<ValidatedDifferential> dphi=
        flow_differential(dphic,dphib,so,to);
    FlowStepTaylorModelType tphi=
        flow_function(dphi,domy,domt,doma,sweeper);

    auto differential_coefficient_mag =
        [](ValidatedDifferential const& d) {
            auto r=mag(d.value());
            for(auto const& term : d.expansion()) {
                r=max(r,mag(term.coefficient()));
            }
            return r;
        };
    auto graded_vector_mag =
        [&](Vector<GradedValidatedDifferential> const& w) {
            auto r=differential_coefficient_mag(w[0u][0u]);
            for(SizeType i=0u; i!=w.size(); ++i) {
                for(SizeType k=0u; k!=w[i].size(); ++k) {
                    r=max(r,differential_coefficient_mag(w[i][k]));
                }
            }
            return r;
        };

    std::cerr << "[AffineProcedureFlowDiagnostic]"
              << " h=" << (domt.upper_bound()-domt.lower_bound())
              << " centre_dphi_coeff_mag=" << graded_vector_mag(dphic)
              << " bounding_dphi_coeff_mag=" << graded_vector_mag(dphib)
              << " tphi_errors=" << tphi.errors()
              << std::endl;

    return tphi;
}


} // namespace


// Flow step using graded differential with fixed degree
FlowStepTaylorModelType
graded_series_flow_step(const Vector<ValidatedProcedure>& f,
                        const ExactBoxType& domx, const ExactIntervalType& domt, const ExactBoxType& doma, const UpperBoxType& bndx,
                        Sweeper<FloatDP> const& sweeper, DegreeType so, DegreeType to)
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN_AT(1,"f="<<f);
    CONCLOG_PRINTLN_AT(1,"domx="<<domx<<", domt="<<domt<<", doma="<<doma<<", bndx="<<bndx);
    CONCLOG_PRINTLN_AT(1,"sweeper="<<sweeper<<", so="<<so<<", to="<<to);

    ARIADNE_PRECONDITION(f.result_size()==domx.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==domx.dimension()+doma.dimension() || f.argument_size()==domx.dimension()+1u+doma.dimension());

    SizeType nx=domx.dimension();

    StepSizeType t=static_cast<StepSizeType>(domt.lower_bound());
    StepSizeType h=static_cast<StepSizeType>(domt.upper_bound())-t;

    ExactIntervalType widt(t-h,t+h);

    Vector<ValidatedNumericType> dx=cast_singleton(domx);
    Scalar<ValidatedNumericType> dt=cast_singleton(domt);
    Scalar<ValidatedNumericType> wdt=cast_singleton(widt);
    Vector<ValidatedNumericType> da=cast_singleton(doma);
    Vector<ValidatedNumericType> bx=cast_singleton(bndx);

    Vector<ValidatedNumericType> mdx=midpoint(dx);
    Scalar<ValidatedNumericType> mdt=midpoint(widt);
    Vector<ValidatedNumericType> mda=midpoint(da);

    Vector<ValidatedNumericType> dc=midpoint(join(dx,wdt,da));
    Vector<ValidatedNumericType> db=join(bx,dt,da);
    ExactBoxType domc=join(domx,domt,doma);

    CONCLOG_PRINTLN_AT(2,"dx="<<dx<<", dt="<<dt<<", da="<<da<<", wdt="<<wdt<<", bx="<<bx);

    ValidatedDifferential dzero(domx.dimension()+doma.dimension(),so,dx.element_characteristics());
    GradedValidatedDifferential null(0u,dzero);
    Vector<GradedValidatedDifferential> dphic(0u,null),fdphic(0u,null),dphib(0u,null),fdphib(0u,null);
    List<GradedValidatedDifferential> tmpdphic,tmpdphib;

    static SizeType graded_internal_diagnostic_count=0u;
    auto differential_coefficient_mag =
        [](ValidatedDifferential const& d) {
            auto r=mag(d.value());
            for(auto const& term : d.expansion()) {
                r=max(r,mag(term.coefficient()));
            }
            return r;
        };
    auto graded_coefficient_mag =
        [&](Vector<GradedValidatedDifferential> const& v) {
            // The fy vectors returned by graded_flow_init contain empty
            // Graded elements until the first graded_flow_iterate evaluates
            // the vector field.  Start from a genuine zero Differential and
            // skip empty graded elements instead of dereferencing [0].
            auto r=mag(dzero.value());
            for(SizeType i=0u; i!=v.size(); ++i) {
                for(SizeType k=0u; k!=v[i].size(); ++k) {
                    r=max(r,differential_coefficient_mag(v[i][k]));
                }
            }
            return r;
        };

    Ariadne::graded_flow_init(f,fdphic,tmpdphic,dphic,mdx,mdt,mda,so,to);
    Ariadne::graded_flow_init(f,fdphib,tmpdphib,dphib,bx,dt,da,so,to);

    if(graded_internal_diagnostic_count<8u) {
        std::cerr << "[GradedIterationDiagnostic]"
                  << " call=" << graded_internal_diagnostic_count
                  << " iteration=init"
                  << " fdphic_coeff_mag=" << graded_coefficient_mag(fdphic)
                  << " fdphib_coeff_mag=" << graded_coefficient_mag(fdphib)
                  << " dphic_coeff_mag=" << graded_coefficient_mag(dphic)
                  << " dphib_coeff_mag=" << graded_coefficient_mag(dphib)
                  << std::endl;
    }

    for(DegreeType i=0; i!=to; ++i) {
        Ariadne::graded_flow_iterate(
            f,fdphic,tmpdphic,dphic,
            graded_internal_diagnostic_count,"centre",i+1u);
        Ariadne::graded_flow_iterate(
            f,fdphib,tmpdphib,dphib,
            graded_internal_diagnostic_count,"bounding",i+1u);
        if(graded_internal_diagnostic_count<8u) {
            std::cerr << "[GradedIterationDiagnostic]"
                      << " call=" << graded_internal_diagnostic_count
                      << " iteration=" << (i+1u)
                      << " fdphic_coeff_mag=" << graded_coefficient_mag(fdphic)
                      << " fdphib_coeff_mag=" << graded_coefficient_mag(fdphib)
                      << " dphic_coeff_mag=" << graded_coefficient_mag(dphic)
                      << " dphib_coeff_mag=" << graded_coefficient_mag(dphib)
                      << std::endl;
        }
    }
    CONCLOG_PRINTLN_AT(3,"dphic="<<dphic);
    CONCLOG_PRINTLN_AT(3,"dphib="<<dphib);

    dphic=project(dphic,range(0,nx));
    dphib=project(dphib,range(0,nx));
    CONCLOG_PRINTLN_AT(3,"dphic="<<dphic);
    CONCLOG_PRINTLN_AT(3,"dphib="<<dphib);

    Vector<ValidatedDifferential> dphi=flow_differential(dphic,dphib,so,to);
    CONCLOG_PRINTLN_AT(2,"dphi="<<dphi);

    FlowStepTaylorModelType tphi=flow_function(dphi,domx,domt,doma,sweeper);

    if(graded_internal_diagnostic_count<8u) {
        auto differential_vector_mag =
            [&](Vector<ValidatedDifferential> const& v) {
                auto r=differential_coefficient_mag(v[0u]);
                for(SizeType i=1u; i!=v.size(); ++i) {
                    r=max(r,differential_coefficient_mag(v[i]));
                }
                return r;
            };

        std::cerr << "[GradedInternalDiagnostic]"
                  << " call=" << graded_internal_diagnostic_count
                  << " domx=" << domx
                  << " domt=" << domt
                  << " bndx=" << bndx
                  << " fdphic_coeff_mag=" << graded_coefficient_mag(fdphic)
                  << " fdphib_coeff_mag=" << graded_coefficient_mag(fdphib)
                  << " dphic_coeff_mag=" << graded_coefficient_mag(dphic)
                  << " dphib_coeff_mag=" << graded_coefficient_mag(dphib)
                  << " dphi_coeff_mag=" << differential_vector_mag(dphi)
                  << " tphi_errors=" << tphi.errors()
                  << std::endl;
        ++graded_internal_diagnostic_count;
    }

    CONCLOG_PRINTLN("phi="<<tphi);

    return tphi;
}

// Flow step using graded differential with varying degree and specified maximum error
FlowStepModelType
graded_series_flow_step(const Vector<ValidatedProcedure>& f,
                        const ExactBoxType& domx, const ExactIntervalType& domt, const ExactBoxType& doma, const UpperBoxType& bndx,
                        ExactDouble max_err, Sweeper<FloatDP> const& sweeper, DegreeType init_so, DegreeType init_to, DegreeType max_so, DegreeType max_to)
{
    CONCLOG_SCOPE_CREATE;
    CONCLOG_PRINTLN_AT(1,"f="<<f);
    CONCLOG_PRINTLN_AT(1,"domx="<<domx<<", domt="<<domt<<", doma="<<doma<<", bndx="<<bndx);
    CONCLOG_PRINTLN_AT(1,"max_err="<<max_err<<", sweeper="<<sweeper<<", "<<
                            "init_so="<<init_so<<", init_to="<<init_to<<", max_so="<<max_so<<", max_to="<<max_to);

    ARIADNE_PRECONDITION(f.result_size()==domx.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==domx.dimension()+doma.dimension() || f.argument_size()==domx.dimension()+1u+doma.dimension());

    static const ExactDouble TRY_SPACIAL_ORDER_INCREASE_FACTOR=4;

    DegreeType so=init_so;
    DegreeType to=init_to;

    FlowStepTaylorModelType phi=graded_series_flow_step(f,domx,domt,doma,bndx, sweeper,so,to);

    CONCLOG_PRINTLN_AT(1,"phi="<<phi);
    SizeType nnz=0; for(SizeType i=0; i!=phi.size(); ++i) { nnz+=phi.model(i).number_of_nonzeros(); }
    CONCLOG_PRINTLN_AT(1,"so="<<so<<" to="<<to<<" nnz="<<nnz<<" err="<<phi.error());

    FloatDPError old_error=phi.error()*FloatDPError(TRY_SPACIAL_ORDER_INCREASE_FACTOR*two,dp);

    while(phi.error().raw()>max_err && (so<max_so || to<max_to) ) {

        // Compare the current error with the error obtained at the previous
        // order.  The old code overwrote old_error immediately before this
        // comparison, making the spatial-order branch unconditionally true
        // for every positive error.
        bool const insufficient_temporal_improvement =
            (phi.error()*FloatDPError(TRY_SPACIAL_ORDER_INCREASE_FACTOR,dp)).raw() > old_error.raw();

        if(so<max_so && (to>=max_to || insufficient_temporal_improvement)) {
            // Try increasing spatial degree and restart the temporal order.
            ++so;
            to=init_to;
        } else if(to<max_to) {
            ++to;
        } else {
            break;
        }

        old_error = phi.error();
        phi=graded_series_flow_step(f,domx,domt,doma,bndx, sweeper,so,to);

        CONCLOG_PRINTLN_AT(2,"so="<<so<<" to="<<to<<" err="<<phi.error());
    }
    CONCLOG_PRINTLN("phi="<<phi);
    return static_cast<ValidatedVectorMultivariateFunctionPatch>(phi);
}


FlowStepModelType
graded_series_flow_step(const ValidatedVectorMultivariateFunction& f,
                        const ExactBoxType& domx, const Interval<StepSizeType>& domt, const ExactBoxType& doma, const UpperBoxType& bndx,
                        ExactDouble max_err, Sweeper<FloatDP> const& sweeper, DegreeType init_so, DegreeType init_to, DegreeType max_so, DegreeType max_to)
{
    Vector<ValidatedProcedure> p(f);
    ExactIntervalType idomt(domt);

    return graded_series_flow_step(p,domx,idomt,doma,bndx,max_err,sweeper,init_so,init_to, max_so, max_to);
}


// FIXME: Should not be necessary, as should be able to construct FloatBounds<DP> from (Float<DP>,DP)
FloatBounds<DoublePrecision> cast_singleton(ExactIntervalType const& ivl, DoublePrecision pr) {
    return FloatBounds<DoublePrecision>(ivl.lower_bound(),ivl.upper_bound()); }


namespace {
// Compute the midpoint of x, and add the error to e
template<class F, class FE> F med(Bounds<F> const& x, Error<FE>& e) {
    e+=x.error(); return x.value(); }
} // namespace


template<class FLT> ValidatedVectorMultivariateTaylorFunctionModel<FLT>
make_taylor_function_model(ExactBoxType domain, Vector<Differential<Bounds<FLT>>> centre_derivatives, Vector<Differential<Bounds<FLT>>> derivative_ranges, Sweeper<FLT> swp) {
    ARIADNE_ASSERT(centre_derivatives.result_size()==derivative_ranges.result_size());
    ARIADNE_ASSERT(centre_derivatives.argument_size()==domain.dimension());
    ARIADNE_ASSERT(derivative_ranges.argument_size()==domain.dimension());
    ARIADNE_ASSERT(centre_derivatives.degree()==derivative_ranges.degree() or centre_derivatives.degree()+1u==derivative_ranges.degree());

    using PR = PrecisionType<FLT>;
    using X = FloatBounds<PR>;
    PR pr = swp.precision();

    const SizeType m=derivative_ranges.result_size();
    const SizeType n=derivative_ranges.argument_size();
    const DegreeType deg = derivative_ranges.degree();
    FloatBounds<PR> z=derivative_ranges.zero_element().zero_coefficient();

    auto scalings = Vector<Differential<FloatBounds<PR>>>(n,[&](SizeType i){return Differential<X>::variable(n,deg,z,i)*rad(domain[i]);});
    auto scaled_centre_derivatives = compose(centre_derivatives,scalings);
    auto scaled_derivative_ranges = compose(derivative_ranges,scalings);


    // Make the models
    ValidatedVectorMultivariateTaylorFunctionModel<FLT> tf(m,domain,swp);

    for(SizeType i=0; i!=m; ++i) {
        Differential<FloatBounds<PR>> const& dc = scaled_centre_derivatives[i];
        Differential<FloatBounds<PR>> const& dr = scaled_derivative_ranges[i];
        ValidatedTaylorModel<FLT>& model=tf.model(i);

        Expansion<MultiIndex,FloatDP>& expansion=model.expansion();
        FloatError<PR>& error=model.error();
        error=0u;
        expansion.reserve(centre_derivatives[i].expansion().number_of_nonzeros());

        auto riter=dr.begin();
        FloatBounds<PR> coef(pr);

        // Since coefficients are stored in increasing total degree, can first do centre and then ranges
        for(auto centre_iter=dc.begin(); centre_iter!=dc.end() && centre_iter->index().degree()<deg; ++centre_iter) {
            expansion.append(centre_iter->index(),med(centre_iter->coefficient(),error));
        }

        if(not (dr.expansion().empty() or dr.expansion().back().index().degree()<deg)) {
            auto range_iter=dr.begin(); while (range_iter->index().degree()<deg) { ++range_iter; }
            for ( ; range_iter!=dr.end(); ++range_iter) {
                expansion.append(range_iter->index(),med(range_iter->coefficient(),error));
            }
        }

        model.cleanup();
    }
    return tf;
}


template<class FLT> ValidatedVectorMultivariateTaylorFunctionModel<FLT>
make_taylor_flow_function_model(ExactBoxType domain, Vector<Differential<Bounds<FLT>>> centre_derivatives, Vector<Differential<Bounds<FLT>>> derivative_ranges, Sweeper<FLT> sweeper) {
    SizeType nx=derivative_ranges.size();
    ExactIntervalType domt=domain[nx];
    domain[nx]=forwards_backwards_time_domain(domt);
    ValidatedVectorMultivariateTaylorFunctionModel<FLT> phi=make_taylor_function_model(domain,centre_derivatives,derivative_ranges,sweeper);
    domain[nx]=domt;
    return restriction(phi,domain);
}

FlowStepModelType
series_flow_step(const ValidatedVectorMultivariateFunction& f,
                 const ExactBoxType& domx,
                 const ExactIntervalType& domt,
                 const ExactBoxType& doma,
                 const UpperBoxType& bndbx,
                 Vector<Differential<Bounds<FloatDP>>> cdphi,
                 DegreeType deg,
                 Sweeper<FloatDP> swp,
                 Nat verbosity=0)
{
    using FLT=FloatDP;
    using X=Bounds<FLT>;
    auto pr=swp.precision();
    bool is_autonomous = (f.argument_size()==domx.size()+doma.size());

    auto fbdomt = forwards_backwards_time_domain(domt);
    ExactBoxType domxta=product(domx,fbdomt,doma);

    Vector<X> bndx=cast_singleton(bndbx);
    Vector<X> rngx=cast_singleton(domx,pr);
    Scalar<X> rngt=cast_singleton(domt,pr);
    Vector<X> rnga=cast_singleton(doma,pr);

    Vector<Differential<X>> rngdf
        = is_autonomous ? f.differential(join(bndx,rnga),deg) : f.differential(join(bndx,rngt,rnga),deg);
    Vector<Differential<X>> rngdphi
        = is_autonomous ? flow(rngdf, bndx,rnga) : flow(rngdf, bndx,rngt,rnga);

    FlowStepModelType phi = make_taylor_function_model(domxta, cdphi, rngdphi, swp);
    domxta[domx.size()]=domt;
    phi=restriction(phi,domxta);
    return phi;
}

// Solve \f$\dt{\phi}(x,t,a)=f(\phi(x,t),t,a)\f$ for x in domx, t in domt, and a in doma, assuming x remains in bndx.
FlowStepModelType
series_flow_step(const ValidatedVectorMultivariateFunction& f,
                 const ExactBoxType& domx,
                 const ExactIntervalType& domt,
                 const ExactBoxType& doma,
                 const UpperBoxType& bndbx,
                 DegreeType deg,
                 Sweeper<FloatDP> swp)
{
    ARIADNE_PRECONDITION(f.result_size()==domx.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==domx.dimension()+doma.dimension()
                            || f.argument_size()==domx.dimension()+domt.dimension()+doma.dimension());
    const bool is_autonomous = f.argument_size()==domx.dimension()+doma.dimension();

    typedef DoublePrecision PR;
    typedef FloatBounds<PR> X;
    PR pr;

    // Extend time domain from [t:t+h] to [t-h:t+h]
    auto wide_domt = forwards_backwards_time_domain(domt);
    ExactBoxType domxta=product(domx,wide_domt,doma);

    Vector<X> cx(midpoint(domx),pr);
    X t0(domt.lower_bound(),pr);
    Vector<X> ca(midpoint(doma),pr);
    Vector<Differential<X>> cdf = is_autonomous ? f.differential(join(cx,ca),deg) : f.differential(join(cx,t0,ca),deg);
    Vector<Differential<X>> centre_flow_derivatives = is_autonomous ? flow(cdf, cx,ca) : flow(cdf, cx,t0,ca);

    Vector<X> bndx=cast_singleton(bndbx);
    Vector<X> rngx=cast_singleton(domx,pr);
    X rngt=cast_singleton(domt,pr);
    Vector<X> rnga=cast_singleton(doma,pr);
    Vector<Differential<X>> rngdf = is_autonomous ? f.differential(join(bndx,rnga),deg) : f.differential(join(bndx,rngt,rnga),deg);
    Vector<Differential<X>> range_flow_derivatives = is_autonomous ? flow(rngdf, bndx,rnga) : flow(rngdf, bndx,rngt,rnga);

    FlowStepModelType forwards_backwards_taylor_function_model = make_taylor_function_model(domxta, centre_flow_derivatives, range_flow_derivatives, swp);
    domxta[domx.size()]=ExactIntervalType(domt);
    FlowStepModelType forwards_taylor_function_model=restriction(forwards_backwards_taylor_function_model,domxta);
    return forwards_taylor_function_model;
}



FlowStepModelType
TaylorSeriesIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& domx, const StepSizeType& h, const UpperBoxType& bndx) const
{
    Interval<StepSizeType> domt(0,h);
    ExactBoxType doma;

    return this->flow_step(f, domx,domt,doma, bndx);
}

FlowStepModelType
TaylorSeriesIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& domx, const Interval<StepSizeType>& rngt, const ExactBoxType& doma, const UpperBoxType& bndx) const
{
    ExactIntervalType domt(rngt);
    FlowStepModelType tphi=Ariadne::series_flow_step(f,domx,domt,doma,bndx,this->order(),this->sweeper());
    return tphi;
}

Void TaylorSeriesIntegrator::_write(OutputStream& os) const {
    os << "TaylorSeriesIntegrator"
       << "( function_factory = " << this->function_factory()
       << ", sweeper = " << this->sweeper()
       << ", order = " << this->order()
       << " )";
}

namespace {
template<class... DS> inline decltype(auto) differential_flow(DS const& ... ds) { return flow(ds...); }
}

FlowStepModelType
TaylorSeriesBounderIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f,
                                          const ExactBoxType& domx, const Suggestion<StepSizeType>& hsug) const
{
    Dyadic max_err=Dyadic(this->_step_maximum_error);
    auto deg=this->order();
    auto swp=this->sweeper();

    StepSizeType h=hsug.suggestion();
    ExactIntervalType domt(0,h);
    ExactBoxType doma;

    ARIADNE_PRECONDITION(f.result_size()==domx.dimension());
    ARIADNE_PRECONDITION(f.argument_size()==domx.dimension()+doma.dimension()
                            || f.argument_size()==domx.dimension()+domt.dimension()+doma.dimension());
    const bool is_autonomous = f.argument_size()==domx.dimension()+doma.dimension();

    typedef DoublePrecision PR;
    typedef FloatBounds<PR> X;
    PR pr;

    Vector<X> cx(midpoint(domx),pr);
    X ct(domt.lower_bound(),pr);
    Vector<X> ca(midpoint(doma),pr);
    Vector<Differential<X>> cdf = is_autonomous ? f.differential(join(cx,ca),deg) : f.differential(join(cx,ct,ca),deg);
    Vector<Differential<X>> cdphi = is_autonomous ? differential_flow(cdf, cx,ca) : differential_flow(cdf, cx,ct,ca);

    // Widen domain by doubling size
    ExactBoxType wdomx=cast_exact_box(domx+(cast_singleton(domx)-domx.midpoint()));
    // Widen time-step to h*3/2
    StepSizeType wh=h+hlf(h);
    ExactIntervalType wdomt(0,wh);

    // Below is not needed, just useful to see minimum possible range
    //   ExactBoxType domxta=join(domx,domt,doma);
    //   FlowStepModelType phi=make_taylor_function_model(domxta,cdphi,cdphi,swp);

    ExactBoxType domxta=join(domx,domt,doma);
    ExactBoxType wdomxta=join(wdomx,wdomt,doma);
    FlowStepModelType phi=make_taylor_function_model(wdomxta,cdphi,cdphi,swp);
    UpperBoxType bndbx=phi.range();

    Bool refined = false;
    Bool accurate = false;
    phi=series_flow_step(f,domx,domt,doma,bndbx,cdphi, deg,swp);
    UpperBoxType xrng=phi.range();

    if (refines(xrng,bndbx)) {
        refined=true;
    } else {
        bndbx=cast_exact_box(xrng+(cast_singleton(domx)-domx.midpoint()));
        // bndbx=hull(bndbx,rphi.range());
    }

    while (not refined or not accurate) {
        phi=series_flow_step(f,domx,domt,doma,bndbx,cdphi, deg,swp);
        if (refines(xrng,bndbx)) {
            refined=true;
            bndbx=xrng;
            if (definitely(phi.error()<max_err)) { accurate=true; return phi; }
        } else {
            bndbx=cast_exact_box(xrng+(cast_singleton(domx)-domx.midpoint()));
            //bndbx=hull(bndbx,xrng);
        }

        StepSizeType nh=hlf(wh);
        wh=h;
        h=nh;
        domt=IntervalDomainType(0,h);
    }

    return phi;
}

Void TaylorSeriesBounderIntegrator::_write(OutputStream& os) const {
    os << "TaylorSeriesBounderIntegrator"
       << "( function_factory = " << this->function_factory()
       << ", step_maximum_error = " << this->_step_maximum_error
       << ", sweeper = " << this->sweeper()
       << ", order = " << this->order()
       << " )";
}

FlowStepModelType
GradedTaylorSeriesIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& domx, const StepSizeType& h, const UpperBoxType& bndx) const
{
    Interval<StepSizeType> domt(0,h);
    ExactBoxType doma;

    return this->flow_step(f, domx,domt,doma, bndx);
}

FlowStepModelType
GradedTaylorSeriesIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& domx, const Interval<StepSizeType>& rngt, const ExactBoxType& doma, const UpperBoxType& bndx) const
{
    ExactIntervalType domt(rngt);
    ExactDouble max_err=this->step_maximum_error();

    DegreeType init_so=this->minimum_spacial_order();
    DegreeType init_to=this->minimum_temporal_order();
    DegreeType max_so=this->maximum_spacial_order();
    DegreeType max_to=this->maximum_temporal_order();

    Vector<ValidatedProcedure> p(f);
    FlowStepModelType tphi=Ariadne::graded_series_flow_step(p,domx,domt,doma,bndx,
        max_err,this->sweeper(), init_so,init_to,max_so,max_to);

    if (possibly(tphi.error()>this->step_maximum_error())) {
        ARIADNE_THROW(FlowTimeStepException,"GradedTaylorSeriesIntegrator::flow_step",
                      "Integration of "<<f<<" over "<<domx<<" for time interval "<<domt<<" has error "<<tphi.errors()<<
                      " using spacial order "<<max_so<<" and temporal order "<<max_to<<
                      ", which exceeds maximum single-step error "<<max_err);
    }

    return tphi;
}

Void GradedTaylorSeriesIntegrator::_write(OutputStream& os) const {
    os << "GradedTaylorSeriesIntegrator"
       << "( function_factory = " << this->function_factory()
       << ", step_maximum_error = " << this->step_maximum_error()
       << ", sweeper = " << this->sweeper()
       << ", minimum_spacial_order = " << this->minimum_spacial_order()
       << ", minimum_temporal_order = " << this->minimum_temporal_order()
       << ", maximum_temporal_order = " << this->maximum_temporal_order()
       << ", maximum_spacial_order = " << this->maximum_spacial_order()
       << " )";
}

Void PreconditionedGradedTaylorSeriesIntegrator::_write(OutputStream& os) const {
    os << "PreconditionedGradedTaylorSeriesIntegrator"
       << "( function_factory = " << this->function_factory()
       << ", step_maximum_error = " << this->step_maximum_error()
       << ", sweeper = " << this->sweeper()
       << ", minimum_spacial_order = " << this->minimum_spacial_order()
       << ", minimum_temporal_order = " << this->minimum_temporal_order()
       << ", maximum_temporal_order = " << this->maximum_temporal_order()
       << ", maximum_spacial_order = " << this->maximum_spacial_order()
       << ", preconditioning = "
       << (this->preconditioning()==TaylorSeriesPreconditioning::QR ? "QR" : "IDENTITY")
       << " )";
}


PreconditionedTaylorSeriesState
PreconditionedGradedTaylorSeriesIntegrator::precondition(
        const ValidatedVectorMultivariateFunctionPatch& state) const
{
    auto const& state_taylor=
        dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(
            state.reference());

    SizeType const n=state_taylor.size();
    auto const& factory=this->function_factory();

    Vector<FloatDP> centre(n,FloatDP(dp));
    ValidatedVectorMultivariateFunctionPatch centred=
        factory.create_zeros(n,state.domain());
    for(SizeType i=0u; i!=n; ++i) {
        centre[i]=state_taylor.model(i).value().raw();
        centred[i]=state[i]-FloatDPBounds(centre[i]);
    }

    // QR preconditioning follows the linear part of the current local initial
    // Taylor map.  Flow* obtains an orthogonal matrix from these linear
    // coefficients; the same idea is used here.  If the parameter dimension
    // does not match the state dimension, retain the identity orientation.
    Matrix<FloatDP> rotation=Matrix<FloatDP>::identity(n,dp);
    if(this->preconditioning()==TaylorSeriesPreconditioning::QR
        && state_taylor.argument_size()==n) {
        // Extract the first-order coefficients directly.  Calling the
        // jacobian_value template here would require a FloatDP instantiation
        // that is not exported by the algebra library.
        Matrix<FloatDPApproximation> approximate_J(n,n,dp);
        MultiIndex a(n);
        for(SizeType i=0u; i!=n; ++i) {
            for(SizeType j=0u; j!=n; ++j) {
                a[j]=1u;
                approximate_J[i][j]=
                    FloatDPApproximation(state_taylor.model(i)[a]);
                a[j]=0u;
            }
        }

        // Flow* sorts the linear-coefficient columns by decreasing
        // Euclidean size before QR.  With unpivoted Gram-Schmidt this makes
        // the first orthogonal direction follow the dominant dependency.
        for(SizeType j=0u; j+1u<n; ++j) {
            SizeType jmax=j;
            FloatDPApproximation max_norm_square(0,dp);
            for(SizeType i=0u; i!=n; ++i) {
                max_norm_square+=sqr(approximate_J[i][j]);
            }
            for(SizeType k=j+1u; k!=n; ++k) {
                FloatDPApproximation norm_square(0,dp);
                for(SizeType i=0u; i!=n; ++i) {
                    norm_square+=sqr(approximate_J[i][k]);
                }
                if(decide(norm_square>max_norm_square)) {
                    max_norm_square=norm_square;
                    jmax=k;
                }
            }
            if(jmax!=j) {
                for(SizeType i=0u; i!=n; ++i) {
                    auto tmp=approximate_J[i][j];
                    approximate_J[i][j]=approximate_J[i][jmax];
                    approximate_J[i][jmax]=tmp;
                }
            }
        }

        auto const approximate_QR=orthogonal_decomposition(approximate_J);
        Matrix<FloatDPApproximation> approximate_Q=std::get<0>(approximate_QR);

        // Ariadne's public orthogonal_decomposition returns orthogonal
        // columns, but they are not required to have unit norm.  Turn them
        // into an actual rotation so that the coordinate change itself has
        // condition number one.
        for(SizeType j=0u; j!=n; ++j) {
            FloatDPApproximation norm_square(0,dp);
            for(SizeType i=0u; i!=n; ++i) {
                norm_square+=sqr(approximate_Q[i][j]);
            }
            FloatDPApproximation const column_norm=sqrt(norm_square);
            if(column_norm.raw()!=FloatDP(0,dp)) {
                for(SizeType i=0u; i!=n; ++i) {
                    approximate_Q[i][j]/=column_norm;
                }
            } else {
                // Degenerate linear direction: use the identity orientation
                // rather than introducing a singular preconditioner.
                approximate_Q=Matrix<FloatDPApproximation>::identity(n,dp);
                break;
            }
        }

        rotation=
            reinterpret_cast<Matrix<FloatDP> const&>(approximate_Q);
    }

    Matrix<FloatDPBounds> const inverse_rotation=inverse(rotation);
    ValidatedVectorMultivariateFunctionPatch rotated=
        factory.create_zeros(n,state.domain());
    for(SizeType i=0u; i!=n; ++i) {
        for(SizeType j=0u; j!=n; ++j) {
            rotated[i]=rotated[i]+centred[j]*inverse_rotation[i][j];
        }
    }

    // Do not explicitly rescale the rotated variables to [-1,1].  A Taylor
    // FunctionPatch in Ariadne is already internally scaled to the unit box
    // of its external domain.  Repeating that scaling here divides the
    // accumulated remainder by the physical set radius at every step and was
    // the source of the rapid error growth seen in the QR diagnostic.
    Matrix<FloatDP> linear_map=rotation;
    ExactBoxType local_domain=
        cast_exact_box(widen(rotated.range()));

    return PreconditionedTaylorSeriesState(
        std::move(centre),std::move(linear_map),std::move(local_domain),
        std::move(rotated));
}

PreconditionedTaylorSeriesStep
PreconditionedGradedTaylorSeriesIntegrator::step(
        const ValidatedVectorMultivariateFunction& f,
        const PreconditionedTaylorSeriesState& state,
        const Suggestion<StepSizeType>& hsug) const
{
    SizeType const n=state.centre().size();
    ARIADNE_PRECONDITION(f.result_size()==n);
    ARIADNE_PRECONDITION(f.argument_size()==n);
    ARIADNE_PRECONDITION(state.linear_map().row_size()==n);
    ARIADNE_PRECONDITION(state.linear_map().column_size()==n);

    auto const& factory=this->function_factory();
    ExactBoxType const& domy=state.local_domain();
    Vector<FloatDP> const& centre=state.centre();
    Matrix<FloatDP> const& A=state.linear_map();
    Matrix<FloatDPBounds> const inverse_A=inverse(A);

    // First compute the flow bound in physical coordinates.  The bounder
    // deliberately enlarges the initial box; constructing the transformed
    // vector field only on domy would therefore make it invalid exactly where
    // the bounder needs to evaluate it.
    ValidatedVectorMultivariateFunctionPatch initial_y=
        factory.create_identity(domy);
    ValidatedVectorMultivariateFunctionPatch initial_x=
        factory.create_zeros(n,domy);
    for(SizeType i=0u; i!=n; ++i) {
        initial_x[i]=factory.create_constant(domy,centre[i]);
        for(SizeType j=0u; j!=n; ++j) {
            initial_x[i]=initial_x[i]+initial_y[j]*FloatDPBounds(A[i][j]);
        }
    }

    ExactBoxType const physical_initial_domain=
        cast_exact_box(widen(initial_x.range()));

    StepSizeType h;
    UpperBoxType physical_bounding_box;
    make_lpair(h,physical_bounding_box)=
        this->flow_bounds(f,physical_initial_domain,hsug);

    // Transform the validated physical flow bound to local coordinates.
    // For a future non-diagonal (QR) A this interval matrix product remains
    // conservative.
    UpperBoxType local_bounding_box(n);
    for(SizeType i=0u; i!=n; ++i) {
        FloatDPBounds yi(0,dp);
        for(SizeType j=0u; j!=n; ++j) {
            FloatDPBounds const xj=cast_singleton(physical_bounding_box[j]);
            yi=yi+inverse_A[i][j]*(xj-centre[j]);
        }
        local_bounding_box[i]=UpperIntervalType(yi.lower(),yi.upper());
    }

    // Build the transformed vector field on the whole validated local flow
    // bound, not merely on the local initial domain.
    ExactBoxType const local_vector_field_domain=
        cast_exact_box(local_bounding_box);
    ValidatedVectorMultivariateFunctionPatch y=
        factory.create_identity(local_vector_field_domain);
    ValidatedVectorMultivariateFunctionPatch x_of_y=
        factory.create_zeros(n,local_vector_field_domain);
    for(SizeType i=0u; i!=n; ++i) {
        x_of_y[i]=factory.create_constant(local_vector_field_domain,centre[i]);
        for(SizeType j=0u; j!=n; ++j) {
            x_of_y[i]=x_of_y[i]+y[j]*FloatDPBounds(A[i][j]);
        }
    }

    // y' = A^{-1} f(c+A*y).  This formulation already supports a full
    // non-diagonal A, so QR preconditioning can later reuse the same core.
    ValidatedVectorMultivariateFunctionPatch physical_vector_field=
        compose(f,x_of_y);
    ValidatedVectorMultivariateFunctionPatch local_vector_field=
        factory.create_zeros(n,local_vector_field_domain);
    for(SizeType i=0u; i!=n; ++i) {
        for(SizeType j=0u; j!=n; ++j) {
            local_vector_field[i]=local_vector_field[i]
                + physical_vector_field[j]*inverse_A[i][j];
        }
    }
    ValidatedVectorMultivariateFunction g=cast_unrestricted(local_vector_field);

    if(this->diagnostics()) {
        std::cerr << "[PreconditionedSecondStepGeometry]"
                  << " mode="
                  << (this->preconditioning()==TaylorSeriesPreconditioning::QR ? "QR" : "IDENTITY")
                  << " domy=" << domy
                  << " A=" << A
                  << " physical_initial_domain=" << physical_initial_domain
                  << " physical_bounding_box=" << physical_bounding_box
                  << " local_bounding_box=" << local_bounding_box
                  << " local_vector_field_range=" << local_vector_field.range()
                  << std::endl;
    }

    ExactBoxType doma;
    Vector<ValidatedProcedure> p(g);

    // Match BoundedIntegratorBase's suggested-step semantics: the flow bound
    // was computed for the initially suggested h and is therefore also valid
    // for every smaller h.  The step tolerance is a tolerance in the physical
    // state coordinates, not in the arbitrarily scaled local coordinates y.
    // Testing local_flow.error() here would make the accepted step depend on
    // the choice of preconditioner A (for example a small diagonal entry
    // multiplies the corresponding normalised error by A when returning to x).
    StepSizeType hprev=h*1.5_dy;
    FlowStepModelType local_flow;
    ValidatedVectorMultivariateFunctionPatch physical_local_flow;
    ValidatedVectorMultivariateFunctionPatch gronwall_physical_local_flow;
    Bool have_gronwall_flow=false;
    ExactIntervalType domt;
    while(true) {
        have_gronwall_flow=false;
        domt=ExactIntervalType(0,h);
        local_flow=Ariadne::graded_series_flow_step(
            p,domy,domt,doma,local_bounding_box,
            this->step_maximum_error(),this->sweeper(),
            this->minimum_spacial_order(),this->minimum_temporal_order(),
            this->maximum_spacial_order(),this->maximum_temporal_order());

        if(this->diagnostics()
            && this->preconditioning()==TaylorSeriesPreconditioning::QR
            && this->minimum_spacial_order()==this->maximum_spacial_order()
            && this->minimum_temporal_order()==this->maximum_temporal_order()) {
            FlowStepTaylorModelType centre_polynomial=
                graded_series_centre_polynomial_step(
                    p,domy,domt,this->sweeper(),
                    this->minimum_spacial_order(),
                    this->minimum_temporal_order());

            // Compute the ODE defect R(y,t)=dP/dt-g(P) of the centre-only
            // Taylor polynomial.  If this is already small, the remaining
            // challenge is to validate a separate remainder around P rather
            // than to propagate a full interval-valued graded recurrence.
            ValidatedVectorMultivariateFunctionPatch field_on_polynomial=
                compose(g,centre_polynomial);
            ValidatedVectorMultivariateFunctionPatch defect=
                factory.create_zeros(n,centre_polynomial.domain());
            SizeType const time_index=centre_polynomial.argument_size()-1u;
            for(SizeType i=0u; i!=n; ++i) {
                ValidatedScalarMultivariateFunctionPatch dpoly =
                    derivative(centre_polynomial.get(i),time_index);
                ValidatedScalarMultivariateFunctionPatch fpoly =
                    field_on_polynomial.get(i);
                defect[i]=dpoly-fpoly;
            }

            ValidatedVectorMultivariateFunctionPatch initial_polynomial=
                partial_evaluate(centre_polynomial,time_index,StepSizeType(0.0));
            ValidatedVectorMultivariateFunctionPatch identity_on_domy=
                factory.create_identity(domy);
            ValidatedVectorMultivariateFunctionPatch initial_defect=
                initial_polynomial-identity_on_domy;

            // Crude infinity-norm Lipschitz bound on the whole validated
            // local flow box.  Together with the defect range this is enough
            // to form the a-posteriori Gronwall estimate
            //   |e(h)| <= exp(Lh)|e(0)| + (exp(Lh)-1)/L * sup|R|.
            // Print the ingredients first; keep this diagnostic independent
            // of any particular scalar-bound implementation.
            Vector<ValidatedNumericType> local_bounding_values=
                cast_singleton(local_bounding_box);
            Vector<ValidatedDifferential> dg=
                g.differential(local_bounding_values,1u);
            auto lipschitz_inf=mag(dg[0u].gradient()[0u]);
            for(SizeType j=1u; j!=n; ++j) {
                lipschitz_inf+=mag(dg[0u].gradient()[j]);
            }
            for(SizeType i=1u; i!=n; ++i) {
                auto row_sum=mag(dg[i].gradient()[0u]);
                for(SizeType j=1u; j!=n; ++j) {
                    row_sum+=mag(dg[i].gradient()[j]);
                }
                lipschitz_inf=max(lipschitz_inf,row_sum);
            }

            std::cerr << "[CentrePolynomialDefectDiagnostic]"
                      << " h=" << h
                      << " polynomial_errors=" << centre_polynomial.errors()
                      << " polynomial_range=" << centre_polynomial.range()
                      << " defect_range=" << defect.range()
                      << " initial_defect_range=" << initial_defect.range()
                      << " lipschitz_inf=" << lipschitz_inf
                      << std::endl;

            // First rigorous polynomial+remainder prototype.  The existing
            // bounder still certifies that the exact local flow stays inside
            // local_bounding_box.  If the centre polynomial also stays in
            // that convex box, the Jacobian bound above applies on every
            // segment joining P(y,t) to the exact solution.
            auto polynomial_range=centre_polynomial.range();
            // Compare against an exact outer box, as elsewhere in the
            // integrator.  Comparing a validated range directly with an
            // UpperBox can remain indeterminate even when the printed
            // endpoints show strict containment.
            ExactBoxType const certification_box=
                cast_exact_box(local_bounding_box);
            Bool const polynomial_in_certification_box=
                definitely(subset(polynomial_range,certification_box));

            if(polynomial_in_certification_box) {
                auto defect_ranges=defect.range();
                auto initial_defect_ranges=initial_defect.range();

                // Use the monotone, fully upper-rounded estimate
                //
                //   |e(t)| <= exp(L h) ( |e(0)| + h sup|R| )
                //
                // for every t in [0,h].  This is slightly looser than the
                // exact scalar Gronwall factor (exp(Lh)-1)/L, but avoids any
                // division by a lower bound for L and is therefore a simple
                // rigorous first prototype.
                FloatDP const hraw(h,dp);
                FloatDP const Lraw=lipschitz_inf.raw();
                FloatDP const amplification_raw=
                    exp(up,mul(up,Lraw,hraw));

                Vector<Error<FloatDP>> gronwall_remainder(
                    n,Error<FloatDP>(0u,dp));
                for(SizeType i=0u; i!=n; ++i) {
                    auto epsilon_i=mag(defect_ranges[i]);
                    auto initial_i=mag(initial_defect_ranges[i]);
                    FloatDP const residual_raw=
                        add(up,initial_i.raw(),
                            mul(up,hraw,epsilon_i.raw()));
                    FloatDP const remainder_raw=
                        mul(up,amplification_raw,residual_raw);
                    gronwall_remainder[i]=Error<FloatDP>(remainder_raw);
                    centre_polynomial.model(i).set_error(
                        centre_polynomial.model(i).error()
                        +gronwall_remainder[i]);
                }

                ValidatedVectorMultivariateFunctionPatch
                    physical_gronwall_polynomial=
                        factory.create_zeros(n,centre_polynomial.domain());
                for(SizeType i=0u; i!=n; ++i) {
                    physical_gronwall_polynomial[i]=
                        factory.create_constant(
                            centre_polynomial.domain(),centre[i]);
                    for(SizeType j=0u; j!=n; ++j) {
                        physical_gronwall_polynomial[i]=
                            physical_gronwall_polynomial[i]
                            +centre_polynomial[j]
                                *FloatDPBounds(A[i][j]);
                    }
                }

                gronwall_physical_local_flow=physical_gronwall_polynomial;
                have_gronwall_flow=true;

                std::cerr << "[GronwallPolynomialPrototype]"
                          << " h=" << h
                          << " amplification=" << amplification_raw
                          << " local_remainder=" << gronwall_remainder
                          << " local_errors=" << centre_polynomial.errors()
                          << " physical_errors="
                          << physical_gronwall_polynomial.errors()
                          << " physical_error="
                          << physical_gronwall_polynomial.error()
                          << " physical_range="
                          << physical_gronwall_polynomial.range()
                          << std::endl;
            } else {
                std::cerr << "[GronwallPolynomialPrototype]"
                          << " h=" << h
                          << " rejected=polynomial_outside_certification_box"
                          << " polynomial_range=" << polynomial_range
                          << " certification_box=" << certification_box
                          << std::endl;
            }
        }

        // Return to physical coordinates before deciding whether the local
        // approximation satisfies StepMaximumError.  This is the quantity
        // corresponding to the flow model returned by ordinary integrators.
        // Controlled comparison for the QR investigation.  Evaluate the same
        // local step with the original sparse physical Procedure and perform
        // c+A*y / A^{-1} only on graded Differential objects.  The returned
        // model is diagnostic only; acceptance below still uses the existing
        // dense transformed Procedure so this experiment does not alter the
        // integration semantics.
        if(this->diagnostics()
            && this->preconditioning()==TaylorSeriesPreconditioning::QR
            && this->minimum_spacial_order()==this->maximum_spacial_order()
            && this->minimum_temporal_order()==this->maximum_temporal_order()) {
            Vector<ValidatedProcedure> physical_p(f);
            FlowStepTaylorModelType affine_procedure_flow=
                graded_series_flow_step_affine_procedure(
                    physical_p,centre,A,inverse_A,
                    domy,domt,local_bounding_box,physical_bounding_box,
                    this->sweeper(),
                    this->minimum_spacial_order(),
                    this->minimum_temporal_order());
            std::cerr << "[AffineProcedureComparison]"
                      << " h=" << h
                      << " dense_local_errors=" << local_flow.errors()
                      << " affine_procedure_errors="
                      << affine_procedure_flow.errors()
                      << std::endl;
        }

        physical_local_flow=
            factory.create_zeros(n,local_flow.domain());
        for(SizeType i=0u; i!=n; ++i) {
            physical_local_flow[i]=
                factory.create_constant(local_flow.domain(),centre[i]);
            for(SizeType j=0u; j!=n; ++j) {
                physical_local_flow[i]=physical_local_flow[i]
                    + local_flow[j]*FloatDPBounds(A[i][j]);
            }
        }

        if(this->diagnostics()) {
            std::cerr << "[PreconditionedSecondStepCandidate]"
                      << " mode="
                      << (this->preconditioning()==TaylorSeriesPreconditioning::QR ? "QR" : "IDENTITY")
                      << " h=" << h
                      << " local_errors=" << local_flow.errors()
                      << " local_error=" << local_flow.error()
                      << " physical_errors=" << physical_local_flow.errors()
                      << " physical_error=" << physical_local_flow.error()
                      << " local_flow_range=" << local_flow.range()
                      << std::endl;
        }

        Bool const gronwall_acceptable=
            have_gronwall_flow
            && definitely(gronwall_physical_local_flow.error()
                          <=this->step_maximum_error());
        Bool const dense_acceptable=
            definitely(physical_local_flow.error()
                       <=this->step_maximum_error());

        if(this->diagnostics() && have_gronwall_flow) {
            std::cerr << "[GronwallAcceptanceComparison]"
                      << " h=" << h
                      << " gronwall_physical_error="
                      << gronwall_physical_local_flow.error()
                      << " dense_physical_error="
                      << physical_local_flow.error()
                      << " gronwall_acceptable=" << gronwall_acceptable
                      << " dense_acceptable=" << dense_acceptable
                      << std::endl;
        }

        if(gronwall_acceptable) {
            // From this point on, propagate the polynomial+remainder flow.
            // The dense graded flow above remains computed only to provide an
            // A/B diagnostic during this experiment.
            physical_local_flow=gronwall_physical_local_flow;
            break;
        }

        // If the certification guard could not build a Gronwall candidate,
        // retain the old validated path as a safety fallback for now.
        if(!have_gronwall_flow && dense_acceptable) {
            break;
        }

        StepSizeType const hnew=hlf(hprev);
        hprev=h;
        h=StepSizeType(hnew.get_d());
        CONCLOG_PRINTLN_AT(1,
            "PreconditionedGradedTaylorSeriesIntegrator reduced h to "<<h);
    }

    // Compose exactly once with the normalised local-initial-set TM y(s).
    // The resulting flowpipe is parameterised by the original enclosure
    // parameters s plus the local time variable.
    ExactBoxType const parameter_domain=state.parameter_domain();
    ExactBoxType const flowpipe_domain=product(parameter_domain,domt);
    ValidatedVectorMultivariateFunctionPatch embedded_mapping=
        embed(state.normalised_mapping(),domt);
    ValidatedScalarMultivariateFunctionPatch time_coordinate=
        factory.create_coordinate(flowpipe_domain,flowpipe_domain.size()-1u);
    ValidatedVectorMultivariateFunctionPatch arguments=
        join(embedded_mapping,time_coordinate);
    ValidatedVectorMultivariateFunctionPatch flowpipe_mapping=
        compose(physical_local_flow,arguments);

    // For the evolved set, evaluate time before composing with the local
    // initial Taylor model.  Composing the complete space-time flowpipe first
    // and only then evaluating t=h introduces unnecessary mixed space/time
    // terms and substantially larger sweep/remainder errors.  This also
    // matches the TM-integration update X_{l+1}=p_l(X_l,delta_l)+I_l.
    ValidatedVectorMultivariateFunctionPatch local_endpoint=
        partial_evaluate(
            physical_local_flow,
            physical_local_flow.argument_size()-1u,h);
    ValidatedVectorMultivariateFunctionPatch evolved_mapping=
        compose(local_endpoint,state.normalised_mapping());

    // Preserve the two-layer TM representation across steps.  Precondition
    // the fresh local endpoint Phi_l(y,h) first, while its remainder is still
    // a single-step remainder, and only then compose the resulting local
    // coordinate map with the accumulated y_l(s).  Re-preconditioning the
    // already-composed physical map rotates its axis-aligned accumulated
    // remainder at every step and causes an artificial wrapping explosion.
    PreconditionedTaylorSeriesState local_transition=
        this->precondition(local_endpoint);
    ValidatedVectorMultivariateFunctionPatch next_normalised_mapping=
        compose(
            local_transition.normalised_mapping(),
            state.normalised_mapping());
    ExactBoxType next_local_domain=
        cast_exact_box(widen(next_normalised_mapping.range()));

    PreconditionedTaylorSeriesState final_state(
        local_transition.centre(),
        local_transition.linear_map(),
        std::move(next_local_domain),
        std::move(next_normalised_mapping));

    return PreconditionedTaylorSeriesStep(
        h,std::move(flowpipe_mapping),std::move(evolved_mapping),
        std::move(final_state));
}

FlowStepModelType
PreconditionedGradedTaylorSeriesIntegrator::flow_step(
        const ValidatedVectorMultivariateFunction& f,
        const ExactBoxType& domx,
        const StepSizeType& h,
        const UpperBoxType& bndx) const
{
    SizeType const n=domx.size();
    ARIADNE_PRECONDITION(f.result_size()==n);
    ARIADNE_PRECONDITION(f.argument_size()==n);

    Vector<FloatDP> centre(n,FloatDP(dp));
    Vector<FloatDP> radius(n,FloatDP(dp));
    for(SizeType i=0u; i!=n; ++i) {
        centre[i]=domx[i].midpoint().raw();
        radius[i]=domx[i].radius().upper().raw();
        if(radius[i]==FloatDP(0,dp)) {
            return GradedTaylorSeriesIntegrator::flow_step(f,domx,h,bndx);
        }
    }

    ExactBoxType unit_domain(n,ExactIntervalType(-1,+1));
    auto const& factory=this->function_factory();

    // x = c + R y, with diagonal R.  Build the transformed vector field
    // y' = R^{-1} f(c+Ry) as a validated function on the normalised domain.
    ValidatedVectorMultivariateFunctionPatch y=factory.create_identity(unit_domain);
    ValidatedVectorMultivariateFunctionPatch x_of_y=factory.create_zeros(n,unit_domain);
    for(SizeType i=0u; i!=n; ++i) {
        x_of_y[i]=factory.create_constant(unit_domain,centre[i])
                + y[i]*FloatDPBounds(radius[i]);
    }

    ValidatedVectorMultivariateFunctionPatch gpatch=compose(f,x_of_y);
    for(SizeType i=0u; i!=n; ++i) {
        gpatch[i]=gpatch[i]/FloatDPBounds(radius[i]);
    }
    ValidatedVectorMultivariateFunction g=cast_unrestricted(gpatch);

    UpperBoxType normalized_bounding_box(n);
    for(SizeType i=0u; i!=n; ++i) {
        FloatDPBounds const bx=cast_singleton(bndx[i]);
        FloatDPBounds const nb=(bx-centre[i])/radius[i];
        normalized_bounding_box[i]=UpperIntervalType(nb.lower(),nb.upper());
    }

    ExactIntervalType domt(0,h);
    ExactBoxType doma;
    Vector<ValidatedProcedure> p(g);
    FlowStepModelType yflow=Ariadne::graded_series_flow_step(
        p,unit_domain,domt,doma,normalized_bounding_box,
        this->step_maximum_error(),this->sweeper(),
        this->minimum_spacial_order(),this->minimum_temporal_order(),
        this->maximum_spacial_order(),this->maximum_temporal_order());

    if(possibly(yflow.error()>this->step_maximum_error())) {
        ARIADNE_THROW(FlowTimeStepException,
                      "PreconditionedGradedTaylorSeriesIntegrator::flow_step",
                      "Integration of "<<f<<" over "<<domx
                      <<" for time interval "<<domt
                      <<" has normalised-flow error "<<yflow.errors()
                      <<", which exceeds maximum single-step error "
                      <<this->step_maximum_error());
    }

    // Convert y(t) back to x(t)=c+Ry while the flow is still expressed over
    // the normalised state coordinates.
    ValidatedVectorMultivariateFunctionPatch physical_yflow=yflow;
    for(SizeType i=0u; i!=n; ++i) {
        physical_yflow[i]=factory.create_constant(yflow.domain(),centre[i])
                         + yflow[i]*FloatDPBounds(radius[i]);
    }

    // Re-express the model on the original physical state domain so callers
    // can use it exactly like any other IntegratorInterface flow model.
    ExactBoxType physical_flow_domain=product(domx,domt);
    ValidatedVectorMultivariateFunctionPatch physical_id=
        factory.create_identity(physical_flow_domain);
    ValidatedVectorMultivariateFunctionPatch normalising_arguments=
        factory.create_zeros(n+1u,physical_flow_domain);
    for(SizeType i=0u; i!=n; ++i) {
        normalising_arguments[i]=
            (physical_id[i]-FloatDPBounds(centre[i]))/FloatDPBounds(radius[i]);
    }
    normalising_arguments[n]=physical_id[n];

    return FlowStepModelType(compose(physical_yflow,normalising_arguments));
}




template<class X> Void truncate(Differential<X>& x, DegreeType spacial_order_, DegreeType temporal_order_) {
    SizeType n=x.argument_size()-1;
    typename Differential<X>::Iterator write_iter=x.begin();
    typename Differential<X>::ConstIterator read_iter=x.begin();
    while(read_iter!=x.end()) {
        UniformConstReference<MultiIndex> index = read_iter->index();
        if(index[n]>temporal_order_ || index[n]+spacial_order_<index.degree()) {
        } else {
            *write_iter=*read_iter;
            ++write_iter;
        }
        ++read_iter;
    }
    x.expansion().resize(static_cast<SizeType>(write_iter-x.begin()));
}

template<class X> Void truncate(Vector< Differential<X> >& x, DegreeType spacial_order_, DegreeType temporal_order_) {
    for(DegreeType i=0; i!=x.size(); ++i) { truncate(x[i],spacial_order_,temporal_order_); }
}

AffineIntegrator::AffineIntegrator(SpacialOrder spacial_order_, TemporalOrder temporal_order_)
    : BoundedIntegratorBase(Sweeper<FloatDP>(), lipschitz_tolerance=DEFAULT_LIPSCHITZ_TOLERANCE), _spacial_order(spacial_order_), _temporal_order(temporal_order_) { }

Vector<ValidatedDifferential>
AffineIntegrator::flow_derivative(const ValidatedVectorMultivariateFunction& f, const Vector<ValidatedNumericType>& dom) const
{
    DegreeType const deg = this->_spacial_order+this->_temporal_order;
    Vector<ValidatedDifferential> dx(f.result_size(),ValidatedDifferential(f.result_size()+1,deg,zero));
    for (SizeType i=0; i<f.result_size(); ++i) {
        dx[i]=ValidatedDifferential::variable(f.result_size()+1,deg,dom[i],i);
    }
    Vector<ValidatedDifferential> dphi = dx;

    for(DegreeType i=0; i!=_temporal_order; ++i) {
        dphi = antiderivative(f.evaluate(dphi),dom.size())+dx;
    }
    truncate(dphi,this->_spacial_order,this->_temporal_order);
    return dphi;
}

FlowStepModelType
AffineIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& dom, const StepSizeType& h, const UpperBoxType& bbox) const
{
    Vector<ValidatedNumericType> dmid = Vector<ValidatedNumericType>(midpoint(dom));

    Vector<ValidatedDifferential> mdphi = this->flow_derivative(f,dmid);
    Vector<ValidatedDifferential> bdphi = this->flow_derivative(f,cast_singleton(bbox));

    const SizeType n=dom.size();
    DoublePrecision prec;
    FloatDPError zero_err(prec);

    Vector<FloatDPError> err(n,zero_err);

    Vector<FloatDPError> rad(n+1,zero_err);
    for(SizeType i=0; i!=n; ++i) {
        rad[i] = cast_positive(max(dom[i].upper_bound()-dmid[i].lower(),dmid[i].upper()-dom[i].lower_bound()));
    }
    rad[n] = abs(h);

    for(SizeType i=0; i!=n; ++i) {
        for(Expansion<MultiIndex,ValidatedNumericType>::ConstIterator iter=bdphi[i].begin(); iter!=bdphi[i].end(); ++iter) {
            UniformConstReference<MultiIndex> a=iter->index();
            if(a[n]==this->_temporal_order && a[n]+this->_spacial_order==a.degree()) {
                UniformConstReference<ValidatedNumericType> rng = iter->coefficient();
                UniformConstReference<ValidatedNumericType> mid = mdphi[i][a];
                ARIADNE_ASSERT(rng.lower().raw()<=mid.lower().raw() && mid.upper().raw()<=rng.upper().raw());
                FloatDPError mag = FloatDPError(max(rng.upper()-mid.lower(),mid.upper()-rng.lower()));
                for(SizeType j=0; j!=n+1u; ++j) { mag *= pow(rad[j],static_cast<DegreeType>(a[j])); }
                err[i] += mag;
            }
        }
    }

    ExactBoxType flow_domain = product(dom,ExactIntervalType(0,h));

    FlowStepModelType id = this->function_factory().create_identity(flow_domain);
    FlowStepModelType res = this->function_factory().create_zeros(n,flow_domain);
    for(SizeType i=0; i!=n; ++i) {
        ValidatedScalarMultivariateFunctionPatch res_model = res[i] + static_cast<ValidatedNumber>(static_cast<ValidatedNumericType>(mdphi[i].expansion()[MultiIndex::zero(n+1)]));
        for(SizeType j=0; j!=mdphi[i].argument_size()-1; ++j) {
            // TODO: Remove casts
            res_model+=static_cast<ValidatedNumber>(static_cast<ValidatedNumericType>(mdphi[i].expansion()[MultiIndex::unit(n+1,j)]))*(id[j]-ValidatedNumericType(midpoint(flow_domain[j])));
        }
        SizeType j=mdphi[i].argument_size()-1u;
        // TODO: Remove casts
        res_model+=static_cast<ValidatedNumber>(static_cast<ValidatedNumericType>(mdphi[i].expansion()[MultiIndex::unit(n+1,j)]))*id[j];
        res_model += FloatDPBounds(-err[i],+err[i]);
        res[i]=res_model;
    }
    return res;
}

FlowStepModelType
AffineIntegrator::flow_step(const ValidatedVectorMultivariateFunction& f, const ExactBoxType& D, const Interval<StepSizeType>& T, const ExactBoxType& A, const UpperBoxType& B) const
{
    ARIADNE_NOT_IMPLEMENTED;
}

Void AffineIntegrator::_write(OutputStream& os) const {
    os << "AffineIntegrator"
       << "( function_factory = " << this->function_factory()
       << ", spacial_order = " << this->spacial_order()
       << ", temporal_order = " << this->temporal_order()
       << " )";
}



} // namespace Ariadne
