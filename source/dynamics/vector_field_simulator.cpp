/***************************************************************************
 *            dynamics/vector_field_simulator.cpp
 *
 *  Copyright  2008-21  Luca Geretti, Mirko Albanese
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

#include "config.hpp"

#include "function/functional.hpp"

#include "algebra/algebra.hpp"

#include "utility/array.hpp"
#include "utility/container.hpp"
#include "utility/tuple.hpp"
#include "utility/stlio.hpp"
#include "symbolic/valuation.hpp"
#include "symbolic/assignment.hpp"
#include "symbolic/space.hpp"
#include "symbolic/expression_set.hpp"

#include "function/function.hpp"
#include "function/formula.hpp"
#include "function/taylor_model.hpp"

#include "solvers/runge_kutta_integrator.hpp"

#include "dynamics/orbit.hpp"
#include "dynamics/vector_field.hpp"
#include "dynamics/vector_field_simulator.hpp"

namespace Ariadne {

template class Orbit<LabelledPoint<Approximation<FloatDP>>>;
template class Orbit<Vector<LabelledPoint<Approximation<FloatDP>>>>;

template<class X> LabelledPoint<X> make_state_auxiliary_point(const Point<X>& spt,
        const RealSpace& sspc, const RealSpace& aspc, const RealSpace& saspc, const EffectiveVectorMultivariateFunction& auxiliary_function) {
    Point<X> sapt(saspc.dimension(),spt.zero_element());
    Point<X> apt = evaluate(auxiliary_function,spt);
    for(SizeType i=0; i!=sapt.size(); ++i) {
        RealVariable var = saspc.variable(i);
        sapt[i]= sspc.contains(var) ? spt[sspc[var]] : apt[aspc[var]];
    }
    return LabelledPoint<X>(saspc,sapt);
}

template<class X> Vector<LabelledPoint<X>> make_state_auxiliary_point(const Vector<Point<X>>& spt,
        const RealSpace& sspc, const RealSpace& aspc, const RealSpace& saspc, const EffectiveVectorMultivariateFunction& auxiliary_function) {
    Vector<LabelledPoint<X>> pointLst(spt.size(), LabelledPoint<X>(saspc, Point<X>()));
    for (SizeType i=0; i<spt.size(); i++){
        Point<X> sapt(saspc.dimension(),spt.at(i).zero_element());
        Point<X> apt = evaluate(auxiliary_function,spt.at(i));
        for(SizeType j=0; j!=sapt.size(); ++j) {
            RealVariable var = saspc.variable(j);
            sapt[j]= sspc.contains(var) ? spt.at(i)[sspc[var]] : apt[aspc[var]];
        }
        LabelledPoint<X> point(saspc, sapt);
        pointLst.at(i) = point;
    }

    return pointLst;
}

template<class X> Vector<LabelledPoint<X>> make_state_auxiliary_point(const Vector<LabelledPoint<X>>& spt,
        const RealSpace& sspc, const RealSpace& aspc, const RealSpace& saspc, const EffectiveVectorMultivariateFunction& auxiliary_function) {
    Vector<LabelledPoint<X>> pointLst(spt.size(), LabelledPoint<X>(saspc, Point<X>()));
    for (SizeType i=0; i<spt.size(); i++){
        Point<X> sapt(saspc.dimension(),spt.at(i).zero_element());
        Point<X> apt = evaluate(auxiliary_function,spt.at(i));
        for(SizeType j=0; j!=sapt.size(); ++j) {
            RealVariable var = saspc.variable(j);
            sapt[j]= sspc.contains(var) ? spt.at(i)[sspc[var]] : apt[aspc[var]];
        }
        LabelledPoint<X> point(saspc, sapt);
        pointLst.at(i) = point;
    }

    return pointLst;
}

GridTreePaving create_paving(UpperBoxType box) {
    Point<FloatDP> tmpPointCenter = Point(box.centre());
    Vector<FloatDP> origin(tmpPointCenter.dimension(), FloatDP(0, dp));
    Vector<FloatDP> lengths(box.widths().size(), FloatDP(0, dp));
    for(SizeType i=0; i<tmpPointCenter.dimension(); i++){
        origin[i] = tmpPointCenter[i];
        lengths[i] = box.widths().at(i).raw();
    }
    Grid grid(origin, lengths);
    GridTreePaving gridPaving(grid);

    return gridPaving;
}

Vector<LabelledPoint<FloatDPApproximation>> create_point_list(GridTreePaving& paving, RealSpace spc, DiscretisationType discretisation_type, Nat mince_dimension){

    if(discretisation_type == DiscretisationType::Mince) paving.mince(mince_dimension);
    else paving.recombine();

    GridTreePaving::ConstIterator iter = paving.begin();
    Vector<LabelledPoint<FloatDPApproximation>> result(paving.size(), LabelledPoint(spc, Point<FloatDPApproximation>()));
    SizeType k(0);
    for( ; iter != paving.end(); ++iter){
        UpperBoxType cell = iter->box();
        auto midpoint = cell.midpoint();
        LabelledPoint<FloatDPApproximation> pt(spc, midpoint);
        result.at(k) = LabelledPoint<FloatDPApproximation>(spc, pt);
        k++;
    }

    return result;
}

VectorFieldSimulator::VectorFieldSimulator(SystemType const& system) : _system(system.clone()), _configuration(new VectorFieldSimulatorConfiguration())
{ }

inline FloatDPApproximation evaluate(const EffectiveScalarMultivariateFunction& f, const Vector<FloatDPApproximation>& x) { return f(x); }
inline Vector<FloatDPApproximation> evaluate(const EffectiveVectorMultivariateFunction& f, const Vector<FloatDPApproximation>& x) { return f(x); }

auto VectorFieldSimulator::orbit(const RealExpressionBoundedConstraintSet& init_set, const TerminationType& termination) const
-> OrbitType
{
    auto spc = _system->state_space();
    UpperBoxType box = init_set.euclidean_set(spc).bounding_box();
    return orbit(box, termination);

}

auto VectorFieldSimulator::orbit(const RealBoxType& init_bx, const TerminationType& termination) const
-> OrbitType
{  
    auto spc = _system->state_space();
    auto box = init_bx.euclidean_set(spc);
    UpperBoxType ubox(box);
    return orbit(ubox, termination);

}

auto VectorFieldSimulator::orbit(UpperBoxType& initial_box, const TerminationType& termination) const
-> OrbitType
{
    auto lengths = initial_box.widths();
    Nat box_width_null = 0;
    for(SizeType i=0; i<lengths.size(); i++){
        if(lengths[i].get_d() > 0) { continue; }
        box_width_null++;
    }
    if(box_width_null == lengths.size()) 
    {
        auto midpoint = initial_box.midpoint();
        ApproximateListPointType pointList(1, LabelledPoint(_system->state_space(), Point<FloatDPApproximation>(midpoint, dp)));
        return orbit(pointList, termination);
    }
    else if (box_width_null > 0)
    {
        FloatDPUpperBound eps(0.0001_q,dp);
        initial_box = widen(initial_box, eps);
    }
    GridTreePaving gridPaving = create_paving(initial_box);
    gridPaving.adjoin_outer_approximation(initial_box, _configuration->num_subdivisions());
    ApproximateListPointType pointList = create_point_list(gridPaving, _system->state_space(),
                                                           _configuration->discretisation_type(), _configuration->num_subdivisions());
    return orbit(pointList, termination);
}

auto VectorFieldSimulator::orbit(const ApproximateListPointType& initial_points, const TerminationType& termination) const -> OrbitType {
    CONCLOG_SCOPE_CREATE;

    CONCLOG_PRINTLN("Simulating from " << initial_points.size() << " initial points")

    auto const& auxiliary_function = _system->auxiliary_function();

    auto state_space=_system->state_space();
    auto auxiliary_space = _system->auxiliary_space();
    auto state_auxiliary_space = _system->state_auxiliary_space();

    auto result = std::make_shared<SynchronisedOrbit>(make_state_auxiliary_point(initial_points, state_space, auxiliary_space, state_auxiliary_space, auxiliary_function));

    WorkloadType workload(std::bind_front(&VectorFieldSimulator::_simulate_from_point,this),termination,result);
    for (SizeType i=0; i<initial_points.size(); i++)
        workload.append({i,initial_points.at(i)});
    workload.process();

    return std::move(*result);
}

void VectorFieldSimulator::_simulate_from_point(Pair<SizeType,ApproximatePointType> indexed_initial, TerminationType const& termination, SharedPointer<SynchronisedOrbit> orbit) const {
    CONCLOG_SCOPE_CREATE

    auto const& curve_number = indexed_initial.first;
    auto const& initial = indexed_initial.second;

    VectorField::TimeType t(0);
    Dyadic h(cast_exact(_configuration->step_size()));
    VectorField::TimeType tmax(termination);
    auto const& dynamic_function = _system->dynamic_function();
    auto const& auxiliary_function = _system->auxiliary_function();

    auto const state_space=_system->state_space();
    auto const auxiliary_space = _system->auxiliary_space();
    auto const state_auxiliary_space = _system->state_auxiliary_space();

    RungeKutta4Integrator integrator(_configuration->step_size().get_d());
    Point<FloatDPApproximation> state_pt = initial;
    while(possibly(t<tmax)) {
        Int old_precision = std::clog.precision();
        CONCLOG_PRINTLN("t=" << std::setw(4) << std::left << t.get(dp).value() << " p=" << state_pt << std::setprecision(old_precision));
        state_pt = integrator.step(dynamic_function, state_pt, _configuration->step_size());

        t += h;
        orbit->insert(t.get(DoublePrecision()).value(), make_state_auxiliary_point(ApproximatePointType(state_space, state_pt), state_space, auxiliary_space, state_auxiliary_space, auxiliary_function), curve_number);
    }
}

VectorFieldSimulatorConfiguration::VectorFieldSimulatorConfiguration() :
        _step_size(0.125_x, dp), _num_subdivisions(0), _discretisation_type(DiscretisationType::Recombine) {
}

OutputStream&
VectorFieldSimulatorConfiguration::_write(OutputStream& os) const
{
    os << "VectorFieldSimulatorConfiguration("
       << "\n  step_size=" << step_size()
       << "\n  discretisation_type=" << discretisation_type()
       << "\n  num_subdivisions=" << num_subdivisions()
       << ")\n";
    return os;
}

}  // namespace Ariadne

