/***************************************************************************
 *            dynamics/vector_field_simulator.cpp
 *
 *  Copyright  2008-21  Luca Geretti
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

template class Orbit<Point<FloatDPApproximation>>;

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

VectorFieldSimulator::VectorFieldSimulator(SystemType const& system) : _system(system.clone()), _configuration(new VectorFieldSimulatorConfiguration())
{ }

inline FloatDPApproximation evaluate(const EffectiveScalarMultivariateFunction& f, const Vector<FloatDPApproximation>& x) { return f(x); }
inline Vector<FloatDPApproximation> evaluate(const EffectiveVectorMultivariateFunction& f, const Vector<FloatDPApproximation>& x) { return f(x); }

auto VectorFieldSimulator::orbit(const RealExpressionBoundedConstraintSet& init_set, const TerminationType& termination) const
-> OrbitType
{
    auto spc = _system->state_space();
    auto midpoint = init_set.euclidean_set(spc).bounding_box().midpoint();
    return orbit(ApproximatePointType(spc,midpoint),termination);
}

auto VectorFieldSimulator::orbit(const RealBoxType& init_bx, const TerminationType& termination) const
-> OrbitType
{
    auto spc = _system->state_space();
    auto midpoint = init_bx.euclidean_set(spc).midpoint();
    return orbit(RealPointType(spc,midpoint),termination);
}

auto VectorFieldSimulator::orbit(const RealPointType& init_pt, const TerminationType& termination) const
    -> OrbitType
{
    return orbit(ApproximatePointType(init_pt,dp),termination);
}

auto VectorFieldSimulator::orbit(const ApproximatePointType& init_pt, const TerminationType& termination) const
    -> OrbitType
{
    CONCLOG_SCOPE_CREATE;

    VectorField::TimeType t;
    Dyadic h(cast_exact(configuration().step_size()));
    VectorField::TimeType tmax(termination);

    auto const& dynamic_function =_system->function();
    auto const& auxiliary_function = _system->auxiliary_function();

    auto state_space=_system->state_space();
    auto auxiliary_space = _system->auxiliary_space();
    auto state_auxiliary_space = _system->state_auxiliary_space();

    RungeKutta4Integrator integrator(configuration().step_size().get_d());

    ApproximatePointType point=init_pt;

    OrbitType orbit(make_state_auxiliary_point(point,state_space,auxiliary_space,state_auxiliary_space,auxiliary_function));

    while(possibly(t<tmax)) {
        Int old_precision = std::clog.precision();
        CONCLOG_PRINTLN("t=" << std::setw(4) << std::left << t.get(dp).value() << " p=" << point << std::setprecision(old_precision));

        Point<FloatDPApproximation> state_pt = integrator.step(dynamic_function,point,configuration().step_size());

        point = ApproximatePointType(state_space,state_pt);
        t += h;

        orbit.insert(t.get(dp).value(),make_state_auxiliary_point(point,state_space,auxiliary_space,state_auxiliary_space,auxiliary_function));
    }

    return orbit;
}

VectorFieldSimulatorConfiguration::VectorFieldSimulatorConfiguration() :
    _step_size(0.125_x, dp) {
}

OutputStream&
VectorFieldSimulatorConfiguration::_write(OutputStream& os) const
{
    os << "VectorFieldSimulatorConfiguration("
       << "\n  step_size=" << step_size()
       << "\n)";
    return os;
}

}  // namespace Ariadne

