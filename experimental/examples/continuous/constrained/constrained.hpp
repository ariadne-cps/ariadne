/***************************************************************************
 *            constrained.hpp
 *
 *  Copyright  2023  Luca Geretti
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

#ifndef ARIADNE_CONSTRAINED_HPP
#define ARIADNE_CONSTRAINED_HPP

#include "ariadne_main.hpp"
#include "helper/stopwatch.hpp"
#include "helper/randomiser.hpp"

using namespace std;
using namespace Ariadne;
using namespace pExplore;
using namespace Helper;

struct SystemSpecification {
    SystemSpecification(Helper::String const& name_, VectorField const& dynamics_, RealExpressionBoundedConstraintSet const& initial_set_, Real const& evolution_time_)
        : name(name_), dynamics(dynamics_), initial_set(initial_set_), evolution_time(evolution_time_) { }

    Helper::String name;
    VectorField dynamics;
    RealExpressionBoundedConstraintSet initial_set;
    Real evolution_time;
};

class UnconstrainedEvolutionResult {
  public:
    UnconstrainedEvolutionResult(Orbit<LabelledEnclosure> const& approximate_orbit, Orbit<LabelledEnclosure> const& rigorous_orbit, ConstraintSatisfaction const& satisfaction) :
        _approximate_orbit(approximate_orbit), _rigorous_orbit(rigorous_orbit), _satisfaction(satisfaction) { }
    Orbit<LabelledEnclosure> const& approximate_orbit() const { return _approximate_orbit; }
    Orbit<LabelledEnclosure> const& rigorous_orbit() const { return _rigorous_orbit; }
    ConstraintSatisfaction const& satisfaction() const { return _satisfaction; }
  private:
    Orbit<LabelledEnclosure> const _approximate_orbit;
    Orbit<LabelledEnclosure> const _rigorous_orbit;
    ConstraintSatisfaction const _satisfaction;
};

LabelledUpperBoxType bounding_box(ListSet<LabelledEnclosure> const& es) {
    HELPER_PRECONDITION(not es.empty())
    UpperBoxType result = es[0].euclidean_set().bounding_box();
    for (auto const& e : es) {
        auto const& bx = e.euclidean_set().bounding_box();
        for (size_t j=0; j<result.dimension(); ++j)
            result[j] = hull(result[j],bx[j]);
    }
    return {es[0].state_space(),result};
}

UnconstrainedEvolutionResult unconstrained_evolution(VectorField const& dynamics, RealExpressionBoundedConstraintSet const& initial_set, Real const& evolution_time,
                                                     Configuration<VectorFieldEvolver> const& configuration, List<RealExpression> const& constraints, double time_budget) {
    CONCLOG_SCOPE_CREATE

    static const double MAXIMUM_ENCLOSURE_FRACTION = 0.25;

    auto satisfaction = ConstraintSatisfaction(constraints,dynamics.state_space(), time_budget);
    auto h = satisfaction.indeterminate_constraints_function();

    Orbit<LabelledEnclosure> approximate_orbit({});
    Orbit<LabelledEnclosure> rigorous_orbit({});

    size_t num_indeterminates = constraints.size();

    Stopwatch<Milliseconds> sw;

    auto analysis_point = configuration.search_space().initial_point();
    CONCLOG_PRINTLN_AT(1,"Using point " << analysis_point << " for approximate analysis.")

    auto singleton_configuration = make_singleton(configuration, analysis_point);
    singleton_configuration.set_enable_clobbering(true);

    VectorFieldEvolver approximate_evolver(dynamics, singleton_configuration);

    CONCLOG_PRINTLN_AT(1,"Computing approximate singleton evolution...")
    approximate_orbit = approximate_evolver.orbit(initial_set, evolution_time, Semantics::LOWER);

    CONCLOG_PRINTLN_AT(1,"Computing rigorous singleton evolution...")
    auto approximate_bounding_box = bounding_box(approximate_orbit);

    for (size_t i=0; ; ++i) {
        CONCLOG_PRINTLN("Round #" << i)

        analysis_point = configuration.search_space().initial_point();
        CONCLOG_PRINTLN_AT(1,"Using point " << analysis_point << " for rigorous analysis.")

        singleton_configuration = make_singleton(configuration, analysis_point);
        singleton_configuration.set_enable_clobbering(false).set_enable_premature_termination(true).set_maximum_enclosure_radius(approximate_bounding_box.radius().get_d()*MAXIMUM_ENCLOSURE_FRACTION);

        VectorFieldEvolver rigorous_evolver(dynamics, singleton_configuration);

        h = satisfaction.indeterminate_constraints_function();
        auto uncontrolled_task_constraints = build_uncontrolled_full_task_constraints(h);
        rigorous_evolver.set_constraints(uncontrolled_task_constraints);
        rigorous_orbit = rigorous_evolver.orbit(initial_set, evolution_time, Semantics::LOWER);

        bool exclude_truth = rigorous_orbit.final().empty();
        if (exclude_truth) {
            CONCLOG_PRINTLN("Prescription TRUE will not be checked on this round.")
        }

        satisfaction.merge_from_uncontrolled(rigorous_evolver.constraining_state(), exclude_truth);

        if (satisfaction.has_expired()) {
            CONCLOG_PRINTLN("Time budget hit, aborting.")
            break;
        }

        num_indeterminates = satisfaction.indeterminate_indexes().size();

        CONCLOG_PRINTLN_VAR_AT(1,satisfaction)
        CONCLOG_PRINTLN_VAR(satisfaction.global_success_ratio())

        if (num_indeterminates == 0) {
            CONCLOG_PRINTLN("All constraints satisfiability determined, terminating.")
            break;
        }
    }

    return {approximate_orbit,rigorous_orbit,satisfaction};
}

Map<SatisfactionPrescriptionKind,double> frequencies(List<ConstraintPrescription> const& constraint_prescriptions) {
    Map<SatisfactionPrescriptionKind,double> result;
    result.insert(SatisfactionPrescriptionKind::TRUE, 0);
    result.insert(SatisfactionPrescriptionKind::FALSE_FOR_ALL, 0);
    result.insert(SatisfactionPrescriptionKind::FALSE_FOR_SOME, 0);

    for (auto const& cp : constraint_prescriptions) {
        result[cp.prescription]++;
    }

    auto num_elements = static_cast<double>(constraint_prescriptions.size());

    result[SatisfactionPrescriptionKind::TRUE]/=num_elements;
    result[SatisfactionPrescriptionKind::FALSE_FOR_ALL]/=num_elements;
    result[SatisfactionPrescriptionKind::FALSE_FOR_SOME]/=num_elements;

    return result;
}

List<RealExpression> constraints(List<ConstraintPrescription> const& constraint_prescriptions) {
    List<RealExpression> result;
    for (auto const& cp : constraint_prescriptions) result.push_back(cp.expression);
    return result;
}

Configuration<VectorFieldEvolver> get_configuration() {
    return Configuration<VectorFieldEvolver>().
            set_integrator(TaylorPicardIntegrator(
                Configuration<TaylorPicardIntegrator>()
                .set_step_maximum_error(1e-6,1e-4)
                .set_sweeper(ThresholdSweeperDP(dp,
                    Configuration<ThresholdSweeperDP>()
                    .set_threshold(1e-8,1e-6)))
                .set_minimum_temporal_order(0,5)
                .set_maximum_temporal_order(5,10)
                .set_bounder(EulerBounder(
                    Configuration<EulerBounder>()
                    .set_lipschitz_tolerance(0.0675,0.5)))));
}

List<ConstraintPrescription> generate_ellipsoidal_hyperbolic_constraints(size_t num_constraints, SystemSpecification const& spec, Configuration<VectorFieldEvolver> const& configuration) {
    auto analysis_point = configuration.search_space().initial_point();
    auto singleton_configuration = make_singleton(configuration, analysis_point);
    singleton_configuration.set_enable_clobbering(true);
    VectorFieldEvolver approximate_evolver(spec.dynamics, singleton_configuration);
    CONCLOG_RUN_MUTED(auto orbit = approximate_evolver.orbit(spec.initial_set, spec.evolution_time, Semantics::LOWER))
    auto bb = bounding_box(orbit.intermediate());

    List<RealExpression> constraints;

    UniformRealRandomiser<double> lower_rnd(0.5,1.0);
    UniformRealRandomiser<double> upper_rnd(1.0,1.5);

    UniformIntRandomiser<size_t> sign_rnd(0,1);

    for (size_t i=0; i<num_constraints; ++i) {
        while (true) {
            List<size_t> signs;
            size_t sum = 0;
            for (size_t j=0; j<spec.dynamics.dimension()+1; ++j) {
                signs.append(sign_rnd.get());
                sum += signs.at(j);
            }
            if (sum != 0 and sum != spec.dynamics.dimension()+1) {
                size_t j=0;
                RealExpression expr = cast_exact(static_cast<double>(signs.at(j++))*2-1);
                bool use_lower = (sign_rnd.get() == 0);
                for (auto const& v : spec.dynamics.state_space().variables())
                    expr = expr + sqr((v - RealConstant(cast_exact(bb[v].midpoint().get_d())))/RealConstant(cast_exact((use_lower ? lower_rnd.get() : upper_rnd.get())*bb[v].radius().get_d())))*cast_exact(static_cast<double>(signs.at(j++))*2-1);
                constraints.push_back(expr);
                break;
            }
        }
    }

    // Discard those with all negative or all positive signs

    List<EffectiveScalarMultivariateFunction> constraint_functions;
    for (auto const& c : constraints)
        constraint_functions.push_back(make_function(spec.dynamics.state_space(),c));

    List<ConstraintPrescription> result;
    for (auto const& c : constraints)
        result.push_back({c, SatisfactionPrescriptionKind::TRUE});

    for (auto const& encl : orbit.intermediate()) {
        for (size_t m=0; m<constraints.size(); ++m) {
            auto const& eval = outer_evaluate_from_function(constraint_functions.at(m),encl);
            auto upper = eval.upper().get_d();
            auto lower = eval.lower().get_d();
            if (result.at(m).prescription != SatisfactionPrescriptionKind::FALSE_FOR_ALL) {
                if (upper < 0) result.at(m).prescription = SatisfactionPrescriptionKind::FALSE_FOR_ALL;
                else if (lower < 0) result.at(m).prescription = SatisfactionPrescriptionKind::FALSE_FOR_SOME;
            }
        }
    }

    return result;
}

List<ConstraintPrescription> generate_ellipsoidal_constraints(size_t num_constraints, SystemSpecification const& spec, Configuration<VectorFieldEvolver> const& configuration) {
    auto analysis_point = configuration.search_space().initial_point();
    auto singleton_configuration = make_singleton(configuration, analysis_point);
    singleton_configuration.set_enable_clobbering(true);
    VectorFieldEvolver approximate_evolver(spec.dynamics, singleton_configuration);
    CONCLOG_RUN_MUTED(auto orbit = approximate_evolver.orbit(spec.initial_set, spec.evolution_time, Semantics::LOWER))
    auto bb = bounding_box(orbit.intermediate());

    List<RealExpression> constraints;

    UniformRealRandomiser<double> lower_rnd(0.5,1.0);
    UniformRealRandomiser<double> upper_rnd(1.0,1.5);

    for (size_t i=0; i<num_constraints/2; ++i) {
        RealExpression expr = -1;
        for (auto const& v : spec.dynamics.state_space().variables())
            expr = expr + sqr((v - RealConstant(cast_exact(bb[v].midpoint().get_d())))/RealConstant(cast_exact(lower_rnd.get()*bb[v].radius().get_d())));
        constraints.push_back(expr);
    }
    for (size_t i=0; i<num_constraints/2; ++i) {
        RealExpression expr = 1;
        for (auto const& v : spec.dynamics.state_space().variables())
            expr = expr - sqr((v - RealConstant(cast_exact(bb[v].midpoint().get_d())))/RealConstant(cast_exact(upper_rnd.get()*bb[v].radius().get_d())));
        constraints.push_back(expr);
    }

    List<EffectiveScalarMultivariateFunction> constraint_functions;
    for (auto const& c : constraints)
        constraint_functions.push_back(make_function(spec.dynamics.state_space(),c));

    List<ConstraintPrescription> result;
    for (auto const& c : constraints)
        result.push_back({c, SatisfactionPrescriptionKind::TRUE});

    for (auto const& encl : orbit.intermediate()) {
        for (size_t m=0; m<constraints.size(); ++m) {
            auto const& eval = outer_evaluate_from_function(constraint_functions.at(m),encl);
            auto upper = eval.upper().get_d();
            auto lower = eval.lower().get_d();
            if (result.at(m).prescription != SatisfactionPrescriptionKind::FALSE_FOR_ALL) {
                if (upper < 0) result.at(m).prescription = SatisfactionPrescriptionKind::FALSE_FOR_ALL;
                else if (lower < 0) result.at(m).prescription = SatisfactionPrescriptionKind::FALSE_FOR_SOME;
            }
        }
    }

    return result;
}

void constrained_execution(SystemSpecification const& spec, Configuration<VectorFieldEvolver> const& configuration, List<RealExpression> const& constraints) {

    CONCLOG_PRINTLN("Dynamics: " << spec.dynamics)
    CONCLOG_PRINTLN_VAR(configuration)
    CONCLOG_PRINTLN_VAR_AT(1,configuration.search_space())
    CONCLOG_PRINTLN_AT(1,"Total points: " << configuration.search_space().total_points())

    auto result = constrained_evolution(spec.dynamics,spec.initial_set,spec.evolution_time,configuration,constraints);

    CONCLOG_PRINTLN("Constraint satisfaction result:" << result.satisfaction())

    auto best_scores = pExplore::TaskManager::instance().best_scores();
    CONCLOG_PRINTLN_AT(2,"Best scores:")
    for (auto const& b : best_scores) {
        CONCLOG_PRINTLN_AT(2,b)
    }

    CONCLOG_PRINTLN_AT(1,"Optimal point: " << pExplore::TaskManager::instance().optimal_point())

    auto variable_names = spec.dynamics.state_space().variable_names();
    auto drawing_box = bounding_box(result.rigorous().reach());

    CONCLOG_PRINTLN("Plotting...")
    for (size_t i=0; i<spec.dynamics.dimension()-1; i++) {
        auto xi = RealVariable(variable_names.at(i));
        for (size_t j=1; j<spec.dynamics.dimension(); j++) {
            auto xj = RealVariable(variable_names.at(j));
            LabelledFigure fig=LabelledFigure({drawing_box[xi].lower_bound().get_d()<=xi<=drawing_box[xi].upper_bound().get_d(),drawing_box[xj].lower_bound().get_d()<=xj<=drawing_box[xj].upper_bound().get_d()});
            fig.clear();
            fig.draw(result.approximate());
            char var_char[64] = "";
            if (spec.dynamics.dimension() > 2) snprintf(var_char,64,"[%s,%s]",xi.name().c_str(),xj.name().c_str());
            CONCLOG_RUN_MUTED(fig.write((spec.name+"_approximate"+var_char).c_str()))
            fig.clear();
            fig.draw(result.rigorous());
            CONCLOG_RUN_MUTED(fig.write((spec.name+"_rigorous"+var_char).c_str()))
            fig.clear();
            fig.draw(result.constrained());
            CONCLOG_RUN_MUTED(fig.write((spec.name+"_constrained"+var_char).c_str()))
        }
    }
}

#endif
