/***************************************************************************
 *            test_constraint_solver.cpp
 *
 *  Copyright  2010-20  Pieter Collins
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

#include <iostream>
#include <fstream>

#include "config.hpp"
#include "../test.hpp"

#include "numeric/numeric.hpp"
#include "foundations/logical.hpp"
#include "algebra/vector.hpp"
#include "function/function.hpp"
#include "function/constraint.hpp"
#include "function/procedure.hpp"
#include "solvers/constraint_solver.hpp"
#include "geometry/box.hpp"
#include "io/command_line_interface.hpp"

using namespace std;
using namespace Ariadne;

class TestConstraintSolver
{
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_empty_reduce_inequality());
        ARIADNE_TEST_CALL(test_empty_reduce_equality());
        ARIADNE_TEST_CALL(test_empty_reduce_mixed());
        ARIADNE_TEST_CALL(test_empty_hull_reduce());
        ARIADNE_TEST_CALL(test_empty_box_reduce());
        ARIADNE_TEST_CALL(test_hull_reduce());
        ARIADNE_TEST_CALL(test_box_reduce());
        ARIADNE_TEST_CALL(test_pruning_well_definedness());
        ARIADNE_TEST_CALL(test_composite_pruning_witness_preservation());
        ARIADNE_TEST_CALL(test_partial_domain_pruning_witness_preservation());
        ARIADNE_TEST_CALL(test_monotone_reduce());
        ARIADNE_TEST_CALL(test_reduce_edge_cases());
        ARIADNE_TEST_CALL(test_check_feasibility());
        ARIADNE_TEST_CALL(test_feasible());
        ARIADNE_TEST_CALL(test_split());
    }

    Void test_empty_reduce_inequality() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0_x,1.0_x},{0.0_x,1.0_x}};
        List<EffectiveConstraint> c = {4<=2*x[0]+x[1]};

        ConstraintSolver propagator;

        ARIADNE_TEST_EXECUTE(propagator.reduce(D,c));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.is_empty());
    }

    Void test_empty_reduce_equality() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0_x,1.0_x},{0.0_x,1.0_x}};
        List<EffectiveConstraint> c = {2*x[0]+x[1]==4};

        ConstraintSolver propagator;

        ARIADNE_TEST_EXECUTE(propagator.reduce(D,c));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.is_empty());
    }

    Void test_empty_reduce_mixed() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0_x,0.25_x},{0.0_x, 2.0_x}};
        List<EffectiveConstraint> c = {x[1]<=1,x[0]+x[1]==2};

        ConstraintSolver propagator;

        ARIADNE_TEST_EXECUTE(propagator.reduce(D,c));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.is_empty());
    }

    Void test_empty_hull_reduce() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0_x,0.25_x},{0.0_x,2.0_x}};
        List<EffectiveConstraint> c = {x[1]<=1, x[0]+x[1]==2};

        ConstraintSolver propagator;

        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[0]));
        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[1]));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.is_empty());
    }

    Void test_empty_box_reduce() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0_x,0.25_x},{0.0_x, 2.0_x}};
        List<EffectiveConstraint> c = {x[1]<=1,x[0]+x[1]==2};

        ConstraintSolver propagator;

        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c[0],0));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c[1],0));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c[0],1));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c[1],1));
        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[0]));
        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[1]));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.is_empty());
    }

    Void test_hull_reduce() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0_x,2.0_x},{0.0_x,2.0_x}};
        List<EffectiveConstraint> c = {-2<=2*x[0]+x[1]<=1};

        ConstraintSolver propagator;

        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[0]));
        ARIADNE_TEST_SAME(D,UpperBoxType({{0.0_x,0.5_x},{0.0_x,1.0_x}}));
    }

    Void test_box_reduce() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0_x,2.0_x},{0.0_x,2.0_x}};
        EffectiveConstraint c = (-2<=2*x[0]+x[1]<=1);

        ConstraintSolver propagator;

        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c,0));
        ARIADNE_TEST_SAME(D,UpperBoxType({{0.0_x,0.75_x},{0.0_x,2.0_x}}));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c,1));
        ARIADNE_TEST_SAME(D,UpperBoxType({{0.0_x,0.75_x},{0.0_x,1.25_x}}));
    }


    Void test_pruning_well_definedness() {
        ConstraintSolver contractor;
        auto x=ValidatedScalarMultivariateFunction::coordinates(2);

        {
            std::cout << "[constraint-prune] hull reduction is reductive and preserves a known solution" << std::endl;
            UpperBoxType before=ExactBoxType({{-2.0_x,2.0_x},{-2.0_x,2.0_x}});
            UpperBoxType after=before;
            ValidatedScalarMultivariateFunction function=x[0]+x[1];
            ExactIntervalType codomain(1.0_x,1.0_x);
            UpperBoxType solution=ExactBoxType({
                {0.25_x,0.25_x},{0.75_x,0.75_x}
            });

            contractor.hull_reduce(after,function,codomain);
            ARIADNE_TEST_ASSERT(refines(after,before));
            ARIADNE_TEST_ASSERT(not definitely(disjoint(solution,after)));
        }

        {
            std::cout << "[constraint-prune] nonlinear hull reduction preserves a known solution" << std::endl;
            UpperBoxType before=ExactBoxType({{-2.0_x,2.0_x},{-2.0_x,2.0_x}});
            UpperBoxType after=before;
            ValidatedScalarMultivariateFunction function=sqr(x[0])+x[1];
            ExactIntervalType codomain(1.0_x,1.0_x);
            UpperBoxType solution=ExactBoxType({
                {0.0_x,0.0_x},{1.0_x,1.0_x}
            });

            contractor.hull_reduce(after,function,codomain);
            ARIADNE_TEST_ASSERT(refines(after,before));
            ARIADNE_TEST_ASSERT(not definitely(disjoint(solution,after)));
        }

        {
            std::cout << "[constraint-prune] box shaving is reductive and preserves a known solution" << std::endl;
            UpperBoxType before=ExactBoxType({{-2.0_x,2.0_x},{-2.0_x,2.0_x}});
            UpperBoxType after=before;
            ValidatedScalarMultivariateFunction function=sqr(x[0])+x[1];
            ExactIntervalType codomain(1.0_x,1.0_x);
            UpperBoxType solution=ExactBoxType({
                {0.0_x,0.0_x},{1.0_x,1.0_x}
            });

            contractor.box_reduce(after,function,codomain,0u);
            contractor.box_reduce(after,function,codomain,1u);
            ARIADNE_TEST_ASSERT(refines(after,before));
            ARIADNE_TEST_ASSERT(not definitely(disjoint(solution,after)));
        }

        {
            std::cout << "[constraint-prune] definitely infeasible hull reduction closes the box" << std::endl;
            auto y=ValidatedScalarMultivariateFunction::coordinates(1);
            UpperBoxType domain=ExactBoxType({{0.0_x,1.0_x}});
            contractor.hull_reduce(
                domain,y[0],ExactIntervalType(2.0_x,3.0_x));
            ARIADNE_TEST_ASSERT(definitely(domain.is_empty()));
        }

        {
            std::cout << "[constraint-prune] definitely infeasible box shaving closes the box" << std::endl;
            auto y=ValidatedScalarMultivariateFunction::coordinates(1);
            UpperBoxType domain=ExactBoxType({{0.0_x,1.0_x}});
            contractor.box_reduce(
                domain,y[0],ExactIntervalType(2.0_x,3.0_x),0u);
            ARIADNE_TEST_ASSERT(definitely(domain.is_empty()));
        }
    }

    Void test_composite_pruning_witness_preservation() {
        ConstraintSolver contractor;
        auto x=ValidatedScalarMultivariateFunction::coordinates(2);

        auto check = [&](String const& label,
                         UpperBoxType const& domain,
                         ValidatedScalarMultivariateFunction const& function,
                         ExactIntervalType const& codomain,
                         UpperBoxType const& witness) {
            std::cout << "[constraint-prune-composite] " << label << std::endl;
            UpperBoxType reduced=domain;
            contractor.hull_reduce(reduced,function,codomain);
            ARIADNE_TEST_ASSERT(refines(reduced,domain));
            ARIADNE_TEST_ASSERT(not definitely(disjoint(witness,reduced)));
        };

        check(
            "negative square branch survives",
            ExactBoxType({{-2.0_x,2.0_x},{0.0_x,0.0_x}}),
            sqr(x[0])+x[1],
            ExactIntervalType(1.0_x,1.0_x),
            ExactBoxType({{-1.0_x,-1.0_x},{0.0_x,0.0_x}}));

        check(
            "abs of repeated nonlinear expression survives",
            ExactBoxType({{-2.0_x,2.0_x},{-1.0_x,1.0_x}}),
            abs(sqr(x[0])-1)+sqr(x[1]),
            ExactIntervalType(0.0_x,0.0_x),
            ExactBoxType({{-1.0_x,-1.0_x},{0.0_x,0.0_x}}));

        check(
            "product and division composition survives",
            ExactBoxType({{-2.0_x,2.0_x},{1.0_x,3.0_x}}),
            (x[0]*x[1])/(x[0]+3),
            ExactIntervalType(-1.0_x,-1.0_x),
            ExactBoxType({{-1.0_x,-1.0_x},{2.0_x,2.0_x}}));

        check(
            "nested elementary composition survives",
            ExactBoxType({{-1.0_x,1.0_x},{-1.0_x,1.0_x}}),
            exp(sqr(x[0]))+abs(x[1]),
            ExactIntervalType(1.0_x,1.0_x),
            ExactBoxType({{0.0_x,0.0_x},{0.0_x,0.0_x}}));

        check(
            "repeated coordinate cancellation does not lose witness",
            ExactBoxType({{-2.0_x,2.0_x},{-2.0_x,2.0_x}}),
            (x[0]-x[0])+x[1],
            ExactIntervalType(0.5_x,0.5_x),
            ExactBoxType({{-1.5_x,-1.5_x},{0.5_x,0.5_x}}));
    }

    Void test_partial_domain_pruning_witness_preservation() {
        ConstraintSolver contractor;
        auto x=ValidatedScalarMultivariateFunction::coordinates(2);

        auto check = [&](String const& label,
                         UpperBoxType const& domain,
                         ValidatedScalarMultivariateFunction const& function,
                         ExactIntervalType const& codomain,
                         UpperBoxType const& witness) {
            std::cout << "[constraint-prune-domain] " << label << std::endl;
            UpperBoxType reduced=domain;
            contractor.hull_reduce(reduced,function,codomain);
            ARIADNE_TEST_ASSERT(refines(reduced,domain));
            ARIADNE_TEST_ASSERT(not definitely(disjoint(witness,reduced)));
        };

        check(
            "reciprocal domain straddles zero",
            ExactBoxType({{-2.0_x,2.0_x},{0.0_x,0.0_x}}),
            rec(x[0])+x[1],
            ExactIntervalType(1.0_x,1.0_x),
            ExactBoxType({{1.0_x,1.0_x},{0.0_x,0.0_x}}));

        check(
            "division denominator straddles zero",
            ExactBoxType({{1.0_x,2.0_x},{-2.0_x,2.0_x}}),
            x[0]/x[1],
            ExactIntervalType(1.0_x,1.0_x),
            ExactBoxType({{1.0_x,1.0_x},{1.0_x,1.0_x}}));

        check(
            "sqrt argument crosses negative side",
            ExactBoxType({{-1.0_x,4.0_x},{0.0_x,0.0_x}}),
            sqrt(x[0])+x[1],
            ExactIntervalType(1.0_x,1.0_x),
            ExactBoxType({{1.0_x,1.0_x},{0.0_x,0.0_x}}));

        check(
            "log argument crosses zero",
            ExactBoxType({{-1.0_x,3.0_x},{0.0_x,0.0_x}}),
            log(x[0])+x[1],
            ExactIntervalType(0.0_x,0.0_x),
            ExactBoxType({{1.0_x,1.0_x},{0.0_x,0.0_x}}));

        check(
            "asin argument crosses both boundaries",
            ExactBoxType({{-2.0_x,2.0_x},{0.0_x,0.0_x}}),
            asin(x[0])+x[1],
            ExactIntervalType(0.0_x,0.0_x),
            ExactBoxType({{0.0_x,0.0_x},{0.0_x,0.0_x}}));

        check(
            "acos argument crosses both boundaries",
            ExactBoxType({{-2.0_x,2.0_x},{-2.0_x,2.0_x}}),
            acos(x[0])+x[1],
            ExactIntervalType(0.0_x,0.0_x),
            ExactBoxType({{1.0_x,1.0_x},{-1.5707963267948966_x,-1.5707963267948966_x}}));
    }

    Void test_monotone_reduce() {
        ConstraintSolver propagator;
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);

        UpperBoxType linear_domain=ExactBoxType{{0.0_x,2.0_x}};
        UpperBoxType linear_original=linear_domain;
        Bool linear_empty=propagator.monotone_reduce(
            linear_domain,
            2*x[0]+1,
            ExactIntervalType(3.0_x,3.0_x),
            0u);
        ARIADNE_TEST_ASSERT(not linear_empty);
        ARIADNE_TEST_ASSERT(refines(linear_domain,linear_original));
        ARIADNE_TEST_ASSERT(
            possibly(contains(linear_domain[0],ExactDouble(1.0_x))));

        UpperBoxType smooth_domain=ExactBoxType{{0.0_x,2.0_x}};
        UpperBoxType smooth_original=smooth_domain;
        auto smooth_function=exp(x[0])+x[0];
        Bool smooth_empty=propagator.monotone_reduce(
            smooth_domain,
            smooth_function,
            smooth_function.derivative(0u),
            ExactIntervalType(3.0_x,3.0_x),
            0u);
        ARIADNE_TEST_ASSERT(not smooth_empty);
        ARIADNE_TEST_ASSERT(refines(smooth_domain,smooth_original));
        ARIADNE_TEST_ASSERT(
            smooth_domain[0].width().raw()<=smooth_original[0].width().raw());

        UpperBoxType singleton_domain=ExactBoxType{{1.0_x,1.0_x}};
        Bool singleton_empty=propagator.monotone_reduce(
            singleton_domain,
            x[0],
            ExactIntervalType(1.0_x,1.0_x),
            0u);
        ARIADNE_TEST_ASSERT(not singleton_empty);
        ARIADNE_TEST_SAME(
            singleton_domain,UpperBoxType(ExactBoxType{{1.0_x,1.0_x}}));

        UpperBoxType lower_boundary_domain=ExactBoxType{{0.0_x,2.0_x}};
        Bool lower_boundary_empty=propagator.monotone_reduce(
            lower_boundary_domain,x[0],ExactIntervalType(0.0_x,0.0_x),0u);
        ARIADNE_TEST_ASSERT(not lower_boundary_empty);
        ARIADNE_TEST_ASSERT(
            possibly(contains(lower_boundary_domain[0],ExactDouble(0.0_x))));

        UpperBoxType upper_boundary_domain=ExactBoxType{{0.0_x,2.0_x}};
        Bool upper_boundary_empty=propagator.monotone_reduce(
            upper_boundary_domain,x[0],ExactIntervalType(2.0_x,2.0_x),0u);
        ARIADNE_TEST_ASSERT(not upper_boundary_empty);
        ARIADNE_TEST_ASSERT(
            possibly(contains(upper_boundary_domain[0],ExactDouble(2.0_x))));

        UpperBoxType constraint_domain=ExactBoxType{{0.0_x,2.0_x}};
        ValidatedConstraint constraint(
            ValidatedNumber(1.0_x),
            x[0],
            ValidatedNumber(1.0_x));
        Bool constraint_empty=propagator.monotone_reduce(
            constraint_domain,constraint,0u);
        ARIADNE_TEST_ASSERT(not constraint_empty);
        ARIADNE_TEST_ASSERT(
            possibly(contains(constraint_domain[0],ExactDouble(1.0_x))));
    }

    Void test_reduce_edge_cases() {
        ConstraintSolver solver;
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);

        std::cout << "[constraint-reduce] already empty vector domain" << std::endl;
        UpperBoxType empty_domain=ExactBoxType({ExactIntervalType::empty_interval()});
        ValidatedVectorMultivariateFunction vector_function({x[0]});
        ExactBoxType vector_codomain({ExactIntervalType(0.0_x,1.0_x)});
        ARIADNE_TEST_ASSERT(
            solver.reduce(empty_domain,vector_function,vector_codomain));

        std::cout << "[constraint-reduce] nonempty vector domain reaches fixed point" << std::endl;
        UpperBoxType stable_domain=ExactBoxType({{0.0_x,1.0_x}});
        ARIADNE_TEST_ASSERT(
            not solver.reduce(stable_domain,vector_function,vector_codomain));
        ARIADNE_TEST_ASSERT(not definitely(stable_domain.is_empty()));
    }

    Void test_check_feasibility() {
        ConstraintSolver solver;
        auto x=ValidatedScalarMultivariateFunction::coordinates(1);
        ValidatedVectorMultivariateFunction function({x[0]});
        ExactBoxType domain({{0.0_x,2.0_x}});
        ExactBoxType codomain({{0.0_x,1.0_x}});

        std::cout << "[constraint-feasible] reject point outside domain" << std::endl;
        ConstraintSolver::ExactPointType outside({FloatDP(-1.0_x,dp)});
        ARIADNE_TEST_ASSERT(
            not possibly(solver.check_feasibility(
                domain,function,codomain,outside)));

        std::cout << "[constraint-feasible] reject image outside codomain" << std::endl;
        ConstraintSolver::ExactPointType image_outside({FloatDP(2.0_x,dp)});
        ARIADNE_TEST_ASSERT(
            not possibly(solver.check_feasibility(
                domain,function,codomain,image_outside)));

        std::cout << "[constraint-feasible] boundary point is indeterminate" << std::endl;
        ConstraintSolver::ExactPointType boundary({FloatDP(0.0_x,dp)});
        ARIADNE_TEST_ASSERT(
            is_indeterminate(solver.check_feasibility(
                domain,function,codomain,boundary)));

        std::cout << "[constraint-feasible] strict interior point is validated" << std::endl;
        ConstraintSolver::ExactPointType interior({FloatDP(0.5_x,dp)});
        ARIADNE_TEST_ASSERT(
            definitely(solver.check_feasibility(
                domain,function,codomain,interior)));
    }

    Void test_split() {
        ConstraintSolver solver;
        auto x=ValidatedScalarMultivariateFunction::coordinates(2);
        UpperBoxType domain=ExactBoxType({{0.0_x,4.0_x},{0.0_x,1.0_x}});
        ValidatedVectorMultivariateFunction function({x[0]+x[1]});
        ExactBoxType codomain({{1.0_x,2.0_x}});

        std::cout << "[constraint-split] delegates to box bisection" << std::endl;
        auto children=solver.split(domain,function,codomain);
        ARIADNE_TEST_ASSERT(refines(children.first,domain));
        ARIADNE_TEST_ASSERT(refines(children.second,domain));
        ARIADNE_TEST_ASSERT(not same(children.first,children.second));
        ARIADNE_TEST_ASSERT(
            children.first[0].upper_bound().raw()
            ==children.second[0].lower_bound().raw());
    }

    Void test_feasible() {

        ConstraintSolver contractor;
        {
            std::cout << "[constraint-feasible] empty constraint system" << std::endl;
            List<ValidatedConstraint> no_constraints;
            ExactBoxType nonempty_domain({{0.0_x,1.0_x}});
            auto nonempty=contractor.feasible(nonempty_domain,no_constraints);
            ARIADNE_TEST_ASSERT(definitely(nonempty.first));
            ARIADNE_TEST_EQUAL(
                nonempty.second.dimension(),nonempty_domain.dimension());

            ExactBoxType empty_domain({ExactIntervalType::empty_interval()});
            auto empty=contractor.feasible(empty_domain,no_constraints);
            ARIADNE_TEST_ASSERT(not possibly(empty.first));
        }

        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(1);
        EffectiveConstraint c = (x[0]-2<=0);

        List<ValidatedConstraint> constraints;
        constraints.append(c);

        ARIADNE_TEST_PRINT(constraints);

        ExactBoxType domain1({{1.9375_x,2.0_x}});
        auto feasible1=contractor.feasible(domain1,constraints);
        ARIADNE_TEST_ASSERT(definitely(feasible1.first));
        ARIADNE_TEST_EQUAL(feasible1.second.dimension(),domain1.dimension());
        ARIADNE_TEST_ASSERT(definitely(contractor.check_feasibility(
            domain1,
            ValidatedVectorMultivariateFunction(
                {constraints[0].function()}),
            ExactBoxType({constraints[0].bounds()}),
            feasible1.second)));

        ExactBoxType domain2({{2.015625_x,2.5_x}});
        ARIADNE_TEST_ASSERT(!possibly(contractor.feasible(domain2,constraints).first));

        ExactBoxType domain3({{2.0_x,2.015625_x}});
        ARIADNE_TEST_ASSERT(is_indeterminate(contractor.feasible(domain3,constraints).first));

        {
            auto xy=ValidatedScalarMultivariateFunction::coordinates(2);
            ValidatedVectorMultivariateFunction function({xy[0]+xy[1]});
            ExactBoxType domain({{0.0_x,1.0_x},{0.0_x,1.0_x}});
            ExactBoxType codomain({{0.19_x,0.21_x}});
            ARIADNE_TEST_ASSERT(not definitely(contractor.check_feasibility(
                domain,function,codomain,domain.midpoint())));
            auto feasibility_result=contractor.feasible(domain,function,codomain);
            ARIADNE_TEST_ASSERT(possibly(feasibility_result.first));
            if(definitely(feasibility_result.first)) {
                ARIADNE_TEST_EQUAL(feasibility_result.second.dimension(),domain.dimension());
                ARIADNE_TEST_ASSERT(definitely(contractor.check_feasibility(
                    domain,function,codomain,feasibility_result.second)));
            }
        }

        {
            auto coordinates=ValidatedScalarMultivariateFunction::coordinates(7);
            ValidatedScalarMultivariateFunction sum=coordinates[0];
            for(SizeType i=1u; i!=7u; ++i) {
                sum=sum+coordinates[i];
            }
            ValidatedVectorMultivariateFunction function({sum});
            ExactBoxType domain({
                {0.0_x,1.0_x},{0.0_x,1.0_x},{0.0_x,1.0_x},
                {0.0_x,1.0_x},{0.0_x,1.0_x},{0.0_x,1.0_x},
                {0.0_x,1.0_x}
            });
            ExactBoxType codomain({{1.175_x,1.425_x}});
            auto feasibility_result=contractor.feasible(domain,function,codomain);
            ARIADNE_TEST_ASSERT(possibly(feasibility_result.first));
            if(definitely(feasibility_result.first)) {
                ARIADNE_TEST_EQUAL(feasibility_result.second.dimension(),domain.dimension());
                ARIADNE_TEST_ASSERT(definitely(contractor.check_feasibility(
                    domain,function,codomain,feasibility_result.second)));
            }
        }
    }
};

Int main(Int argc, const char* argv[]) {
    if (not CommandLineInterface::instance().acquire(argc,argv)) return -1;
    TestConstraintSolver().test();
    return ARIADNE_TEST_FAILURES;
}

