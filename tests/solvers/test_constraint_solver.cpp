/***************************************************************************
 *            test_constraint_solver.cpp
 *
 *  Copyright  2010  Pieter Collins
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

#include <chrono>

#include "config.hpp"
#include "../test.hpp"

#include "numeric/numeric.hpp"
#include "numeric/logical.hpp"
#include "algebra/vector.hpp"
#include "function/function.hpp"
#include "function/constraint.hpp"
#include "function/procedure.hpp"
#include "solvers/constraint_solver.hpp"
#include "geometry/box.hpp"

#include "function/formula.hpp"
#include "RandomFunction.hpp"

#include "solvers/nonlinear_programming.hpp"


using namespace std;
using namespace Ariadne;


Pair<ValidatedKleenean,ExactPoint> __feasible__(const ExactBoxType& domain, const List<ValidatedConstraint>& constraints)
{
  if(constraints.empty()) { return make_pair(!domain.is_empty(),domain.midpoint()); }

  ValidatedVectorMultivariateFunction function(constraints.size(),constraints[0].function().domain());
  ExactBoxType codomain(constraints.size());
  ValidatedScalarMultivariateFunction barrier_function(domain.size());
  EffectiveVectorMultivariateFunction empty_function(0u, domain.size());
  ExactBoxType empty_box = ExactBoxType{};
  NonlinearSQPOptimiser nlsqp;
  NonlinearInteriorPointOptimiser nlipm;
  nlsqp.verbosity=6;

  Ariadne::EffectiveScalarMultivariateFunction mu(
    Ariadne::EuclideanDomain(domain.size()),
    Ariadne::simplify(Ariadne::EffectiveFormula::constant(
      Real(0.5)
    )));

  for(Nat i=0; i!=constraints.size(); ++i) {
      function[i]=constraints[i].function();
      codomain[i]=constraints[i].bounds();
      // std::cerr<<"l: "<<codomain[i].lower()<<", u: "<<codomain[i].upper()<<"\n\n";
      if(codomain[i].lower()==codomain[i].upper())
      {
        barrier_function = barrier_function + 1/(2*mu)*pow(function[i],2);
        continue;
      }
      // if(!is_inf(cast_raw(codomain[i].upper())))
      //   barrier_function = barrier_function - log(function[i]-Real(codomain[i].lower()));
      // if(!is_inf(cast_raw(codomain[i].lower())))
      //   barrier_function = barrier_function - log(Real(codomain[i].upper())-function[i]);
      barrier_function = barrier_function - 1/(function[i]);
  }

  UpperBoxType image=apply(function,domain);
  for(Nat i=0; i!=image.size(); ++i) {
      if(definitely(disjoint(image[i],codomain[i]))) {
          return make_pair(false,ExactPoint());
      }
  }

  if(decide(nlsqp.check_feasibility(domain,function,codomain,midpoint(domain))))
  {
    return make_pair(true,midpoint(domain));
  }


  auto optimal_x = nlsqp.minimise(barrier_function,domain,empty_function,empty_box);
  // std::cerr<<"Minimum of nlsqp: "<<optimal_x<<"\n";
  if(decide(nlsqp.check_feasibility(domain,function,codomain,cast_exact(optimal_x))))
  {
    return make_pair(true,cast_exact(optimal_x));
  }
  return make_pair(indeterminate,cast_exact(optimal_x));
}

class TestConstraintSolver
{
    Nat verbosity;
  public:
    TestConstraintSolver(Nat v) : verbosity(v) { }

    Void test() {
        ARIADNE_TEST_CALL(test_empty_reduce_inequality());
        ARIADNE_TEST_CALL(test_empty_reduce_equality());
        ARIADNE_TEST_CALL(test_empty_reduce_mixed());
        ARIADNE_TEST_CALL(test_empty_hull_reduce());
        ARIADNE_TEST_CALL(test_empty_box_reduce());
        ARIADNE_TEST_CALL(test_hull_reduce());
        ARIADNE_TEST_CALL(test_box_reduce());
        ARIADNE_TEST_CALL(test_monotone_reduce());
        ARIADNE_TEST_CALL(test_feasible());
        ARIADNE_TEST_CALL(test_split());
        ARIADNE_TEST_CALL(test_speed());
    }

    Void test_empty_reduce_inequality() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0,1.0},{0.0,1.0}};
        List<EffectiveConstraint> c = {4<=2*x[0]+x[1]};

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.reduce(D,c));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.is_empty());
    }

    Void test_empty_reduce_equality() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0,1.0},{0.0,1.0}};
        List<EffectiveConstraint> c = {2*x[0]+x[1]==4};

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.reduce(D,c));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.is_empty());
    }

    Void test_empty_reduce_mixed() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0,0.25},{0.0, 2.0}};
        List<EffectiveConstraint> c = {x[1]<=1,x[0]+x[1]==2};

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.reduce(D,c));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.is_empty());
    }

    Void test_empty_hull_reduce() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0,0.25},{0.0,2.0}};
        List<EffectiveConstraint> c = {x[1]<=1, x[0]+x[1]==2};

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[0]));
        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[1]));
        ARIADNE_TEST_PRINT(D);
        ARIADNE_TEST_ASSERT(D.is_empty());
    }

    Void test_empty_box_reduce() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0,0.25},{0.0, 2.0}};
        List<EffectiveConstraint> c = {x[1]<=1,x[0]+x[1]==2};

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

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
        UpperBoxType D = ExactBoxType{{0.0,2.0},{0.0,2.0}};
        List<EffectiveConstraint> c = {-2<=2*x[0]+x[1]<=1};

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.hull_reduce(D,c[0]));
        ARIADNE_TEST_SAME(D,UpperBoxType({{0.0,0.5},{0.0,1.0}}));
    }

    Void test_box_reduce() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0,2.0},{0.0,2.0}};
        EffectiveConstraint c = (-2<=2*x[0]+x[1]<=1);

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c,0));
        ARIADNE_TEST_SAME(D,UpperBoxType({{0.0,0.75},{0.0,2.0}}));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c,1));
        ARIADNE_TEST_SAME(D,UpperBoxType({{0.0,0.75},{0.0,1.25}}));
    }


    Void test_monotone_reduce() {
        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(2);
        UpperBoxType D = ExactBoxType{{0.0,2.0},{0.0,2.0}};
        EffectiveConstraint c = (-2<=2*x[0]+x[1]<=1);

        ConstraintSolver propagator;
        propagator.verbosity=this->verbosity;

        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c,0));
        ARIADNE_TEST_SAME(D,UpperBoxType({{0.0,0.75},{0.0,2.0}}));
        ARIADNE_TEST_EXECUTE(propagator.box_reduce(D,c,1));
        ARIADNE_TEST_SAME(D,UpperBoxType({{0.0,0.75},{0.0,1.25}}));
    }

    Void test_split() {
        ARIADNE_TEST_WARN("test_split: Not implemented");
    }

    Void test_feasible() {

        List<EffectiveScalarMultivariateFunction> x=EffectiveScalarMultivariateFunction::coordinates(1);
        EffectiveConstraint c = (x[0]-2<=0);
        EffectiveScalarMultivariateFunction f(x[0]+1);
        // EffectiveNumber bound_l = EffectiveNumber::constant(-inf);
        EffectiveNumber bound_u = 0;
        EffectiveConstraint d = (-10<=f<=bound_u);

        List<ValidatedConstraint> constraints;
        constraints.append(c);
        // constraints.append(d);

        ConstraintSolver contractor;

        ARIADNE_TEST_PRINT(constraints);

        ExactBoxType domain1({{1.9375,2.0}});
        ARIADNE_TEST_ASSERT(definitely(contractor.feasible(domain1,constraints).first));

        ExactBoxType domain2({{2.015625,2.5}});
        ARIADNE_TEST_ASSERT(!possibly(contractor.feasible(domain2,constraints).first));

        ExactBoxType domain3({{2.0,2.015625}});
        std::cerr<<"contractor.fesible: "<<contractor.feasible(domain3,constraints)<<"\n";
        ARIADNE_TEST_ASSERT(is_indeterminate(contractor.feasible(domain3,constraints).first));
    }
    Void test_speed()
    {
      const unsigned nVars = 2;
      const unsigned nCons = 2;
      const unsigned maxExp = 2;
      ConstraintSolver contractor;

      RandomPolynomialEngine engine(0, maxExp);
      FunctionDistribution<EffectiveScalarMultivariateFunction, F_TEST_6>
          generator(nVars);

      ExactBoxType domain(nVars,ExactIntervalType(1,10));
      ExactBoxType codomain(nCons,ExactIntervalType(-1,100));

      auto f = generator(engine);
      EffectiveConstraint c1 = (-1<=f<=10);
      f = generator(engine);
      EffectiveConstraint c2 = (-1<=f<=10);
      List<ValidatedConstraint> constraints;
      constraints.append(c1);
      ARIADNE_TEST_PRINT(constraints);

      // std::cerr<<"sol: "<<f(contractor.feasible(domain,constraints).second)<<"\n";

      clock_t s_time = clock();
      ARIADNE_TEST_PRINT(contractor.feasible(domain,constraints));
      clock_t e_time = clock();

      float elapsed_time = static_cast<float>(e_time-s_time) / CLOCKS_PER_SEC;

      std::cerr<<"contractor.feasible() -> time elapsed: "<<elapsed_time<<"\n";

      s_time = clock();
      ARIADNE_TEST_PRINT(__feasible__(domain,constraints));
      e_time = clock();

      elapsed_time = static_cast<float>(e_time-s_time) / CLOCKS_PER_SEC;
      std::cerr<<"__feasible()__ -> time elapsed: "<<elapsed_time<<"\n";
    }
};

Int main(Int argc, const char* argv[]) {
    TestConstraintSolver(get_verbosity(argc,argv)).test();
    return ARIADNE_TEST_FAILURES;
}
