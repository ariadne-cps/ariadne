/***************************************************************************
 *            test_function_generator_nonlinear_programming.cpp
 *
 *          Copyright  2018  Nicola Dessì
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

#include <chrono>
#include <fstream>
#include <iostream>

#include "config.hpp"

#include "../test.hpp"

#include "solvers/solver.hpp"

#include "algebra/algebra.hpp"
#include "algebra/vector.hpp"
#include "function/formula.hpp"
#include "function/function.hpp"
#include "function/function_model.hpp"
#include "function/taylor_model.hpp"
#include "geometry/box.hpp"
#include "numeric/numeric.hpp"
#include "solvers/nonlinear_programming.hpp"

#include "function/symbolic_function.hpp"

#include "DataStream.hpp"
#include "RandomFunction.hpp"

using namespace std;
using namespace Ariadne;

#define M_PI 3.14159265358979323846 /* pi */
#define RED "\033[1;31m\n"
#define GREEN "\033[1;32m\n"
#define WHITE "\033[0m\n"

#define ci(a) (static_cast<int>(a))
#define cd(a) (static_cast<double>(a))
#define cfdp(a) (static_cast<FloatDP>(a))

inline char activity_symbol(SizeType step) {
  switch (step % 4) {
  case 0:
    return '\\';
  case 1:
    return '|';
  case 2:
    return '/';
  default:
    return '-';
  }
}

class TestOptimiser {
private:
  std::unique_ptr<OptimiserInterface> optimiser;
  // DoublePrecision pr;

public:
  TestOptimiser(const OptimiserInterface &opt) : optimiser(opt.clone()) {}

  Void test(const std::string &name, int seed) {
    optimiser->problem_name = "h_v";
    optimiser->complexity_order = 3u;
    variables(name + "_h_v", seed, 3, 80u);
    optimiser->complexity_order = 8u;
    optimiser->problem_name = "h_v_h_s";
    variables(name + "_h_v_h_s", seed, 8, 40u);
    optimiser->complexity_order = 14u;
    optimiser->problem_name = "h_s";
    variables(name + "_h_s", seed, 14, 20u);

    optimiser->complexity_order = 3u;
    optimiser->problem_name = "c_h_c";
    constraints_low_active(name + "_c_h_c", seed, 3, 20u, 20u);
    optimiser->complexity_order = 10u;
    optimiser->problem_name = "c_h_s";
    constraints_low_active(name + "_c_h_s", seed, 10, 20u, 10u);
    optimiser->complexity_order = 5u;
    optimiser->problem_name = "c_h_c_h_s";
    constraints_low_active(name + "_c_h_c_h_s", seed, 5, 20u, 20u);

    optimiser->complexity_order = 3u;
    optimiser->problem_name = "c_h_c_h_a";
    optimiser->active_set_high = true;
    constraints_high_active(name + "_c_h_c_h_a", seed, 3, 20u, 20u);
    optimiser->complexity_order = 10u;
    optimiser->problem_name = "c_h_s_h_a";
    optimiser->active_set_high = true;
    constraints_high_active(name + "_c_h_s_h_a", seed, 10, 20u, 10u);
    optimiser->complexity_order = 5u;
    optimiser->problem_name = "c_h_c_h_s_h_a";
    optimiser->active_set_high = true;
    constraints_high_active(name + "_c_h_c_h_s_h_a", seed, 5, 20u, 20u);
  }

  // VARIABLES
  Void variables(const std::string &name, int seed, int order,
                 unsigned maxNVars) {

    std::cerr << "Testing " << name << "...\n";
    DataStream ds = csv::make_datastream();
    ds.setHeader({"#sequence_" + name, "nVars_" + name, "order_" + name,
                  "correct_" + name, "f_obj_" + name, "time_" + name,
                  "status_" + name});
    ds.open(optimiser->problem_name + "/" + name + ".csv");

    EffectiveScalarMultivariateFunction f;
    FloatBoundsVector x_optimal;

    float elapsed_time;
    FloatDP f_obj;
    int index = 0;
    for (unsigned nVars = 1; nVars <= maxNVars; ++nVars, ++index) {
      std::cerr << "\r[" << activity_symbol(static_cast<unsigned>(index))
                << "] " << static_cast<int>((index * 100) / (ci(maxNVars)))
                << "% " << std::flush;
      try {
        // std::cerr<<nVars<<", ";
        std::vector<float> bounds = {1, 1, 1, 1, 1, 1};
        single_variables(index, nVars, order, f_obj, elapsed_time, f, x_optimal,
                         bounds);
        ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order), cfdp(1),
                            f_obj, cfdp(elapsed_time), 1});
      } catch (DomainException &de) {
        ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order), cfdp(-9999),
                            f_obj, cfdp(elapsed_time), -9999});
      } catch (std::runtime_error &e) {
        std::cerr << e.what() << "\n";
        ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order), cfdp(-15),
                            f_obj, cfdp(elapsed_time), -15});
      }
    }

    ds.close();
    std::cerr << "\r[" << activity_symbol(static_cast<unsigned>(index)) << "] "
              << 100 << "% " << std::flush;
    std::cerr << "\n";
  }
  Void single_variables(int seed, unsigned nVars, int _order, FloatDP &f_obj,
                        float &elapsed_time,
                        EffectiveScalarMultivariateFunction &f,
                        FloatBoundsVector &x_optimal,
                        std::vector<float> _bounds) {
    RandomPolynomialEngine engine(seed, _order);
    FunctionDistribution<EffectiveScalarMultivariateFunction, F_TEST_1>
        generator(nVars, _bounds);
    FunctionDistribution<EffectiveScalarMultivariateFunction, F_TEST_1>
        generator_g(nVars, _bounds);
    RandomPolynomialEngine engine_g(seed, 1);

    ExactBoxType D(nVars + 1, ExactIntervalType(-10, 10));
    ExactBoxType test(nVars + 1, ExactIntervalType(-10.01, 10.01));
    ExactBoxType C = {};
    ExactBoxType testC = {};
    EffectiveVectorMultivariateFunction g(0u, nVars + 1);
    // ExactBoxType C = {{-10,10},{-20,1},{-20,20},{-40,40}};
    // ExactBoxType testC =
    // {{-10.01,10.01},{-20.01,1.01},{-20.01,20.01},{-40.01,40.01}};

    //
    f = generator(engine);

    clock_t s_time = clock();
    // run code
    x_optimal = optimiser->minimise(f, D, g, C);
    // End time
    clock_t e_time = clock();
    f_obj = cfdp(f(x_optimal));

    if (!decide(element(x_optimal, test)) ||
        !decide(element(g(x_optimal), testC))) {
      throw DomainException("Unfeasible");
    }

    elapsed_time = static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
  }

  // CONSTRAINTS
  Void constraints_high_active(const std::string &name, int seed, int order,
                               unsigned maxNVars, unsigned maxNCons) {
    std::cerr << "Testing " << name << "...\n";
    DataStream ds = csv::make_datastream();
    ds.setHeader({"#sequence_" + name, "nVars_" + name, "order_" + name,
                  "correct_" + name, "f_obj_" + name, "time_" + name,
                  "status_" + name});
    ds.open(optimiser->problem_name + "/" + name + ".csv");

    EffectiveScalarMultivariateFunction f;
    FloatBoundsVector x_optimal;

    float elapsed_time;
    FloatDP f_obj;
    int index = 0;
    for (unsigned nVars = 1; nVars <= maxNVars; ++nVars) {
      for (unsigned nCons = 5; nCons <= maxNCons; ++nCons, ++index) {
        std::cerr << "\r[" << activity_symbol(static_cast<unsigned>(index))
                  << "] "
                  << static_cast<int>((index * 100) / (ci(maxNVars * maxNCons)))
                  << "% " << std::flush;
        try {
          std::vector<float> bounds(nCons, 1);
          single_constraints_high_active(index, nVars, order, f_obj,
                                         elapsed_time, f, x_optimal, bounds,
                                         nCons);
          ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order), cfdp(1),
                              f_obj, cfdp(elapsed_time), 1});
        } catch (DomainException &de) {
          ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order),
                              cfdp(-9999), f_obj, cfdp(elapsed_time), -9999});
        } catch (std::runtime_error &e) {
          // std::cerr<<e.what()<<"\n";
          ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order), cfdp(-15),
                              f_obj, cfdp(elapsed_time), -15});
        }
      }
    }

    ds.close();
    std::cerr << "\r[" << activity_symbol(static_cast<unsigned>(index)) << "] "
              << 100 << "% " << std::flush;
    std::cerr << "\n";
  }
  Void single_constraints_high_active(int seed, unsigned nVars, int _order,
                                      FloatDP &f_obj, float &elapsed_time,
                                      EffectiveScalarMultivariateFunction &f,
                                      FloatBoundsVector &x_optimal,
                                      std::vector<float> _bounds, int nCons) {
    RandomPolynomialEngine engine(seed, _order);
    FunctionDistribution<EffectiveScalarMultivariateFunction, F_TEST_1>
        generator(nVars, _bounds);
    FunctionDistribution<EffectiveScalarMultivariateFunction, F_TEST_1>
        generator_g(nVars, _bounds);
    RandomPolynomialEngine engine_g(seed, _order);

    ExactBoxType D(nVars + 1, ExactIntervalType(-10, 10));
    ExactBoxType test(nVars + 1, ExactIntervalType(-10.01, 10.01));
    ExactBoxType C(nCons, ExactIntervalType(-5, +10));
    ExactBoxType testC(nCons, ExactIntervalType(-5.01, +10.01));

    //
    f = generator(engine);
    EffectiveVectorMultivariateFunction g(nCons, nVars + 1);
    for (unsigned i = 0; i < nCons; ++i) {
      g[i] = generator_g(engine_g, seed * nCons + i);
    }

    // std::cerr<<"g: "<<g<<"\n";
    clock_t s_time = clock();
    // run code
    x_optimal = optimiser->minimise(f, D, g, C);
    // End time
    clock_t e_time = clock();
    f_obj = cfdp(f(x_optimal));

    if (!decide(element(x_optimal, test)) ||
        !decide(element(g(x_optimal), testC))) {
      throw DomainException("Unfeasible");
    }

    elapsed_time = static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
    // std::cerr<<x_optimal<<"\n\n";
  }
  // ACTIVE SET
  Void constraints_low_active(const std::string &name, int seed, int order,
                              unsigned maxNVars, unsigned maxNCons) {
    std::cerr << "Testing " << name << "...\n";
    DataStream ds = csv::make_datastream();
    ds.setHeader({"#sequence_" + name, "nVars_" + name, "order_" + name,
                  "correct_" + name, "f_obj_" + name, "time_" + name,
                  "status_" + name});
    ds.open(optimiser->problem_name + "/" + name + ".csv");

    EffectiveScalarMultivariateFunction f;
    FloatBoundsVector x_optimal;

    float elapsed_time;
    FloatDP f_obj;
    int index = 0;
    for (unsigned nVars = 1; nVars <= maxNVars; ++nVars) {
      for (unsigned nCons = 5; nCons <= maxNCons; ++nCons, ++index) {
        std::cerr << "\r[" << activity_symbol(static_cast<unsigned>(index))
                  << "] "
                  << static_cast<int>((index * 100) / (ci(maxNVars * maxNCons)))
                  << "% " << std::flush;
        try {
          std::vector<float> bounds(nCons, 1);
          single_constraints_low_active(index, nVars, order, f_obj,
                                        elapsed_time, f, x_optimal, bounds,
                                        nCons);
          ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order), cfdp(1),
                              f_obj, cfdp(elapsed_time), 1});
        } catch (DomainException &de) {
          ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order),
                              cfdp(-9999), f_obj, cfdp(elapsed_time), -9999});
        } catch (std::runtime_error &e) {
          // std::cerr<<e.what()<<"\n";
          ds.append<FloatDP>({cfdp(index), cfdp(nVars), cfdp(order), cfdp(-15),
                              f_obj, cfdp(elapsed_time), -15});
        }
      }
    }

    ds.close();
    std::cerr << "\r[" << activity_symbol(static_cast<unsigned>(index)) << "] "
              << 100 << "% " << std::flush;
    std::cerr << "\n";
  }
  Void single_constraints_low_active(int seed, unsigned nVars, int _order,
                                     FloatDP &f_obj, float &elapsed_time,
                                     EffectiveScalarMultivariateFunction &f,
                                     FloatBoundsVector &x_optimal,
                                     std::vector<float> _bounds, int nCons) {
    RandomPolynomialEngine engine(seed, _order);
    FunctionDistribution<EffectiveScalarMultivariateFunction, F_TEST_1>
        generator(nVars, _bounds);
    FunctionDistribution<EffectiveScalarMultivariateFunction, F_TEST_1>
        generator_g(nVars, _bounds);
    RandomPolynomialEngine engine_g(seed, _order);

    ExactBoxType D(nVars + 1, ExactIntervalType(-10, 10));
    ExactBoxType test(nVars + 1, ExactIntervalType(-10.01, 10.01));
    ExactBoxType C(nCons, ExactIntervalType(-99999999999999, +99999999999999));
    ExactBoxType testC(
        nCons, ExactIntervalType(-99999999999999.01, +99999999999999.01));

    //
    f = generator(engine);
    EffectiveVectorMultivariateFunction g(nCons, nVars + 1);
    for (unsigned i = 0; i < nCons; ++i) {
      g[i] = generator_g(engine_g, seed * nCons + i);
    }

    // std::cerr<<"g: "<<g<<"\n";
    clock_t s_time = clock();
    // run code
    x_optimal = optimiser->minimise(f, D, g, C);
    // End time
    clock_t e_time = clock();
    f_obj = cfdp(f(x_optimal));

    if (!decide(element(x_optimal, test)) ||
        !decide(element(g(x_optimal), testC))) {
      throw DomainException("Unfeasible");
    }

    elapsed_time = static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
    // std::cerr<<x_optimal<<"\n\n";
  }
};

<<<<<<< HEAD
Int main(Int argc, const char *argv[]) {
  Nat optimiser_verbosity = get_verbosity(argc, argv);
  int seed = 0;

  NonlinearSQPOptimiser nlsqp;
  nlsqp.verbosity = optimiser_verbosity;
  TestOptimiser(nlsqp).test("nlsqp", seed);

  NonlinearInteriorPointOptimiser nlipm;
  nlipm.verbosity = optimiser_verbosity;
  TestOptimiser(nlipm).test("nlipm", seed);

  NonlinearMixedOptimiser nlmop;
  nlmop.verbosity = optimiser_verbosity;
  TestOptimiser(nlmop).test("nlmop", seed);

  return ARIADNE_TEST_FAILURES;
=======
template<class X>
void dfs(const Formula<X>& a)
{
  switch(a.kind()) {
      case OperatorKind::UNARY:
        std::cerr<<"Compute unary: "<<compute(a.op(),Real(1)).get_d()<<"\n";
      case OperatorKind::SCALAR:
      case OperatorKind::GRADED:
      case OperatorKind::COORDINATE:
      case OperatorKind::NULLARY:
        std::cerr<<"Found operator "<<a.kind()<<" is "<<a.op()<<"\n";

        return;
      case OperatorKind::BINARY:
        dfs(a.arg1());
        dfs(a.arg2());
        return;
      default: ARIADNE_FAIL_MSG("WTF\n");
    }
}
void bfs();

ValidatedVectorMultivariateFunctionModelDP
forward_backwardCtr(const EffectiveVectorMultivariateFunction& f,
  const ValidatedVectorMultivariateFunctionModelDP& id,
  const ValidatedVectorMultivariateFunctionModelDP& h)
  {
    auto tmp = make_formula(f);

    std::cerr<<"tmp: "<<tmp<<"\n";

    auto formula = tmp[0];

    std::cerr<<"formula: "<<formula<<"\n";

    Vector<Interval<Value<FloatDP> > > domain = f.domain();
    Vector<Real> a = {1,1};

    auto evaluation = evaluate(formula,a);

    std::cerr<<evaluate(formula,a)<<"\n";

    dfs(formula);

    return h;
  }

FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory();

Int main(Int argc, const char *argv[])
{
    // Nat optimiser_verbosity = get_verbosity(argc, argv);
    //
    // // unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    // unsigned seed = 0;
    //
    //
    // NonlinearSQPOptimiser nlsqp;
    // nlsqp.verbosity = optimiser_verbosity;
    // TestOptimiser(nlsqp).test("nlsqp", seed);
    // // return ARIADNE_TEST_FAILURES;
    //
    // NonlinearInteriorPointOptimiser nlipm;
    // nlipm.verbosity = optimiser_verbosity;
    // TestOptimiser(nlipm).test("nlipm", seed);
    // // return ARIADNE_TEST_FAILURES;
    //
    //
    // NonlinearInfeasibleInteriorPointOptimiser nliipm;
    // nliipm.verbosity = optimiser_verbosity;
    // TestOptimiser(nliipm).test("nliipm", seed);
    // // return ARIADNE_TEST_FAILURES;
    //
    //
    // NonlinearMixedOptimiser nlmop;
    // nlmop.verbosity = optimiser_verbosity;
    // TestOptimiser(nlmop).test("nlmop",seed);

    std::shared_ptr< SolverBase > solver_ptr=std::shared_ptr<SolverBase>(new IntervalNewtonSolver(1e-8,12));
    List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(2);

    ExactBoxType ip = {{1,1},{1,1}};
    ExactBoxType ix = {{1,1},{1,1}};

    ExactBoxType C(1, ExactIntervalType(-1, 100));
    ExactBoxType testC(1, ExactIntervalType(-1.01, 100.01));

    EffectiveVectorMultivariateFunction f={x[0] + 2 * (sin(x[1]) + x[0]) - x[1]};
    ValidatedVectorMultivariateFunctionModelDP id(solver_ptr->function_factory().create_identity(ip));
    ValidatedVectorMultivariateFunctionModelDP h(solver_ptr->function_factory().create_constants(ip,cast_singleton(ix)));
    forward_backwardCtr(f,id,h);
    return ARIADNE_TEST_FAILURES;
>>>>>>> Small fixes. Implemented temporary __feasible__ function to test barrier method.
}
