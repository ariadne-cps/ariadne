/***************************************************************************
 *            test_nonlinear_programming.cpp
 *
 *         Copyright  2010  Pieter Collins
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

#include <fstream>
#include <iostream>

#include "config.hpp"

#include "../test.hpp"

#include "algebra/algebra.hpp"
#include "algebra/vector.hpp"
#include "function/formula.hpp"
#include "function/function.hpp"
#include "function/taylor_model.hpp"
#include "geometry/box.hpp"
#include "numeric/numeric.hpp"
#include "solvers/nonlinear_programming.hpp"

using namespace std;
using namespace Ariadne;

#define M_PI 3.14159265358979323846 /* pi */
#define RED "\033[1;31m\n"
#define GREEN "\033[1;32m\n"
#define WHITE "\033[0m\n"

class TestOptimiser {
private:
  std::unique_ptr<OptimiserInterface> optimiser;
  DoublePrecision pr;

public:
  TestOptimiser(const OptimiserInterface &opt) : optimiser(opt.clone()) {}

  Void test() {
    // ARIADNE_TEST_CALL(test_feasibility_check());
    // ARIADNE_TEST_CALL(test_unconstrained_optimisation());
    // ARIADNE_TEST_CALL(test_constrained_optimisation()); //example of failure of Active set methods for quadratic progamming used in sequential quadratic programming
    // ARIADNE_TEST_CALL(test_equality_constrained_optimisation());
    // ARIADNE_TEST_CALL(test_linear_feasibility());
    // ARIADNE_TEST_CALL(test_nonlinear_feasibility());
    // ARIADNE_TEST_CALL(test_nonlinear_equality_feasibility());
    //
    // ARIADNE_TEST_CALL(test_1());
    // ARIADNE_TEST_CALL(test_2());
    // ARIADNE_TEST_CALL(test_3());
    // ARIADNE_TEST_CALL(test_4());
    // ARIADNE_TEST_CALL(test_5());
    // ARIADNE_TEST_CALL(test_6());
    // ARIADNE_TEST_CALL(test_7());
    // ARIADNE_TEST_CALL(test_8());
    // ARIADNE_TEST_CALL(test_9());
    // ARIADNE_TEST_CALL(test_10());
    ARIADNE_TEST_CALL(temp_test());
  }

  Void test_unconstrained_optimisation() {
    // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
    List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveScalarMultivariateFunction x0s = sqr(x[0]);
    EffectiveScalarMultivariateFunction f(x0s * (12 + x0s * (Decimal(6.3) + x0s)) +
                              6 * x[1] * (x[1] - x[0]));
    ARIADNE_TEST_PRINT(f);
    EffectiveVectorMultivariateFunction g(0u, 2u);
    ARIADNE_TEST_PRINT(g);
    ExactBoxType D = ExactBoxType{{-1.0, 2.0}, {-3.0, 5.0}};
    ExactBoxType C = ExactBoxType{};
    ARIADNE_TEST_PRINT(Ariadne::make_tuple(f, D, g, C));

    FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);
    ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, D);
    ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), C);
    FloatDPValue required_accuracy(1e-8);
    // ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
  }

  Void test_equality_constrained_optimisation() {
    List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveScalarMultivariateFunction f = (sqr(x[0]) + sqr(x[1]));
    ARIADNE_TEST_PRINT(f);
    Real a(1.5);
    Real b(0.25);
    EffectiveVectorMultivariateFunction g = {a + x[0] + 2 * x[1] + b * x[0] * x[1]};
    ARIADNE_TEST_PRINT(g);
    ExactIntervalVectorType C = {{0.0, 0.0}};
    ExactBoxType D = ExactBoxType{{-1.0, 2.0}, {-3.0, 5.0}};
    ARIADNE_TEST_PRINT(Ariadne::make_tuple(f, D, g, C));

    // FloatDPValue required_accuracy(1e-7);
    FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);
    ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, D);
    // ARIADNE_TEST_LESS(norm(g(x_optimal)), required_accuracy);
  }

  Void test_constrained_optimisation() {
    List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(3);
    EffectiveScalarMultivariateFunction x0s = sqr(x[0]);
    EffectiveScalarMultivariateFunction f = x0s * (12 + x0s * (Decimal(6.3) + x0s)) +
                                6 * x[1] * (x[1] - x[0]) + x[2];
    ARIADNE_TEST_PRINT(f);
    // EffectiveVectorMultivariateFunction g( (x[0]-1, x[0]+x[1]*x[1], x[1]*x[1]) );
    ExactBoxType D = ExactBoxType{{-1.0, 2.0}, {-3.0, 5.0}, {-3.0, 5.0}};
    ARIADNE_TEST_PRINT(D);
    EffectiveVectorMultivariateFunction g = {2 * x[1] + x[0],
                                 x[0] + x[1] * x[1] - Real(0.875)};
    ARIADNE_TEST_PRINT(g);
    ExactBoxType C = ExactBoxType{{0.0, inf}, {0.0, inf}};
    ARIADNE_TEST_PRINT(C);

    FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);
    ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, D);
    ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), C);
    FloatDPValue required_accuracy(1e-6);
    // ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
  }

  Void test_mixed_constrained_optimisation() {
    List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(3);
    EffectiveScalarMultivariateFunction f(+(sqr(x[0]) + sqr(x[1]) + x[1] * x[2]));
    ARIADNE_TEST_PRINT(f);
    ExactBoxType D = ExactBoxType{{-1.0, 2.0}, {-3.0, 5.0}, {1.25, 2.25}};
    ARIADNE_TEST_PRINT(D);
    EffectiveScalarMultivariateFunction g = x[0] * x[1] - x[0] * Real(1.25);
    EffectiveVectorMultivariateFunction h = {Real(1.5) + x[0] + 2 * x[1] +
                                 Real(0.25) * x[0] * x[1]};
    EffectiveVectorMultivariateFunction gh = join(g, h);
    ARIADNE_TEST_PRINT(gh);
    ExactBoxType C = ExactBoxType{{-1.0, -0.5}, {0.0, 0.0}};
    ARIADNE_TEST_PRINT(C);

    FloatBoundsVector x_optimal = optimiser->minimise(f, D, gh, C);
    FloatDPValue required_accuracy(1e-8);
    ARIADNE_TEST_LESS(norm(h(x_optimal)), required_accuracy);
  }

  Void test_linear_feasibility() {
    // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
    List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveVectorMultivariateFunction g = EffectiveVectorMultivariateFunction(1u, 2 * x[0] + x[1]);
    ARIADNE_TEST_PRINT(g);
    ExactBoxType D = ExactBoxType{{0.0, 2.0}, {0.0, 2.0}};
    ExactBoxType C = ExactBoxType{{-2.0, 1.0}};

    ARIADNE_TEST_ASSERT(optimiser->feasible(D, g, C));
    C = ExactBoxType{{1.0, 1.5}};
    ARIADNE_TEST_ASSERT(optimiser->feasible(D, g, C));
    D = ExactBoxType{{1.0, 1.5}, {0.5, 1.0}};
    ARIADNE_TEST_ASSERT(!optimiser->feasible(D, g, C));
  }

  Void test_nonlinear_feasibility() {
    // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
    List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveVectorMultivariateFunction g = {2 * x[0] + x[1] + x[0] * x[1] / 8};
    ARIADNE_TEST_PRINT(g);
    ExactBoxType D = ExactBoxType{{0.0, 2.0}, {0.0, 2.0}};
    ExactBoxType C = ExactBoxType{{-2.0, 1.0}};

    ARIADNE_TEST_ASSERT(optimiser->feasible(D, g, C));
    return;
    C = ExactBoxType{{1.0, 1.5}};
    ARIADNE_TEST_ASSERT(optimiser->feasible(D, g, C));
    D = ExactBoxType{{1.0, 1.5}, {0.5, 1.0}};
    ARIADNE_TEST_ASSERT(!optimiser->feasible(D, g, C));
  }

  Void test_nonlinear_equality_feasibility() {
    // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
    List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveVectorMultivariateFunction h = {2 * x[0] - x[1] + x[0] * x[1] / 8};
    ARIADNE_TEST_PRINT(h);
    ExactBoxType D = ExactBoxType{{0.0, 2.0}, {0.0, 2.0}};
    ExactBoxType C = ExactBoxType{{0.0, 0.0}};

    ARIADNE_TEST_ASSERT(optimiser->feasible(D, h, C));
  }

  Void test_feasibility_check() {
    EffectiveVectorMultivariateFunction x = EffectiveVectorMultivariateFunction::identity(2);
    ARIADNE_TEST_CONSTRUCT(EffectiveVectorMultivariateFunction, g,
                           ({sqr(x[0]) + 2 * sqr(x[1]) - 1}));
    ARIADNE_TEST_CONSTRUCT(ExactIntervalVectorType, D,
                           ({{-1.0, 1.0}, {-1.0, 1.0}}));
    ARIADNE_TEST_CONSTRUCT(ExactIntervalVectorType, C, ({{0.0, 0.0}}));

    ARIADNE_TEST_CONSTRUCT(FloatBoundsVector, X1,
                           ({{0.30, 0.40}, {0.60, 0.70}}, pr));
    ARIADNE_TEST_ASSERT(
        definitely(optimiser->contains_feasible_point(D, g, C, X1)));

    // The following test fails since it is difficult to find the feasible
    // point in the box.
    ARIADNE_TEST_CONSTRUCT(FloatBoundsVector, X2,
                           ({{0.30, 0.40}, {0.65, 0.65}}, pr));
    ARIADNE_TEST_ASSERT(optimiser->contains_feasible_point(D, g, C, X2));

    ARIADNE_TEST_CONSTRUCT(FloatBoundsVector, X3,
                           ({{0.30, 0.40}, {0.65, 0.68}}, pr));
    ARIADNE_TEST_ASSERT(
        definitely(optimiser->contains_feasible_point(D, g, C, X3)));

    ARIADNE_TEST_CONSTRUCT(ExactFloatVector, x2, ({0.35, 0.655}, pr));
    ARIADNE_TEST_ASSERT(optimiser->validate_feasibility(D, g, C, x2));
  }

  // Optimum at [ 5.10619449, -0.39380551]
      // Hard constraint
      Void test_1()
      {
          List<EffectiveScalarMultivariateFunction> x =
              EffectiveScalarMultivariateFunction::coordinates(2);
          EffectiveScalarMultivariateFunction f(exp(x[0] + x[1]));
          ExactBoxType D = ExactBoxType{{1.0, 10.0}, {-10.0, 10.0}};
          EffectiveVectorMultivariateFunction g = {log(x[0]) + Decimal(1), cos(x[0] + x[1])};
          ExactBoxType C = ExactBoxType{{0.0, inf}, {0.0, inf}};

          // Start time
          clock_t s_time = clock();
          // run code
          FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);
          // End time
          clock_t e_time = clock();

          ExactBoxType test_D = ExactBoxType{{0.99, 10.01}, {-10.1, 10.1}};
          ExactBoxType test_C = ExactBoxType{{-0.01, 10.1}, {-0.01, 10.01}};

          ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, test_D);
          ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), test_C);
          FloatDPValue required_accuracy(1e-6);
          // ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
          // end code

          float elapsed_time =
              static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
          std::cout << "Elapsed time: " << elapsed_time << " sec\n";
          std::cerr<<"f(x_opt)="<<f(x_optimal)<<", g(x_opt)="<<g(x_optimal)<<"\n";
      }

      // Optimum at [ 1]
      // epxp2
      Void test_2()
      {
          List<EffectiveScalarMultivariateFunction> x =
              EffectiveScalarMultivariateFunction::coordinates(1);
          EffectiveScalarMultivariateFunction f(exp(x[0]));
          ExactBoxType D = ExactBoxType{{1, 10}};
          EffectiveVectorMultivariateFunction g = {log(x[0])};
          ExactBoxType C = ExactBoxType{{0.0, +inf}};

          // Start time
          clock_t s_time = clock();
          // run code
          FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);
          // End time
          clock_t e_time = clock();

          ExactBoxType test_D = ExactBoxType{{0.99, 1.01}};
          ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, test_D);
          ExactBoxType test_C = ExactBoxType{{-0.01,+inf}};
          ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), test_C);
          FloatDPValue required_accuracy(1e-6);
          // ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
          // end code

          float elapsed_time =
              static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
          std::cout << "Elapsed time: " << elapsed_time << " sec\n";
          std::cerr<<"f(x_opt)="<<f(x_optimal)<<", g(x_opt)="<<g(x_optimal)<<"\n";
      }

      // Optimum at [ 0]
      // exp1
      Void test_3()
      {
          List<EffectiveScalarMultivariateFunction> x =
              EffectiveScalarMultivariateFunction::coordinates(1);
          EffectiveScalarMultivariateFunction f(exp(x[0]));
          ExactBoxType D = ExactBoxType{{0, 10}};
          EffectiveVectorMultivariateFunction g = {x[0]};
          ExactBoxType C = ExactBoxType{{0, 100}};

          // Start time
          clock_t s_time = clock();
          // run code
          FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);
          // End time
          clock_t e_time = clock();

          ExactBoxType test_D = ExactBoxType{{-0.01, 0.01}};
          ExactBoxType test_C = ExactBoxType{{-0.01, 100.01}};

          ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, test_D);
          ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), test_C);
          FloatDPValue required_accuracy(1e-6);
          // ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
          // end code

          float elapsed_time =
              static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
          std::cout << "Elapsed time: " << elapsed_time << " sec\n";
          std::cerr<<"f(x_opt)="<<f(x_optimal)<<", g(x_opt)="<<g(x_optimal)<<"\n";
      }

      // branin function
      Void test_4(Real a = Real(1.0),
                  Real b = Real(5.1 / (4.0 * (M_PI * M_PI))), // 0.129184509144
                  Real c = Real(5.0 / M_PI),                  // 1.59154943092
                  Real r = Real(6.0), Real s = Real(10.0),
                  Real t = Real(1.0 / (8.0 * M_PI))) // 0.039788735773
      {
          List<EffectiveScalarMultivariateFunction> x =
              EffectiveScalarMultivariateFunction::coordinates(2);
          EffectiveScalarMultivariateFunction f(
              a * ((x[1] - b * (x[0] * x[0]) + c * x[0] - r) *
                   (x[1] - b * (x[0] * x[0]) + c * x[0] - r)) +
              s * (Decimal(1) - t) * cos(x[0]) + s);
          ExactBoxType D = ExactBoxType{{-5.0, 10.0}, {0.0, 15.0}};
          EffectiveVectorMultivariateFunction g(0u, 2u);
          ExactBoxType C = ExactBoxType{};

          // Start time
          clock_t s_time = clock();
          // run code
          FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);
          // End time
          clock_t e_time = clock();

          ExactBoxType test_D = ExactBoxType{{-5.01, 10.01}, {-0.01, 15.01}};
          ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, test_D);
          // ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
          // end code

          float elapsed_time =
              static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
          std::cout << "Elapsed time: " << elapsed_time << " sec\n";
          std::cerr<<"f(x_opt)="<<f(x_optimal)<<", g(x_opt)="<<g(x_optimal)<<"\n";
      }

      Void test_5()
  {
      List<EffectiveScalarMultivariateFunction> x =
          EffectiveScalarMultivariateFunction::coordinates(2);
      EffectiveScalarMultivariateFunction f((x[0] * x[0]) + (x[1] * x[1]));
      ExactBoxType D = ExactBoxType{{-10, 10}, {-1, 10}};
      EffectiveVectorMultivariateFunction g = {x[0] + x[1] - 1, exp(x[0]),
                                               x[1]*x[0],x[1]*x[1]+x[0],sin(x[1])};
      ExactBoxType C = ExactBoxType{{0, 100}, {0, 100}, {0, 100},{0,10},{0,10}};

      // Start time
      clock_t s_time = clock();
      // run code
      FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);
      // End time
      clock_t e_time = clock();

      ExactBoxType test_D = ExactBoxType{{-0.01,10.01},{-0.01,10.01}};
      ExactBoxType test_C = ExactBoxType{{-0.01, 100.01}, {-0.01, 100.01}, {-0.01, 100.01},{-0.01,10.01},{-0.01,10.01}};
      ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, test_D);
      ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), test_C);
      FloatDPValue required_accuracy(1e-6);
      // ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
      // end code

      float elapsed_time =
          static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
      std::cout << "Elapsed time: " << elapsed_time << " sec\n";
      std::cerr<<"f(x_opt)="<<f(x_optimal)<<", g(x_opt)="<<g(x_optimal)<<"\n";

  }


  Void test_6()
{
  List<EffectiveScalarMultivariateFunction> x =
      EffectiveScalarMultivariateFunction::coordinates(4);
  EffectiveScalarMultivariateFunction f(exp(x[0]*x[1])-(x[1]*x[2])+10*(x[3]*x[3]));
  ExactBoxType D = ExactBoxType{{0, 6}, {0, 6}, {0,6}, {0,6}};
  EffectiveVectorMultivariateFunction g = {x[0],x[1]};
  ExactBoxType C = ExactBoxType{{0, 100},{0,100}};

  // Start time
  clock_t s_time = clock();
  // run code
  FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);
  // End time
  clock_t e_time = clock();

  ExactBoxType test_D = ExactBoxType{{-0.01, 6.01}, {-0.01, 6.01}, {-0.01,6.01}, {-0.01,6.01}};
  ExactBoxType test_C = ExactBoxType{{-0.01, 100.01}, {-0.01, 100.01}};
  ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, test_D);
  ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), test_C);
  FloatDPValue required_accuracy(1e-6);
  // ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
  // end code

  float elapsed_time =
      static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
  std::cout << "Elapsed time: " << elapsed_time << " sec\n";
  std::cerr<<"f(x_opt)="<<f(x_optimal)<<", g(x_opt)="<<g(x_optimal)<<"\n";

}

Void test_7()
{
List<EffectiveScalarMultivariateFunction> x =
    EffectiveScalarMultivariateFunction::coordinates(2);
EffectiveScalarMultivariateFunction f(exp(-x[0])*log(x[1])+x[0]);
ExactBoxType D = ExactBoxType{{-10, 10}, {1, 10}};
EffectiveVectorMultivariateFunction g = {x[0]*x[1],x[1]*x[1],log(x[1]),Decimal(3)*x[0]+x[1],sin(x[1]),sqrt(x[1])};
ExactBoxType C = ExactBoxType{{-10, 100},{0,100},{-1,2},{-10,10},{-3,3},{0,10}};

// Start time
clock_t s_time = clock();
// run code
FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);
// End time
clock_t e_time = clock();

ExactBoxType test_D = ExactBoxType{{-10.01, 10.01}, {0.99, 10.01}};
ExactBoxType test_C = ExactBoxType{{-10.01, 100.01},{-0.01,100.01},{-1.01,2.01},{-10.01,10.01},{-3.01,3.01},{-0.01,10.01}};
ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, test_D);
ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), test_C);
FloatDPValue required_accuracy(1e-6);
// ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
// end code

float elapsed_time =
    static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
std::cout << "Elapsed time: " << elapsed_time << " sec\n";
std::cerr<<"f(x_opt)="<<f(x_optimal)<<", g(x_opt)="<<g(x_optimal)<<"\n";

}

Void test_8()
{
List<EffectiveScalarMultivariateFunction> x =
    EffectiveScalarMultivariateFunction::coordinates(2);
EffectiveScalarMultivariateFunction f(x[0]*(sin(x[0])+cos(x[1]))-exp(-x[1]));
ExactBoxType D = ExactBoxType{{-10, 10}, {0, 10}};
EffectiveVectorMultivariateFunction g = {x[0]*x[0],x[1]*x[1],x[0]*x[1]};
ExactBoxType C = ExactBoxType{{0,10},{0,10},{0,8}};

// Start time
clock_t s_time = clock();
// run code
FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);
// End time
clock_t e_time = clock();

ExactBoxType test_D = ExactBoxType{{-10.01, 10.01}, {-0.01, 10.01}};
ExactBoxType test_C = ExactBoxType{{-0.01, 10.01},{-0.01,100.01},{-10.01,8.01}};
ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, test_D);
ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), test_C);
FloatDPValue required_accuracy(1e-6);
// ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
// end code

float elapsed_time =
    static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
std::cout << "Elapsed time: " << elapsed_time << " sec\n";
std::cerr<<"f(x_opt)="<<f(x_optimal)<<", g(x_opt)="<<g(x_optimal)<<"\n";

}

Void test_9()
{
List<EffectiveScalarMultivariateFunction> x =
    EffectiveScalarMultivariateFunction::coordinates(2);
EffectiveScalarMultivariateFunction f(log(x[1])*exp(x[0]));
ExactBoxType D = ExactBoxType{{-10, 10}, {1, 10}};
EffectiveVectorMultivariateFunction g = {x[1]*x[0]};
ExactBoxType C = ExactBoxType{{-1, 100}};

// Start time
clock_t s_time = clock();
// run code
FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);
// End time
clock_t e_time = clock();

ExactBoxType test_D = ExactBoxType{{-10.01, 10.01}, {0.99, 10.01}};
ExactBoxType test_C = ExactBoxType{{-1.01, 100.01}};
ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, test_D);
ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), test_C);
FloatDPValue required_accuracy(1e-6);
// ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
// end code

float elapsed_time =
    static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
std::cout << "Elapsed time: " << elapsed_time << " sec\n";
std::cerr<<"f(x_opt)="<<f(x_optimal)<<", g(x_opt)="<<g(x_optimal)<<"\n";

}

Void test_10()
{
List<EffectiveScalarMultivariateFunction> x =
    EffectiveScalarMultivariateFunction::coordinates(3);
EffectiveScalarMultivariateFunction f(x[0]+(x[1])*x[2]);
ExactBoxType D = ExactBoxType{{0, 10}, {0, 10},{0,10}};
EffectiveVectorMultivariateFunction g = {x[0]+x[1]+x[2]};
ExactBoxType C = ExactBoxType{{-10, 10}};

// Start time
clock_t s_time = clock();
// run code
FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);
// End time
clock_t e_time = clock();

ExactBoxType test_D = ExactBoxType{{-0.01, 10.01}, {-0.01, 10.01},{-0.01,10.01}};
ExactBoxType test_C = ExactBoxType{{-10.01,10.01}};
ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, test_D);
ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), test_C);
FloatDPValue required_accuracy(1e-6);
// ARIADNE_TEST_LESS(norm(x_optimal),required_accuracy);
// end code

float elapsed_time =
    static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
std::cout << "Elapsed time: " << elapsed_time << " sec\n";
std::cerr<<"f(x_opt)="<<f(x_optimal)<<", g(x_opt)="<<g(x_optimal)<<"\n";

}

Void temp_test() {
  // Test the feasibility of x0>0, x1>0, 2x1+x2<1 using box [0,2]x[0,2]
  unsigned nVars = 2;
  List<EffectiveScalarMultivariateFunction> x = EffectiveScalarMultivariateFunction::coordinates(nVars);

  EffectiveScalarMultivariateFunction x0=x[0];
  EffectiveScalarMultivariateFunction x1=x[1];
  // EffectiveScalarMultivariateFunction x2=x[2];
  // EffectiveScalarMultivariateFunction x3=x[3];
  // EffectiveScalarMultivariateFunction x4=x[4];
  // EffectiveScalarMultivariateFunction x5=x[5];
  // EffectiveScalarMultivariateFunction x6=x[6];

  EffectiveScalarMultivariateFunction f(
    pow(x0+1,2)+1
  );
  ARIADNE_TEST_PRINT(f);
  EffectiveVectorMultivariateFunction g = {
    13*x0*x1+76*x1*x0+46*x1*x1
  };
  ARIADNE_TEST_PRINT(g);
  ExactBoxType D(nVars, ExactIntervalType(-10, 11));
  ExactBoxType testD(nVars, ExactIntervalType(-10.01, 11.01));
  ExactBoxType C(1, ExactIntervalType(-1, +inf));
  ExactBoxType testC(1, ExactIntervalType(-1.01, 100.01));
  ARIADNE_TEST_PRINT(Ariadne::make_tuple(f, D, g, C));

  FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);
  ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, testD);
  ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), testC);

  std::cerr<<"f(opt) = "<<f(x_optimal)<<"\n";
}
};

Int main(Int argc, const char *argv[]) {
  Nat optimiser_verbosity = get_verbosity(argc, argv);

  NonlinearSQPOptimiser nlsqp;
  nlsqp.verbosity = optimiser_verbosity;
  TestOptimiser(nlsqp).test();
  // return ARIADNE_TEST_FAILURES;
  // // //
  // NonlinearInfeasibleInteriorPointOptimiser nliipo;
  // nliipo.verbosity = optimiser_verbosity;
  // TestOptimiser(nliipo).test();
  //
  // return ARIADNE_TEST_FAILURES;

  NonlinearMixedOptimiser nlhop;
  nlhop.verbosity = optimiser_verbosity;
  TestOptimiser(nlhop).test();
  return ARIADNE_TEST_FAILURES;

  NonlinearInteriorPointOptimiser nlo;
  nlo.verbosity = optimiser_verbosity;
  TestOptimiser(nlo).test();
  return ARIADNE_TEST_FAILURES;
  ApproximateOptimiser appo;
  appo.verbosity = optimiser_verbosity;
  TestOptimiser(appo).test_nonlinear_equality_feasibility();

  IntervalOptimiser ivlo;
  ivlo.verbosity = optimiser_verbosity;
  // TestOptimiser(ivlo).test_nonlinear_equality_feasibility();
  return ARIADNE_TEST_FAILURES;
}
