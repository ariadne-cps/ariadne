/***************************************************************************
 *            test_sequential_quadratic_programming.cpp
 *
 *         Copyright  2010  Nicola Dessì
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

class TestOptimiser
{
private:
  std::unique_ptr<OptimiserInterface> optimiser;
  DoublePrecision                     pr;

public:
  TestOptimiser(const OptimiserInterface &opt) : optimiser(opt.clone()) {}
  Void test()
  {
    // ARIADNE_TEST_CALL(benchmark_allinit());
    // ARIADNE_TEST_CALL(benchmark_extrasim());
    // ARIADNE_TEST_CALL(benchmark_hs54());
    ARIADNE_TEST_CALL(benchmark_booth());
    ARIADNE_TEST_CALL(benchmark_himmelbc());
    // ARIADNE_TEST_CALL(benchmark_hs44());
  }

  /*
  Results from AmplOnline with LOQO:

  LOQO 7.03: optimal solution (11 iterations, 11 evaluations)
  primal objective 1.000000008
  dual objective 1
  f = 1

  x = 7.97704e-09
  y = 1
  */
  Void benchmark_extrasim()
  {

    List<EffectiveScalarMultivariateFunction> x =
        EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveScalarMultivariateFunction f(x[0] + 1);
    ExactBoxType D = ExactBoxType{{0, +inf}, {-inf, +inf}};
    EffectiveVectorMultivariateFunction g = {x[0] + 2 * x[1] - 2};
    ExactBoxType                        C = ExactBoxType{{0, 0}};
    FloatBoundsVector x_optimal           = optimiser->minimise(f, D, g, C);

    ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, D);
    ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), C);
  }

  /*
  Results from AmplOnline with LOQO:

    LOQO 7.03: optimal solution (11 iterations, 11 evaluations)
    primal objective 24.46390782
      dual objective 24.46390771
    f = 24.4639

    x [*] :=
    1  -0.635594
    2   1
    3   7.78044e-09
    4   2
  */
  Void benchmark_allinit()
  {
    List<EffectiveScalarMultivariateFunction> x =
        EffectiveScalarMultivariateFunction::coordinates(4);
    EffectiveScalarMultivariateFunction f(
        x[2] - 1 + pow(x[0], 2) + pow(x[1], 2) + pow((x[2] + x[3]), 2) +
        pow(sin(x[2]), 2) + pow(x[0], 2) * pow(x[1], 2) + x[3] - 3 +
        pow(sin(x[2]), 2) + pow((x[3] - 1), 2) + pow((pow(x[1], 2)), 2) +
        pow((pow(x[2], 2) + pow((x[3] + x[0]), 2)), 2) +
        pow((x[0] - 4 + pow(sin(x[3]), 2) + pow(x[1], 2) * pow(x[2], 2)), 2) +
        pow(sin(x[3]), 4));
    ExactBoxType D = ExactBoxType{{-inf, +inf}, {1, +inf}, {0, 1}, {2, 2}};
    EffectiveVectorMultivariateFunction g(0u, 4u);
    ExactBoxType                        C = {};
    FloatBoundsVector x_optimal           = optimiser->minimise(f, D, g, C);

    ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, D);
    ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), C);
  }

  /*
    Precision test-case:
      smaller rtol is in SQP framework, better is the result. Asyntotically it
    is correct.
    With rtol = 1.49012e-10 the result  is
    [-0.00000000000000:0.00000000000364] which is close to 0 (error of order
    2^12)

    LOQO 7.03: optimal solution (33 iterations, 44 evaluations)
      primal objective -0.9957643067
        dual objective -0.9957643126
      x [*] :=
      1     14007.3
      2         0.898165
      3   4611680
      4        10
      5         0.001
      6  56006500
      ;

      Obj = -0.995764
  */
  Void benchmark_hs54()
  {
    List<EffectiveScalarMultivariateFunction> x =
        EffectiveScalarMultivariateFunction::coordinates(6);
    EffectiveScalarMultivariateFunction h(

        (pow(x[0] - Decimal(1e6), 2) / Decimal(6.4e13) +
         (x[0] - Decimal(1e4)) * (x[1] - 1) / Decimal(2e4) + pow(x[1] - 1, 2)) *
            pow(x[2] - Decimal(2e6), 2) / (Decimal(0.96) * Decimal(4.9e13)) +
        pow(x[3] - 10, 2) / Decimal(2.5e3) +
        pow(x[4] - Decimal(1e-3), 2) / Decimal(2.5e-3) +
        pow(x[5] - Decimal(1e8), 2) / Decimal(2.5e17)

    );
    EffectiveScalarMultivariateFunction f(-exp(-h / 2));

    ExactBoxType D = ExactBoxType{{0, 2E+4}, {-10, +10}, {0, 1E+7},
                                  {0, 20},   {-1, +1},   {0, 2E+8}};
    EffectiveVectorMultivariateFunction g = {x[0] + Decimal(4E+3) * x[1] -
                                             Decimal(1.76E+4)};
    ExactBoxType                        C = {{0, 0}};
    FloatBoundsVector x_optimal           = optimiser->minimise(f, D, g, C);

    std::cout << "f(x_optimal): " << f(x_optimal) << "\n";
    ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, D);
    ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), C);
  }

  /*
  Empty bound on variables test:
    If there aren't bounds on variables this will generate an empty vector which
  generate an error on line nonlinear_programming.cpp:3442.
  Solution: if vector is empty generate an empty vector, or implement on Vector
  class level the operation with empty sets.

  LOQO 7.03: optimal solution (2 QP iterations, 3 evaluations)
  primal objective 7.222223417e-21
    dual objective 7.222223417e-21
  f = 7.22222e-21

  x [*] :=
  1  1
  2  3
  ;
  */
  Void benchmark_booth()
  {
    List<EffectiveScalarMultivariateFunction> x =
        EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveScalarMultivariateFunction f(pow(x[0] + 2 * x[1] - 7, 2) +
                                          pow(2 * x[0] + x[1] - 5, 2));
    ExactBoxType D = ExactBoxType{{-inf, +inf}, {-inf, +inf}};
    EffectiveVectorMultivariateFunction g(0u, 2u);
    ExactBoxType                        C = {};
    FloatBoundsVector x_optimal           = optimiser->minimise(f, D, g, C);

    ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, D);
    ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), C);
  }

  /*

  If there aren't bounds on variables this will generate an empty vector which
  generate an error on line nonlinear_programming.cpp:3442.
  Solution: if vector is empty generate an empty vector, or implement on Vector
  class level the operation with empty sets.

  LOQO 7.03: optimal solution (10 iterations, 11 evaluations)
  primal objective 1.009641594e-15
    dual objective 1.009641594e-15
  f = 1.00964e-15

  x [*] :=
  1  3
  2  2
  ;
  */
  Void benchmark_himmelbc()
  {

    List<EffectiveScalarMultivariateFunction> x =
        EffectiveScalarMultivariateFunction::coordinates(2);
    EffectiveScalarMultivariateFunction f(pow(x[1] - 11 + pow(x[0], 2), 2) +
                                          pow(x[0] - 7 + pow(x[1], 2), 2));
    ExactBoxType D = ExactBoxType{{-inf, +inf}, {-inf, +inf}};
    EffectiveVectorMultivariateFunction g(0u, 2u);
    ExactBoxType                        C = {};
    FloatBoundsVector x_optimal           = optimiser->minimise(f, D, g, C);

    ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, D);
    ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), C);
  }

  /*
    
    LOQO 7.03: optimal solution (16 QP iterations, 16 evaluations)
      primal objective -14.99999999
        dual objective -15.00000003
      x [*] :=
      1  1.55544e-10
      2  3
      3  1.17188e-09
      4  4
      ;
      
      f = -15
  */
  Void benchmark_hs44()
  {
    List<EffectiveScalarMultivariateFunction> x =
        EffectiveScalarMultivariateFunction::coordinates(4);
    EffectiveScalarMultivariateFunction f(x[0] - x[1] - x[2] - x[0] * x[2] +
                                              x[0] * x[3] + x[1] * x[2] -
                                              x[1] * x[3]);

    ExactBoxType D = ExactBoxType{{0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}};
    EffectiveVectorMultivariateFunction g = {
        8 - x[0] - 2 * x[1], 12 - 4 * x[0] - x[1], 12 - 3 * x[0] - 4 * x[1],
        8 - 2 * x[2] - x[3], 8 - x[2] - 2 * x[3],  5 - x[2] - x[3]};
    ExactBoxType      C         = {{0, +inf}, {0, +inf}, {0, +inf},
                      {0, +inf}, {0, +inf}, {0, +inf}};
    FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);

    std::cout << "f(x_optimal): " << f(x_optimal) << "\n";
    ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, D);
    ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), C);
  }
};

Int main(Int argc, const char *argv[])
{
  Nat optimiser_verbosity = get_verbosity(argc, argv);

  NonlinearSQPOptimiser nlsqp;
  nlsqp.verbosity = optimiser_verbosity;
  TestOptimiser(nlsqp).test();
  // return ARIADNE_TEST_FAILURES;

  NonlinearInteriorPointOptimiser nlo;
  nlo.verbosity = optimiser_verbosity;
  TestOptimiser(nlo).test();
  return ARIADNE_TEST_FAILURES;

  NonlinearMixedOptimiser nlhop;
  nlhop.verbosity = optimiser_verbosity;
  TestOptimiser(nlhop).test();
  return ARIADNE_TEST_FAILURES;
}
