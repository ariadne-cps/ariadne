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

#include "../../test.hpp"

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
  // DoublePrecision                     pr;
  /*
   * -> near boundary feasible region exception
   ** -> x not contained in D
   */
public:
  TestOptimiser(const OptimiserInterface& opt)
    : optimiser(opt.clone())
  {}
  Void test()
  {
    ARIADNE_TEST_CALL(benchmark_allinit());
    // ARIADNE_TEST_CALL(benchmark_extrasim());
    // ARIADNE_TEST_CALL(benchmark_hs54()); //**
    // ARIADNE_TEST_CALL(benchmark_booth());
    // ARIADNE_TEST_CALL(benchmark_himmelbc());
    // ARIADNE_TEST_CALL(benchmark_hs44());
    // ARIADNE_TEST_CALL(benchmark_s394());
    // ARIADNE_TEST_CALL(benchmark_dualc2());
    // ARIADNE_TEST_CALL(benchmark_loadbal());  //*
    // ARIADNE_TEST_CALL(benchmark_mistake());  //*
    // ARIADNE_TEST_CALL(benchmark_ssnlbeam()); //*
    // ARIADNE_TEST_CALL(benchmark_optprloc()); //**
    // ARIADNE_TEST_CALL(benchmark_eigminc());  //*
    // ARIADNE_TEST_CALL(benchmark_smmpsf());   //*
    // ARIADNE_TEST_CALL(benchmark_reading3()); //*
    // ARIADNE_TEST_CALL(benchmark_dittert());  //*
    // ARIADNE_TEST_CALL(benchmark_avion2());
    // ARIADNE_TEST_CALL(benchmark_degenlpa());
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
  Void benchmark_allinit();

  /*
  Results from AmplOnline with LOQO:

  LOQO 7.03: optimal solution (11 iterations, 11 evaluations)
  primal objective 1.000000008
  dual objective 1
  f = 1

  x = 7.97704e-09
  y = 1
  */
  Void benchmark_extrasim();

  /*
    Precision test-case:
      smaller rtol is in SQP framework, better is the result. Asyntotically it
    is correct.
    With rtol = 1.49012e-10 the result  is
    [-0.00000000000000:0.00000000000364] which is close to 0 (error of order
    2^-12)

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
  Void benchmark_hs54();

  /*
    Empty bound on variables test:
      If there aren't bounds on variables this will generate an empty vector
    which generate an error on line nonlinear_programming.cpp:3442. Solution: if
    vector is empty generate an empty vector, or implement on Vector class level
    the operation with empty sets.

    LOQO 7.03: optimal solution (2 QP iterations, 3 evaluations)
    primal objective 7.222223417e-21
      dual objective 7.222223417e-21
    f = 7.22222e-21

    x [*] :=
    1  1
    2  3
    ;
  */
  Void benchmark_booth();

  /*
    If there aren't bounds on variables this will generate an empty vector which
    generate an error on line nonlinear_programming.cpp:3442.
    Solution: if vector is empty generate an empty vector, or implement on
    Vector class level the operation with empty sets.

    LOQO 7.03: optimal solution (10 iterations, 11 evaluations)
    primal objective 1.009641594e-15
      dual objective 1.009641594e-15
    f = 1.00964e-15

    x [*] :=
    1  3
    2  2
    ;
  */
  Void benchmark_himmelbc();

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
  Void benchmark_hs44();

  /*
    Hard problem where without correct starting point (2.0) algorithm keeps
    going! LOQO 7.03: iteration limit (500 iterations, 5301 evaluations) primal
    objective 11.12319128 dual objective 11.10868168 x [*] := 1  -1.69479 2
    -0.0170978 3  -0.00454764 4  -0.002012 5  -0.00110918 6  -0.00069215 7
    -0.000467859 8  -0.000334474 9  -0.000249263 10  -0.000191824 11 -0.00015145
      12  -0.000122102
      13  -0.000100173
      14  -8.34052e-05
      15  -7.03299e-05
      16  -5.99608e-05
      17  -5.16164e-05
      18  -4.48142e-05
      19  -3.92055e-05
      20  -3.45337e-05
      ;

      Using starting point correctly LOQO gives this:
      LOQO 7.03: optimal solution (17 iterations, 18 evaluations)
        primal objective 1.916673223
          dual objective 1.916673223
        x [*] :=
         1   0.912872
         2   0.40825
         3  -2.76502e-08
         4  -3.4376e-08
         5   1.15118e-08
         6   1.91728e-08
         7   2.17833e-08
         8   2.29519e-08
         9   2.35228e-08
        10   2.38022e-08
        11   2.39273e-08
        12   2.3967e-08
        13   2.39582e-08
        14   2.3922e-08
        15   2.38705e-08
        16   2.38109e-08
        17   2.37476e-08
        18   2.36833e-08
        19   2.36197e-08
        20   2.35577e-08
        ;

      SQP give the same f value of LOQO if it starts from < 2, 2, ...., 2 >
  */
  Void benchmark_s394();

  /*
    Constraints fail also in LOQO, MINOS and others...
      LOQO 7.03: optimal solution (21 QP iterations, 21 evaluations)
      primal objective 3551.306378
        dual objective 3551.306356
      x1 = 0.717116

        x2 = 2.5277e-10

        x3 = 1.51067e-11

        x4 = 0.0968997

        x5 = 8.89946e-11

        x6 = 2.29969e-12

        x7 = 0.185984

      obj = 3551.31
  */
  Void benchmark_dualc2();

  /*
    Decrease rtol to have correct solution, also LOQO returns variables that
    does not satisfy constraints (especially equalities) LOQO 7.03: optimal
    solution (23 iterations, 23 evaluations) primal objective 0.4528510407 dual
    objective 0.4528510345 x4_1 = 1.42687e-07
  */
  Void benchmark_loadbal();

  /*
     SQP is not so accurate here, max error of +/-0.001 (really big) if starting
     from given point. Starting from 0, it remains in 0 and it is feasible. (the
     same for LOQO). LOQO 7.03: optimal solution (15 iterations, 15 evaluations)
      primal objective -1.000000017
        dual objective -1.000000036
      x1 = 0.491772
      x2 = 1.03702
      x3 = 0.0144265
      x4 = 0.957078
      x5 = 0.88515
      x6 = 0.465305
      x7 = 0.0264672
      x8 = 1.92217
      x9 = 0.166293
      obj = -1
  */
  Void benchmark_mistake();

  /*
    NB: this need higher accuracy (epsilon=10e-16) and starting point as
    specified The same result of LOQO. LOQO 7.03: optimal solution (60
    iterations, 60 evaluations) primal objective 337.7724788 dual objective
    337.772476 t0 = 0.827038 t1 = 0.26724 t2 = -0.26724 t3 = -0.52151 t4 =
    -0.378766 t5 = 2.71011e-09 t6 = 0.378766 t7 = 0.52151 t8 = 0.26724 t9 =
    -0.26724 t10 = -0.827038 x0 = 0 x1 = 0.05 x2 = 0.05 x3 = 0.011887 x4 =
    -0.0315113 x5 = -0.05 x6 = -0.0315113 x7 = 0.011887 x8 = 0.05 x9 = 0.05 x10
    = 0 u0 = -5.2353 u1 = -5.96067 u2 = -4.72892 u3 = -0.356498 u4 = 3.21138 u5
    = 4.36395 u6 = 3.21138 u7 = -0.356498 u8 = -4.72892 u9 = -5.96067 u10 =
    -5.2353 obj = 337.772
  */
  Void benchmark_ssnlbeam();

  /*
      LOQO 7.03: optimal solution (23 iterations, 23 evaluations)
        primal objective -16.4197738
          dual objective -16.4197739
        x1 = 2
        x2 = 8
        x3 = 7.32852
        x4 = 3.52381
        x5 = 4
        y1 = 0.931507
        y2 = 0.709683
        y3 = 0.675482
        y4 = 0.50181
        y5 = 0.775392
        y6 = 1
        y7 = 0.781821
        y8 = 1
        y9 = 0.829286
        y10 = 0.111667
        y11 = 0.817812
        y12 = 0.743768
        y13 = 0.938481
        y14 = 0.613581
        y15 = 1
        y16 = 0.691209
        y17 = 1
        y18 = 0.919566
        y19 = 0.830827
        y20 = 0.97454
        y21 = 0.933849
        y22 = 0.571531
        y23 = 0.498603
        y24 = 0.910953
        y25 = 1

        obj = -16.4198

        SQP gives wrong results on constraints (really hard to solve!)
  */
  Void benchmark_optprloc();

  /*
    LOQO 7.03: optimal solution (13 iterations, 13 evaluations)
      primal objective 1
        dual objective 0.9999999927
      d = 1
      q1 = -2.46398e-06
      q2 = 2.2176e-05
      q3 = -0.000174944
      q4 = 0.00120243
      q5 = -0.00703963
      q6 = 0.0339957
      q7 = -0.128943
      q8 = 0.352834
      q9 = -0.576725
      q10 = 0.223891
      q11 = 0.576725
      q12 = 0.352834
      q13 = 0.128943
      q14 = 0.0339957
      q15 = 0.00703963
      q16 = 0.00120243
      q17 = 0.000174944
      q18 = 2.21796e-05
      q19 = 2.49235e-06
      q20 = 2.51533e-07
      q21 = 2.28776e-08

      obj = 1
  */
  Void benchmark_eigminc();

  /*
    Cannot run with n > 300 variables on ampl demo version
  */
  Void benchmark_smmpsf();

  /*
    LOQO 7.03: optimal solution (49 iterations, 52 evaluations)
      primal objective -2.243213634e-14
        dual objective -1.1197009e-08
      obj = -2.24321e-14
  */
  Void benchmark_reading3();

  /*
    Cannot run with n > 300 variables on ampl demo version
  */
  Void benchmark_dittert();

  Void benchmark_degenlpa();

  Void benchmark_avion2();
};

#if defined HAVE_EIGEN3_H && defined HAVE_GLPK_H
#include "benchmark/benchmark.hpp"
#endif

Int
main(Int argc, const char* argv[])
{
  Nat optimiser_verbosity = get_verbosity(argc, argv);

  // std::cout << "NonlinearInteriorPointOptimiser\n";
  // NonlinearInteriorPointOptimiser nlo;
  // nlo.verbosity = optimiser_verbosity;
  // TestOptimiser(nlo).test();
  // return ARIADNE_TEST_FAILURES;

  // NonlinearInfeasibleInteriorPointOptimiser nlio;
  // nlio.verbosity = optimiser_verbosity;
  // TestOptimiser(nlio).test();
  // return ARIADNE_TEST_FAILURES;

#if defined HAVE_EIGEN3_H && defined HAVE_GLPK_H

  std::cout << "NonlinearSQPOptimiser\n";
  NonlinearSQPOptimiser nlsqp;
  nlsqp.verbosity = optimiser_verbosity;
  TestOptimiser(nlsqp).test();
  return ARIADNE_TEST_FAILURES;

  std::cout << "NonlinearMixedOptimiser\n";
  NonlinearMixedOptimiser nlhop;
  nlhop.verbosity = optimiser_verbosity;
  TestOptimiser(nlhop).test();
  // return ARIADNE_TEST_FAILURES;
#endif
}
