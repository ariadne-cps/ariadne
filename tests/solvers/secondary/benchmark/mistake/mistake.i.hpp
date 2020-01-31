#pragma once

Void TestOptimiser::benchmark_mistake()
{
  List<EffectiveScalarMultivariateFunction> x =
      EffectiveScalarMultivariateFunction::coordinates(9);
  EffectiveScalarMultivariateFunction f(

      -Real(0.5) * x[0] * x[3] + Real(0.5) * x[1] * x[2] -
      Real(0.5) * x[2] * x[8] + Real(0.5) * x[4] * x[8] -
      Real(0.5) * x[4] * x[7] + Real(0.5) * x[5] * x[6]

  );

  ExactBoxType D = ExactBoxType{
      {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000},
      {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000}, {0.0, +1000000},
  };
  EffectiveVectorMultivariateFunction g = {
      x[2] * x[2] + x[3] * x[3] - Real(1.0),
      x[4] * x[4] + x[5] * x[5] - Real(1.0),
      x[8] * x[8] - Real(1.0),
      x[0] * x[0] + (x[1] - x[8]) * (x[1] - x[8]) - Real(1.0),
      (x[0] - x[4]) * (x[0] - x[4]) + (x[1] - x[5]) * (x[1] - x[5]) - Real(1.0),
      (x[0] - x[6]) * (x[0] - x[6]) + (x[1] - x[7]) * (x[1] - x[7]) - Real(1.0),
      (x[2] - x[4]) * (x[2] - x[4]) + (x[3] - x[5]) * (x[3] - x[5]) - Real(1.0),
      (x[2] - x[6]) * (x[2] - x[6]) + (x[3] - x[7]) * (x[3] - x[7]) - Real(1.0),
      x[6] * x[6] + x[7] * x[8] - Real(1.0),
      x[7] * x[8],
      x[4] * x[7] - x[5] * x[6],
      x[0] * x[3] - x[1] * x[2],
      -x[4] * x[8]};

  ExactBoxType C = {
      {-1000000, 0.0}, {-1000000, 0.0}, {-1000000, 0.0}, {-1000000, 0.0}, {-1000000, 0.0},
      {-1000000, 0.0}, {-1000000, 0.0}, {-1000000, 0.0}, {-1000000, 0.0}, {0.0, +1000000},
      {0.0, +1000000}, {0.0, +1000000}, {-1000000, 0.0}
  };

    ExactBoxType testC = {
            {-1000000.01, 0.01}, {-1000000.01, 0.01}, {-1000000.01, 0.01}, {-1000000.01, 0.01}, {-1000000.01, 0.01},
            {-1000000.01, 0.01}, {-1000000.01, 0.01}, {-1000000.01, 0.01}, {-1000000.01, 0.01}, {-0.01, +1000000.01},
            {-0.01, +1000000.01}, {-0.01, +1000000.01}, {-1000000.01, 0.01}
    };

  optimiser->initial_guess     = Vector<FloatDP>(x.size(), 1.0);
  optimiser->use_initial_guess = true;

  float             elapsed_time = 0;
  clock_t           s_time       = clock();
  FloatBoundsVector x_optimal    = optimiser->minimise(f, D, g, C);
  clock_t           e_time       = clock();

  elapsed_time = static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
  std::cout << "Elapsed time: " << elapsed_time << " sec\n";

  std::cout << "f(x_optimal): " << f(x_optimal) << "\n";
  ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, D);
  ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), testC);
}