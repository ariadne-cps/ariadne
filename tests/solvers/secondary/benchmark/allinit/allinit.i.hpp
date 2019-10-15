#pragma once

Void TestOptimiser::benchmark_allinit()
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
  ExactBoxType D = ExactBoxType{{-inf, +inf}, {1, +inf}, {0, 1}, {2, 2.001}};
  EffectiveVectorMultivariateFunction g(0u, 4u);
  ExactBoxType                        C            = {};
  float                               elapsed_time = 0;
  clock_t                             s_time       = clock();
  FloatBoundsVector x_optimal = optimiser->minimise(f, D, g, C);
  clock_t           e_time    = clock();

  elapsed_time = static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
  std::cout << "Elapsed time: " << elapsed_time << " sec\n";

  std::cout << "f(x_optimal): " << f(x_optimal) << "\n";
  ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, D);
  ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), C);
}