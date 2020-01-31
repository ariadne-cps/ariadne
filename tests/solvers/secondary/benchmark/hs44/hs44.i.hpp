#pragma once

Void TestOptimiser::benchmark_hs44()
{
  List<EffectiveScalarMultivariateFunction> x =
      EffectiveScalarMultivariateFunction::coordinates(4);
  EffectiveScalarMultivariateFunction f(x[0] - x[1] - x[2] - x[0] * x[2] +
                                        x[0] * x[3] + x[1] * x[2] -
                                        x[1] * x[3]);

  ExactBoxType D = ExactBoxType{{0, +1000000}, {0, +1000000}, {0, +1000000}, {0, +1000000}};
  EffectiveVectorMultivariateFunction g = {
      8 - x[0] - 2 * x[1], 12 - 4 * x[0] - x[1], 12 - 3 * x[0] - 4 * x[1],
      8 - 2 * x[2] - x[3], 8 - x[2] - 2 * x[3],  5 - x[2] - x[3]};
  ExactBoxType      C            = {{0, +1000000}, {0, +1000000}, {0, +1000000},
                    {0, +1000000}, {0, +1000000}, {0, +1000000}};
  float             elapsed_time = 0;
  clock_t           s_time       = clock();
  FloatBoundsVector x_optimal    = optimiser->minimise(f, D, g, C);
  clock_t           e_time       = clock();

  elapsed_time = static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
  std::cout << "Elapsed time: " << elapsed_time << " sec\n";

  std::cout << "f(x_optimal): " << f(x_optimal) << "\n";
  ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, D);
  ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), C);
}