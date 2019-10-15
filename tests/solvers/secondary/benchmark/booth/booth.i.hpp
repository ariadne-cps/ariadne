#pragma once

Void TestOptimiser::benchmark_booth()
{
  List<EffectiveScalarMultivariateFunction> x =
      EffectiveScalarMultivariateFunction::coordinates(2);
  EffectiveScalarMultivariateFunction f(pow(x[0] + 2 * x[1] - 7, 2) +
                                        pow(2 * x[0] + x[1] - 5, 2));
  ExactBoxType D = ExactBoxType{{-inf, +inf}, {-inf, +inf}};
  EffectiveVectorMultivariateFunction g(0u, 2u);
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