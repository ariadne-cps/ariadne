#pragma once

Void TestOptimiser::benchmark_hs54()
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
  ExactBoxType                        C = {{-1, 1}};
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