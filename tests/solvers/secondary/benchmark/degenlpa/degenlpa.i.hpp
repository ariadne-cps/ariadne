#pragma once

Void TestOptimiser::benchmark_degenlpa()
{

  List<EffectiveScalarMultivariateFunction> x =
      EffectiveScalarMultivariateFunction::coordinates(20);
  EffectiveScalarMultivariateFunction f(
      Real(0.01) * x[1] + Real(33.333) * x[2] + Real(100.0) * x[3] +
      Real(0.01) * x[4] + Real(33.343) * x[5] + Real(100.01) * x[6] +
      Real(33.333) * x[7] + Real(133.33) * x[8] + Real(100.0) * x[9]);
  ExactBoxType D =
      ExactBoxType{{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0},
                   {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0},
                   {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0},
                   {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}};
  EffectiveVectorMultivariateFunction g = {
      -Real(0.70785) + x[0] + 2 * x[1] + 2 * x[2] + 2 * x[3] + x[4] + 2 * x[5] +
          2 * x[6] + x[7] + 2 * x[8] + x[9],
      Real(0.326) * x[0] - 101 * x[1] + 200 * x[4] + Real(0.06) * x[5] +
          Real(0.02) * x[6],
      Real(0.0066667) * x[0] - Real(1.03) * x[2] + 200 * x[5] +
          Real(0.06) * x[7] + Real(0.02) * x[8],
      Real(0.00066667) * x[0] - Real(1.01) * x[3] + 200 * x[6] +
          Real(0.06) * x[8] + Real(0.02) * x[9],
      Real(0.978) * x[1] - 201 * x[4] + 100 * x[10] + Real(0.03) * x[11] +
          Real(0.01) * x[12],
      Real(0.01) * x[1] + Real(0.489) * x[2] - Real(101.03) * x[5] +
          100 * x[11] + Real(0.03) * x[13] + Real(0.01) * x[14],
      Real(0.001) * x[1] + Real(0.489) * x[3] - Real(101.03) * x[6] +
          100 * x[12] + Real(0.03) * x[14] + Real(0.01) * x[15],
      Real(0.001) * x[2] + Real(0.01) * x[3] - Real(1.04) * x[8] + 100 * x[14] +
          Real(0.03) * x[17] + Real(0.01) * x[18],
      Real(0.02) * x[2] - Real(1.06) * x[7] + 100 * x[13] + Real(0.03) * x[16] +
          Real(0.01) * x[18],
      Real(0.002) * x[3] - Real(1.02) * x[9] + 100 * x[15] +
          Real(0.03) * x[18] + Real(0.01) * x[19],
      -Real(2.5742e-6) * x[10] + Real(0.00252) * x[12] - Real(0.61975) * x[15] +
          Real(1.03) * x[19],
      -Real(0.00257) * x[10] + Real(0.25221) * x[11] - Real(6.2) * x[13] +
          Real(1.09) * x[16],
      Real(0.00629) * x[10] - Real(0.20555) * x[11] - Real(4.1106) * x[12] +
          Real(101.04) * x[14] + Real(505.1) * x[15] - Real(256.72) * x[18],
      Real(0.00841) * x[11] - Real(0.08406) * x[12] - Real(0.20667) * x[13] +
          Real(20.658) * x[15] + Real(1.07) * x[17] - Real(10.5) * x[18]};
  ExactBoxType C =
      ExactBoxType{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};

  optimiser->initial_guess     = Vector<FloatDP>(20u, 1.0);
  optimiser->use_initial_guess = true;

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