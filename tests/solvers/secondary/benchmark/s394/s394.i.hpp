#pragma once

Void TestOptimiser::benchmark_s394()
{
  List<EffectiveScalarMultivariateFunction> x =
      EffectiveScalarMultivariateFunction::coordinates(20);
  EffectiveScalarMultivariateFunction f(
      1 * (pow(x[0], 2) + pow(x[0], 4)) + 2 * (pow(x[1], 2) + pow(x[1], 4)) +
      3 * (pow(x[2], 2) + pow(x[2], 4)) + 4 * (pow(x[3], 2) + pow(x[3], 4)) +
      5 * (pow(x[4], 2) + pow(x[4], 4)) + 6 * (pow(x[5], 2) + pow(x[5], 4)) +
      7 * (pow(x[6], 2) + pow(x[6], 4)) + 8 * (pow(x[7], 2) + pow(x[7], 4)) +
      9 * (pow(x[8], 2) + pow(x[8], 4)) + 10 * (pow(x[9], 2) + pow(x[9], 4)) +
      11 * (pow(x[10], 2) + pow(x[10], 4)) +
      12 * (pow(x[11], 2) + pow(x[11], 4)) +
      13 * (pow(x[12], 2) + pow(x[12], 4)) +
      14 * (pow(x[13], 2) + pow(x[13], 4)) +
      15 * (pow(x[14], 2) + pow(x[14], 4)) +
      16 * (pow(x[15], 2) + pow(x[15], 4)) +
      17 * (pow(x[16], 2) + pow(x[16], 4)) +
      18 * (pow(x[17], 2) + pow(x[17], 4)) +
      19 * (pow(x[18], 2) + pow(x[18], 4)) +
      20 * (pow(x[19], 2) + pow(x[19], 4)));

  ExactBoxType D = ExactBoxType{
      {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000},
      {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000},
      {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000},
      {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000}, {-1000000, +1000000}};
  EffectiveVectorMultivariateFunction g = {
      pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2) + pow(x[3], 2) + pow(x[4], 2) +
      pow(x[5], 2) + pow(x[6], 2) + pow(x[7], 2) + pow(x[8], 2) + pow(x[9], 2) +
      pow(x[10], 2) + pow(x[11], 2) + pow(x[12], 2) + pow(x[13], 2) +
      pow(x[14], 2) + pow(x[15], 2) + pow(x[16], 2) + pow(x[17], 2) +
      pow(x[18], 2) + pow(x[19], 2)};

  ExactBoxType C = {{1, 1.001}};

  optimiser->initial_guess     = Vector<FloatDP>(x.size(), 2.0);
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
