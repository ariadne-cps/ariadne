#pragma once

Void TestOptimiser::benchmark_eigminc()
{
  List<EffectiveScalarMultivariateFunction> x =
      EffectiveScalarMultivariateFunction::coordinates(22);
  EffectiveScalarMultivariateFunction f(x[0]);

  ExactBoxType D = ExactBoxType{
      {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0},
      {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0},
      {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0},
      {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0},
      {-1.0, 1.0}, {-1.0, 1.0}};

  ExactBoxType testD = ExactBoxType{
      {-1.0, 1.0001}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0},
      {-1.0, 1.0},    {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0},
      {-1.0, 1.0},    {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0},
      {-1.0, 1.0},    {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0},
      {-1.0, 1.0},    {-1.0, 1.0}};

  EffectiveVectorMultivariateFunction g = {
      x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + x[4] * x[4] + x[5] * x[5] +
          x[6] * x[6] + x[7] * x[7] + x[8] * x[8] + x[9] * x[9] +
          x[10] * x[10] + x[11] * x[11] + x[12] * x[12] + x[13] * x[13] +
          x[14] * x[14] + x[15] * x[15] + x[16] * x[16] + x[17] * x[17] +
          x[18] * x[18] + x[19] * x[19] + x[20] * x[20] + x[21] * x[21] -
          Real(1.0),
      x[1] * x[0] - Real(10.0) * x[1] - x[2],
      x[2] * x[0] - x[1] - Real(9.0) * x[2] - x[3],
      x[3] * x[0] - x[2] - Real(8.0) * x[3] - x[4],
      x[4] * x[0] - x[3] - Real(7.0) * x[4] - x[5],
      x[5] * x[0] - x[4] - Real(6.0) * x[5] - x[6],
      x[6] * x[0] - x[5] - Real(5.0) * x[6] - x[7],
      x[7] * x[0] - x[6] - Real(4.0) * x[7] - x[8],
      x[8] * x[0] - x[7] - Real(3.0) * x[8] - x[9],
      x[9] * x[0] - x[8] - Real(2.0) * x[9] - x[10],
      x[10] * x[0] - x[9] - x[10] - x[11],
      x[11] * x[0] - x[10] - x[12],
      x[12] * x[0] - x[11] + x[12] - x[13],
      x[13] * x[0] - x[12] + Real(2.0) * x[13] - x[14],
      x[14] * x[0] - x[13] + Real(3.0) * x[14] - x[15],
      x[15] * x[0] - x[14] + Real(4.0) * x[15] - x[16],
      x[16] * x[0] - x[15] + Real(5.0) * x[16] - x[17],
      x[17] * x[0] - x[16] + Real(6.0) * x[17] - x[18],
      x[18] * x[0] - x[17] + Real(7.0) * x[18] - x[19],
      x[19] * x[0] - x[18] + Real(8.0) * x[19] - x[20],
      x[20] * x[0] - x[19] + Real(9.0) * x[20] - x[21],
      x[21] * x[0] - x[20] + Real(10.0) * x[21]

  };

  ExactBoxType C = {
      {-1, +1}, {-1, +1}, {-1, +1}, {-1, +1}, {-1, +1}, {-1, +1},
      {-1, +1}, {-1, +1}, {-1, +1}, {-1, +1}, {-1, +1}, {-1, +1},
      {-1, +1}, {-1, +1}, {-1, +1}, {-1, +1}, {-1, +1}, {-1, +1},
      {-1, +1}, {-1, +1}, {-1, +1}, {-1, +1},
  };

  ExactBoxType testC = {
      {-1.01, +1.01}, {-1.01, +1.01}, {-1.01, +1.01},
      {-1.01, +1.01}, {-1.01, +1.01}, {-1.01, +1.01},
      {-1.01, +1.01}, {-1.01, +1.01}, {-1.01, +1.01},
      {-1.01, +1.01}, {-1.01, +1.01}, {-1.01, +1.01},
      {-1.01, +1.01}, {-1.01, +1.01}, {-1.01, +1.01},
      {-1.01, +1.01}, {-1.01, +1.01}, {-1.01, +1.01},
      {-1.01, +1.01}, {-1.01, +1.01}, {-1.01, +1.01},
      {-1.01, +1.01},
  };

  optimiser->initial_guess     = {1.0,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924,
                              0.2182178902359924};
  optimiser->use_initial_guess = true;

  float             elapsed_time = 0;
  clock_t           s_time       = clock();
  FloatBoundsVector x_optimal    = optimiser->minimise(f, D, g, C);
  clock_t           e_time       = clock();

  elapsed_time = static_cast<float>(e_time - s_time) / CLOCKS_PER_SEC;
  std::cout << "Elapsed time: " << elapsed_time << " sec\n";

  std::cout << "f(x_optimal): " << f(x_optimal) << "\n";
  ARIADNE_TEST_BINARY_PREDICATE(element, x_optimal, testD);
  ARIADNE_TEST_BINARY_PREDICATE(element, g(x_optimal), testC);
}