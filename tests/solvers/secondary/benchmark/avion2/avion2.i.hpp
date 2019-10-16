#pragma once

Void TestOptimiser::benchmark_avion2()
{

  List<EffectiveScalarMultivariateFunction> x =
      EffectiveScalarMultivariateFunction::coordinates(49);
  EffectiveScalarMultivariateFunction f(
      pow((x[7] - Real(0.01) * x[2] * x[0]), 2) +
      pow((x[12] - (x[15] - x[14] - x[13] * x[10]) / (pow(x[10], 2))), 2) +
      pow((-2 * x[11] + x[14] + x[15] + Real(0.01) * x[3] / x[10]), 2) +
      pow((x[11] - Real(0.025) * x[14] * pow(x[13], 2) / x[12]), 2) +
      pow((x[16] - Real(27.5) * x[6] - Real(1.3) * pow(x[6], 2)), 2) +
      pow((x[17] - 70 * x[7] + Real(8.6) * pow(x[7], 2)), 2) +
      pow((x[19] - 1000 + pow(x[21], 2) / 24000), 2) +
      pow((1000 * x[20] - x[23] * x[24]), 2) +
      pow((x[31] + x[26] + x[33] / 790 + 2 - x[30] / x[27] + x[29] * x[20]),
          2) +
      pow((x[35] - 1000 * x[20] / (x[28] + 20) - 12 * sqrt(x[20])), 2) +
      pow((x[8] - Real(1.25) * x[0] * x[46]), 2) +
      pow((x[0] - x[23] / x[41]), 2) +
      pow((x[32] - Real(2.4) * x[4] * sqrt(x[4]) * x[42] / sqrt(x[5])), 2) +
      pow((x[14] - Real(0.785) * pow(x[43], 2) * x[20]), 2) +
      pow((x[15] - Real(0.785) * pow(x[44], 2) * x[20]), 2) +
      pow((x[13] - 2 * (x[31] - x[12] * pow(x[10], 3)) /
                       (pow(x[10], 2) * (3 - x[14] * x[10]))),
          2) +
      pow((x[45] - Real(1.15) * x[4] * (15 + Real(0.15) * x[4]) *
                       (8 + pow(sqrt(x[22] * x[5] / (50 * x[0] * x[42])), 3))),
          2));
  ExactBoxType D = ExactBoxType{
      {10, 150},    {0, 10},       {0, 10},       {0, 5},        {7, 120},
      {1.5, 8},     {2, 20},       {2, 30},       {30, 500},     {20, 200},
      {0.01, 20},   {0, 10},       {0.2, -0.001}, {0.1, 2},      {0, 1},
      {0, 2},       {100, 1000},   {500, 5000},   {500, 5000},   {1000, 20000},
      {2, 30},      {2000, 20000}, {3000, 30000}, {5000, 50000}, {0.2, 0.8},
      {1, 5},       {0, 20},       {100, 400},    {4, 15},       {0, 10},
      {500, 10000}, {10, 50},      {250, 5000},   {750, 15000},  {250, 3000},
      {10, 5000},   {35, 70},      {100, 3000},   {200, 400},    {120, 240},
      {700, 1900},  {100, 1000},   {2, 20},       {0, 1},        {0, 2},
      {500, 5000},  {1, 2},        {1, 2},        {1, 2},
  };
  EffectiveVectorMultivariateFunction g = {
      x[6] - Real(0.13) * x[0],
      x[4] - Real(0.7) * x[0],
      x[5] - x[1],
      x[9] - x[8] - 2 * x[6] - 2 * x[4] - 2 * x[7],
      x[18] - 20 * x[9],
      x[23] - 2 * x[21],
      x[33] - x[19] - x[32],
      x[34] - Real(0.137) * x[21],
      x[36] - 35 * x[46],
      x[37] - Real(0.043) * x[19],
      x[38] - 200 * x[47],
      x[39] - 120 * x[48],
      x[40] - 300 * x[25] - 400,
      x[22] - x[21] + 95 * x[47] + 70 * x[48] + 660 * x[46] +
          Real(0.5) * x[19] - 380,
      x[30] - x[34] + x[36] + x[37] + x[38] + x[39] + x[40] + 290};
  ExactBoxType C =
      ExactBoxType{{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
                   {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};

  optimiser->initial_guess = {
      2.7452e+01, 1.5000e+00, 1.0000e+01,  0.0000e+00, 1.9217e+01,
      1.5000e+00, 3.5688e+00, 4.0696e+00,  3.4315e+01, 8.8025e+01,
      5.1306e+00, 0.0000e+00, -1.4809e-01, 7.5980e-01, 0.0000e+00,
      0.0000e+00, 1.1470e+02, 5.0000e+02,  1.7605e+03, 2.3256e+03,
      5.6788e+00, 1.4197e+04, 1.2589e+04,  2.8394e+04, 2.0000e-01,
      1.0000e+00, 0.0000e+00, 1.0000e+02,  1.5000e+01, 0.0000e+00,
      5.0000e+02, 1.0000e+01, 8.1490e+02,  3.1405e+03, 1.9450e+03,
      1.9085e+02, 3.5000e+01, 1.0000e+02,  2.0000e+02, 1.2000e+02,
      7.0000e+02, 1.0000e+03, 4.9367e+00,  0.0000e+00, 0.0000e+00,
      5.0000e+03, 1,          1,           1};
  optimiser->use_initial_guess   = false;
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
