#pragma once

Void TestOptimiser::benchmark_ssnlbeam()
{
  List<EffectiveScalarMultivariateFunction> x =
      EffectiveScalarMultivariateFunction::coordinates(33);
  EffectiveScalarMultivariateFunction f(
      Real(0.05) * x[23] * x[23] + Real(0.05) * x[22] * x[22] +
      Real(17.5) * (cos(x[1])) + Real(17.5) * (cos(x[0])) +
      Real(0.05) * x[24] * x[24] + Real(0.05) * x[23] * x[23] +
      Real(17.5) * (cos(x[2])) + Real(17.5) * (cos(x[1])) +
      Real(0.05) * x[25] * x[25] + Real(0.05) * x[24] * x[24] +
      Real(17.5) * (cos(x[3])) + Real(17.5) * (cos(x[2])) +
      Real(0.05) * x[26] * x[26] + Real(0.05) * x[25] * x[25] +
      Real(17.5) * (cos(x[4])) + Real(17.5) * (cos(x[3])) +
      Real(0.05) * x[27] * x[27] + Real(0.05) * x[26] * x[26] +
      Real(17.5) * (cos(x[5])) + Real(17.5) * (cos(x[4])) +
      Real(0.05) * x[28] * x[28] + Real(0.05) * x[27] * x[27] +
      Real(17.5) * (cos(x[6])) + Real(17.5) * (cos(x[5])) +
      Real(0.05) * x[29] * x[29] + Real(0.05) * x[28] * x[28] +
      Real(17.5) * (cos(x[7])) + Real(17.5) * (cos(x[6])) +
      Real(0.05) * x[30] * x[30] + Real(0.05) * x[29] * x[29] +
      Real(17.5) * (cos(x[8])) + Real(17.5) * (cos(x[7])) +
      Real(0.05) * x[31] * x[31] + Real(0.05) * x[30] * x[30] +
      Real(17.5) * (cos(x[9])) + Real(17.5) * (cos(x[8])) +
      Real(0.05) * x[32] * x[32] + Real(0.05) * x[31] * x[31] +
      Real(17.5) * (cos(x[10])) + Real(17.5) * (cos(x[9])));

  ExactBoxType D = ExactBoxType{
      {-1.0, 1.0},   {-1.0, 1.0},   {-1.0, 1.0},   {-1.0, 1.0},   {-1.0, 1.0},
      {-1.0, 1.0},   {-1.0, 1.0},   {-1.0, 1.0},   {-1.0, 1.0},   {-1.0, 1.0},
      {-1.0, 1.0},   {0.0, 0.0},    {-0.05, 0.05}, {-0.05, 0.05}, {-0.05, 0.05},
      {-0.05, 0.05}, {-0.05, 0.05}, {-0.05, 0.05}, {-0.05, 0.05}, {-0.05, 0.05},
      {-0.05, 0.05}, {0.0, 0.0},    {-inf, +inf},  {-inf, +inf},  {-inf, +inf},
      {-inf, +inf},  {-inf, +inf},  {-inf, +inf},  {-inf, +inf},  {-inf, +inf},
      {-inf, +inf},  {-inf, +inf},  {-inf, +inf}};

  EffectiveVectorMultivariateFunction g = {
      -Real(0.05) * (sin(x[1])) - Real(0.05) * (sin(x[0])) + x[12] - x[11],
      x[1] - x[0] - Real(0.05) * x[23] - Real(0.05) * x[22],
      -Real(0.05) * (sin(x[2])) - Real(0.05) * (sin(x[1])) + x[13] - x[12],
      x[2] - x[1] - Real(0.05) * x[24] - Real(0.05) * x[23],
      -Real(0.05) * (sin(x[3])) - Real(0.05) * (sin(x[2])) + x[14] - x[13],
      x[3] - x[2] - Real(0.05) * x[25] - Real(0.05) * x[24],
      -Real(0.05) * (sin(x[4])) - Real(0.05) * (sin(x[3])) + x[15] - x[14],
      x[4] - x[3] - Real(0.05) * x[26] - Real(0.05) * x[25],
      -Real(0.05) * (sin(x[5])) - Real(0.05) * (sin(x[4])) + x[16] - x[15],
      x[5] - x[4] - Real(0.05) * x[27] - Real(0.05) * x[26],
      -Real(0.05) * (sin(x[6])) - Real(0.05) * (sin(x[5])) + x[17] - x[16],
      x[6] - x[5] - Real(0.05) * x[28] - Real(0.05) * x[27],
      -Real(0.05) * (sin(x[7])) - Real(0.05) * (sin(x[6])) + x[18] - x[17],
      x[7] - x[6] - Real(0.05) * x[29] - Real(0.05) * x[28],
      -Real(0.05) * (sin(x[8])) - Real(0.05) * (sin(x[7])) + x[19] - x[18],
      x[8] - x[7] - Real(0.05) * x[30] - Real(0.05) * x[29],
      -Real(0.05) * (sin(x[9])) - Real(0.05) * (sin(x[8])) + x[20] - x[19],
      x[9] - x[8] - Real(0.05) * x[31] - Real(0.05) * x[30],
      -Real(0.05) * (sin(x[10])) - Real(0.05) * (sin(x[9])) + x[21] - x[20],
      x[10] - x[9] - Real(0.05) * x[32] - Real(0.05) * x[31]

  };

  ExactBoxType C = {
      {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
      {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
      {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
      {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
  };

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