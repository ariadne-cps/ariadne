#pragma once

Void TestOptimiser::benchmark_optprloc()
{
  List<EffectiveScalarMultivariateFunction> x =
      EffectiveScalarMultivariateFunction::coordinates(30);
  EffectiveScalarMultivariateFunction f(
      Real(0.6) * pow((x[0] - Real(0.0)), 2) +
      Real(0.1) * pow((x[3] - Real(0.0)), 2) - x[5] - Real(0.2) * x[6] - x[7] -
      Real(0.2) * x[8] - Real(0.9) * x[9] - Real(0.9) * x[10] -
      Real(0.1) * x[11] - Real(0.8) * x[12] - x[13] - Real(0.4) * x[14] -
      x[15] - Real(0.3) * x[16] - Real(0.1) * x[17] - Real(0.3) * x[18] -
      Real(0.5) * x[19] - Real(0.9) * x[20] - Real(0.8) * x[21] -
      Real(0.1) * x[22] - Real(0.9) * x[23] - x[24] - x[25] - x[26] -
      Real(0.2) * x[27] - Real(0.7) * x[28] - Real(0.7) * x[29] -
      Real(0.9) * x[1] - Real(0.5) * x[2] + x[4]);

  ExactBoxType D = ExactBoxType{
      {2.0, 4.5}, {0.0, 8.0}, {3.0, 9.0}, {0.0, 5.0}, {4.0, 10.0}, {0.0, 1.0},
      {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0},  {0.0, 1.0},
      {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0},  {0.0, 1.0},
      {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0},  {0.0, 1.0},
      {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0},  {0.0, 1.0}};

  ExactBoxType testD = ExactBoxType{
      {2.0, 4.5},   {0.0, 8.0},   {3.0, 9.0},   {-0.000001, 5.0}, {4.0, 10.0},
      {0.0, 1.001}, {0.0, 1.001}, {0.0, 1.001}, {0.0, 1.001},     {0.0, 1.001},
      {0.0, 1.001}, {0.0, 1.001}, {0.0, 1.001}, {0.0, 1.001},     {0.0, 1.001},
      {0.0, 1.001}, {0.0, 1.001}, {0.0, 1.001}, {0.0, 1.001},     {0.0, 1.001},
      {0.0, 1.001}, {0.0, 1.001}, {0.0, 1.001}, {0.0, 1.001},     {0.0, 1.001},
      {0.0, 1.001}, {0.0, 1.001}, {0.0, 1.001}, {0.0, 1.001},     {0.0, 1.001}};

  EffectiveVectorMultivariateFunction g = {
      Real(9.57) * pow((x[0] - Real(2.26)), 2) +
          Real(2.74) * pow((x[1] - Real(5.15)), 2) +
          Real(9.75) * pow((x[2] - Real(4.03)), 2) +
          Real(3.96) * pow((x[3] - Real(1.74)), 2) +
          Real(8.67) * pow((x[4] - Real(4.74)), 2) + Real(1000.0) * x[5] -
          Real(1077.83985),
      Real(8.38) * pow((x[0] - Real(5.51)), 2) +
          Real(3.93) * pow((x[1] - Real(9.01)), 2) +
          Real(5.18) * pow((x[2] - Real(3.84)), 2) +
          Real(5.2) * pow((x[3] - Real(1.47)), 2) +
          Real(7.82) * pow((x[4] - Real(9.92)), 2) + Real(1000.0) * x[6] -
          Real(1175.971),
      Real(9.81) * pow((x[0] - Real(4.06)), 2) +
          Real(0.04) * pow((x[1] - Real(1.8)), 2) +
          Real(4.21) * pow((x[2] - Real(0.71)), 2) +
          Real(7.38) * pow((x[3] - Real(9.09)), 2) +
          Real(4.11) * pow((x[4] - Real(8.13)), 2) + Real(1000.0) * x[7] -
          Real(1201.8226),
      Real(7.41) * pow((x[0] - Real(6.3)), 2) +
          Real(6.08) * pow((x[1] - Real(0.11)), 2) +
          Real(5.46) * pow((x[2] - Real(4.08)), 2) +
          Real(4.86) * pow((x[3] - Real(7.29)), 2) +
          Real(1.48) * pow((x[4] - Real(4.24)), 2) + Real(1000.0) * x[8] -
          Real(1143.9533000000001),
      Real(9.96) * pow((x[0] - Real(2.81)), 2) +
          Real(9.13) * pow((x[1] - Real(1.65)), 2) +
          Real(2.95) * pow((x[2] - Real(8.08)), 2) +
          Real(8.25) * pow((x[3] - Real(3.99)), 2) +
          Real(3.58) * pow((x[4] - Real(3.51)), 2) + Real(1000.0) * x[9] -
          Real(1154.3895),
      Real(9.39) * pow((x[0] - Real(4.29)), 2) +
          Real(4.27) * pow((x[1] - Real(9.49)), 2) +
          Real(5.09) * pow((x[2] - Real(2.24)), 2) +
          Real(1.81) * pow((x[3] - Real(9.78)), 2) +
          Real(7.58) * pow((x[4] - Real(1.52)), 2) + Real(1000.0) * x[10] -
          Real(1433.3177),
      Real(1.88) * pow((x[0] - Real(9.76)), 2) +
          Real(7.2) * pow((x[1] - Real(3.64)), 2) +
          Real(6.65) * pow((x[2] - Real(6.62)), 2) +
          Real(1.74) * pow((x[3] - Real(3.66)), 2) +
          Real(2.86) * pow((x[4] - Real(9.08)), 2) + Real(1000.0) * x[11] -
          Real(1109.0764),
      Real(4.01) * pow((x[0] - Real(1.37)), 2) +
          Real(2.67) * pow((x[1] - Real(6.99)), 2) +
          Real(4.86) * pow((x[2] - Real(7.19)), 2) +
          Real(2.55) * pow((x[3] - Real(3.03)), 2) +
          Real(6.91) * pow((x[4] - Real(3.39)), 2) + Real(1000.0) * x[12] -
          Real(1041.59592),
      Real(4.18) * pow((x[0] - Real(8.89)), 2) +
          Real(1.92) * pow((x[1] - Real(8.29)), 2) +
          Real(2.6) * pow((x[2] - Real(6.05)), 2) +
          Real(7.15) * pow((x[3] - Real(7.48)), 2) +
          Real(2.86) * pow((x[4] - Real(4.09)), 2) + Real(1000.0) * x[13] -
          Real(1144.0623),
      Real(7.81) * pow((x[0] - Real(7.42)), 2) +
          Real(2.14) * pow((x[1] - Real(4.6)), 2) +
          Real(9.63) * pow((x[2] - Real(0.3)), 2) +
          Real(7.61) * pow((x[3] - Real(0.97)), 2) +
          Real(9.17) * pow((x[4] - Real(8.77)), 2) + Real(1000.0) * x[14] -
          Real(1099.8341599999999),
      Real(8.96) * pow((x[0] - Real(1.54)), 2) +
          Real(3.47) * pow((x[1] - Real(7.06)), 2) +
          Real(5.49) * pow((x[2] - Real(0.01)), 2) +
          Real(4.73) * pow((x[3] - Real(1.23)), 2) +
          Real(9.43) * pow((x[4] - Real(3.11)), 2) + Real(1000.0) * x[15] -
          Real(1149.1791),
      Real(9.94) * pow((x[0] - Real(7.74)), 2) +
          Real(1.63) * pow((x[1] - Real(4.4)), 2) +
          Real(1.23) * pow((x[2] - Real(7.93)), 2) +
          Real(4.33) * pow((x[3] - Real(5.95)), 2) +
          Real(7.08) * pow((x[4] - Real(4.88)), 2) + Real(1000.0) * x[16] -
          Real(1123.8074),
      Real(0.31) * pow((x[0] - Real(9.94)), 2) +
          Real(5.0) * pow((x[1] - Real(5.21)), 2) +
          Real(0.16) * pow((x[2] - Real(8.58)), 2) +
          Real(2.52) * pow((x[3] - Real(0.13)), 2) +
          Real(3.08) * pow((x[4] - Real(4.57)), 2) + Real(1000.0) * x[17] -
          Real(1027.22197),
      Real(6.02) * pow((x[0] - Real(9.54)), 2) +
          Real(0.92) * pow((x[1] - Real(1.57)), 2) +
          Real(7.47) * pow((x[2] - Real(9.66)), 2) +
          Real(9.74) * pow((x[3] - Real(5.24)), 2) +
          Real(1.76) * pow((x[4] - Real(7.9)), 2) + Real(1000.0) * x[18] -
          Real(1089.9268299999999),
      Real(5.06) * pow((x[0] - Real(7.46)), 2) +
          Real(4.52) * pow((x[1] - Real(8.81)), 2) +
          Real(1.89) * pow((x[2] - Real(1.67)), 2) +
          Real(1.22) * pow((x[3] - Real(6.47)), 2) +
          Real(9.05) * pow((x[4] - Real(1.81)), 2) + Real(1000.0) * x[19] -
          Real(1293.0765999999999),
      Real(5.92) * pow((x[0] - Real(0.56)), 2) +
          Real(2.56) * pow((x[1] - Real(8.1)), 2) +
          Real(7.74) * pow((x[2] - Real(0.19)), 2) +
          Real(6.96) * pow((x[3] - Real(6.11)), 2) +
          Real(5.18) * pow((x[4] - Real(6.4)), 2) + Real(1000.0) * x[20] -
          Real(1174.317),
      Real(6.45) * pow((x[0] - Real(3.86)), 2) +
          Real(1.52) * pow((x[1] - Real(6.68)), 2) +
          Real(0.06) * pow((x[2] - Real(6.42)), 2) +
          Real(5.34) * pow((x[3] - Real(7.29)), 2) +
          Real(8.47) * pow((x[4] - Real(4.66)), 2) + Real(1000.0) * x[21] -
          Real(1125.1028000000001),
      Real(1.04) * pow((x[0] - Real(2.98)), 2) +
          Real(1.36) * pow((x[1] - Real(2.98)), 2) +
          Real(5.99) * pow((x[2] - Real(3.03)), 2) +
          Real(8.1) * pow((x[3] - Real(0.02)), 2) +
          Real(5.22) * pow((x[4] - Real(0.67)), 2) + Real(1000.0) * x[22] -
          Real(1222.8417),
      Real(1.4) * pow((x[0] - Real(3.61)), 2) +
          Real(1.35) * pow((x[1] - Real(7.62)), 2) +
          Real(0.59) * pow((x[2] - Real(1.79)), 2) +
          Real(8.58) * pow((x[3] - Real(7.8)), 2) +
          Real(1.21) * pow((x[4] - Real(9.81)), 2) + Real(1000.0) * x[23] -
          Real(1050.48593),
      Real(6.68) * pow((x[0] - Real(5.68)), 2) +
          Real(9.48) * pow((x[1] - Real(4.24)), 2) +
          Real(1.6) * pow((x[2] - Real(4.17)), 2) +
          Real(6.74) * pow((x[3] - Real(6.75)), 2) +
          Real(8.92) * pow((x[4] - Real(1.08)), 2) + Real(1000.0) * x[24] -
          Real(1361.1973),
      Real(1.95) * pow((x[0] - Real(5.48)), 2) +
          Real(0.46) * pow((x[1] - Real(3.74)), 2) +
          Real(2.9) * pow((x[2] - Real(3.34)), 2) +
          Real(1.79) * pow((x[3] - Real(6.22)), 2) +
          Real(0.99) * pow((x[4] - Real(7.94)), 2) + Real(1000.0) * x[25] -
          Real(1040.32642),
      Real(5.18) * pow((x[0] - Real(8.13)), 2) +
          Real(5.1) * pow((x[1] - Real(8.72)), 2) +
          Real(8.81) * pow((x[2] - Real(3.93)), 2) +
          Real(3.27) * pow((x[3] - Real(8.8)), 2) +
          Real(9.63) * pow((x[4] - Real(8.56)), 2) + Real(1000.0) * x[26] -
          Real(1161.8518),
      Real(1.47) * pow((x[0] - Real(1.37)), 2) +
          Real(5.71) * pow((x[1] - Real(0.54)), 2) +
          Real(6.95) * pow((x[2] - Real(1.55)), 2) +
          Real(1.42) * pow((x[3] - Real(5.56)), 2) +
          Real(3.49) * pow((x[4] - Real(5.85)), 2) + Real(1000.0) * x[27] -
          Real(1066.85827),
      Real(5.4) * pow((x[0] - Real(8.79)), 2) +
          Real(3.12) * pow((x[1] - Real(5.04)), 2) +
          Real(5.37) * pow((x[2] - Real(4.83)), 2) +
          Real(6.1) * pow((x[3] - Real(6.94)), 2) +
          Real(3.71) * pow((x[4] - Real(0.38)), 2) + Real(1000.0) * x[28] -
          Real(1340.5807),
      Real(6.32) * pow((x[0] - Real(2.66)), 2) +
          Real(0.81) * pow((x[1] - Real(4.19)), 2) +
          Real(6.12) * pow((x[2] - Real(6.49)), 2) +
          Real(6.73) * pow((x[3] - Real(8.04)), 2) +
          Real(7.93) * pow((x[4] - Real(1.66)), 2) + Real(1000.0) * x[29] -
          Real(1407.52),
      x[0] - x[1] + x[2] + x[3] + x[4] - Real(10.0),
      Real(0.6) * x[0] - Real(0.9) * x[1] - Real(0.5) * x[2] +
          Real(0.1) * x[3] + x[4] + Real(0.64),
      x[0] - x[1] + x[2] - x[3] + x[4] - Real(0.69),
      Real(0.157) * x[0] + Real(0.05) * x[1] - Real(1.5),
      Real(0.25) * x[1] + Real(1.05) * x[3] - Real(0.3) * x[4] - Real(4.5)};

  ExactBoxType C = {
      {0.0, +inf}, {0.0, +inf}, {0.0, +inf}, {0.0, +inf}, {0.0, +inf},
      {0.0, +inf}, {0.0, +inf}, {0.0, +inf}, {0.0, +inf}, {0.0, +inf},
      {0.0, +inf}, {0.0, +inf}, {0.0, +inf}, {0.0, +inf}, {0.0, +inf},
      {0.0, +inf}, {0.0, +inf}, {0.0, +inf}, {0.0, +inf}, {0.0, +inf},
      {0.0, +inf}, {0.0, +inf}, {0.0, +inf}, {0.0, +inf}, {0.0, +inf},
      {0.0, +inf}, {-inf, 0.0}, {0.0, +inf}, {-inf, 0.0},
  };

  ExactBoxType testC = {
      {-0.0001, +inf}, {-0.0001, +inf}, {-0.0001, +inf}, {-0.0001, +inf},
      {-0.0001, +inf}, {-0.0001, +inf}, {-0.0001, +inf}, {-0.0001, +inf},
      {-0.0001, +inf}, {-0.0001, +inf}, {-0.0001, +inf}, {-0.0001, +inf},
      {-0.0001, +inf}, {-0.0001, +inf}, {-0.0001, +inf}, {-0.0001, +inf},
      {-0.0001, +inf}, {-0.0001, +inf}, {-0.0001, +inf}, {-0.0001, +inf},
      {-0.0001, +inf}, {-0.0001, +inf}, {-0.0001, +inf}, {-0.0001, +inf},
      {-0.0001, +inf}, {-0.0001, +inf}, {-inf, 0.0001},  {-0.0001, +inf},
      {-inf, 0.0001},
  };

  optimiser->initial_guess     = Vector<FloatDP>(x.size(), 0.0);
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