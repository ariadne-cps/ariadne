#pragma once

Void TestOptimiser::benchmark_dualc2()
{

  List<EffectiveScalarMultivariateFunction> x =
      EffectiveScalarMultivariateFunction::coordinates(7);

  EffectiveScalarMultivariateFunction f(
      Real(12465.0) * Real(0.5) * x[0] * x[0] + Real(14562.0) * x[0] * x[1] +
      Real(41412.0) * x[0] * x[2] - Real(3984.0) * x[0] * x[3] +
      Real(8262.0) * x[0] * x[4] + Real(62178.0) * x[0] * x[5] -
      Real(20577.0) * x[0] * x[6] + Real(17708.0) * Real(0.5) * x[1] * x[1] +
      Real(44808.0) * x[1] * x[2] - Real(1278.0) * x[1] * x[3] +
      Real(10878.0) * x[1] * x[4] + Real(74380.0) * x[1] * x[5] -
      Real(22456.0) * x[1] * x[6] + Real(164408.0) * Real(0.5) * x[2] * x[2] -
      Real(41278.0) * x[2] * x[3] + Real(4578.0) * x[2] * x[4] +
      Real(158680.0) * x[2] * x[5] - Real(70106.0) * x[2] * x[6] +
      Real(31161.0) * Real(0.5) * x[3] * x[3] + Real(24199.0) * x[3] * x[4] +
      Real(37662.0) * x[3] * x[5] + Real(6221.0) * x[3] * x[6] +
      Real(39937.0) * Real(0.5) * x[4] * x[4] + Real(120170.0) * x[4] * x[5] -
      Real(23267.0) * x[4] * x[6] + Real(492812.0) * Real(0.5) * x[5] * x[5] -
      Real(127852.0) * x[5] * x[6] + Real(42338.0) * Real(0.5) * x[6] * x[6] +
      Real(2.7078336115) * x[1] + Real(41008.309175) * x[2] +
      Real(3406.2993615) * x[3] + Real(11325.133357) * x[4] +
      Real(239354.78723) * x[5] + Real(11004.893252) * x[6]);

  ExactBoxType D = ExactBoxType{{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0},
                                {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}};

  EffectiveVectorMultivariateFunction g = {
      x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6],
      Decimal(282.0) * x[0] + Decimal(282.0) * x[1] + Decimal(282.0) * x[2] +
          Decimal(282.0) * x[3] + Decimal(282.0) * x[4] +
          Decimal(353.0) * x[5] + Decimal(330.0) * x[6],
      Decimal(401.0) * x[0] + Decimal(401.0) * x[1] + Decimal(385.0) * x[2] +
          Decimal(385.0) * x[3] + Decimal(385.0) * x[4] +
          Decimal(383.0) * x[5] + Decimal(362.0) * x[6],
      Decimal(743.0) * x[0] + Decimal(743.0) * x[1] + Decimal(743.0) * x[2] +
          Decimal(743.0) * x[3] + Decimal(746.0) * x[4] +
          Decimal(743.0) * x[5] + Decimal(797.0) * x[6],
      Decimal(746.0) * x[0] + Decimal(746.0) * x[1] + Decimal(746.0) * x[2] +
          Decimal(746.0) * x[3] + Decimal(746.0) * x[4] +
          Decimal(754.0) * x[5] + Decimal(763.0) * x[6],
      Decimal(1449.0) * x[0] + Decimal(1449.0) * x[1] + Decimal(1445.0) * x[2] +
          Decimal(1459.0) * x[3] + Decimal(1449.0) * x[4] +
          Decimal(1441.0) * x[5] + Decimal(1451.0) * x[6],
      Decimal(411.0) * x[0] + Decimal(411.0) * x[1] + Decimal(411.0) * x[2] +
          Decimal(411.0) * x[3] + Decimal(410.0) * x[4] +
          Decimal(407.0) * x[5] + Decimal(397.0) * x[6],
      Decimal(96.0) * x[0] + Decimal(96.0) * x[1] + Decimal(116.0) * x[2] +
          Decimal(112.0) * x[3] + Decimal(112.0) * x[4] +
          Decimal(113.0) * x[5] + Decimal(111.0) * x[6],
      Decimal(834.0) * x[0] + Decimal(834.0) * x[1] + Decimal(834.0) * x[2] +
          Decimal(834.0) * x[3] + Decimal(834.0) * x[4] +
          Decimal(818.0) * x[5] + Decimal(776.0) * x[6],
      Decimal(976.0) * x[0] + Decimal(976.0) * x[1] + Decimal(976.0) * x[2] +
          Decimal(976.0) * x[3] + Decimal(966.0) * x[4] +
          Decimal(918.0) * x[5] + Decimal(994.0) * x[6],
      Decimal(1161.0) * x[0] + Decimal(1161.0) * x[1] + Decimal(1161.0) * x[2] +
          Decimal(1161.0) * x[3] + Decimal(1169.0) * x[4] +
          Decimal(1166.0) * x[5] + Decimal(1184.0) * x[6],
      Decimal(1461.0) * x[0] + Decimal(1461.0) * x[1] + Decimal(1461.0) * x[2] +
          Decimal(1425.0) * x[3] + Decimal(1461.0) * x[4] +
          Decimal(1462.0) * x[5] + Decimal(1377.0) * x[6],
      Decimal(1541.0) * x[0] + Decimal(1541.0) * x[1] + Decimal(1490.0) * x[2] +
          Decimal(1541.0) * x[3] + Decimal(1483.0) * x[4] +
          Decimal(1498.0) * x[5] + Decimal(1521.0) * x[6],
      Decimal(576.0) * x[0] + Decimal(576.0) * x[1] + Decimal(576.0) * x[2] +
          Decimal(576.0) * x[3] + Decimal(576.0) * x[4] +
          Decimal(606.0) * x[5] + Decimal(621.0) * x[6],
      Decimal(2048.0) * x[0] + Decimal(2048.0) * x[1] + Decimal(2048.0) * x[2] +
          Decimal(2048.0) * x[3] + Decimal(2048.0) * x[4] +
          Decimal(2048.0) * x[5] + Decimal(2090.0) * x[6],
      Decimal(1911.0) * x[0] + Decimal(1911.0) * x[1] + Decimal(1906.0) * x[2] +
          Decimal(1911.0) * x[3] + Decimal(1911.0) * x[4] +
          Decimal(1928.0) * x[5] + Decimal(1926.0) * x[6],
      Decimal(1445.0) * x[0] + Decimal(1445.0) * x[1] + Decimal(1449.0) * x[2] +
          Decimal(1481.0) * x[3] + Decimal(1445.0) * x[4] +
          Decimal(1418.0) * x[5] + Decimal(1514.0) * x[6],
      Decimal(1347.0) * x[0] + Decimal(1347.0) * x[1] + Decimal(1358.0) * x[2] +
          Decimal(1347.0) * x[3] + Decimal(1355.0) * x[4] +
          Decimal(1341.0) * x[5] + Decimal(1331.0) * x[6],
      Decimal(1716.0) * x[0] + Decimal(1716.0) * x[1] + Decimal(1717.0) * x[2] +
          Decimal(1716.0) * x[3] + Decimal(1716.0) * x[4] +
          Decimal(1717.0) * x[5] + Decimal(1681.0) * x[6],
      Decimal(988.0) * x[0] + Decimal(988.0) * x[1] + Decimal(988.0) * x[2] +
          Decimal(987.0) * x[3] + Decimal(970.0) * x[4] +
          Decimal(970.0) * x[5] + Decimal(1006.0) * x[6],
      Decimal(342.0) * x[0] + Decimal(342.0) * x[1] + Decimal(342.0) * x[2] +
          Decimal(342.0) * x[3] + Decimal(342.0) * x[4] +
          Decimal(344.0) * x[5] + Decimal(333.0) * x[6],
      Decimal(748.0) * x[0] + Decimal(748.0) * x[1] + Decimal(744.0) * x[2] +
          Decimal(748.0) * x[3] + Decimal(740.0) * x[4] +
          Decimal(738.0) * x[5] + Decimal(744.0) * x[6],
      Decimal(826.0) * x[0] + Decimal(826.0) * x[1] + Decimal(826.0) * x[2] +
          Decimal(826.0) * x[3] + Decimal(826.0) * x[4] +
          Decimal(782.0) * x[5] + Decimal(862.0) * x[6],
      Decimal(586.0) * x[0] + Decimal(586.0) * x[1] + Decimal(586.0) * x[2] +
          Decimal(586.0) * x[3] + Decimal(615.0) * x[4] +
          Decimal(615.0) * x[5] + Decimal(564.0) * x[6],
      Decimal(1898.0) * x[0] + Decimal(1898.0) * x[1] + Decimal(1898.0) * x[2] +
          Decimal(1898.0) * x[3] + Decimal(1898.0) * x[4] +
          Decimal(1899.0) * x[5] + Decimal(1896.0) * x[6],
      Decimal(1560.0) * x[0] + Decimal(1560.0) * x[1] + Decimal(1560.0) * x[2] +
          Decimal(1560.0) * x[3] + Decimal(1559.0) * x[4] +
          Decimal(1559.0) * x[5] + Decimal(1567.0) * x[6],
      Decimal(160.0) * x[0] + Decimal(160.0) * x[1] + Decimal(160.0) * x[2] +
          Decimal(154.0) * x[3] + Decimal(154.0) * x[4] +
          Decimal(111.0) * x[5] + Decimal(159.0) * x[6],
      Decimal(1138.0) * x[0] + Decimal(1138.0) * x[1] + Decimal(1138.0) * x[2] +
          Decimal(1134.0) * x[3] + Decimal(1138.0) * x[4] +
          Decimal(1097.0) * x[5] + Decimal(1133.0) * x[6],
      Decimal(1584.0) * x[0] + Decimal(1584.0) * x[1] + Decimal(1606.0) * x[2] +
          Decimal(1584.0) * x[3] + Decimal(1584.0) * x[4] +
          Decimal(1627.0) * x[5] + Decimal(1623.0) * x[6],
      Decimal(911.0) * x[0] + Decimal(911.0) * x[1] + Decimal(911.0) * x[2] +
          Decimal(911.0) * x[3] + Decimal(911.0) * x[4] +
          Decimal(905.0) * x[5] + Decimal(915.0) * x[6],
      Decimal(1386.0) * x[0] + Decimal(1386.0) * x[1] + Decimal(1386.0) * x[2] +
          Decimal(1386.0) * x[3] + Decimal(1386.0) * x[4] +
          Decimal(1386.0) * x[5] + Decimal(1432.0) * x[6],
      Decimal(1747.0) * x[0] + Decimal(1747.0) * x[1] + Decimal(1759.0) * x[2] +
          Decimal(1753.0) * x[3] + Decimal(1753.0) * x[4] +
          Decimal(1746.0) * x[5] + Decimal(1769.0) * x[6],
      Decimal(1686.0) * x[0] + Decimal(1686.0) * x[1] + Decimal(1686.0) * x[2] +
          Decimal(1686.0) * x[3] + Decimal(1686.0) * x[4] +
          Decimal(1686.0) * x[5] + Decimal(1690.0) * x[6],
      Decimal(136.0) * x[0] + Decimal(140.0) * x[1] + Decimal(140.0) * x[2] +
          Decimal(140.0) * x[3] + Decimal(136.0) * x[4] + Decimal(55.0) * x[5] +
          Decimal(174.0) * x[6],
      Decimal(257.0) * x[0] + Decimal(258.0) * x[1] + Decimal(262.0) * x[2] +
          Decimal(258.0) * x[3] + Decimal(257.0) * x[4] +
          Decimal(261.0) * x[5] + Decimal(245.0) * x[6],
      Decimal(1887.0) * x[0] + Decimal(1887.0) * x[1] + Decimal(1887.0) * x[2] +
          Decimal(1830.0) * x[3] + Decimal(1860.0) * x[4] +
          Decimal(1916.0) * x[5] + Decimal(1854.0) * x[6],
      Decimal(1797.0) * x[0] + Decimal(1797.0) * x[1] + Decimal(1797.0) * x[2] +
          Decimal(1797.0) * x[3] + Decimal(1797.0) * x[4] +
          Decimal(1797.0) * x[5] + Decimal(1802.0) * x[6],
      Decimal(404.0) * x[0] + Decimal(404.0) * x[1] + Decimal(416.0) * x[2] +
          Decimal(377.0) * x[3] + Decimal(367.0) * x[4] +
          Decimal(305.0) * x[5] + Decimal(415.0) * x[6],
      Decimal(1501.0) * x[0] + Decimal(1501.0) * x[1] + Decimal(1501.0) * x[2] +
          Decimal(1501.0) * x[3] + Decimal(1501.0) * x[4] +
          Decimal(1514.0) * x[5] + Decimal(1529.0) * x[6],
      Decimal(1732.0) * x[0] + Decimal(1732.0) * x[1] + Decimal(1697.0) * x[2] +
          Decimal(1698.0) * x[3] + Decimal(1698.0) * x[4] +
          Decimal(1677.0) * x[5] + Decimal(1749.0) * x[6],
      Decimal(222.0) * x[0] + Decimal(222.0) * x[1] + Decimal(222.0) * x[2] +
          Decimal(221.0) * x[3] + Decimal(221.0) * x[4] +
          Decimal(304.0) * x[5] + Decimal(216.0) * x[6],
      Decimal(368.0) * x[0] + Decimal(369.0) * x[1] + Decimal(367.0) * x[2] +
          Decimal(369.0) * x[3] + Decimal(368.0) * x[4] +
          Decimal(384.0) * x[5] + Decimal(325.0) * x[6],
      Decimal(1483.0) * x[0] + Decimal(1483.0) * x[1] + Decimal(1478.0) * x[2] +
          Decimal(1490.0) * x[3] + Decimal(1468.0) * x[4] +
          Decimal(1458.0) * x[5] + Decimal(1483.0) * x[6],
      Decimal(824.0) * x[0] + Decimal(824.0) * x[1] + Decimal(824.0) * x[2] +
          Decimal(824.0) * x[3] + Decimal(847.0) * x[4] +
          Decimal(828.0) * x[5] + Decimal(854.0) * x[6],
      Decimal(1032.0) * x[0] + Decimal(1044.0) * x[1] + Decimal(1055.0) * x[2] +
          Decimal(1044.0) * x[3] + Decimal(1044.0) * x[4] +
          Decimal(1076.0) * x[5] + Decimal(1042.0) * x[6],
      Decimal(630.0) * x[0] + Decimal(630.0) * x[1] + Decimal(630.0) * x[2] +
          Decimal(630.0) * x[3] + Decimal(630.0) * x[4] +
          Decimal(596.0) * x[5] + Decimal(610.0) * x[6],
      Decimal(1371.0) * x[0] + Decimal(1371.0) * x[1] + Decimal(1355.0) * x[2] +
          Decimal(1371.0) * x[3] + Decimal(1371.0) * x[4] +
          Decimal(1364.0) * x[5] + Decimal(1371.0) * x[6],
      Decimal(1599.0) * x[0] + Decimal(1599.0) * x[1] + Decimal(1599.0) * x[2] +
          Decimal(1599.0) * x[3] + Decimal(1599.0) * x[4] +
          Decimal(1599.0) * x[5] + Decimal(1586.0) * x[6],
      Decimal(1800.0) * x[0] + Decimal(1800.0) * x[1] + Decimal(1800.0) * x[2] +
          Decimal(1800.0) * x[3] + Decimal(1800.0) * x[4] +
          Decimal(1770.0) * x[5] + Decimal(1837.0) * x[6],
      Decimal(448.0) * x[0] + Decimal(450.0) * x[1] + Decimal(450.0) * x[2] +
          Decimal(450.0) * x[3] + Decimal(450.0) * x[4] +
          Decimal(455.0) * x[5] + Decimal(445.0) * x[6],
      Decimal(1101.0) * x[0] + Decimal(1100.0) * x[1] + Decimal(1100.0) * x[2] +
          Decimal(1100.0) * x[3] + Decimal(1101.0) * x[4] +
          Decimal(1101.0) * x[5] + Decimal(1096.0) * x[6],
      Decimal(1831.0) * x[0] + Decimal(1831.0) * x[1] + Decimal(1837.0) * x[2] +
          Decimal(1831.0) * x[3] + Decimal(1829.0) * x[4] +
          Decimal(1873.0) * x[5] + Decimal(1854.0) * x[6],
      Decimal(1679.0) * x[0] + Decimal(1682.0) * x[1] + Decimal(1682.0) * x[2] +
          Decimal(1682.0) * x[3] + Decimal(1682.0) * x[4] +
          Decimal(1736.0) * x[5] + Decimal(1681.0) * x[6],
      Decimal(1040.0) * x[0] + Decimal(1040.0) * x[1] + Decimal(1040.0) * x[2] +
          Decimal(1040.0) * x[3] + Decimal(1040.0) * x[4] +
          Decimal(1083.0) * x[5] + Decimal(1033.0) * x[6],
      Decimal(921.0) * x[0] + Decimal(921.0) * x[1] + Decimal(921.0) * x[2] +
          Decimal(921.0) * x[3] + Decimal(996.0) * x[4] +
          Decimal(965.0) * x[5] + Decimal(1001.0) * x[6],
      Decimal(576.0) * x[0] + Decimal(576.0) * x[1] + Decimal(576.0) * x[2] +
          Decimal(576.0) * x[3] + Decimal(549.0) * x[4] +
          Decimal(588.0) * x[5] + Decimal(536.0) * x[6],
      Decimal(120.0) * x[0] + Decimal(120.0) * x[1] + Decimal(104.0) * x[2] +
          Decimal(104.0) * x[3] + Decimal(104.0) * x[4] + Decimal(82.0) * x[5] +
          Decimal(102.0) * x[6],
      Decimal(665.0) * x[0] + Decimal(665.0) * x[1] + Decimal(665.0) * x[2] +
          Decimal(663.0) * x[3] + Decimal(663.0) * x[4] +
          Decimal(694.0) * x[5] + Decimal(647.0) * x[6],
      Decimal(537.0) * x[0] + Decimal(537.0) * x[1] + Decimal(537.0) * x[2] +
          Decimal(537.0) * x[3] + Decimal(537.0) * x[4] +
          Decimal(564.0) * x[5] + Decimal(528.0) * x[6],
      Decimal(1489.0) * x[0] + Decimal(1489.0) * x[1] + Decimal(1489.0) * x[2] +
          Decimal(1530.0) * x[3] + Decimal(1530.0) * x[4] +
          Decimal(1521.0) * x[5] + Decimal(1524.0) * x[6],
      Decimal(674.0) * x[0] + Decimal(674.0) * x[1] + Decimal(674.0) * x[2] +
          Decimal(674.0) * x[3] + Decimal(674.0) * x[4] +
          Decimal(684.0) * x[5] + Decimal(684.0) * x[6],
      Decimal(631.0) * x[0] + Decimal(631.0) * x[1] + Decimal(631.0) * x[2] +
          Decimal(592.0) * x[3] + Decimal(619.0) * x[4] +
          Decimal(639.0) * x[5] + Decimal(587.0) * x[6],
      Decimal(2021.0) * x[0] + Decimal(2021.0) * x[1] + Decimal(2021.0) * x[2] +
          Decimal(2021.0) * x[3] + Decimal(2021.0) * x[4] +
          Decimal(2027.0) * x[5] + Decimal(2038.0) * x[6],
      Decimal(1618.0) * x[0] + Decimal(1618.0) * x[1] + Decimal(1618.0) * x[2] +
          Decimal(1618.0) * x[3] + Decimal(1618.0) * x[4] +
          Decimal(1642.0) * x[5] + Decimal(1626.0) * x[6],
      Decimal(1419.0) * x[0] + Decimal(1419.0) * x[1] + Decimal(1413.0) * x[2] +
          Decimal(1435.0) * x[3] + Decimal(1435.0) * x[4] +
          Decimal(1418.0) * x[5] + Decimal(1437.0) * x[6],
      Decimal(1358.0) * x[0] + Decimal(1355.0) * x[1] + Decimal(1331.0) * x[2] +
          Decimal(1319.0) * x[3] + Decimal(1343.0) * x[4] +
          Decimal(1380.0) * x[5] + Decimal(1297.0) * x[6],
      Decimal(635.0) * x[0] + Decimal(635.0) * x[1] + Decimal(635.0) * x[2] +
          Decimal(635.0) * x[3] + Decimal(647.0) * x[4] +
          Decimal(634.0) * x[5] + Decimal(683.0) * x[6],
      Decimal(339.0) * x[0] + Decimal(339.0) * x[1] + Decimal(330.0) * x[2] +
          Decimal(339.0) * x[3] + Decimal(339.0) * x[4] +
          Decimal(342.0) * x[5] + Decimal(346.0) * x[6],
      Decimal(1969.0) * x[0] + Decimal(1969.0) * x[1] + Decimal(1973.0) * x[2] +
          Decimal(1948.0) * x[3] + Decimal(1985.0) * x[4] +
          Decimal(1990.0) * x[5] + Decimal(1998.0) * x[6],
      Decimal(1049.0) * x[0] + Decimal(1049.0) * x[1] + Decimal(1049.0) * x[2] +
          Decimal(1050.0) * x[3] + Decimal(1050.0) * x[4] +
          Decimal(1078.0) * x[5] + Decimal(1081.0) * x[6],
      Decimal(345.0) * x[0] + Decimal(345.0) * x[1] + Decimal(345.0) * x[2] +
          Decimal(339.0) * x[3] + Decimal(339.0) * x[4] +
          Decimal(314.0) * x[5] + Decimal(328.0) * x[6],
      Decimal(1194.0) * x[0] + Decimal(1194.0) * x[1] + Decimal(1194.0) * x[2] +
          Decimal(1194.0) * x[3] + Decimal(1194.0) * x[4] +
          Decimal(1140.0) * x[5] + Decimal(1205.0) * x[6],
      Decimal(668.0) * x[0] + Decimal(666.0) * x[1] + Decimal(666.0) * x[2] +
          Decimal(666.0) * x[3] + Decimal(666.0) * x[4] +
          Decimal(677.0) * x[5] + Decimal(668.0) * x[6],
      Decimal(1460.0) * x[0] + Decimal(1448.0) * x[1] + Decimal(1478.0) * x[2] +
          Decimal(1482.0) * x[3] + Decimal(1482.0) * x[4] +
          Decimal(1504.0) * x[5] + Decimal(1446.0) * x[6],
      Decimal(1346.0) * x[0] + Decimal(1346.0) * x[1] + Decimal(1348.0) * x[2] +
          Decimal(1348.0) * x[3] + Decimal(1348.0) * x[4] +
          Decimal(1430.0) * x[5] + Decimal(1345.0) * x[6],
      Decimal(925.0) * x[0] + Decimal(925.0) * x[1] + Decimal(925.0) * x[2] +
          Decimal(907.0) * x[3] + Decimal(907.0) * x[4] +
          Decimal(943.0) * x[5] + Decimal(953.0) * x[6],
      Decimal(1354.0) * x[0] + Decimal(1368.0) * x[1] + Decimal(1334.0) * x[2] +
          Decimal(1334.0) * x[3] + Decimal(1334.0) * x[4] +
          Decimal(1346.0) * x[5] + Decimal(1350.0) * x[6],
      Decimal(974.0) * x[0] + Decimal(974.0) * x[1] + Decimal(974.0) * x[2] +
          Decimal(974.0) * x[3] + Decimal(974.0) * x[4] +
          Decimal(939.0) * x[5] + Decimal(974.0) * x[6],
      Decimal(930.0) * x[0] + Decimal(930.0) * x[1] + Decimal(930.0) * x[2] +
          Decimal(930.0) * x[3] + Decimal(930.0) * x[4] +
          Decimal(919.0) * x[5] + Decimal(944.0) * x[6],
      Decimal(1132.0) * x[0] + Decimal(1132.0) * x[1] + Decimal(1093.0) * x[2] +
          Decimal(1132.0) * x[3] + Decimal(1136.0) * x[4] +
          Decimal(1099.0) * x[5] + Decimal(1148.0) * x[6],
      Decimal(1213.0) * x[0] + Decimal(1213.0) * x[1] + Decimal(1213.0) * x[2] +
          Decimal(1213.0) * x[3] + Decimal(1213.0) * x[4] +
          Decimal(1210.0) * x[5] + Decimal(1197.0) * x[6],
      Decimal(715.0) * x[0] + Decimal(715.0) * x[1] + Decimal(715.0) * x[2] +
          Decimal(715.0) * x[3] + Decimal(715.0) * x[4] +
          Decimal(733.0) * x[5] + Decimal(734.0) * x[6],
      Decimal(1338.0) * x[0] + Decimal(1338.0) * x[1] + Decimal(1338.0) * x[2] +
          Decimal(1345.0) * x[3] + Decimal(1347.0) * x[4] +
          Decimal(1297.0) * x[5] + Decimal(1351.0) * x[6],
      Decimal(1186.0) * x[0] + Decimal(1186.0) * x[1] + Decimal(1188.0) * x[2] +
          Decimal(1188.0) * x[3] + Decimal(1186.0) * x[4] +
          Decimal(1186.0) * x[5] + Decimal(1175.0) * x[6],
      Decimal(1623.0) * x[0] + Decimal(1623.0) * x[1] + Decimal(1623.0) * x[2] +
          Decimal(1613.0) * x[3] + Decimal(1623.0) * x[4] +
          Decimal(1556.0) * x[5] + Decimal(1629.0) * x[6],
      Decimal(670.0) * x[0] + Decimal(670.0) * x[1] + Decimal(670.0) * x[2] +
          Decimal(663.0) * x[3] + Decimal(676.0) * x[4] +
          Decimal(719.0) * x[5] + Decimal(688.0) * x[6],
      Decimal(164.0) * x[0] + Decimal(174.0) * x[1] + Decimal(138.0) * x[2] +
          Decimal(174.0) * x[3] + Decimal(164.0) * x[4] +
          Decimal(123.0) * x[5] + Decimal(213.0) * x[6],
      Decimal(890.0) * x[0] + Decimal(890.0) * x[1] + Decimal(890.0) * x[2] +
          Decimal(890.0) * x[3] + Decimal(890.0) * x[4] +
          Decimal(888.0) * x[5] + Decimal(898.0) * x[6],
      Decimal(1501.0) * x[0] + Decimal(1501.0) * x[1] + Decimal(1501.0) * x[2] +
          Decimal(1501.0) * x[3] + Decimal(1501.0) * x[4] +
          Decimal(1501.0) * x[5] + Decimal(1506.0) * x[6],
      Decimal(508.0) * x[0] + Decimal(508.0) * x[1] + Decimal(508.0) * x[2] +
          Decimal(508.0) * x[3] + Decimal(500.0) * x[4] +
          Decimal(490.0) * x[5] + Decimal(519.0) * x[6],
      Decimal(1607.0) * x[0] + Decimal(1607.0) * x[1] + Decimal(1607.0) * x[2] +
          Decimal(1607.0) * x[3] + Decimal(1606.0) * x[4] +
          Decimal(1606.0) * x[5] + Decimal(1620.0) * x[6],
      Decimal(184.0) * x[0] + Decimal(184.0) * x[1] + Decimal(184.0) * x[2] +
          Decimal(184.0) * x[3] + Decimal(184.0) * x[4] +
          Decimal(189.0) * x[5] + Decimal(184.0) * x[6],
      Decimal(1483.0) * x[0] + Decimal(1483.0) * x[1] + Decimal(1487.0) * x[2] +
          Decimal(1483.0) * x[3] + Decimal(1483.0) * x[4] +
          Decimal(1535.0) * x[5] + Decimal(1395.0) * x[6],
      Decimal(788.0) * x[0] + Decimal(790.0) * x[1] + Decimal(796.0) * x[2] +
          Decimal(783.0) * x[3] + Decimal(782.0) * x[4] +
          Decimal(781.0) * x[5] + Decimal(847.0) * x[6],
      Decimal(768.0) * x[0] + Decimal(778.0) * x[1] + Decimal(790.0) * x[2] +
          Decimal(768.0) * x[3] + Decimal(778.0) * x[4] +
          Decimal(750.0) * x[5] + Decimal(851.0) * x[6],
      Decimal(721.0) * x[0] + Decimal(721.0) * x[1] + Decimal(721.0) * x[2] +
          Decimal(721.0) * x[3] + Decimal(721.0) * x[4] +
          Decimal(719.0) * x[5] + Decimal(700.0) * x[6],
      Decimal(1568.0) * x[0] + Decimal(1558.0) * x[1] + Decimal(1515.0) * x[2] +
          Decimal(1575.0) * x[3] + Decimal(1559.0) * x[4] +
          Decimal(1498.0) * x[5] + Decimal(1608.0) * x[6],
      Decimal(385.0) * x[0] + Decimal(385.0) * x[1] + Decimal(385.0) * x[2] +
          Decimal(385.0) * x[3] + Decimal(385.0) * x[4] +
          Decimal(383.0) * x[5] + Decimal(331.0) * x[6],
      Decimal(1226.0) * x[0] + Decimal(1226.0) * x[1] + Decimal(1226.0) * x[2] +
          Decimal(1226.0) * x[3] + Decimal(1224.0) * x[4] +
          Decimal(1224.0) * x[5] + Decimal(1224.0) * x[6],
      Decimal(155.0) * x[0] + Decimal(155.0) * x[1] + Decimal(116.0) * x[2] +
          Decimal(155.0) * x[3] + Decimal(155.0) * x[4] + Decimal(70.0) * x[5] +
          Decimal(152.0) * x[6],
      Decimal(441.0) * x[0] + Decimal(441.0) * x[1] + Decimal(378.0) * x[2] +
          Decimal(449.0) * x[3] + Decimal(447.0) * x[4] +
          Decimal(303.0) * x[5] + Decimal(459.0) * x[6],
      Decimal(712.0) * x[0] + Decimal(712.0) * x[1] + Decimal(706.0) * x[2] +
          Decimal(719.0) * x[3] + Decimal(720.0) * x[4] +
          Decimal(693.0) * x[5] + Decimal(743.0) * x[6],
      Decimal(1535.0) * x[0] + Decimal(1535.0) * x[1] + Decimal(1535.0) * x[2] +
          Decimal(1559.0) * x[3] + Decimal(1530.0) * x[4] +
          Decimal(1491.0) * x[5] + Decimal(1611.0) * x[6],
      Decimal(1196.0) * x[0] + Decimal(1196.0) * x[1] + Decimal(1140.0) * x[2] +
          Decimal(1220.0) * x[3] + Decimal(1194.0) * x[4] +
          Decimal(977.0) * x[5] + Decimal(1239.0) * x[6],
      Decimal(786.0) * x[0] + Decimal(786.0) * x[1] + Decimal(719.0) * x[2] +
          Decimal(930.0) * x[3] + Decimal(730.0) * x[4] +
          Decimal(569.0) * x[5] + Decimal(934.0) * x[6],
      Decimal(84.0) * x[0] + Decimal(76.0) * x[1] + Decimal(76.0) * x[2] +
          Decimal(77.0) * x[3] + Decimal(77.0) * x[4] + Decimal(48.0) * x[5] +
          Decimal(99.0) * x[6],
      Decimal(539.0) * x[0] + Decimal(539.0) * x[1] + Decimal(539.0) * x[2] +
          Decimal(538.0) * x[3] + Decimal(537.0) * x[4] +
          Decimal(523.0) * x[5] + Decimal(674.0) * x[6],
      Decimal(1105.0) * x[0] + Decimal(1105.0) * x[1] + Decimal(1066.0) * x[2] +
          Decimal(1098.0) * x[3] + Decimal(1087.0) * x[4] +
          Decimal(998.0) * x[5] + Decimal(1068.0) * x[6],
      Decimal(1524.0) * x[0] + Decimal(1506.0) * x[1] + Decimal(1506.0) * x[2] +
          Decimal(1524.0) * x[3] + Decimal(1528.0) * x[4] +
          Decimal(1479.0) * x[5] + Decimal(1508.0) * x[6],
      Decimal(167.0) * x[0] + Decimal(167.0) * x[1] + Decimal(167.0) * x[2] +
          Decimal(167.0) * x[3] + Decimal(147.0) * x[4] +
          Decimal(107.0) * x[5] + Decimal(168.0) * x[6],
      Decimal(1530.0) * x[0] + Decimal(1527.0) * x[1] + Decimal(1527.0) * x[2] +
          Decimal(1527.0) * x[3] + Decimal(1527.0) * x[4] +
          Decimal(1534.0) * x[5] + Decimal(1509.0) * x[6],
      Decimal(1993.0) * x[0] + Decimal(1993.0) * x[1] + Decimal(1993.0) * x[2] +
          Decimal(1993.0) * x[3] + Decimal(2013.0) * x[4] +
          Decimal(2002.0) * x[5] + Decimal(1967.0) * x[6],
      Decimal(1830.0) * x[0] + Decimal(1830.0) * x[1] + Decimal(1830.0) * x[2] +
          Decimal(1819.0) * x[3] + Decimal(1829.0) * x[4] +
          Decimal(1830.0) * x[5] + Decimal(1822.0) * x[6],
      Decimal(1518.0) * x[0] + Decimal(1518.0) * x[1] + Decimal(1518.0) * x[2] +
          Decimal(1518.0) * x[3] + Decimal(1505.0) * x[4] +
          Decimal(1508.0) * x[5] + Decimal(1494.0) * x[6],
      Decimal(790.0) * x[0] + Decimal(790.0) * x[1] + Decimal(790.0) * x[2] +
          Decimal(790.0) * x[3] + Decimal(790.0) * x[4] +
          Decimal(794.0) * x[5] + Decimal(805.0) * x[6],
      Decimal(1065.0) * x[0] + Decimal(1065.0) * x[1] + Decimal(1065.0) * x[2] +
          Decimal(1065.0) * x[3] + Decimal(1055.0) * x[4] +
          Decimal(1041.0) * x[5] + Decimal(1068.0) * x[6],
      Decimal(692.0) * x[0] + Decimal(692.0) * x[1] + Decimal(692.0) * x[2] +
          Decimal(692.0) * x[3] + Decimal(688.0) * x[4] +
          Decimal(642.0) * x[5] + Decimal(717.0) * x[6],
      Decimal(1128.0) * x[0] + Decimal(1128.0) * x[1] + Decimal(1116.0) * x[2] +
          Decimal(1128.0) * x[3] + Decimal(1128.0) * x[4] +
          Decimal(1118.0) * x[5] + Decimal(1109.0) * x[6],
      Decimal(1328.0) * x[0] + Decimal(1328.0) * x[1] + Decimal(1328.0) * x[2] +
          Decimal(1330.0) * x[3] + Decimal(1330.0) * x[4] +
          Decimal(1313.0) * x[5] + Decimal(1339.0) * x[6],
      Decimal(836.0) * x[0] + Decimal(836.0) * x[1] + Decimal(833.0) * x[2] +
          Decimal(836.0) * x[3] + Decimal(836.0) * x[4] +
          Decimal(851.0) * x[5] + Decimal(917.0) * x[6],
      Decimal(823.0) * x[0] + Decimal(823.0) * x[1] + Decimal(823.0) * x[2] +
          Decimal(823.0) * x[3] + Decimal(823.0) * x[4] +
          Decimal(823.0) * x[5] + Decimal(857.0) * x[6],
      Decimal(672.0) * x[0] + Decimal(672.0) * x[1] + Decimal(672.0) * x[2] +
          Decimal(680.0) * x[3] + Decimal(670.0) * x[4] +
          Decimal(517.0) * x[5] + Decimal(740.0) * x[6],
      Decimal(519.0) * x[0] + Decimal(519.0) * x[1] + Decimal(519.0) * x[2] +
          Decimal(492.0) * x[3] + Decimal(480.0) * x[4] +
          Decimal(479.0) * x[5] + Decimal(510.0) * x[6],
      Decimal(1690.0) * x[0] + Decimal(1693.0) * x[1] + Decimal(1693.0) * x[2] +
          Decimal(1693.0) * x[3] + Decimal(1693.0) * x[4] +
          Decimal(1696.0) * x[5] + Decimal(1677.0) * x[6],
      Decimal(410.0) * x[0] + Decimal(410.0) * x[1] + Decimal(410.0) * x[2] +
          Decimal(410.0) * x[3] + Decimal(418.0) * x[4] +
          Decimal(405.0) * x[5] + Decimal(421.0) * x[6],
      Decimal(2102.0) * x[0] + Decimal(2102.0) * x[1] + Decimal(2102.0) * x[2] +
          Decimal(2102.0) * x[3] + Decimal(2094.0) * x[4] +
          Decimal(2094.0) * x[5] + Decimal(2073.0) * x[6],
      Decimal(82.0) * x[0] + Decimal(100.0) * x[1] + Decimal(100.0) * x[2] +
          Decimal(100.0) * x[3] + Decimal(100.0) * x[4] +
          Decimal(146.0) * x[5] + Decimal(115.0) * x[6],
      Decimal(303.0) * x[0] + Decimal(303.0) * x[1] + Decimal(303.0) * x[2] +
          Decimal(303.0) * x[3] + Decimal(303.0) * x[4] +
          Decimal(303.0) * x[5] + Decimal(319.0) * x[6],
      Decimal(481.0) * x[0] + Decimal(481.0) * x[1] + Decimal(481.0) * x[2] +
          Decimal(481.0) * x[3] + Decimal(481.0) * x[4] +
          Decimal(484.0) * x[5] + Decimal(498.0) * x[6],
      Decimal(1546.0) * x[0] + Decimal(1546.0) * x[1] + Decimal(1546.0) * x[2] +
          Decimal(1546.0) * x[3] + Decimal(1535.0) * x[4] +
          Decimal(1532.0) * x[5] + Decimal(1505.0) * x[6],
      Decimal(1384.0) * x[0] + Decimal(1384.0) * x[1] + Decimal(1384.0) * x[2] +
          Decimal(1384.0) * x[3] + Decimal(1384.0) * x[4] +
          Decimal(1382.0) * x[5] + Decimal(1385.0) * x[6],
      Decimal(824.0) * x[0] + Decimal(824.0) * x[1] + Decimal(824.0) * x[2] +
          Decimal(824.0) * x[3] + Decimal(824.0) * x[4] +
          Decimal(767.0) * x[5] + Decimal(859.0) * x[6],
      Decimal(1778.0) * x[0] + Decimal(1778.0) * x[1] + Decimal(1778.0) * x[2] +
          Decimal(1778.0) * x[3] + Decimal(1758.0) * x[4] +
          Decimal(1752.0) * x[5] + Decimal(1854.0) * x[6],
      Decimal(192.0) * x[0] + Decimal(192.0) * x[1] + Decimal(192.0) * x[2] +
          Decimal(193.0) * x[3] + Decimal(185.0) * x[4] +
          Decimal(236.0) * x[5] + Decimal(169.0) * x[6],
      Decimal(1476.0) * x[0] + Decimal(1476.0) * x[1] + Decimal(1478.0) * x[2] +
          Decimal(1478.0) * x[3] + Decimal(1459.0) * x[4] +
          Decimal(1457.0) * x[5] + Decimal(1574.0) * x[6],
      Decimal(429.0) * x[0] + Decimal(429.0) * x[1] + Decimal(429.0) * x[2] +
          Decimal(429.0) * x[3] + Decimal(430.0) * x[4] +
          Decimal(430.0) * x[5] + Decimal(430.0) * x[6],
      Decimal(1076.0) * x[0] + Decimal(1076.0) * x[1] + Decimal(1074.0) * x[2] +
          Decimal(1068.0) * x[3] + Decimal(1087.0) * x[4] +
          Decimal(1079.0) * x[5] + Decimal(1105.0) * x[6],
      Decimal(691.0) * x[0] + Decimal(691.0) * x[1] + Decimal(691.0) * x[2] +
          Decimal(693.0) * x[3] + Decimal(693.0) * x[4] +
          Decimal(695.0) * x[5] + Decimal(747.0) * x[6],
      Decimal(1569.0) * x[0] + Decimal(1569.0) * x[1] + Decimal(1569.0) * x[2] +
          Decimal(1569.0) * x[3] + Decimal(1559.0) * x[4] +
          Decimal(1562.0) * x[5] + Decimal(1582.0) * x[6],
      Decimal(1454.0) * x[0] + Decimal(1454.0) * x[1] + Decimal(1453.0) * x[2] +
          Decimal(1454.0) * x[3] + Decimal(1454.0) * x[4] +
          Decimal(1388.0) * x[5] + Decimal(1478.0) * x[6],
      Decimal(220.0) * x[0] + Decimal(220.0) * x[1] + Decimal(220.0) * x[2] +
          Decimal(220.0) * x[3] + Decimal(209.0) * x[4] +
          Decimal(208.0) * x[5] + Decimal(175.0) * x[6],
      Decimal(2237.0) * x[0] + Decimal(2237.0) * x[1] + Decimal(2203.0) * x[2] +
          Decimal(2237.0) * x[3] + Decimal(2231.0) * x[4] +
          Decimal(2157.0) * x[5] + Decimal(2198.0) * x[6],
      Decimal(756.0) * x[0] + Decimal(760.0) * x[1] + Decimal(713.0) * x[2] +
          Decimal(760.0) * x[3] + Decimal(756.0) * x[4] +
          Decimal(709.0) * x[5] + Decimal(731.0) * x[6],
      Decimal(256.0) * x[0] + Decimal(256.0) * x[1] + Decimal(256.0) * x[2] +
          Decimal(256.0) * x[3] + Decimal(328.0) * x[4] +
          Decimal(334.0) * x[5] + Decimal(251.0) * x[6],
      Decimal(21.0) * x[0] + Decimal(21.0) * x[1] + Decimal(21.0) * x[2] +
          Decimal(21.0) * x[3] + Decimal(21.0) * x[4] - Decimal(22.0) * x[5] +
          Decimal(30.0) * x[6],
      Decimal(1079.0) * x[0] + Decimal(1079.0) * x[1] + Decimal(1079.0) * x[2] +
          Decimal(1079.0) * x[3] + Decimal(1079.0) * x[4] +
          Decimal(1087.0) * x[5] + Decimal(1066.0) * x[6],
      Decimal(1517.0) * x[0] + Decimal(1517.0) * x[1] + Decimal(1510.0) * x[2] +
          Decimal(1511.0) * x[3] + Decimal(1511.0) * x[4] +
          Decimal(1510.0) * x[5] + Decimal(1495.0) * x[6],
      Decimal(408.0) * x[0] + Decimal(408.0) * x[1] + Decimal(408.0) * x[2] +
          Decimal(408.0) * x[3] + Decimal(432.0) * x[4] +
          Decimal(434.0) * x[5] + Decimal(455.0) * x[6],
      Decimal(2150.0) * x[0] + Decimal(2150.0) * x[1] + Decimal(2150.0) * x[2] +
          Decimal(2149.0) * x[3] + Decimal(2149.0) * x[4] +
          Decimal(2182.0) * x[5] + Decimal(2141.0) * x[6],
      Decimal(1615.0) * x[0] + Decimal(1615.0) * x[1] + Decimal(1615.0) * x[2] +
          Decimal(1615.0) * x[3] + Decimal(1615.0) * x[4] +
          Decimal(1568.0) * x[5] + Decimal(1624.0) * x[6],
      Decimal(855.0) * x[0] + Decimal(855.0) * x[1] + Decimal(855.0) * x[2] +
          Decimal(855.0) * x[3] + Decimal(855.0) * x[4] +
          Decimal(821.0) * x[5] + Decimal(870.0) * x[6],
      Decimal(1006.0) * x[0] + Decimal(1006.0) * x[1] + Decimal(1007.0) * x[2] +
          Decimal(1006.0) * x[3] + Decimal(950.0) * x[4] +
          Decimal(948.0) * x[5] + Decimal(1023.0) * x[6],
      Decimal(232.0) * x[0] + Decimal(232.0) * x[1] + Decimal(198.0) * x[2] +
          Decimal(198.0) * x[3] + Decimal(198.0) * x[4] +
          Decimal(164.0) * x[5] + Decimal(220.0) * x[6],
      Decimal(430.0) * x[0] + Decimal(430.0) * x[1] + Decimal(430.0) * x[2] +
          Decimal(430.0) * x[3] + Decimal(430.0) * x[4] +
          Decimal(424.0) * x[5] + Decimal(419.0) * x[6],
      Decimal(766.0) * x[0] + Decimal(766.0) * x[1] + Decimal(766.0) * x[2] +
          Decimal(766.0) * x[3] + Decimal(766.0) * x[4] +
          Decimal(766.0) * x[5] + Decimal(816.0) * x[6],
      Decimal(885.0) * x[0] + Decimal(885.0) * x[1] + Decimal(885.0) * x[2] +
          Decimal(885.0) * x[3] + Decimal(866.0) * x[4] +
          Decimal(917.0) * x[5] + Decimal(883.0) * x[6],
      Decimal(1617.0) * x[0] + Decimal(1617.0) * x[1] + Decimal(1617.0) * x[2] +
          Decimal(1615.0) * x[3] + Decimal(1615.0) * x[4] +
          Decimal(1601.0) * x[5] + Decimal(1607.0) * x[6],
      Decimal(789.0) * x[0] + Decimal(789.0) * x[1] + Decimal(785.0) * x[2] +
          Decimal(813.0) * x[3] + Decimal(789.0) * x[4] +
          Decimal(719.0) * x[5] + Decimal(867.0) * x[6],
      Decimal(154.0) * x[0] + Decimal(154.0) * x[1] + Decimal(154.0) * x[2] +
          Decimal(154.0) * x[3] + Decimal(154.0) * x[4] +
          Decimal(173.0) * x[5] + Decimal(238.0) * x[6],
      Decimal(1442.0) * x[0] + Decimal(1442.0) * x[1] + Decimal(1426.0) * x[2] +
          Decimal(1426.0) * x[3] + Decimal(1426.0) * x[4] +
          Decimal(1455.0) * x[5] + Decimal(1425.0) * x[6],
      Decimal(983.0) * x[0] + Decimal(983.0) * x[1] + Decimal(983.0) * x[2] +
          Decimal(983.0) * x[3] + Decimal(1002.0) * x[4] +
          Decimal(909.0) * x[5] + Decimal(1009.0) * x[6],
      Decimal(192.0) * x[0] + Decimal(192.0) * x[1] + Decimal(192.0) * x[2] +
          Decimal(192.0) * x[3] + Decimal(192.0) * x[4] +
          Decimal(219.0) * x[5] + Decimal(183.0) * x[6],
      Decimal(192.0) * x[0] + Decimal(192.0) * x[1] + Decimal(192.0) * x[2] +
          Decimal(184.0) * x[3] + Decimal(176.0) * x[4] +
          Decimal(160.0) * x[5] + Decimal(168.0) * x[6],
      Decimal(1479.0) * x[0] + Decimal(1479.0) * x[1] + Decimal(1486.0) * x[2] +
          Decimal(1477.0) * x[3] + Decimal(1495.0) * x[4] +
          Decimal(1542.0) * x[5] + Decimal(1527.0) * x[6],
      Decimal(1844.0) * x[0] + Decimal(1844.0) * x[1] + Decimal(1857.0) * x[2] +
          Decimal(1810.0) * x[3] + Decimal(1844.0) * x[4] +
          Decimal(1870.0) * x[5] + Decimal(1753.0) * x[6],
      Decimal(214.0) * x[0] + Decimal(211.0) * x[1] + Decimal(211.0) * x[2] +
          Decimal(193.0) * x[3] + Decimal(184.0) * x[4] +
          Decimal(162.0) * x[5] + Decimal(183.0) * x[6],
      Decimal(1387.0) * x[0] + Decimal(1387.0) * x[1] + Decimal(1387.0) * x[2] +
          Decimal(1363.0) * x[3] + Decimal(1387.0) * x[4] +
          Decimal(1365.0) * x[5] + Decimal(1365.0) * x[6],
      Decimal(436.0) * x[0] + Decimal(436.0) * x[1] + Decimal(436.0) * x[2] +
          Decimal(436.0) * x[3] + Decimal(428.0) * x[4] +
          Decimal(428.0) * x[5] + Decimal(452.0) * x[6],
      Decimal(334.0) * x[0] + Decimal(334.0) * x[1] + Decimal(332.0) * x[2] +
          Decimal(332.0) * x[3] + Decimal(344.0) * x[4] +
          Decimal(344.0) * x[5] + Decimal(316.0) * x[6],
      Decimal(1897.0) * x[0] + Decimal(1896.0) * x[1] + Decimal(1896.0) * x[2] +
          Decimal(1893.0) * x[3] + Decimal(1894.0) * x[4] +
          Decimal(1895.0) * x[5] + Decimal(1898.0) * x[6],
      Decimal(1937.0) * x[0] + Decimal(1937.0) * x[1] + Decimal(1940.0) * x[2] +
          Decimal(1937.0) * x[3] + Decimal(1937.0) * x[4] +
          Decimal(1950.0) * x[5] + Decimal(1927.0) * x[6],
      Decimal(1256.0) * x[0] + Decimal(1256.0) * x[1] + Decimal(1256.0) * x[2] +
          Decimal(1256.0) * x[3] + Decimal(1256.0) * x[4] +
          Decimal(1250.0) * x[5] + Decimal(1261.0) * x[6],
      Decimal(1254.0) * x[0] + Decimal(1254.0) * x[1] + Decimal(1254.0) * x[2] +
          Decimal(1254.0) * x[3] + Decimal(1255.0) * x[4] +
          Decimal(1252.0) * x[5] + Decimal(1252.0) * x[6],
      Decimal(896.0) * x[0] + Decimal(896.0) * x[1] + Decimal(896.0) * x[2] +
          Decimal(896.0) * x[3] + Decimal(907.0) * x[4] +
          Decimal(898.0) * x[5] + Decimal(891.0) * x[6],
      Decimal(617.0) * x[0] + Decimal(620.0) * x[1] + Decimal(620.0) * x[2] +
          Decimal(610.0) * x[3] + Decimal(620.0) * x[4] +
          Decimal(658.0) * x[5] + Decimal(630.0) * x[6],
      Decimal(53.0) * x[0] + Decimal(53.0) * x[1] + Decimal(49.0) * x[2] +
          Decimal(88.0) * x[3] + Decimal(41.0) * x[4] + Decimal(43.0) * x[5] +
          Decimal(92.0) * x[6],
      Decimal(608.0) * x[0] + Decimal(608.0) * x[1] + Decimal(608.0) * x[2] +
          Decimal(607.0) * x[3] + Decimal(605.0) * x[4] +
          Decimal(584.0) * x[5] + Decimal(575.0) * x[6],
      Decimal(1246.0) * x[0] + Decimal(1246.0) * x[1] + Decimal(1246.0) * x[2] +
          Decimal(1246.0) * x[3] + Decimal(1246.0) * x[4] +
          Decimal(1246.0) * x[5] + Decimal(1246.0) * x[6],
      Decimal(807.0) * x[0] + Decimal(807.0) * x[1] + Decimal(807.0) * x[2] +
          Decimal(807.0) * x[3] + Decimal(807.0) * x[4] +
          Decimal(774.0) * x[5] + Decimal(833.0) * x[6],
      Decimal(546.0) * x[0] + Decimal(558.0) * x[1] + Decimal(519.0) * x[2] +
          Decimal(558.0) * x[3] + Decimal(557.0) * x[4] +
          Decimal(547.0) * x[5] + Decimal(565.0) * x[6],
      Decimal(202.0) * x[0] + Decimal(202.0) * x[1] + Decimal(143.0) * x[2] +
          Decimal(202.0) * x[3] + Decimal(131.0) * x[4] - x[5] +
          Decimal(244.0) * x[6],
      Decimal(1757.0) * x[0] + Decimal(1757.0) * x[1] + Decimal(1757.0) * x[2] +
          Decimal(1757.0) * x[3] + Decimal(1737.0) * x[4] +
          Decimal(1754.0) * x[5] + Decimal(1715.0) * x[6],
      Decimal(1152.0) * x[0] + Decimal(1152.0) * x[1] + Decimal(1149.0) * x[2] +
          Decimal(1152.0) * x[3] + Decimal(1152.0) * x[4] +
          Decimal(1111.0) * x[5] + Decimal(1163.0) * x[6],
      Decimal(829.0) * x[0] + Decimal(829.0) * x[1] + Decimal(829.0) * x[2] +
          Decimal(829.0) * x[3] + Decimal(829.0) * x[4] +
          Decimal(829.0) * x[5] + Decimal(832.0) * x[6],
      Decimal(1380.0) * x[0] + Decimal(1380.0) * x[1] + Decimal(1345.0) * x[2] +
          Decimal(1380.0) * x[3] + Decimal(1380.0) * x[4] +
          Decimal(1366.0) * x[5] + Decimal(1431.0) * x[6],
      Decimal(759.0) * x[0] + Decimal(759.0) * x[1] + Decimal(759.0) * x[2] +
          Decimal(759.0) * x[3] + Decimal(759.0) * x[4] +
          Decimal(784.0) * x[5] + Decimal(748.0) * x[6],
      Decimal(127.0) * x[0] + Decimal(127.0) * x[1] + Decimal(127.0) * x[2] +
          Decimal(127.0) * x[3] + Decimal(127.0) * x[4] +
          Decimal(113.0) * x[5] + Decimal(126.0) * x[6],
      Decimal(1628.0) * x[0] + Decimal(1628.0) * x[1] + Decimal(1628.0) * x[2] +
          Decimal(1589.0) * x[3] + Decimal(1589.0) * x[4] +
          Decimal(1626.0) * x[5] + Decimal(1642.0) * x[6],
      Decimal(2052.0) * x[0] + Decimal(2052.0) * x[1] + Decimal(2052.0) * x[2] +
          Decimal(2052.0) * x[3] + Decimal(2052.0) * x[4] +
          Decimal(2054.0) * x[5] + Decimal(1991.0) * x[6],
      Decimal(1130.0) * x[0] + Decimal(1139.0) * x[1] + Decimal(1112.0) * x[2] +
          Decimal(1139.0) * x[3] + Decimal(1110.0) * x[4] +
          Decimal(1110.0) * x[5] + Decimal(1125.0) * x[6],
      Decimal(590.0) * x[0] + Decimal(590.0) * x[1] + Decimal(590.0) * x[2] +
          Decimal(590.0) * x[3] + Decimal(590.0) * x[4] +
          Decimal(613.0) * x[5] + Decimal(651.0) * x[6],
      Decimal(1449.0) * x[0] + Decimal(1449.0) * x[1] + Decimal(1449.0) * x[2] +
          Decimal(1449.0) * x[3] + Decimal(1449.0) * x[4] +
          Decimal(1401.0) * x[5] + Decimal(1544.0) * x[6],
      Decimal(201.0) * x[0] + Decimal(201.0) * x[1] + Decimal(201.0) * x[2] +
          Decimal(201.0) * x[3] + Decimal(201.0) * x[4] +
          Decimal(215.0) * x[5] + Decimal(274.0) * x[6],
      Decimal(859.0) * x[0] + Decimal(847.0) * x[1] + Decimal(847.0) * x[2] +
          Decimal(847.0) * x[3] + Decimal(835.0) * x[4] +
          Decimal(838.0) * x[5] + Decimal(815.0) * x[6],
      Decimal(1799.0) * x[0] + Decimal(1799.0) * x[1] + Decimal(1799.0) * x[2] +
          Decimal(1792.0) * x[3] + Decimal(1799.0) * x[4] +
          Decimal(1799.0) * x[5] + Decimal(1843.0) * x[6],
      Decimal(1033.0) * x[0] + Decimal(1033.0) * x[1] + Decimal(1038.0) * x[2] +
          Decimal(1031.0) * x[3] + Decimal(1025.0) * x[4] +
          Decimal(1027.0) * x[5] + Decimal(1030.0) * x[6],
      Decimal(183.0) * x[0] + Decimal(183.0) * x[1] + Decimal(183.0) * x[2] +
          Decimal(190.0) * x[3] + Decimal(183.0) * x[4] +
          Decimal(188.0) * x[5] + Decimal(152.0) * x[6],
      Decimal(654.0) * x[0] + Decimal(654.0) * x[1] + Decimal(648.0) * x[2] +
          Decimal(654.0) * x[3] + Decimal(665.0) * x[4] +
          Decimal(664.0) * x[5] + Decimal(632.0) * x[6],
      Decimal(1345.0) * x[0] + Decimal(1345.0) * x[1] + Decimal(1345.0) * x[2] +
          Decimal(1345.0) * x[3] + Decimal(1345.0) * x[4] +
          Decimal(1360.0) * x[5] + Decimal(1361.0) * x[6],
      Decimal(1361.0) * x[0] + Decimal(1361.0) * x[1] + Decimal(1370.0) * x[2] +
          Decimal(1351.0) * x[3] + Decimal(1361.0) * x[4] +
          Decimal(1370.0) * x[5] + Decimal(1343.0) * x[6],
      Decimal(1068.0) * x[0] + Decimal(1068.0) * x[1] + Decimal(1068.0) * x[2] +
          Decimal(1068.0) * x[3] + Decimal(1068.0) * x[4] +
          Decimal(1148.0) * x[5] + Decimal(1050.0) * x[6],
      Decimal(457.0) * x[0] + Decimal(457.0) * x[1] + Decimal(461.0) * x[2] +
          Decimal(457.0) * x[3] + Decimal(444.0) * x[4] +
          Decimal(475.0) * x[5] + Decimal(451.0) * x[6],
      Decimal(244.0) * x[0] + Decimal(244.0) * x[1] + Decimal(235.0) * x[2] +
          Decimal(256.0) * x[3] + Decimal(248.0) * x[4] +
          Decimal(131.0) * x[5] + Decimal(306.0) * x[6],
      Decimal(1819.0) * x[0] + Decimal(1819.0) * x[1] + Decimal(1819.0) * x[2] +
          Decimal(1809.0) * x[3] + Decimal(1819.0) * x[4] +
          Decimal(1830.0) * x[5] + Decimal(1802.0) * x[6],
      Decimal(484.0) * x[0] + Decimal(484.0) * x[1] + Decimal(484.0) * x[2] +
          Decimal(484.0) * x[3] + Decimal(484.0) * x[4] +
          Decimal(484.0) * x[5] + Decimal(467.0) * x[6],
      Decimal(334.0) * x[0] + Decimal(334.0) * x[1] + Decimal(334.0) * x[2] +
          Decimal(334.0) * x[3] + Decimal(334.0) * x[4] +
          Decimal(323.0) * x[5] + Decimal(335.0) * x[6],
      Decimal(1781.0) * x[0] + Decimal(1781.0) * x[1] + Decimal(1781.0) * x[2] +
          Decimal(1781.0) * x[3] + Decimal(1777.0) * x[4] +
          Decimal(1777.0) * x[5] + Decimal(1829.0) * x[6],
      Decimal(654.0) * x[0] + Decimal(654.0) * x[1] + Decimal(658.0) * x[2] +
          Decimal(654.0) * x[3] + Decimal(653.0) * x[4] +
          Decimal(666.0) * x[5] + Decimal(662.0) * x[6],
      Decimal(45.0) * x[0] + Decimal(45.0) * x[1] + Decimal(45.0) * x[2] +
          Decimal(45.0) * x[3] + Decimal(46.0) * x[4] + Decimal(46.0) * x[5] +
          Decimal(46.0) * x[6],
      Decimal(137.0) * x[0] + Decimal(137.0) * x[1] + Decimal(137.0) * x[2] +
          Decimal(137.0) * x[3] + Decimal(137.0) * x[4] +
          Decimal(113.0) * x[5] + Decimal(123.0) * x[6],
      Decimal(1660.0) * x[0] + Decimal(1660.0) * x[1] + Decimal(1660.0) * x[2] +
          Decimal(1660.0) * x[3] + Decimal(1735.0) * x[4] +
          Decimal(1732.0) * x[5] + Decimal(1655.0) * x[6],
      Decimal(1134.0) * x[0] + Decimal(1134.0) * x[1] + Decimal(1134.0) * x[2] +
          Decimal(1134.0) * x[3] + Decimal(1134.0) * x[4] +
          Decimal(1164.0) * x[5] + Decimal(1150.0) * x[6],
      Decimal(1664.0) * x[0] + Decimal(1664.0) * x[1] + Decimal(1660.0) * x[2] +
          Decimal(1664.0) * x[3] + Decimal(1664.0) * x[4] +
          Decimal(1621.0) * x[5] + Decimal(1677.0) * x[6],
      Decimal(1473.0) * x[0] + Decimal(1473.0) * x[1] + Decimal(1473.0) * x[2] +
          Decimal(1473.0) * x[3] + Decimal(1473.0) * x[4] +
          Decimal(1473.0) * x[5] + Decimal(1457.0) * x[6],
      Decimal(54.0) * x[0] + Decimal(54.0) * x[1] + Decimal(56.0) * x[2] +
          Decimal(56.0) * x[3] + Decimal(54.0) * x[4] + Decimal(54.0) * x[5] +
          Decimal(58.0) * x[6],
      Decimal(936.0) * x[0] + Decimal(936.0) * x[1] + Decimal(938.0) * x[2] +
          Decimal(938.0) * x[3] + Decimal(938.0) * x[4] +
          Decimal(944.0) * x[5] + Decimal(923.0) * x[6],
      Decimal(1131.0) * x[0] + Decimal(1131.0) * x[1] + Decimal(1131.0) * x[2] +
          Decimal(1131.0) * x[3] + Decimal(1131.0) * x[4] +
          Decimal(1160.0) * x[5] + Decimal(1136.0) * x[6],
      Decimal(1273.0) * x[0] + Decimal(1273.0) * x[1] + Decimal(1273.0) * x[2] +
          Decimal(1273.0) * x[3] + Decimal(1281.0) * x[4] +
          Decimal(1275.0) * x[5] + Decimal(1276.0) * x[6],
      Decimal(486.0) * x[0] + Decimal(486.0) * x[1] + Decimal(486.0) * x[2] +
          Decimal(486.0) * x[3] + Decimal(473.0) * x[4] +
          Decimal(486.0) * x[5] + Decimal(518.0) * x[6],
      Decimal(781.0) * x[0] + Decimal(781.0) * x[1] + Decimal(782.0) * x[2] +
          Decimal(781.0) * x[3] + Decimal(769.0) * x[4] +
          Decimal(762.0) * x[5] + Decimal(776.0) * x[6],
      Decimal(1323.0) * x[0] + Decimal(1323.0) * x[1] + Decimal(1323.0) * x[2] +
          Decimal(1323.0) * x[3] + Decimal(1323.0) * x[4] +
          Decimal(1338.0) * x[5] + Decimal(1302.0) * x[6],
      Decimal(50.0) * x[0] + Decimal(50.0) * x[1] + Decimal(48.0) * x[2] +
          Decimal(48.0) * x[3] + Decimal(48.0) * x[4] + Decimal(51.0) * x[5] +
          Decimal(47.0) * x[6],
      Decimal(1920.0) * x[0] + Decimal(1920.0) * x[1] + Decimal(1881.0) * x[2] +
          Decimal(1920.0) * x[3] + Decimal(1920.0) * x[4] +
          Decimal(1851.0) * x[5] + Decimal(1934.0) * x[6],
      Decimal(166.0) * x[0] + Decimal(166.0) * x[1] + Decimal(166.0) * x[2] +
          Decimal(165.0) * x[3] + Decimal(234.0) * x[4] +
          Decimal(234.0) * x[5] - Decimal(34.0) * x[6],
      Decimal(90.0) * x[0] + Decimal(93.0) * x[1] + Decimal(132.0) * x[2] +
          Decimal(132.0) * x[3] + Decimal(163.0) * x[4] +
          Decimal(180.0) * x[5] - Decimal(3.0) * x[6],
      Decimal(38.0) * x[0] + Decimal(38.0) * x[1] - Decimal(20.0) * x[2] +
          Decimal(36.0) * x[3] + Decimal(154.0) * x[4] + Decimal(516.0) * x[5] -
          Decimal(66.0) * x[6],
      Decimal(119.0) * x[0] + Decimal(127.0) * x[1] + Decimal(127.0) * x[2] +
          Decimal(117.0) * x[3] + Decimal(138.0) * x[4] +
          Decimal(384.0) * x[5] - Decimal(68.0) * x[6],
      Decimal(3.0) * x[0] - Decimal(2.0) * x[1] - Decimal(52.0) * x[2] +
          Decimal(86.0) * x[3] + Decimal(168.0) * x[4] + Decimal(414.0) * x[5] -
          Decimal(85.0) * x[6],
      Decimal(106.0) * x[0] + Decimal(106.0) * x[1] + Decimal(156.0) * x[2] +
          Decimal(156.0) * x[3] + Decimal(156.0) * x[4] +
          Decimal(417.0) * x[5] - Decimal(101.0) * x[6]

  };
  ExactBoxType C = ExactBoxType{
      {1, 1},    {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf}, {0, +inf},
      {0, +inf}};
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