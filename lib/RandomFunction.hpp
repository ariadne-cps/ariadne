#pragma once
#include <random>

enum F_TYPE {
  F_TEST_1 = 0,
  G_TEST_1 = 1,
  F_TEST_2 = 2,
  G_TEST_2 = 3,
  G_TEST_3 = 4,
  F_TEST_4 = 5,
  G_TEST_4 = 6,
  F_TEST_5 = 7,
  F_TEST_6 = 8,
  F_TEST_7 = 9
};

template <class X, F_TYPE Y> class FunctionDistribution;

//! Class to create a polynomial random
class RandomPolynomialEngine {
  template <class, F_TYPE> friend class FunctionDistribution;

public:
  //! Create an engine using a seed (to generate random coefficients) and a
  //! max order
  RandomPolynomialEngine(int _seed, int _max_order = 0) {
    this->seed = _seed;
    this->max_order = static_cast<unsigned>(_max_order);
  }

private:
  int seed;
  unsigned max_order;
};

//! Class to create a random distribution, using a random engine
template <class X, F_TYPE Y> class FunctionDistribution {
public:
  //! Default constructor
  FunctionDistribution() = default;
  //! Overloading of function operator to create a random distribution using
  //! an engine and a seed. If no seed are given, engine.seed is used.
  X operator()(RandomPolynomialEngine &_engine, int _seed);
};

#include "impl/RandomFunction.i.hpp"
