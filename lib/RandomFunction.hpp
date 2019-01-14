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
    F_TEST_5 = 7
};

template <class X, F_TYPE Y> class FunctionDistribution;

//! Class to create a polynomial random
class RandomPolynomialEngine
{
    template <class, F_TYPE> friend class FunctionDistribution;

  public:
    //! Create an engine using a seed (to generate random coefficients) and a
    //! max order
    RandomPolynomialEngine(unsigned seed, int max_order = 0)
    {
        this->seed = seed;
        this->max_order = static_cast<unsigned>(max_order) + 1;
        int_distribution = std::uniform_int_distribution<int>(0, max_order);
    }

  private:
    unsigned seed;
    std::uniform_int_distribution<int> int_distribution;
    unsigned max_order;
};

//! Class to create a random distribution, using a random engine
template <class X, F_TYPE Y> class FunctionDistribution
{
  public:
    //! Default constructor
    FunctionDistribution() = default;
    //! Overloading of funciton operator to create a random distribution using
    //! an engine and a seed. If no seed are given, engine.seed is used.
    X operator()(RandomPolynomialEngine &engine, int seed = -1);

  private:
    size_t size;
    unsigned seed;
};

#include "impl/RandomFunction.i.hpp"
