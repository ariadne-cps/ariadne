<<<<<<< HEAD
template <>
class FunctionDistribution<Ariadne::EffectiveScalarMultivariateFunction,
<<<<<<< HEAD
                           F_TEST_1> {
  using VType = Ariadne::EffectiveFormula;
  using X = Ariadne::EffectiveScalarMultivariateFunction;

public:
  FunctionDistribution(size_t _size, std::vector<float> &_v) {
    this->size = _size;
    generator_bounds = _v;
    x = std::vector<VType>(size + 1);
    unsigned i = 0;
    for (; i < size; ++i)
      x[i] = VType::coordinate(i);
    x[size] = VType::coordinate(size);
  }
  X operator()(RandomPolynomialEngine &_engine, int _seed = -1) {
    if (_seed == -1) {
      seed = _engine.seed;
    } else {
      seed = _seed;
    }

    std::default_random_engine generator(seed);

    int seed_tmp = seed;

    auto guess = [&seed_tmp, &generator]() -> bool {
      std::uniform_int_distribution<int> coin(0, 1);
      return coin(generator);
    };

    Ariadne::EffectiveFormula polynomial;

    // for (int j = static_cast<int>(_engine.max_order); j > 0; --j) {
    //   for (unsigned i = 0; i < size; ++i) {
    //     if (guess())
    //       continue;
    //     std::uniform_real_distribution<float> alpha_distribution(
    //         -generator_bounds[static_cast<unsigned>(j) - 1u],
    //         generator_bounds[static_cast<unsigned>(j) - 1u]);
    //     polynomial = Ariadne::EffectiveFormula::binary(
    //         Ariadne::OperatorCode::ADD, polynomial,
    //         Ariadne::EffectiveFormula::binary(
    //             Ariadne::OperatorCode::MUL,
    //             Ariadne::EffectiveFormula::constant(Ariadne::Real(
    //                 generator_bounds[static_cast<unsigned>(j) - 1u])),
    //             Ariadne::EffectiveFormula::binary(
    //                 Ariadne::OperatorCode::MUL,
    //                 Ariadne::EffectiveFormula::graded(
    //                     Ariadne::OperatorCode::POW, x[size], j),
    //                 x[i])));
    //   }
    // }

    // for (unsigned i = 0; i < size; ++i) {
    //   if (guess())
    //     continue;
    //   polynomial = Ariadne::EffectiveFormula::binary(Ariadne::OperatorCode::ADD,
    //                                                  polynomial, x[i]);
    // }

    return Ariadne::EffectiveScalarMultivariateFunction(
        Ariadne::EuclideanDomain(size + 1), Ariadne::simplify(polynomial));
  }

private:
  size_t size;
  int seed;
  std::vector<VType> x;
  std::vector<float> generator_bounds;
=======
=======
//----------------------------------------------------------------------------//
//                F_1(X) = sum_i sum_j rand(-100,100)*x_i^j
//----------------------------------------------------------------------------//
template <>
class FunctionDistribution<Ariadne::EffectiveScalarMultivariateFunction,
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
                           F_TEST_1>
{
    using VType = Ariadne::EffectiveFormula;
    using X = Ariadne::EffectiveScalarMultivariateFunction;

  public:
    FunctionDistribution(size_t size)
    {
        this->size = size;
        x = std::vector<VType>(size);
        unsigned i = 0;
        for (; i < size; ++i)
            x[i] = VType::coordinate(i);
    }
    X operator()(RandomPolynomialEngine &engine, int _seed = -1)
    {
        if (_seed == -1)
            seed = engine.seed;

        Ariadne::Vector<Ariadne::EffectiveFormula> monomials(size *
                                                             engine.max_order);
        Ariadne::Int expi = 2;
        std::default_random_engine generator(seed);
        std::uniform_int_distribution<int> distribution(1, 1);

        for (unsigned i = 0; i < size; ++i) {
            for (unsigned j = 0; j < engine.max_order; ++j, expi++)
            {
                monomials[engine.max_order * i + j] =
                    Ariadne::EffectiveFormula::binary(
                        Ariadne::OperatorCode::MUL,
                        Ariadne::EffectiveFormula::constant(
                            distribution(generator)),
                        Ariadne::EffectiveFormula::graded(
                            Ariadne::OperatorCode::POW, x[i]-distribution(generator), expi));
            }
            expi = 0;
        }

        Ariadne::EffectiveFormula polynomial;
        for (unsigned k = 0; k < size * engine.max_order; ++k)
            polynomial = Ariadne::EffectiveFormula::binary(
                Ariadne::OperatorCode::ADD, polynomial, monomials[k]);

        return Ariadne::EffectiveScalarMultivariateFunction(
            Ariadne::EuclideanDomain(size), Ariadne::simplify(polynomial));
    }

  private:
    size_t size;
    unsigned seed;
    std::vector<VType> x;
};
//----------------------------------------------------------------------------//
//                    G_1(X) = sum_i sum_j (x_i-x_j)^2
//----------------------------------------------------------------------------//
template <>
class FunctionDistribution<Ariadne::EffectiveVectorMultivariateFunction,
                           G_TEST_1>
{
    using VType = Ariadne::EffectiveFormula;
    using X = Ariadne::EffectiveVectorMultivariateFunction;

  public:
    FunctionDistribution(size_t nVars, size_t nCons = 1)
    {
        this->nVars = nVars;
        this->nCons = nCons;
        x = std::vector<VType>(nVars);
        unsigned i = 0;
        for (; i < nVars; ++i)
            x[i] = VType::coordinate(i);
    }
    X operator()(RandomPolynomialEngine &engine, int _seed = -1)
    {
        if (_seed == -1)
            seed = engine.seed;

        Ariadne::Vector<Ariadne::EffectiveFormula> monomials(nVars * nVars);
        for (unsigned i = 0; i < nVars; ++i) {
            for (unsigned j = i + 1; j < nVars; ++j) {
                if (i != j)
                    monomials[nVars * i + j] =
                        Ariadne::EffectiveFormula::graded(
                            Ariadne::OperatorCode::POW, x[i] - x[j], 2);
                else
                    monomials[nVars * i + j] =
                        Ariadne::EffectiveFormula::constant(0);
            }
        }

        Ariadne::EffectiveFormula polynomial;
        for (unsigned k = 0; k < nVars * nVars; ++k)
            polynomial = Ariadne::EffectiveFormula::binary(
                Ariadne::OperatorCode::ADD, polynomial, monomials[k]);
        Ariadne::Vector<Ariadne::EffectiveFormula> polynomials(nCons);
        polynomials[0] = Ariadne::simplify(polynomial);

        return Ariadne::EffectiveVectorMultivariateFunction(
            Ariadne::EuclideanDomain(nVars), polynomials);
    }

  private:
    size_t nVars;
    size_t nCons;
    unsigned seed;
    std::vector<VType> x;
};
//----------------------------------------------------------------------------//
//                    F_2(X) = sum_i sum_j x_i * x_j
//----------------------------------------------------------------------------//
template <>
class FunctionDistribution<Ariadne::EffectiveScalarMultivariateFunction,
                           F_TEST_2>
{
    using VType = Ariadne::EffectiveFormula;
    using X = Ariadne::EffectiveScalarMultivariateFunction;

  public:
    FunctionDistribution(size_t size)
    {
        this->size = size;
        x = std::vector<VType>(size);
        unsigned i = 0;
        for (; i < size; ++i)
            x[i] = VType::coordinate(i);
    }
    X operator()(RandomPolynomialEngine &engine, int _seed = -1)
    {
        if (_seed == -1)
            seed = engine.seed;

        Ariadne::Vector<Ariadne::EffectiveFormula> monomials(size * size);
        std::default_random_engine generator(seed);
        std::uniform_int_distribution<int> distribution(0, 100);

        for (unsigned i = 0; i < size; ++i) {
            for (unsigned j = 0; j < size; ++j)
                monomials[size * i + j] = Ariadne::EffectiveFormula::binary(
                    Ariadne::OperatorCode::MUL,
                    Ariadne::EffectiveFormula::constant(
                        distribution(generator)),
                    Ariadne::EffectiveFormula::binary(
                        Ariadne::OperatorCode::MUL, x[i], x[j]));
        }

        Ariadne::EffectiveFormula polynomial;
        for (unsigned k = 0; k < size * size; ++k)
            polynomial = Ariadne::EffectiveFormula::binary(
                Ariadne::OperatorCode::ADD, polynomial, monomials[k]);

        return Ariadne::EffectiveScalarMultivariateFunction(
            Ariadne::EuclideanDomain(size), Ariadne::simplify(polynomial));
    }

  private:
    size_t size;
    unsigned seed;
    std::vector<VType> x;
};
//----------------------------------------------------------------------------//
//                       G_2(X) = sum_i x_i
//----------------------------------------------------------------------------//
template <>
class FunctionDistribution<Ariadne::EffectiveVectorMultivariateFunction,
                           G_TEST_2>
{
    using VType = Ariadne::EffectiveFormula;
    using X = Ariadne::EffectiveVectorMultivariateFunction;

  public:
    FunctionDistribution(size_t nVars)
    {
        this->nVars = nVars;
        this->nCons = 2;
        x = std::vector<VType>(nVars);
        unsigned i = 0;
        for (; i < nVars; ++i)
            x[i] = VType::coordinate(i);
    }
    X operator()(RandomPolynomialEngine &engine, int _seed = -1)
    {
        if (_seed == -1)
            seed = engine.seed;

        std::default_random_engine generator(seed);
        std::uniform_real_distribution<float> distribution(0, 1);

        Ariadne::Vector<Ariadne::EffectiveFormula> budget(nVars);
        Ariadne::Vector<Ariadne::EffectiveFormula> risk(nVars);
        for (unsigned i = 0; i < nVars; ++i) {
            risk[i] = Ariadne::EffectiveFormula::binary(
                Ariadne::OperatorCode::MUL,
                Ariadne::EffectiveFormula::constant(1), x[i]);
            budget[i] = Ariadne::EffectiveFormula::binary(
                Ariadne::OperatorCode::MUL,
                Ariadne::EffectiveFormula::constant(1), x[i]);
        }

        Ariadne::EffectiveFormula budgets, risks;
        for (unsigned i = 0; i < nVars; ++i) {
            budgets = Ariadne::EffectiveFormula::binary(
                Ariadne::OperatorCode::ADD, budgets, budget[i]);
            risks = Ariadne::EffectiveFormula::binary(
                Ariadne::OperatorCode::ADD, risks, risk[i]);
        }
        Ariadne::Vector<Ariadne::EffectiveFormula> polynomials(nCons);
        polynomials[0] = Ariadne::simplify(budgets);
        polynomials[1] = Ariadne::simplify(risks);

        return Ariadne::EffectiveVectorMultivariateFunction(
            Ariadne::EuclideanDomain(nVars), polynomials);
    }

  private:
    size_t nVars;
    size_t nCons;
    unsigned seed;
    std::vector<VType> x;
};
//----------------------------------------------------------------------------//
//                  G_3(X) = sum_i sum_j rand(0,100)* x_i*x_j
//----------------------------------------------------------------------------//
template <>
class FunctionDistribution<Ariadne::EffectiveVectorMultivariateFunction,
                           G_TEST_3>
{
    using VType = Ariadne::EffectiveFormula;
    using X = Ariadne::EffectiveVectorMultivariateFunction;

  public:
    FunctionDistribution(size_t nVars)
    {
        this->nVars = nVars;
        this->nCons = 2;
        x = std::vector<VType>(nVars);
        unsigned i = 0;
        for (; i < nVars; ++i)
            x[i] = VType::coordinate(i);
    }
    X operator()(RandomPolynomialEngine &engine, int _seed = -1)
    {
      if (_seed == -1)
          seed = engine.seed;

      Ariadne::Vector<Ariadne::EffectiveFormula> monomials(nVars * nVars);
      std::default_random_engine generator(seed);
      std::uniform_int_distribution<int> distribution(0, 100);

      for (unsigned i = 0; i < nVars; ++i) {
          for (unsigned j = 0; j < nVars; ++j)
              monomials[nVars * i + j] = Ariadne::EffectiveFormula::binary(
                  Ariadne::OperatorCode::MUL,
                  Ariadne::EffectiveFormula::constant(
                      distribution(generator)),
                  Ariadne::EffectiveFormula::binary(
                      Ariadne::OperatorCode::MUL, x[i], x[j]));
      }

      Ariadne::EffectiveFormula polynomial;
      for (unsigned k = 0; k < nVars * nVars; ++k)
          polynomial = Ariadne::EffectiveFormula::binary(
              Ariadne::OperatorCode::ADD, polynomial, monomials[k]);

      Ariadne::Vector<Ariadne::EffectiveFormula> polynomials(1);
      polynomials[0] = Ariadne::simplify(polynomial);

      return Ariadne::EffectiveVectorMultivariateFunction(
          Ariadne::EuclideanDomain(nVars), polynomials);
    }

  private:
    size_t nVars;
    size_t nCons;
    unsigned seed;
    std::vector<VType> x;
};
//----------------------------------------------------------------------------//
//              F_4(X) = sum_i x_i * [rand(-40,-1)+ln(x_i/(sum_j x_i))]
//----------------------------------------------------------------------------//
template<>
class FunctionDistribution<Ariadne::EffectiveScalarMultivariateFunction,
                           F_TEST_4>
{
    using VType = Ariadne::EffectiveFormula;
    using X = Ariadne::EffectiveScalarMultivariateFunction;

  public:
    FunctionDistribution(size_t size)
    {
        this->size = size;
        x = std::vector<VType>(size);
        unsigned i = 0;
        for (; i < size; ++i)
            x[i] = VType::coordinate(i);
    }
    X operator()(RandomPolynomialEngine &engine, int _seed = -1)
    {
        if (_seed == -1)
            seed = engine.seed;

        Ariadne::EffectiveFormula                   sum;
        Ariadne::Vector<Ariadne::EffectiveFormula>  monomials(size * size);
        std::default_random_engine                  generator(seed);
        std::uniform_int_distribution<int>          distribution(-40, -1);

        for(unsigned i=0;i<size;++i)
        {
          sum = Ariadne::EffectiveFormula::binary(
            Ariadne::OperatorCode::ADD, sum, x[i]);
        }

        for (unsigned i = 0; i < size; ++i) {
                monomials[i] = Ariadne::EffectiveFormula::binary(
                    Ariadne::OperatorCode::MUL, x[i],
                    Ariadne::EffectiveFormula::binary(
                        Ariadne::OperatorCode::ADD,
                        Ariadne::EffectiveFormula::constant(distribution(generator)),
                        log(x[i]/sum)));
        }

        Ariadne::EffectiveFormula polynomial;
        for (unsigned k = 0; k < size * size; ++k)
            polynomial = Ariadne::EffectiveFormula::binary(
                Ariadne::OperatorCode::ADD, polynomial, monomials[k]);

        return Ariadne::EffectiveScalarMultivariateFunction(
            Ariadne::EuclideanDomain(size), Ariadne::simplify(polynomial));
    }

  private:
    size_t size;
    unsigned seed;
    std::vector<VType> x;
};
//----------------------------------------------------------------------------//
//              G_4(X) = sum_i rand(0,10)*x_i
//----------------------------------------------------------------------------//
template <>
class FunctionDistribution<Ariadne::EffectiveVectorMultivariateFunction,
                           G_TEST_4>
{
    using VType = Ariadne::EffectiveFormula;
    using X = Ariadne::EffectiveVectorMultivariateFunction;

  public:
    FunctionDistribution(size_t nVars)
    {
        this->nVars = nVars;
        this->nCons = nVars;
        x = std::vector<VType>(nVars);
        unsigned i = 0;
        for (; i < nVars; ++i)
            x[i] = VType::coordinate(i);
    }
    X operator()(RandomPolynomialEngine &engine, int _seed = -1)
    {
        if (_seed == -1)
            seed = engine.seed;

        std::default_random_engine                  generator(seed);
        std::uniform_real_distribution<float>       distribution(0, 10);
        Ariadne::Vector<Ariadne::EffectiveFormula>  balances(nCons);


        for (unsigned j = 0; j < nCons; ++j) {
          for(unsigned i=0; i<nVars;++i)
          {
            balances[j] = Ariadne::EffectiveFormula::binary(
              Ariadne::OperatorCode::ADD, balances[j],
              Ariadne::EffectiveFormula::binary(
                Ariadne::OperatorCode::MUL,
                Ariadne::EffectiveFormula::constant(distribution(generator)),
                x[i]
              )
            );
          }
        }

        return Ariadne::EffectiveVectorMultivariateFunction(
            Ariadne::EuclideanDomain(nVars), Ariadne::simplify(balances));
    }

  private:
    size_t nVars;
    size_t nCons;
    unsigned seed;
    std::vector<VType> x;
};

//----------------------------------------------------------------------------//
//              F_5(X) = sum_i rand(-10,10) * e^(rand(-1,1)*x_i)
//----------------------------------------------------------------------------//
template<>
class FunctionDistribution<Ariadne::EffectiveScalarMultivariateFunction,
                           F_TEST_5>
{
    using VType = Ariadne::EffectiveFormula;
    using X = Ariadne::EffectiveScalarMultivariateFunction;

  public:
    FunctionDistribution(size_t size)
    {
        this->size = size;
        x = std::vector<VType>(size);
        unsigned i = 0;
        for (; i < size; ++i)
            x[i] = VType::coordinate(i);
    }
    X operator()(RandomPolynomialEngine &engine, int _seed = -1)
    {
        if (_seed == -1)
            seed = engine.seed;

        Ariadne::EffectiveFormula                   polynomial;
        Ariadne::Vector<Ariadne::EffectiveFormula>  monomials(size * size);
        std::default_random_engine                  generator(seed);
        std::uniform_int_distribution<int>          a_i(-10, 10);
        std::uniform_int_distribution<int>          l_i(-1, 1);

        for(unsigned i=0;i<size;++i)
        {
          polynomial = Ariadne::EffectiveFormula::binary(
            Ariadne::OperatorCode::ADD, polynomial,
            Ariadne::EffectiveFormula::binary(Ariadne::OperatorCode::MUL,
              Ariadne::EffectiveFormula::constant(a_i(generator)),
              exp(Ariadne::EffectiveFormula::binary(
                Ariadne::OperatorCode::MUL,
                Ariadne::EffectiveFormula::constant(l_i(generator)),
                x[i]
              ))
            )
            );
        }

        return Ariadne::EffectiveScalarMultivariateFunction(
            Ariadne::EuclideanDomain(size), Ariadne::simplify(polynomial));
    }

  private:
    size_t size;
    unsigned seed;
    std::vector<VType> x;
<<<<<<< HEAD
>>>>>>> Small fixes. Implemented temporary __feasible__ function to test barrier method.
=======
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
};

//----------------------------------------------------------------------------//
//              F_5(X) = sum_i rand(-10,10) * e^(rand(-1,1)*x_i)
//----------------------------------------------------------------------------//

template <>
class FunctionDistribution<Ariadne::EffectiveScalarMultivariateFunction,
                           F_TEST_6>
{
    using VType = Ariadne::EffectiveFormula;
    using X = Ariadne::EffectiveScalarMultivariateFunction;

  public:
    FunctionDistribution(size_t size)
    {
        this->size = size;
        x = std::vector<VType>(size);
        unsigned i = 0;
        for (; i < size; ++i)
            x[i] = VType::coordinate(i);
    }
    X operator()(RandomPolynomialEngine &engine, int _seed = -1)
    {
        if (_seed == -1)
            seed = engine.seed;

        Ariadne::Vector<Ariadne::EffectiveFormula> monomials(size *
                                                             engine.max_order);
        Ariadne::Int expi = 0;
        std::default_random_engine generator(seed);
        std::uniform_int_distribution<int> distribution(1, 1);

        for (unsigned i = 0; i < size; ++i) {
            for (unsigned j = 0; j < engine.max_order; ++j, expi++)
            {
                monomials[engine.max_order * i + j] =
                    Ariadne::EffectiveFormula::binary(
                        Ariadne::OperatorCode::MUL,
                        Ariadne::EffectiveFormula::constant(
                            distribution(generator)),
                        Ariadne::EffectiveFormula::graded(
                            Ariadne::OperatorCode::POW, x[i], expi));
            }
            expi = 0;
        }

        Ariadne::EffectiveFormula polynomial;
        for (unsigned k = 0; k < size * engine.max_order; ++k)
            polynomial = Ariadne::EffectiveFormula::binary(
                Ariadne::OperatorCode::ADD, polynomial, monomials[k]);

        return Ariadne::EffectiveScalarMultivariateFunction(
            Ariadne::EuclideanDomain(size), Ariadne::simplify(polynomial));
    }

  private:
    size_t size;
    unsigned seed;
    std::vector<VType> x;
};
