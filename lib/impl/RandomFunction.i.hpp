template <>
class FunctionDistribution<Ariadne::EffectiveScalarMultivariateFunction,
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
};
