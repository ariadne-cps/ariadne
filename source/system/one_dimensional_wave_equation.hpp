// Parameters of the physical one-dimensional vibrating-string model.
#ifndef ARIADNE_SYSTEM_ONE_DIMENSIONAL_WAVE_EQUATION_HPP
#define ARIADNE_SYSTEM_ONE_DIMENSIONAL_WAVE_EQUATION_HPP
#include "numeric/numeric.hpp"
namespace Ariadne {
    template<class PR>
    class WaveEquation1D
        {
            public:
                WaveEquation1D(PR pr) :    length(cast_exact(ApproximateDouble(0.0)), pr),
                                        tension(cast_exact(ApproximateDouble(0.0)), pr),
                                        mass(cast_exact(ApproximateDouble(0.0)), pr),
                                        damping(cast_exact(ApproximateDouble(0.0)), pr),
                                        x0(cast_exact(ApproximateDouble(0.0)), pr),
                                        amp(cast_exact(ApproximateDouble(0.0)), pr),
                                        k(cast_exact(ApproximateDouble(0.0)), pr),
                                        omega(cast_exact(ApproximateDouble(0.0)), pr),
                                        is_triangular(false)
                {}
                Float<PR> length;
                Float<PR> tension;
                Float<PR> mass;
                Float<PR> damping;
                Float<PR> x0;
                Float<PR> amp;
                Float<PR> k;
                Float<PR> omega;
                Bool is_triangular;

        };
} // namespace Ariadne
#endif
