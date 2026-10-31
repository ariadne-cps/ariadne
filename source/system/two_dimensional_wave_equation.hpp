// Physical geometry and damping of one axis of a two-dimensional wave model.
#ifndef ARIADNE_SYSTEM_TWO_DIMENSIONAL_WAVE_EQUATION_HPP
#define ARIADNE_SYSTEM_TWO_DIMENSIONAL_WAVE_EQUATION_HPP
#include "numeric/numeric.hpp"
namespace Ariadne {
    template<class PR>
    class Parameter2D
    {
        public:
            Parameter2D(PR pr)  :   length(0, pr),
                                    damping(0, pr),
                                    x0(0, pr)
            {}
            Float<PR> length;     //Length
            Float<PR> damping;    //damping
            Float<PR> x0;
    };
} // namespace Ariadne
#endif
