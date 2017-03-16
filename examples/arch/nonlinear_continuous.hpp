#include "ariadne.hpp"

namespace Ariadne {

template<class IS, class SYS, class SS> struct Problem {
    IS initial_set; SYS system; SS safe_set;
};
template<class IS, class SYS, class SS> Problem<IS,SYS,SS> make_problem(IS initial_set, SYS system, SS safe_set) {
    return Problem<IS,SYS,SS>{initial_set,system,safe_set};
}

RealVariable x("x"), y("y");

decltype(auto) make_dai_problem() {
    VectorField system({dot(x)=2*x-x*y,dot(y)=2*x*x-y});
    RealExpressionBoundedConstraintSet initial_set({-1<=x<=+1,-3<=y<=-1},{sqr(x)+sqr(y+2)<=1});
    RealExpressionBoundedConstraintSet safe_set = {{-5<=x<=+5,-4<=y<=+6},{sqr(x)+sqr(y-1)>=0.09_dec}};
    return make_problem(initial_set,system,safe_set);
}


decltype(auto) make_fitzhugh_nagumo_problem() {
    VectorField system({dot(x)=-x*x*x/3+x-y+0.875_dec, dot(y)=0.04_dec*(x-0.8_dec*y+0.7_dec)});

    RealVariablesBox initial_box = {-1<=x<=0.5_dec,1<=y<=1.5_dec};
    RealVariablesBox unsafe_box = {-2.5_dec<=x<=-2,-2_dec<=y<=-1.5_dec};

    return make_problem(initial_box,system,unsafe_box);
}

} // namespace Ariadne
