#include "ariadne.hpp"

namespace Ariadne {

template<class IS, class SYS, class US> struct Problem {
    IS initial_set; SYS system; US unsafe_set;
};
template<class IS, class SYS, class US> Problem<IS,SYS,US> make_problem(IS initial_set, SYS system, US unsafe_set) {
    return Problem<IS,SYS,US>{initial_set,system,unsafe_set};
}

RealVariable x("x"), y("y");

decltype(auto) make_dai_problem() {
    VectorField system({dot(x)=2*x-x*y,dot(y)=2*x*x-y});
//    RealExpressionBoundedConstraintSet initial_patch({x.in(-0.6_dec,-0.4_dec),y.in(-2.1_dec,-1.9_dec)},{});
//    BoundedConstraintSet initial_patch_set=initial_patch.euclidean_set({x,y});
    RealExpressionBoundedConstraintSet initial_points({-1<=x<=+1,-3<=y<=-1},{sqr(x)+sqr(y+2)<=1});
    BoundedConstraintSet initial_set=initial_points.euclidean_set({x,y});
    RealVariablesBox unsafe_box = {-2.5_dec<=x<=-2,-2_dec<=y<=-1.5_dec};
    return make_problem(initial_set,system,unsafe_box);
}


decltype(auto) make_fitzhugh_nagumo_problem() {
    VectorField system({dot(x)=-x*x*x/3+x-y+0.875_dec, dot(y)=0.04_dec*(x-0.8_dec*y+0.7_dec)});

    RealVariablesBox initial_box = {-1<=x<=0.5_dec,1<=y<=1.5_dec};
    RealVariablesBox unsafe_box = {-2.5_dec<=x<=-2,-2_dec<=y<=-1.5_dec};

    return make_problem(initial_box,system,unsafe_box);
}

} // namespace Ariadne
