/***************************************************************************
 *            profile_taylor_calculus.cc
 *
 *  Copyright 2009  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <cassert>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <string>
using std::string;

#include <boost/timer.hpp>
#include <boost/progress.hpp>

#include "numeric.h"
#include "real.h"
#include "taylor_set.h"
#include "taylor_function.h"
#include "taylor_calculus.h"
#include "user_function.h"
#include "box.h"

namespace Ariadne {

typedef Vector<Real> RealVector;
typedef Vector<Interval> IntervalVector;

ScalarFunction f(uint as, uint j) { return ScalarFunction::variable(as,j); }
RealVector e(int rs, int i) { return RealVector::unit(rs,i); }

ScalarFunction x=f(2,0); ScalarFunction y=f(2,1);
ScalarFunction x0=f(3,0); ScalarFunction y0=f(3,1); ScalarFunction t=f(3,2);
RealVector e0=e(2,0); RealVector e1=e(2,1);


struct ProfileReset {
    ProfileReset(TaylorCalculus* c_, VectorFunctionInterface& f_, TaylorSet& s_, int n_=1) : c(c_), f(f_), s(s_), n(n_) { }
    TaylorCalculus* c; VectorFunctionInterface& f; TaylorSet s; int n; typedef TaylorSet Result;
    TaylorSet operator()() const { TaylorSet r=c->reset_step(f,s); for(int i=1; i<n; ++i) { r=c->reset_step(f,s); } return r; }
};

struct ProfileBounds {
    ProfileBounds(TaylorCalculus* c_, VectorFunctionInterface& f_, IntervalVector d_, Float h_)
        : c(c_), f(f_), d(d_), h(h_) { }
    TaylorCalculus* c; VectorFunctionInterface& f; IntervalVector d; Float h; typedef VectorTaylorFunction Result;
    Result operator()() const { c->flow_bounds(f,d,h,std::numeric_limits<double>::max()); VectorTaylorFunction r; return r; }
};

struct ProfileFlow {
    ProfileFlow(TaylorCalculus* c_, VectorFunctionInterface& f_, IntervalVector d_, Float h_, IntervalVector b_)
        : c(c_), f(f_), d(d_), h(h_), b(b_) { }
    TaylorCalculus* c; VectorFunctionInterface& f; IntervalVector d; Float h; IntervalVector b;
    typedef VectorTaylorFunction Result;
    VectorTaylorFunction operator()() const { return c->flow_model(f,d,h,b); }
};

struct ProfileCrossing {
    ProfileCrossing(TaylorCalculus* c_, ScalarFunctionInterface& g_, VectorTaylorFunction& f_, TaylorSet s_) : c(c_), g(g_), f(f_), s(s_) { }
    TaylorCalculus* c; ScalarFunctionInterface& g; VectorTaylorFunction& f; TaylorSet s; typedef TaylorModel Result;
    TaylorModel operator()() const {
        return c->crossing_time(g,f,s); }
};


template<class T>
double error(const T& f) {
    Float max_error=0.0;
    for(uint i=0; i!=f.size(); ++i) {
        max_error=std::max(max_error,f.models()[i].error());
    }
    return max_error;
}

template<class T>
unsigned int number_of_nonzeros(const T& f) {
    unsigned int nnz=0;
    for(uint i=0; i!=f.size(); ++i) {
        nnz+=f.models()[i].number_of_nonzeros();
    }
    return nnz;
}

template<> double error(const ScalarTaylorFunction& t) { return t.model().error(); }
template<> unsigned int number_of_nonzeros(const ScalarTaylorFunction& t) { return t.model().number_of_nonzeros(); }
template<> double error(const TaylorModel& t) { return t.error(); }
template<> unsigned int number_of_nonzeros(const TaylorModel& t) { return t.number_of_nonzeros(); }

template<class Test>
void profile(const char* name, const Test& test, unsigned int tries)
{

    boost::timer tm; double t=0;

    typename Test::Result res=test();
    //VectorTaylorFunction res=static_cast<VectorTaylorFunction>(test(calc,args1,args2));
    unsigned int nnz=number_of_nonzeros(res);
    double err=error(res);

    tm.restart();
    for(uint i=0; i!=tries; ++i) {
        test();
    }
    t=tm.elapsed();

    std::cout << std::fixed
              << std::setw(19) << std::left << name << std::right
              << std::setw(10) << std::setprecision(1) << 1000000*(t/tries)
              << std::setw(14) << std::setprecision(4) << std::scientific << err
              << std::setw(6) << nnz
              << std::setw(9) << std::setprecision(4) << std::fixed << t
              << "/"  << tries
              << std::endl;
}

struct ForcedVanDerPol : VectorFunctionData<3,3,3> {
    template<class R, class A, class P>
    static void compute(R& r, const A& x, const P& p) {
        r[0]=x[1];
        std::cerr<<p[0]*x[0]<<"\n";
        std::cerr<<(1.0-x[0]*x[0])<<"\n";
        std::cerr<<(1.0-x[0]*x[0])*x[1]-x[0]<<"\n";
        std::cerr<<sin(p[2]*x[2])<<"\n";
        r[1]=p[0]*(1.0-x[0]*x[0])*x[1]-x[0]+p[1]*sin(p[2]*x[2]);
        r[2]=1.0;
    }
};

int main(int argc, const char* argv[]) {

    int n=100;

    double maximum_step_size=0.125;
    double maximum_domain_extent=128;

    double a=1.5; double b=0.375;
    ScalarFunction x=ScalarFunction::variable(2,0);
    ScalarFunction y=ScalarFunction::variable(2,1);
    VectorFunction henon_map ( (a+x*x+b*y)*e0+x*e1 );
    TaylorSet henon_initial_set = Box(2, 0.875,1.125, 0.125,0.250);

    a=-0.25; b=0.75; double c=1.0;
    VectorFunction spiral_vector_field = (c+a*x-b*y)*e0+(b*x+a*y)*e1;
    IntervalVector spiral_domain = Box(2, 1.25,1.5, 0.5,0.75);
    Float spiral_step_size; IntervalVector spiral_bounding_box;
    make_lpair(spiral_step_size,spiral_bounding_box) =
        TaylorCalculus().flow_bounds(spiral_vector_field,spiral_domain,maximum_step_size,maximum_domain_extent);

    Float g=9.8;
    VectorFunction c1=y*e0;
    ScalarFunction c2=(g+0.0*x);
    VectorFunction c3=c2*e1;
    VectorFunction c4=c1+c3;
    VectorFunction c5=(g+0.0*x)*e1;
    VectorFunction c6=y*e0+(g+0.0*x)*e1;
    VectorFunction c7=y*e0-(g+0.0*x)*e1;
    //VectorFunction ball_vector_field(2,2);
    //VectorFunction ball_flow(2,3);
    VectorFunction ball_vector_field =  y*e0 - (g + 0.0*x)*e1 ;
    VectorFunction ball_flow = (ScalarFunction(x0+(0.5*y0)*t)*t)*e0 + (y0-g*t)*e1 ;
    IntervalVector ball_domain = Box(2, 1.25,1.75, 0.5,0.5);
    Float ball_step_size; IntervalVector ball_bounding_box;
    make_lpair(ball_step_size,ball_bounding_box) =
        TaylorCalculus().flow_bounds(ball_vector_field,ball_domain,maximum_step_size,maximum_domain_extent);

    VectorTaylorFunction ball_flow_model(join(ball_domain,Interval(0,5.0)),ball_flow);


    double mu=0.5; a=1.0; double omega=1.0;
    VectorUserFunction<ForcedVanDerPol> vdp_vf(Vector<Float>(3u,mu,a,omega));
    Vector<Interval> vdp_dom = Box(3, 1.25,1.5, 0.5,0.75, 0.0,0.0);


    TaylorCalculus calculus;

    std::cout<<"name                   time(us)       error  size  time(s)/tries\n";

    profile("apply-henon",ProfileReset(&calculus,henon_map,henon_initial_set),n*1000);
    profile("bounds-spiral",ProfileBounds(&calculus,spiral_vector_field,spiral_domain,maximum_step_size),n*8);
    profile("bounds-ball",ProfileBounds(&calculus,ball_vector_field,ball_domain,maximum_step_size),n*8);
    profile("flow-spiral",ProfileFlow(&calculus,spiral_vector_field,spiral_domain,spiral_step_size,spiral_bounding_box),n*1);
    profile("flow-ball",ProfileFlow(&calculus,ball_vector_field,ball_domain,ball_step_size,ball_bounding_box),n*4);
    //profile("flow-vdp",ProfileFlow(&calculus,vdp_vf,vdp_dom,h),n*1);

    return 0;
}

} // namespace Ariadne


int main(int argc, const char* argv[]) {
    Ariadne::main(argc,argv);
}
