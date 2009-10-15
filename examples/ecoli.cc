/***************************************************************************
 *            tutorial.cc
 *
 *  Copyright  2008  Pieter Collins
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


#include "ariadne.h"
using namespace Ariadne;

// The sigmoid function 1/(1+exp(-a*(x-t)))
template<class X, class C>
X sigmoid(const X& x, const C& t, const C& a)
{
    X r=x;
    r-=t;
    r*=(-a);
    r=Ariadne::exp(r);
    r+=1;
    r=Ariadne::rec(r);
    return r;
}


// Example from:
//   Samuel Drulhe, Giancarlo Ferrari-Trecate, Hidde de Jong
//   "Reconstruction of Switching Thresholds in
//     Piecewise-Affine Models of Genetic Regulatory Networks"
//   INRIA Report No 0322, 2006
//   HSCC06_TR.pdf
struct EColi : VectorFunctionData<4,4,17> {

    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p)
    {
        typedef typename A::value_type X;
        typedef Float C;

        X& dx_CRP=r[0]; X& dx_Fis=r[1]; X& dx_GyrAB=r[2]; X& dx_rrn=r[3];
        const X& x_CRP=x[0]; const X& x_Fis=x[1]; const X& x_GyrAB=x[2]; const X& x_rrn=x[3];

        const C k0_CRP  = 0.25;
        const C k1_CRP  = 0.40;
        const C k1_Fis  = 0.60;
        const C k2_Fis  = 1.15;
        const C k_GyrAB = 0.75;
        const C k_rrn   = 1.12;

        const C theta1_CRP = 0.33;
        //const C theta2_CRP = 0.67;
        const C theta1_Fis = 0.10;
        const C theta2_Fis = 0.50;
        const C theta3_Fis = 0.75;
        //const C theta_rrn  = 0.50;
        const C theta_S    = 0.50;

        const C gamma_CRP   = 0.5;
        const C gamma_Fis   = 2.0;
        const C gamma_GyrAB = 1.0;
        const C gamma_rrn   = 1.5;

        const C x_S         = 1.0;
        const C scal = 2.0;

        X sigma1_CRP  = sigmoid(x_CRP,theta1_CRP,scal);
        //X sigma2_CRP  = sigmoid(x_CRP,theta2_CRP,scal);
        X sigma1_Fis  = sigmoid(x_Fis,theta1_Fis,scal);
        X sigma2_Fis  = sigmoid(x_Fis,theta2_Fis,scal);
        X sigma3_Fis  = sigmoid(x_Fis,theta3_Fis,scal);
        X sigma_GyrAB = sigmoid(x_Fis,theta3_Fis,scal);
        C sigma_S     = sigmoid(x_S,theta_S,scal);

        dx_CRP   = k0_CRP  +  k1_CRP * (1-sigma1_Fis) * sigma1_CRP * sigma_S  -  gamma_CRP * x_CRP;
        dx_Fis   = k1_Fis * (1 - sigma1_CRP * sigma_S)  +  k2_Fis * (1-sigma_GyrAB) * (1-sigma1_CRP) * sigma_S  -  gamma_Fis*x_Fis;
        dx_GyrAB = k_GyrAB * (1-sigma3_Fis) - gamma_GyrAB * x_GyrAB;
        dx_rrn   = k_rrn * sigma2_Fis - gamma_rrn  *x_rrn;
    }
};

int main()
{
    VectorUserFunction<EColi> ecoli_function(Vector<Interval>(17));

    HybridAutomaton ecoli_system;
    DiscreteState starvation_mode(1);
    ecoli_system.new_mode(starvation_mode,ecoli_function);

    PlanarProjectionMap projection(4,0,1);
    Box bounding_box(2, 0.0,1.0, 0.0,1.0);
    Colour turquoise(0.0,0.5,1.0);
    Colour green(0.0,1.0,0.25);

    /// Create a HybridEvolver object
    /// Set the evolution parameters
    HybridEvolver evolver;
    evolver.parameters().maximum_enclosure_radius = 0.25;
    evolver.parameters().maximum_step_size = 0.25;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the types to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;

    Point initial_state(Vector<Float>(4,0.33,0.1,0.1,0.1));
    Box initial_box=Vector<Float>(initial_state)+Vector<Interval>(4,Interval(-1,1)*0.025);
    HybridEnclosureType initial_enclosure(starvation_mode,initial_box);

    double lock_to_grid_time=5.0;
    HybridTime evolution_time(5.0,1);
    HybridTime reachability_time(5.0,1);
    HybridTime zero_time(0.0,1);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(ecoli_system,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done.\n" << std::endl;
    plot("ecoli-orbit.png", projection, bounding_box, turquoise, orbit);


    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().lock_to_grid_time=lock_to_grid_time;
    analyser.parameters().initial_grid_depth=4;
    analyser.parameters().maximum_grid_depth=4;
    std::cout<<"Discrete evolution parameters="<<analyser.parameters()<<"\n";

    HybridImageSet initial_set;
    initial_box=Vector<Float>(initial_state)+Vector<Interval>(4,Interval(-1,1)*0.025);
    initial_set[starvation_mode]=initial_box;
    HybridGridTreeSet upper_initial_set = analyser.upper_reach(ecoli_system,initial_set,zero_time);
    std::cout << "Computing reachable and evolved set... " << std::flush;
    HybridGridTreeSet upper_evolve_set, upper_reach_set;
    make_lpair(upper_reach_set,upper_evolve_set)= analyser.upper_reach_evolve(ecoli_system,initial_set,reachability_time);
    std::cout << "done.\n" << std::endl;

    plot("ecoli-initial.png", projection, bounding_box, turquoise, upper_initial_set);
    plot("ecoli-evolve.png", projection, bounding_box, turquoise, upper_evolve_set);
    plot("ecoli-reach.png", projection, bounding_box, turquoise, upper_reach_set);


}

