/*****************************************************************************************************
 *            cmos_inverter_analoginput.cc
 *  
 *            by Luca Geretti
 *
 * Provides the behavior of a CMOS inverter fed by a sinusoidal (thus analog) input.
 *
 *****************************************************************************************************/

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;

/// Variables:
/*
 t: absolute time;
 Vi: input voltage, consequently gate-source voltage for the nMOS, while Vdd-Vi is the gate-source voltage for the pMOS;
 Vo: output voltage, consequently drain-source voltage for the nMOS, while Vdd-Vo is the drain-source voltage for the 
 
*/

/// Function for the behavior of the system in the nMOS linear mode, pMOS subthreshold mode (Vi >= Vth, Vo <= Vi-Vth)
/// (t' = 1; Vi' = 2*pi*f*Vdd*cos(2*pi*f*t); Vo' = -beta_n*Sn/Cl*((Vi-Vth)*Vo-Vo^2/2) + Id0/Cl*e^((-Vi-Vth+Vdd)/(nVT)) )
struct nl_pt_df : FunctionData<3,3,11> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
        r[1] = p[4]*2*pi<Float>()*p[10]*Ariadne::cos(2*pi<Float>()*p[10]*x[0]);
        r[2] = -p[5]*p[7]/p[3] * ((x[1]-p[1])*x[2] - x[2]*x[2]/2) + p[0]/p[3] * Ariadne::exp((-x[1]-p[1]+p[4])/p[2]);
    }
};

/// Function for the behavior of the system in the nMOS saturation mode, pMOS subthreshold mode (Vi >= Vdd-Vth, Vo >= Vi-Vth)
/// (t' = 1; Vi' = 2*pi*f*Vdd*cos(2*pi*f*t); Vo' = -beta_n*Sn/Cl/2*(Vi-Vth)^2 * (1+lambda*Vo) + Id0/Cl*e^((-Vi-Vth+Vdd)/(nVT)) )
struct ns_pt_df : FunctionData<3,3,11> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
        r[1] = p[4]*2*pi<Float>()*p[10]*Ariadne::cos(2*pi<Float>()*p[10]*x[0]);
        r[2] = -p[5]*p[7]/p[3]/2 * (x[1]-p[1]) * (x[1]-p[1]) * (1.0 + p[9]*x[2]) + p[0]/p[3] * Ariadne::exp((-x[1]-p[1]+p[4])/p[2]);
    }
};

/// Function for the behavior of the system in the nMOS subthreshold mode, pMOS linear mode (Vi <= Vth, Vo >= Vi+Vth)
/// (t' = 1; Vi' = 2*pi*f*Vdd*cos(2*pi*f*t); Vo' = -Id0/Cl*e^((Vi-Vth)/(nVT)) + beta_p*Sp/Cl*((Vi-Vdd+Vth)*(Vo-Vdd)-(Vo-Vdd)^2/2) )
struct nt_pl_df : FunctionData<3,3,11> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
        r[1] = p[4]*2*pi<Float>()*p[10]*Ariadne::cos(2*pi<Float>()*p[10]*x[0]);
        r[2] = -p[0]/p[3] * Ariadne::exp((x[1]-p[1])/p[2]) + p[6]*p[8]/p[3] * ((x[1]-p[4]+p[1])*(x[2]-p[4]) - (x[2]-p[4])*(x[2]-p[4])/2);
    }
};

/// Function for the behavior of the system in the nMOS subthreshold mode, pMOS saturation mode (Vi <= Vth, Vo <= Vi+Vth)
/// (t' = 1; Vi' = 2*pi*f*Vdd*cos(2*pi*f*t); Vo' = -Id0/Cl*e^((Vi-Vth)/(nVT)) + beta_p*Sp/Cl/2*(Vi-Vdd+Vth)^2 * (1-lambda*(Vo-Vdd)) )
struct nt_ps_df : FunctionData<3,3,11> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
        r[1] = p[4]*2*pi<Float>()*p[10]*Ariadne::cos(2*pi<Float>()*p[10]*x[0]);
        r[2] = -p[0]/p[3] * Ariadne::exp((x[1]-p[1])/p[2]) + p[6]*p[8]/p[3]/2 * (x[1]-p[4]+p[1]) * (x[1]-p[4]+p[1]) * (1.0 - p[9]*(x[2]-p[4]));
    }
};

/// Function for the behavior of the system in the nMOS linear mode, pMOS saturation mode (Vth <= Vi <= Vdd-Vth, Vo <= Vi-Vth)
/// (t' = 1; Vi' = 2*pi*f*Vdd*cos(2*pi*f*t); Vo' = -beta_n*Sn/Cl*((Vi-Vth)*Vo-Vo^2/2) + beta_p*Sp/Cl/2*(Vi-Vdd+Vth)^2 * (1-lambda*(Vo-Vdd)) )
struct nl_ps_df : FunctionData<3,3,11> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
        r[1] = p[4]*2*pi<Float>()*p[10]*Ariadne::cos(2*pi<Float>()*p[10]*x[0]);
        r[2] = -p[5]*p[7]/p[3] * ((x[1]-p[1])*x[2] - x[2]*x[2]/2) + p[6]*p[8]/p[3]/2 * (x[1]-p[4]+p[1]) * (x[1]-p[4]+p[1]) * (1.0 - p[9]*(x[2]-p[4]));
    }
};

/// Function for the behavior of the system in the nMOS saturation mode, pMOS linear mode (Vth <= Vi <= Vdd-Vth, Vo >= Vi+Vth)
/// (t' = 1; Vi' = 2*pi*f*Vdd*cos(2*pi*f*t); Vo' = -beta_n*Sn/Cl/2*(Vi-Vth)^2 * (1+lambda*Vo) + beta_p*Sp/Cl*((Vi-Vdd+Vth)*(Vo-Vdd)-(Vo-Vdd)^2/2) )
struct ns_pl_df : FunctionData<3,3,11> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
        r[1] = p[4]*2*pi<Float>()*p[10]*Ariadne::cos(2*pi<Float>()*p[10]*x[0]);
        r[2] = -p[5]*p[7]/p[3]/2 * (x[1]-p[1]) * (x[1]-p[1]) * (1.0 + p[9]*x[2]) + p[6]*p[8]/p[3] * ((x[1]-p[4]+p[1])*(x[2]-p[4]) - (x[2]-p[4])*(x[2]-p[4])/2);
    }
};

/// Function for the behavior of the system in the nMOS saturation mode, pMOS saturation mode (Vth <= Vi <= Vdd-Vth, Vi-Vth <= Vo <= Vi+Vth)
/// (t' = 1; Vi' = 2*pi*f*Vdd*cos(2*pi*f*t); Vo' = -beta_n*Sn/Cl/2*(Vi-Vth)^2 * (1+lambda*Vo) + beta_p*Sp/Cl/2*(Vi-Vdd+Vth)^2 * (1-lambda*(Vo-Vdd)) )
struct ns_ps_df : FunctionData<3,3,11> {
    template<class R, class A, class P> static void 
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
        r[1] = p[4]*2*pi<Float>()*p[10]*Ariadne::cos(2*pi<Float>()*p[10]*x[0]);
        r[2] = -p[5]*p[7]/p[3]/2 * (x[1]-p[1]) * (x[1]-p[1]) * (1.0 + p[9]*x[2]) + p[6]*p[8]/p[3]/2 * (x[1]-p[4]+p[1]) * (x[1]-p[4]+p[1]) * (1.0 - p[9]*(x[2]-p[4]));
    }
};

/// Function for the behavior of the system in the "rising" or "falling" support locations
/// (t' = 0; Vi' = 0; Vo'=0 )
struct support_df : FunctionData<3,3,11> {
    template<class R, class A, class P> static void  
    compute(R& r, const A& x, const P& p) {
	r[0] = 0.0;
        r[1] = 0.0;
        r[2] = 0.0;
    }
};

/// Function for plotting the orbit and reachability set
template<class SET> void plot(const char* filename, const int& xaxis, const int& yaxis, const int& numVariables, const Box& bbox, const Colour& fc, const SET& set, const int& MAX_GRID_DEPTH) { 
    // Assigns local variables
    Figure fig; 
    array<uint> xy(2,xaxis,yaxis);

    fig.set_projection_map(ProjectionFunction(xy,numVariables)); 
    fig.set_bounding_box(bbox); 

    // If the grid must be shown
    if (MAX_GRID_DEPTH >= 0)
    {
	// The rectangle to be drawn
	Box rect = Box(numVariables);
	// Chooses the fill colour
        fig << fill_colour(Colour(1.0,1.0,1.0));

	// Gets the number of times each variable interval would be divided by 2
        int numDivisions = MAX_GRID_DEPTH / numVariables;
	// Gets the step in the x direction, by 1/2^(numDivisions+h), where h is 1 if the step is to be further divided by 2, 0 otherwise
	double step_x = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > xaxis) ? 1 : 0)));
	// Initiates the x position to the bounding box left bound
        double pos_x = bbox[0].lower();
        // Sets the rectangle 2-nd interval to the corresponding bounding box interval (while the >2 intervals are kept at [0,0])
	rect[yaxis] = bbox[1];
        // While between the interval
        while (pos_x < bbox[0].upper())
        {
	    rect[xaxis] = Interval(pos_x,pos_x+step_x); // Sets the rectangle x coordinate
	    pos_x += step_x; // Shifts the x position
	    fig << rect; // Appends the rectangle
        }

	// Repeats for the rectangles in the y direction
	double step_y = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > yaxis) ? 1 : 0)));  
        double pos_y = bbox[1].lower();
	rect[xaxis] = bbox[0];
        while (pos_y < bbox[1].upper())
        {
	    rect[yaxis] = Interval(pos_y,pos_y+step_y);
   	    fig << rect;
	    pos_y += step_y;
        }
    }
    // Draws and creates file
    fig.set_fill_colour(fc); 
    fig << set; 
    fig.write(filename); 
}


int main() 
{    
    /// Dynamics parameters
    Vector<Float> dp(11);

    dp[0] = 1e-6; /// Subthreshold current, Id0
    dp[1] = 0.15; /// Threshold voltage, Vth
    dp[2] = 0.035; /// n-Thermal voltage, n*VT
    dp[3] = 2e-4; /// Load capacitance, Cl
    dp[4] = 1.0; /// Operating voltage, Vdd
    dp[5] = 1.0e-2; /// Beta of the nMOS, beta_n = mu_n * Cox
    dp[6] = 0.5e-2; /// Beta of the pMOS, beta_p = mu_p * Cox
    dp[7] = 1.0; /// Shape factor of the nMOS, Sn
    dp[8] = 2.0; /// Shape factor of the pMOS, Sp
    dp[9] = 0.01; /// Early effect constant, lambda
    dp[10] = 1.0; /// Frequency, f

    /// Constants
    float EVOL_TIME = 2.0/dp[10]; /// Evolution time
    int EVOL_TRANS = 4; /// Evolution transitions
    float MAX_ENCL_RADIUS = 0.005; /// Maximum enclosure radius
    float MAX_STEP_SIZE = 0.005; /// Maximum step size

    // std::cout << "Enter Maximum number of discrete transitions:";
    // std::cin >> EVOL_TRANS;
    // std::cout << std::endl << "Maximum discrete transitions = "<< EVOL_TRANS << std::endl;
    
    /// Build the Hybrid System
  
    /// Create a HybridAutomaton object
    HybridAutomaton inverter;
  
    /// Create the discrete states
    DiscreteState nt_pl(1);
    DiscreteState nt_ps(2);
    DiscreteState nl_pt(3);
    DiscreteState ns_pt(4);
    DiscreteState nl_ps(5);
    DiscreteState ns_pl(6);
    DiscreteState ns_ps(7);
    DiscreteState falling_nt(8);
    DiscreteState falling_pt(9);
    DiscreteState rising_ns(10);
    DiscreteState rising_ps(11);

    /// Create the discrete events
    DiscreteEvent to_nt_pl(1);
    DiscreteEvent to_nt_ps(2);
    DiscreteEvent to_nl_pt(3);
    DiscreteEvent to_ns_pt(4);
    DiscreteEvent to_nl_ps(5);
    DiscreteEvent to_ns_pl(6);
    DiscreteEvent to_ns_ps(7);
    DiscreteEvent to_falling_nt(8);
    DiscreteEvent to_falling_pt(9);
    DiscreteEvent to_rising_ns(10);
    DiscreteEvent to_rising_ps(11);
 
    /// Create the dynamics
    Function<nt_pl_df> nt_pl_d(dp);
    Function<nt_ps_df> nt_ps_d(dp);
    Function<nl_pt_df> nl_pt_d(dp);
    Function<ns_pt_df> ns_pt_d(dp);
    Function<nl_ps_df> nl_ps_d(dp);
    Function<ns_pl_df> ns_pl_d(dp);    
    Function<ns_ps_df> ns_ps_d(dp);
    Function<support_df> support_d(dp);
  
    /// Build the automaton
    inverter.new_mode(nt_pl,nt_pl_d);
    inverter.new_mode(nt_ps,nt_ps_d);
    inverter.new_mode(nl_pt,nl_pt_d);
    inverter.new_mode(ns_pt,ns_pt_d);
    inverter.new_mode(nl_ps,nl_ps_d);
    inverter.new_mode(ns_pl,ns_pl_d);
    inverter.new_mode(ns_ps,ns_ps_d);
    inverter.new_mode(falling_nt,support_d);
    inverter.new_mode(falling_pt,support_d);
    inverter.new_mode(rising_ns,support_d);
    inverter.new_mode(rising_ps,support_d);

    /// Create the resets
   
    /// Reset for the transitions between modes
    IdentityFunction noop_r(3);

    /// Create the guards

    /// Guard for the switch-on of the nMOS transistor (Vi >= Vth)
    AffineFunction non_g(Matrix<Float>(1,3, 0.0,1.0,0.0),Vector<Float>(1,-dp[1]));
    /// Guard for the switch-on of the pMOS transistor (Vi <= Vdd-Vth)
    AffineFunction pon_g(Matrix<Float>(1,3, 0.0,-1.0,0.0),Vector<Float>(1,dp[4]-dp[1]));
    /// Guard for the switch-off of the nMOS transistor (Vi <= Vth)
    AffineFunction noff_g(Matrix<Float>(1,3, 0.0,-1.0,0.0),Vector<Float>(1,dp[1]));
    /// Guard for the switch-off of the pMOS transistor (Vi >= Vdd-Vth)
    AffineFunction poff_g(Matrix<Float>(1,3, 0.0,1.0,0.0),Vector<Float>(1,-dp[4]+dp[1]));
    /// Guard for the linear region of the nMOS transistor (Vo <= Vi - Vth)
    AffineFunction nl_g(Matrix<Float>(1,3, 0.0,1.0,-1.0),Vector<Float>(1,-dp[1]));
    /// Guard for the saturation region of the nMOS transistor (Vo >= Vi - Vth)
    AffineFunction ns_g(Matrix<Float>(1,3, 0.0,-1.0,1.0),Vector<Float>(1,dp[1]));
    /// Guard for the linear region of the pMOS transistor (Vo >= Vi + Vth)
    AffineFunction pl_g(Matrix<Float>(1,3, 0.0,-1.0,1.0),Vector<Float>(1,-dp[1]));
    /// Guard for the saturation region of the pMOS transistor (Vo <= Vi + Vth)
    AffineFunction ps_g(Matrix<Float>(1,3, 0.0,1.0,-1.0),Vector<Float>(1,dp[1]));

    /// Create the transitions

    /// From subthreshold nMOS, linear pMOS
	    /// To subthreshold nMOS, saturation pMOS
	    inverter.new_forced_transition(to_nt_ps,nt_pl,nt_ps,noop_r,ps_g);
	    /// To saturation nMOS, linear pMOS
	    inverter.new_forced_transition(to_ns_pl,nt_pl,ns_pl,noop_r,non_g);
    /// From subthreshold nMOS, saturation pMOS
	    /// To subthreshold nMOS, linear pMOS
	    inverter.new_forced_transition(to_nt_pl,nt_ps,nt_pl,noop_r,pl_g);
 	    /// To rising towards saturation pMOS
            inverter.new_forced_transition(to_rising_ps,nt_ps,rising_ps,noop_r,non_g);
    /// From rising towards saturation pMOS
            /// To linear nMOS, saturation pMOS
	    inverter.new_forced_transition(to_nl_ps,rising_ps,nl_ps,noop_r,nl_g);
            /// To saturation nMOS, saturation pMOS
	    inverter.new_forced_transition(to_ns_ps,rising_ps,ns_ps,noop_r,ns_g);
    /// From falling towards subthreshold nMOS
            /// To subthreshold nMOS, linear pMOS
	    inverter.new_forced_transition(to_nt_pl,falling_nt,nt_pl,noop_r,pl_g);
            /// To subthreshold nMOS, saturation pMOS
	    inverter.new_forced_transition(to_nt_ps,falling_nt,nt_ps,noop_r,ps_g);
    /// From linear nMOS, saturation pMOS
	    /// To falling towards subthreshold nMOS
	    inverter.new_forced_transition(to_falling_nt,nl_ps,falling_nt,noop_r,noff_g);
 	    /// To saturation nMOS, saturation pMOS
            inverter.new_forced_transition(to_ns_ps,nl_ps,ns_ps,noop_r,ns_g);
 	    /// To falling towards subthreshold pMOS
            inverter.new_forced_transition(to_falling_pt,nl_ps,falling_pt,noop_r,poff_g);
    /// From saturation nMOS, linear pMOS
	    /// To falling towards subthreshold nMOS
	    inverter.new_forced_transition(to_falling_nt,ns_pl,falling_nt,noop_r,noff_g);
 	    /// To saturation nMOS, saturation pMOS
            inverter.new_forced_transition(to_ns_ps,ns_pl,ns_ps,noop_r,ps_g);
 	    /// To falling towards subthreshold pMOS
            inverter.new_forced_transition(to_falling_pt,ns_pl,falling_pt,noop_r,poff_g);
    /// From saturation nMOS, saturation pMOS
	    /// To falling towards subthreshold nMOS
	    inverter.new_forced_transition(to_falling_nt,ns_ps,falling_nt,noop_r,noff_g);
	    /// To linear nMOS, saturation pMOS
	    inverter.new_forced_transition(to_nl_ps,ns_ps,nl_ps,noop_r,nl_g);
 	    /// To falling towards subthreshold pMOS
            inverter.new_forced_transition(to_falling_pt,ns_ps,falling_pt,noop_r,poff_g);
    /// From falling towards subthreshold pMOS
            /// To linear nMOS, subthreshold pMOS
	    inverter.new_forced_transition(to_nl_pt,falling_pt,nl_pt,noop_r,nl_g);
            /// To saturation nMOS, subthreshold pMOS
	    inverter.new_forced_transition(to_ns_pt,falling_pt,ns_pt,noop_r,ns_g);
    /// From rising towards saturation nMOS
            /// To saturation nMOS, linear pMOS
	    inverter.new_forced_transition(to_ns_pl,rising_ns,ns_pl,noop_r,pl_g);
            /// To saturation nMOS, saturation pMOS
	    inverter.new_forced_transition(to_ns_ps,rising_ns,ns_ps,noop_r,ps_g);
    /// From linear nMOS, subthreshold pMOS
            /// To linear nMOS, saturation pMOS
	    inverter.new_forced_transition(to_nl_ps,nl_pt,nl_ps,noop_r,pon_g);
            /// To saturation nMOS, subthreshold pMOS
	    inverter.new_forced_transition(to_ns_pt,nl_pt,ns_pt,noop_r,ns_g);
    /// From saturation nMOS, subthreshold pMOS
            /// To rising towards saturation nMOS
	    inverter.new_forced_transition(to_rising_ns,ns_pt,rising_ns,noop_r,pon_g);
            /// To linear nMOS, subthreshold pMOS
	    inverter.new_forced_transition(to_nl_pt,ns_pt,nl_pt,noop_r,nl_g);

    /// Finished building the automaton

    // cout << "Automaton = " << inverter << endl << endl;

    /// Compute the system evolution

    /// Create a HybridEvolver object
    StableHybridEvolver evolver;
    evolver.verbosity = 1;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = MAX_ENCL_RADIUS;
    evolver.parameters().maximum_step_size = MAX_STEP_SIZE;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    Box initial_box(3, 0.0,0.0, 0.0,0.0, 0.0,0.0);
    HybridEnclosureType initial_enclosure(nt_ps,initial_box);

    HybridTime evolution_time(EVOL_TIME,EVOL_TRANS); 

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(inverter,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.final="<<orbit.final()<<std::endl;

    Box graphic_box_vo(2, 0.0, EVOL_TIME, 0.0, dp[4]);
    Box graphic_box_vi(2, 0.0, EVOL_TIME, -dp[4], dp[4]);

    plot("cmos_inverter_analoginput_orbit_t_vi", 0, 1, 3, graphic_box_vi, Colour(0.0,0.5,1.0), orbit, -1);
    plot("cmos_inverter_analoginput_orbit_t_vo", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit, -1);
    if (orbit.reach().find(nt_pl)!=orbit.reach().locations_end())
        plot("cmos_inverter_analoginput_orbit_t_vo_nsubplin", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit.reach()[nt_pl], -1);
    if (orbit.reach().find(nt_ps)!=orbit.reach().locations_end())
        plot("cmos_inverter_analoginput_orbit_t_vo_nsubpsat", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit.reach()[nt_ps], -1);
    if (orbit.reach().find(nl_pt)!=orbit.reach().locations_end())
    	plot("cmos_inverter_analoginput_orbit_t_vo_nlinpsub", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit.reach()[nl_pt], -1);
    if (orbit.reach().find(ns_pt)!=orbit.reach().locations_end())
    	plot("cmos_inverter_analoginput_orbit_t_vo_nsatpsub", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit.reach()[ns_pt], -1);
    if (orbit.reach().find(nl_ps)!=orbit.reach().locations_end())
    	plot("cmos_inverter_analoginput_orbit_t_vo_nlinpsat", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit.reach()[nl_ps], -1);
    if (orbit.reach().find(ns_pl)!=orbit.reach().locations_end())
    	plot("cmos_inverter_analoginput_orbit_t_vo_nsatplin", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit.reach()[ns_pl], -1);
    if (orbit.reach().find(ns_ps)!=orbit.reach().locations_end())
    	plot("cmos_inverter_analoginput_orbit_t_vo_nsatpsat", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit.reach()[ns_ps], -1); 
}
