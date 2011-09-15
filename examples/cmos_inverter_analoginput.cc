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

/// Function for plotting the orbit and reachability set
template<class SET> void plot(const char* filename, const int& xaxis, const int& yaxis, const int& numVariables, const Box& bbox, const Colour& fc, const SET& set, const int& MAX_GRID_DEPTH) {
    // Assigns local variables
    Figure fig;
    fig.set_projection(numVariables,xaxis,yaxis);
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
        Float step_x = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > xaxis) ? 1 : 0)));
        // Initiates the x position to the bounding box left bound
        Float pos_x = bbox[0].lower();
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
        Float step_y = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > yaxis) ? 1 : 0)));
        Float pos_y = bbox[1].lower();
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
    Vector<double> dp(11);

    double Id0 = 1e-6; /// Subthreshold current
    double Vth = 0.15; /// Threshold voltage
    double nVT = 0.035; /// n-Thermal voltage, n*VT
    double Cl = 1e-3; /// Load capacitance
    double Vdd = 1.0; /// Operating voltage
    double beta_n = 1.0e-2; /// Beta of the nMOS, beta_n = mu_n * Cox
    double beta_p = 0.5e-2; /// Beta of the pMOS, beta_p = mu_p * Cox
    double Sn = 1.0; /// Shape factor of the nMOS
    double Sp = 2.0; /// Shape factor of the pMOS
    double lambda = 0.01; /// Early effect constant
    double freq = 0.25; /// Frequency, f

    /// Constants
    float EVOL_TIME = 2.0/freq;   /// Evolution time
    int EVOL_TRANS = 22;            /// Evolution transitions
    float MAX_ENCL_RADIUS = 0.1;   /// Maximum enclosure radius
    float MAX_STEP_SIZE = 1e-3;     /// Maximum step size
    int VERBOSITY = 1;              /// Verbosity of the GeneralHybridEvolver




    /// Build the Hybrid System

    /// Coordinates
    RealScalarFunction t=RealScalarFunction::coordinate(3,0); // Time
    RealScalarFunction vi=RealScalarFunction::coordinate(3,1); // Input voltage, consequently gate-source voltage for the nMOS, while Vdd-vi is the gate-source voltage for the pMOS
    RealScalarFunction vo=RealScalarFunction::coordinate(3,2); // Output voltage, consequently drain-source voltage for the nMOS, while Vdd-vo is the drain-source voltage for the pMOS

    RealScalarFunction zero=RealScalarFunction::constant(3,0.0);
    RealScalarFunction one=RealScalarFunction::constant(3,1.0);

    /// Create a MonolithicHybridAutomaton object
    MonolithicHybridAutomaton inverter;

    /// Create the discrete states
    DiscreteLocation nt_pl(1);
    DiscreteLocation nt_ps(2);
    DiscreteLocation nl_pt(3);
    DiscreteLocation ns_pt(4);
    DiscreteLocation nl_ps(5);
    DiscreteLocation ns_pl(6);
    DiscreteLocation ns_ps(7);
    DiscreteLocation falling_nt(8);
    DiscreteLocation falling_pt(9);
    DiscreteLocation rising_ns(10);
    DiscreteLocation rising_ps(11);

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

    /// Function for the behavior of the system in the nMOS linear mode, pMOS subthreshold mode (Vi >= Vth, Vo <= Vi-Vth)
    /// (t' = 1; Vi' = 2*pi*f*Vdd*cos(2*pi*f*t); Vo' = -beta_n*Sn/Cl*((Vi-Vth)*Vo-Vo^2/2) + Id0/Cl*e^((-Vi-Vth+Vdd)/(nVT)) )
    RealVectorFunction nl_pt_d((one,2.0*pi*freq*Vdd*Ariadne::cos(2.0*pi*freq*t),-beta_n*Sn/Cl*((vi-Vth)*vo-vo*vo/2)+Id0/Cl*Ariadne::exp((-vi-Vth+Vdd)/nVT)));

    /// Function for the behavior of the system in the nMOS saturation mode, pMOS subthreshold mode (Vi >= Vdd-Vth, Vo >= Vi-Vth)
    /// (t' = 1; Vi' = 2*pi*f*Vdd*cos(2*pi*f*t); Vo' = -beta_n*Sn/Cl/2*(Vi-Vth)^2 * (1+lambda*Vo) + Id0/Cl*e^((-Vi-Vth+Vdd)/(nVT)) )
    RealVectorFunction ns_pt_d((one,2.0*pi*freq*Vdd*Ariadne::cos(2.0*pi*freq*t),-beta_n*Sn/Cl/2*((vi-Vth)*(vi-Vth))+Id0/Cl*Ariadne::exp((-vi-Vth+Vdd)/nVT)));

    /// Function for the behavior of the system in the nMOS subthreshold mode, pMOS linear mode (Vi <= Vth, Vo >= Vi+Vth)
    /// (t' = 1; Vi' = 2*pi*f*Vdd*cos(2*pi*f*t); Vo' = -Id0/Cl*e^((Vi-Vth)/(nVT)) + beta_p*Sp/Cl*((Vi-Vdd+Vth)*(Vo-Vdd)-(Vo-Vdd)^2/2) )
    RealVectorFunction nt_pl_d((one,2.0*pi*freq*Vdd*Ariadne::cos(2.0*pi*freq*t),-Id0/Cl*Ariadne::exp((vi-Vth)/nVT)+beta_p*Sp/Cl*((vi-Vdd+Vth)*(vo-Vdd)-(vo-Vdd)*(vo-Vdd)/2)));

    /// Function for the behavior of the system in the nMOS subthreshold mode, pMOS saturation mode (Vi <= Vth, Vo <= Vi+Vth)
    /// (t' = 1; Vi' = 2*pi*f*Vdd*cos(2*pi*f*t); Vo' = -Id0/Cl*e^((Vi-Vth)/(nVT)) + beta_p*Sp/Cl/2*(Vi-Vdd+Vth)^2 * (1-lambda*(Vo-Vdd)) )
    RealVectorFunction nt_ps_d((one,2.0*pi*freq*Vdd*Ariadne::cos(2.0*pi*freq*t),-Id0/Cl*Ariadne::exp((vi-Vth)/nVT)+beta_p*Sp/Cl/2*(vi-Vdd+Vth)*(vi-Vdd+Vth)*(1-lambda*(vo-Vdd))));

    /// Function for the behavior of the system in the nMOS linear mode, pMOS saturation mode (Vth <= Vi <= Vdd-Vth, Vo <= Vi-Vth)
    /// (t' = 1; Vi' = 2*pi*f*Vdd*cos(2*pi*f*t); Vo' = -beta_n*Sn/Cl*((Vi-Vth)*Vo-Vo^2/2) + beta_p*Sp/Cl/2*(Vi-Vdd+Vth)^2 * (1-lambda*(Vo-Vdd)) )
    RealVectorFunction nl_ps_d((one,2.0*pi*freq*Vdd*Ariadne::cos(2.0*pi*freq*t),-beta_n*Sn/Cl*((vi-Vth)*vo-vo*vo/2)+beta_p*Sp/Cl/2*(vi-Vdd+Vth)*(vi-Vdd+Vth)*(1-lambda*(vo-Vdd))));

    /// Function for the behavior of the system in the nMOS saturation mode, pMOS linear mode (Vth <= Vi <= Vdd-Vth, Vo >= Vi+Vth)
    /// (t' = 1; Vi' = 2*pi*f*Vdd*cos(2*pi*f*t); Vo' = -beta_n*Sn/Cl/2*(Vi-Vth)^2 * (1+lambda*Vo) + beta_p*Sp/Cl*((Vi-Vdd+Vth)*(Vo-Vdd)-(Vo-Vdd)^2/2) )
    RealVectorFunction ns_pl_d((one,2.0*pi*freq*Vdd*Ariadne::cos(2.0*pi*freq*t),-beta_n*Sn/Cl/2*(vi-Vth)*(vi-Vth)*(1+lambda*vo)+beta_p*Sp/Cl/2*(vi-Vdd+Vth)*(vi-Vdd+Vth)*(1-lambda*(vo-Vdd))));

    /// Function for the behavior of the system in the nMOS saturation mode, pMOS saturation mode (Vth <= Vi <= Vdd-Vth, Vi-Vth <= Vo <= Vi+Vth)
    /// (t' = 1; Vi' = 2*pi*f*Vdd*cos(2*pi*f*t); Vo' = -beta_n*Sn/Cl/2*(Vi-Vth)^2 * (1+lambda*Vo) + beta_p*Sp/Cl/2*(Vi-Vdd+Vth)^2 * (1-lambda*(Vo-Vdd)) )
    RealVectorFunction ns_ps_d((one,2.0*pi*freq*Vdd*Ariadne::cos(2.0*pi*freq*t),-beta_n*Sn/Cl/2*(vi-Vth)*(vi-Vth)*(1+lambda*vo)+beta_p*Sp/Cl/2*(vi-Vdd+Vth)*(vi-Vdd+Vth)*(1-lambda*(vo-Vdd))));

    /// Function for the behavior of the system in the "rising" or "falling" support locations
    /// (t' = 0; Vi' = 0; Vo'=0 )
    RealVectorFunction support_d((zero,zero,zero));

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
    ScalarAffineFunction non_g(Vector<Real>(3, 0.0,1.0,0.0),-dp[1]);
    /// Guard for the switch-on of the pMOS transistor (Vi <= Vdd-Vth)
    ScalarAffineFunction pon_g(Vector<Real>(3, 0.0,-1.0,0.0),dp[4]-dp[1]);
    /// Guard for the switch-off of the nMOS transistor (Vi <= Vth)
    ScalarAffineFunction noff_g(Vector<Real>(3, 0.0,-1.0,0.0),dp[1]);
    /// Guard for the switch-off of the pMOS transistor (Vi >= Vdd-Vth)
    ScalarAffineFunction poff_g(Vector<Real>(3, 0.0,1.0,0.0),-dp[4]+dp[1]);
    /// Guard for the linear region of the nMOS transistor (Vo <= Vi - Vth)
    ScalarAffineFunction nl_g(Vector<Real>(3, 0.0,1.0,-1.0),-dp[1]);
    /// Guard for the saturation region of the nMOS transistor (Vo >= Vi - Vth)
    ScalarAffineFunction ns_g(Vector<Real>(3, 0.0,-1.0,1.0),dp[1]);
    /// Guard for the linear region of the pMOS transistor (Vo >= Vi + Vth)
    ScalarAffineFunction pl_g(Vector<Real>(3, 0.0,-1.0,1.0),-dp[1]);
    /// Guard for the saturation region of the pMOS transistor (Vo <= Vi + Vth)
    ScalarAffineFunction ps_g(Vector<Real>(3, 0.0,1.0,-1.0),dp[1]);

    /// Create the transitions

    /// From subthreshold nMOS, linear pMOS
            /// To subthreshold nMOS, saturation pMOS
            inverter.new_transition(nt_pl,to_nt_ps,nt_ps,noop_r,ps_g,urgent);
            /// To saturation nMOS, linear pMOS
            inverter.new_transition(nt_pl,to_ns_pl,ns_pl,noop_r,non_g,urgent);
    /// From subthreshold nMOS, saturation pMOS
            /// To subthreshold nMOS, linear pMOS
            inverter.new_transition(nt_ps,to_nt_pl,nt_pl,noop_r,pl_g,urgent);
             /// To rising towards saturation pMOS
            inverter.new_transition(nt_ps,to_rising_ps,rising_ps,noop_r,non_g,urgent);
    /// From rising towards saturation pMOS
            /// To linear nMOS, saturation pMOS
            inverter.new_transition(rising_ps,to_nl_ps,nl_ps,noop_r,nl_g,urgent);
            /// To saturation nMOS, saturation pMOS
            inverter.new_transition(rising_ps,to_ns_ps,ns_ps,noop_r,ns_g,urgent);
    /// From falling towards subthreshold nMOS
            /// To subthreshold nMOS, linear pMOS
            inverter.new_transition(falling_nt,to_nt_pl,nt_pl,noop_r,pl_g,urgent);
            /// To subthreshold nMOS, saturation pMOS
            inverter.new_transition(falling_nt,to_nt_ps,nt_ps,noop_r,ps_g,urgent);
    /// From linear nMOS, saturation pMOS
            /// To falling towards subthreshold nMOS
            inverter.new_transition(nl_ps,to_falling_nt,falling_nt,noop_r,noff_g,urgent);
             /// To saturation nMOS, saturation pMOS
            inverter.new_transition(nl_ps,to_ns_ps,ns_ps,noop_r,ns_g,urgent);
             /// To falling towards subthreshold pMOS
            inverter.new_transition(nl_ps,to_falling_pt,falling_pt,noop_r,poff_g,urgent);
    /// From saturation nMOS, linear pMOS
            /// To falling towards subthreshold nMOS
            inverter.new_transition(ns_pl,to_falling_nt,falling_nt,noop_r,noff_g,urgent);
             /// To saturation nMOS, saturation pMOS
            inverter.new_transition(ns_pl,to_ns_ps,ns_ps,noop_r,ps_g,urgent);
             /// To falling towards subthreshold pMOS
            inverter.new_transition(ns_pl,to_falling_pt,falling_pt,noop_r,poff_g,urgent);
    /// From saturation nMOS, saturation pMOS
            /// To falling towards subthreshold nMOS
            inverter.new_transition(ns_ps,to_falling_nt,falling_nt,noop_r,noff_g,urgent);
            /// To linear nMOS, saturation pMOS
            inverter.new_transition(ns_ps,to_nl_ps,nl_ps,noop_r,nl_g,urgent);
             /// To falling towards subthreshold pMOS
            inverter.new_transition(ns_ps,to_falling_pt,falling_pt,noop_r,poff_g,urgent);
    /// From falling towards subthreshold pMOS
            /// To linear nMOS, subthreshold pMOS
            inverter.new_transition(falling_pt,to_nl_pt,nl_pt,noop_r,nl_g,urgent);
            /// To saturation nMOS, subthreshold pMOS
            inverter.new_transition(falling_pt,to_ns_pt,ns_pt,noop_r,ns_g,urgent);
    /// From rising towards saturation nMOS
            /// To saturation nMOS, linear pMOS
            inverter.new_transition(rising_ns,to_ns_pl,ns_pl,noop_r,pl_g,urgent);
            /// To saturation nMOS, saturation pMOS
            inverter.new_transition(rising_ns,to_ns_ps,ns_ps,noop_r,ps_g,urgent);
    /// From linear nMOS, subthreshold pMOS
            /// To linear nMOS, saturation pMOS
            inverter.new_transition(nl_pt,to_nl_ps,nl_ps,noop_r,pon_g,urgent);
            /// To saturation nMOS, subthreshold pMOS
            inverter.new_transition(nl_pt,to_ns_pt,ns_pt,noop_r,ns_g,urgent);
    /// From saturation nMOS, subthreshold pMOS
            /// To rising towards saturation nMOS
            inverter.new_transition(ns_pt,to_rising_ns,rising_ns,noop_r,pon_g,urgent);
            /// To linear nMOS, subthreshold pMOS
            inverter.new_transition(ns_pt,to_nl_pt,nl_pt,noop_r,nl_g,urgent);

    /// Finished building the automaton

    // cout << "Automaton = " << inverter << endl << endl;

    /// Compute the system evolution

    /// Create a GeneralHybridEvolver object
    GeneralHybridEvolver evolver;
    evolver.verbosity = VERBOSITY;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = MAX_ENCL_RADIUS;
    evolver.parameters().maximum_step_size = MAX_STEP_SIZE;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef GeneralHybridEvolver::EnclosureType HybridEnclosureType;
    typedef GeneralHybridEvolver::OrbitType OrbitType;
    typedef GeneralHybridEvolver::EnclosureListType EnclosureListType;

    Box initial_box(3, 0.0,0.0, 0.0,0.0, 0.0,0.0);
    RealSpace initial_space=inverter.continuous_state_space(nt_ps);
    HybridSet initial_set(nt_ps,RealVariableBox(initial_space,RealBox(initial_box)));

    HybridTime evolution_time(EVOL_TIME,EVOL_TRANS);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(inverter,initial_set,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.final="<<orbit.final().size()<<std::endl;

    Box graphic_box_vo(2, 0.0, EVOL_TIME, 0.0, dp[4]);
    Box graphic_box_vi(2, 0.0, EVOL_TIME, -dp[4], dp[4]);

    plot("cmos_inverter_analoginput_orbit_t_vi", 0, 1, 3, graphic_box_vi, Colour(0.0,0.5,1.0), orbit, -1);
    plot("cmos_inverter_analoginput_orbit_t_vo", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit, -1);

    plot("cmos_inverter_analoginput_orbit_t_vo_nsubplin", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit.reach()[nt_pl], -1);
    plot("cmos_inverter_analoginput_orbit_t_vo_nsubpsat", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit.reach()[nt_ps], -1);
    plot("cmos_inverter_analoginput_orbit_t_vo_nlinpsub", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit.reach()[nl_pt], -1);
    plot("cmos_inverter_analoginput_orbit_t_vo_nsatpsub", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit.reach()[ns_pt], -1);
    plot("cmos_inverter_analoginput_orbit_t_vo_nlinpsat", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit.reach()[nl_ps], -1);
    plot("cmos_inverter_analoginput_orbit_t_vo_nsatplin", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit.reach()[ns_pl], -1);
    plot("cmos_inverter_analoginput_orbit_t_vo_nsatpsat", 0, 2, 3, graphic_box_vo, Colour(0.0,0.5,1.0), orbit.reach()[ns_ps], -1);
}
