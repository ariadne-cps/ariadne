/***************************************************************************
 *            rectifier.cc
 ****************************************************************************/

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;

/// Control variables:
/// t: time
/// vi: input (sinusoidal) voltage
/// vo: output (rectified) voltage

/// Non-affine dynamics

/// Dynamics for the case of both diodes being off
/// t'=1, vi'= A*cos(2*pi*f*t), vo'=-vo/(Rl*Cl)
struct offoff_df : VectorFunctionData<3,3,5> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
        r[1] = p[0]*2.0*pi<Float>()*p[1]*Ariadne::cos(2.0*pi<Float>()*p[1]*x[0]);
	    r[2] = -x[2]/(p[4]*p[3]);
//        r[2] = -x[2]/(p[4]*p[3]) - 2e-12/p[3];
    }
};

/// Dynamics for the case of the first diode being on, the second being off
/// t'=1, vi'= A*cos(2*pi*f*t), vo'=-vo/(Rl*Cl)+(vi-vo)/(Ron*Cl)
struct onoff_df : VectorFunctionData<3,3,5> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
        r[1] = p[0]*2.0*pi<Real>()*p[1]*Ariadne::cos(2.0*pi<Real>()*p[1]*x[0]);
	r[2] = -x[2]/(p[4]*p[3]) + (x[1]-x[2])/(p[2]*p[3]);
	// r[2] = -x[2]/(p[4]*p[3]) + 1e-12/p[3]*(Ariadne::exp((x[1]-x[2])/0.035)-2);
    }
};

/// Dynamics for the case of the first diode being off, the second being on
/// t'=1, vi'= A*cos(2*pi*f*t), vo'=-vo/(Rl*Cl)-(vi+vo)/(Ron*Cl)
struct offon_df : VectorFunctionData<3,3,5> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
        r[1] = p[0]*2.0*pi<Real>()*p[1]*Ariadne::cos(2.0*pi<Real>()*p[1]*x[0]);
	r[2] = -x[2]/(p[4]*p[3]) - (x[1]+x[2])/(p[2]*p[3]);
	// r[2] = -x[2]/(p[4]*p[3]) + 1e-12/p[3]*(Ariadne::exp((-x[1]-x[2])/0.035)-2);
    }
};

/// Dynamics for the case of both diodes being on
/// t'=1, vi'= A*cos(2*pi*f*t), vo'=-vo/(Rl*Cl)-2*vo/(Ron*Cl)
struct onon_df : VectorFunctionData<3,3,5> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
        std::cerr<<p[0]*2.0*pi<Real>()<<"\n";
        r[1] = p[0]*2.0*pi<Real>()*p[1]*Ariadne::cos(2.0*pi<Real>()*p[1]*x[0]);
	r[2] = -x[2]/(p[4]*p[3]) -2*x[2]/(p[2]*p[3]);
	//r[2] = -x[2]/(p[4]*p[3]) + 1e-12/p[3]*(Ariadne::exp((x[1]-x[2])/0.035) + Ariadne::exp((-x[1]-x[2])/0.035)-2);
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
    /// Introduces the dynamics parameters
    Vector<Float> dp(5);
    dp[0] = 4.0; /// Amplitude of the input voltage, Vi
    dp[1] = 50.0; /// Sinusoid frequency, f
    dp[2] = 10.0; /// Diode resistance when on, Ron
    dp[3] = 0.0001; /// Load capacitance, Cl
    dp[4] = 1000.0; /// Load resistance, Rl

    /// Introduces the global parameters
    float TIME_LIMIT = 1.0/dp[1];
//    float TIME_LIMIT = 0.0042;
    float TRAN_LIMIT = 1;
    float MAX_ENCL_RADIUS = 1.0;
    float MAX_STEP_SIZE = 1e-5/dp[1];
//    float LOCK_TOGRID_TIME = 2.0/dp[1];
    float LOCK_TOGRID_TIME = 0.25/dp[1];
    float MAX_GRID_DEPTH = 7;
    int VERBOSITY=3;
    bool ENABLE_SUBDIV=false;

    /// Build the Hybrid System

    /// Create a MonolithicHybridAutomaton object
    MonolithicHybridAutomaton rectifier;

    /// Create the discrete states
    AtomicDiscreteLocation offoff(1);
    AtomicDiscreteLocation onoff(2);
    AtomicDiscreteLocation offon(3);
    AtomicDiscreteLocation onon(4);

    /// Create the discrete events
    DiscreteEvent resettime(1);
    DiscreteEvent jump1(2), jump2(3), jump3(4);

    /// Create the resets

    /// Reset the time (t^=0,vi^=vi,vo^=vo)
    VectorAffineFunction resettime_r(Matrix<Float>(3,3,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0),
                               Vector<Float>(3,0.0,0.0,0.0));
    /// Do nothing (t^=t,vi^=vi,vo^=vo)
    IdentityFunction noop_r(3);

    /// Create the guards

    /// Guard for the reset of time (t>=1/f)
    VectorAffineFunction resettime_g(Matrix<Float>(1,3,1.0,0.0,0.0),Vector<Float>(1,-1/dp[1]));
    /// Guard for the jump from onoff to offoff (vi-vo<=0)
    VectorAffineFunction onoff_offoff_g(Matrix<Float>(1,3,0.0,-1.0,1.0),Vector<Float>(1,0.0));
    /// Guard for the jump from offon to offoff (-vi-vo<=0)
    VectorAffineFunction offon_offoff_g(Matrix<Float>(1,3,0.0,1.0,1.0),Vector<Float>(1,0.0));
    /// Guard for the jump from offoff to onoff (vi-vo>=0)
    VectorAffineFunction offoff_onoff_g(Matrix<Float>(1,3,0.0,1.0,-1.0),Vector<Float>(1,0.0));
    /// Guard for the jump from onon to onoff (-vi-vo<=0)
    VectorAffineFunction onon_onoff_g(Matrix<Float>(1,3,0.0,1.0,1.0),Vector<Float>(1,0.0));
    /// Guard for the jump from offoff to offon (-vi-vo>=0)
    VectorAffineFunction offoff_offon_g(Matrix<Float>(1,3,0.0,-1.0,-1.0),Vector<Float>(1,0.0));
    /// Guard for the jump from onon to offon (vi-vo<=0)
    VectorAffineFunction onon_offon_g(Matrix<Float>(1,3,0.0,-1.0,1.0),Vector<Float>(1,0.0));
    /// Guard for the jump from offon to onon (vi-vo>=0)
    VectorAffineFunction offon_onon_g(Matrix<Float>(1,3,0.0,1.0,-1.0),Vector<Float>(1,0.0));
    /// Guard for the jump from onoff to onon (-vi-vo>=0)
    VectorAffineFunction onoff_onon_g(Matrix<Float>(1,3,0.0,-1.0,-1.0),Vector<Float>(1,0.0));

    /// Create the dynamics
    VectorUserFunction<offoff_df> offoff_d(dp);
    VectorUserFunction<onoff_df> onoff_d(dp);
    VectorUserFunction<offon_df> offon_d(dp);
    VectorUserFunction<onon_df> onon_d(dp);

    /// Build the automaton

    /// Locations
    rectifier.new_mode(offoff,offoff_d);
    rectifier.new_mode(onoff,onoff_d);
    rectifier.new_mode(offon,offon_d);
    rectifier.new_mode(onon,onon_d);
    /// OffOff events
    rectifier.new_forced_transition(resettime,offoff,offoff,resettime_r,resettime_g);
    rectifier.new_forced_transition(jump1,offoff,onoff,noop_r,offoff_onoff_g);
    rectifier.new_forced_transition(jump2,offoff,offon,noop_r,offoff_offon_g);
    /// OnOff events
    rectifier.new_forced_transition(resettime,onoff,onoff,resettime_r,resettime_g);
    rectifier.new_forced_transition(jump1,onoff,offoff,noop_r,onoff_offoff_g);
    rectifier.new_forced_transition(jump3,onoff,onon,noop_r,onoff_onon_g);
    /// OffOn events
    rectifier.new_forced_transition(resettime,offon,offon,resettime_r,resettime_g);
    rectifier.new_forced_transition(jump1,offon,offoff,noop_r,offon_offoff_g);
    rectifier.new_forced_transition(jump3,offon,onon,noop_r,offon_onon_g);
    /// OnOn events
    rectifier.new_forced_transition(resettime,onon,onon,resettime_r,resettime_g);
    rectifier.new_forced_transition(jump2,onon,onoff,noop_r,onon_onoff_g);
    rectifier.new_forced_transition(jump3,onon,offon,noop_r,onon_offon_g);


    /// Finished building the automaton

    cout << "Automaton = " << rectifier << endl << endl;

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver;
    evolver.verbosity = 0;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = MAX_ENCL_RADIUS;
    evolver.parameters().maximum_step_size = MAX_STEP_SIZE;
    evolver.parameters().enable_subdivisions = ENABLE_SUBDIV;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    std::cout << "Computing evolution..." << std::endl;

    Box initial_box(3, 0.0, 0.0, 0.0,0.0, 0.8*dp[0],0.8*dp[0]);
    HybridEnclosureType initial_enclosure(offoff,initial_box);

//    Box initial_box(3, 0.002836,0.002836, 3.110529,3.110529, 3.110529,3.110529);
//    HybridEnclosureType initial_enclosure(onoff,initial_box);

    std::cout << "Initial set=" << initial_enclosure << std::endl;

    HybridTime evolution_time(TIME_LIMIT,TRAN_LIMIT);

/*
    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(rectifier,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.final size="<<orbit.final().size()<<std::endl;
*/
    Box graphic_box(2,0.0,1.0/dp[1],-dp[0],dp[0]);
    Box graphic_box2(2,-dp[0],dp[0],2.0,dp[0]);

    std::cout << "Plotting results..." << std::flush;
/*
    textplot("rectifier_orbit.txt", orbit);

    plot("rectifier_orbit_t_vin", 0, 1, 3, graphic_box, Colour(0.0,0.5,1.0), orbit, -1);
    plot("rectifier_orbit_t_vout", 0, 2, 3, graphic_box, Colour(0.0,0.5,1.0), orbit, -1);
    plot("rectifier_orbit_vin_vout", 1, 2, 3, graphic_box2, Colour(0.0,0.5,1.0), orbit, -1);

    std::cout << "done." << std::endl;
*/
    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().lock_to_grid_time = LOCK_TOGRID_TIME;
    analyser.parameters().maximum_grid_depth= MAX_GRID_DEPTH;
    rectifier.set_grid(Grid(Vector<Float>(3, 0.25/dp[1], 1.0, 0.5)));
    std::cout <<  analyser.parameters() << std::endl;

    analyser.verbosity=VERBOSITY;

    HybridImageSet initial_set;
    initial_set[offoff]=initial_box;
    HybridTime reach_time(TIME_LIMIT,TRAN_LIMIT);

    std::cout << "Computing upper reach set... " << std::endl << std::flush;
    HybridGridTreeSet reach = analyser.upper_reach(rectifier,initial_set,reach_time);
    std::cout << "done." << std::endl;

    plot("rectifier_reach_t_vin", 0, 1, 3, graphic_box, Colour(0.0,0.5,1.0), reach, -1);
    plot("rectifier_reach_t_vout", 0, 2, 3, graphic_box, Colour(0.0,0.5,1.0), reach, -1);
    plot("rectifier_reach_vin_vout", 1, 2, 3, graphic_box2, Colour(0.0,0.5,1.0), reach, -1);

}
