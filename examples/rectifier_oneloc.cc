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

/// Dynamics
struct running_df : VectorFunctionData<3,3,5> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
	r[0] = 1.0;
        r[1] = p[0]*2.0*pi<Float>()*p[1]*Ariadne::cos(2.0*pi<Float>()*p[1]*x[0]);
	r[2] = -x[2]/(p[4]*p[3]) + 1e-12/p[3]*(Ariadne::exp((x[1]-x[2])/0.035) + Ariadne::exp((-x[1]-x[2])/0.035)-2.0);
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
    dp[2] = 0.1; /// Diode resistance when on, Ron
    dp[3] = 0.0001; /// Load capacitance, Cl
    dp[4] = 1000.0; /// Load resistance, Rl

    /// Introduces the global parameters
    float TIME_LIMIT = 1.0/dp[1];
    float TRAN_LIMIT = 12;
    float MAX_ENCL_RADIUS = 0.001/dp[1];
    float MAX_STEP_SIZE = 0.001/dp[1];
    float LOCK_TOGRID_TIME = 1.0/dp[1];
    float MAX_GRID_DEPTH = 4;

    /// Build the Hybrid System

    /// Create a MonolithicHybridAutomaton object
    MonolithicHybridAutomaton rectifier;

    /// Create the discrete state
    AtomicDiscreteLocation running(1);

    /// Create the discrete events
    DiscreteEvent resettime(1);

    /// Create the resets

    /// Reset the time (t^=0,vi^=vi,vo^=vo)
    VectorAffineFunction resettime_r(Matrix<Float>(3,3,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0),
                               Vector<Float>(3,0.0,0.0,0.0));

    /// Create the guards

    /// Guard for the reset of time (t>=1/f)
    VectorAffineFunction resettime_g(Matrix<Float>(1,3,1.0,0.0,0.0),Vector<Float>(1,-1/dp[1]));

    /// Create the dynamics
    VectorUserFunction<running_df> running_d(dp);

    /// Build the automaton

    /// Locations
    rectifier.new_mode(running,running_d);
    /// Transitions
    rectifier.new_forced_transition(resettime,running,running,resettime_r,resettime_g);

    /// Finished building the automaton

    cout << "Automaton = " << rectifier << endl << endl;

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver;
    evolver.verbosity = 0;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = MAX_ENCL_RADIUS;
    evolver.parameters().maximum_step_size = MAX_STEP_SIZE;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    std::cout << "Computing evolution..." << std::endl;

    Box initial_box(3, 0.0, 0.0, 0.0,0.0, 0.75*dp[0],0.75*dp[0]);
    HybridEnclosureType initial_enclosure(running,initial_box);

    std::cout << "Initial set=" << initial_enclosure << std::endl;

    HybridTime evolution_time(TIME_LIMIT,TRAN_LIMIT);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(rectifier,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.initial="<<orbit.initial()<<std::endl;
    std::cout << "Orbit.final="<<orbit.final()<<std::endl;

    Box graphic_box(2,-0.1/dp[1],1.0/dp[1]*1.1,-dp[0]*1.1,dp[0]*1.1);
    Box graphic_box2(2,-dp[0],dp[0],0.0,dp[0]);

    plot("rectifier_orbit_t_vin", 0, 1, 3, graphic_box, Colour(0.0,0.5,1.0), orbit, -1);
    plot("rectifier_orbit_t_vout", 0, 2, 3, graphic_box, Colour(0.0,0.5,1.0), orbit, -1);
    plot("rectifier_orbit_vin_vout", 1, 2, 3, graphic_box2, Colour(0.0,0.5,1.0), orbit, -1);
/*
    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().lock_to_grid_time = LOCK_TOGRID_TIME;
    analyser.parameters().maximum_grid_depth= MAX_GRID_DEPTH;
    analyser.parameters().grid = Grid(Vector<Float>(5,1.0/dp[1],1.0,1.0),Vector<Float>(5, 0.0, 0.0, 0.0));
    std::cout <<  analyser.parameters() << std::endl;

    HybridImageSet initial_set;
    initial_set[running]=initial_box;
    HybridTime reach_time(TIME_LIMIT,TRAN_LIMIT);

    std::cout << "Computing upper reach set... " << std::endl << std::flush;
    HybridGridTreeSet* reach = analyser.upper_reach(rectifier,initial_set,reach_time);
    //HybridGridTreeSet* reach = analyser.chain_reach(rectifier,initial_set);
    std::cout << "done." << std::endl;

    plot("rectifier_reach_t_vin", 0, 1, 3, graphic_box, Colour(0.0,0.5,1.0), *reach, MAX_GRID_DEPTH);
    plot("rectifier_reach_t_vout", 0, 2, 3, graphic_box, Colour(0.0,0.5,1.0), *reach, MAX_GRID_DEPTH);
    plot("rectifier_reach_vin_vout", 1, 2, 3, graphic_box2, Colour(0.0,0.5,1.0), *reach, MAX_GRID_DEPTH);
*/
}
