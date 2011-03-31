/***************************************************************************
 *            sinusoid_cos_timed_infinite.cc
 ****************************************************************************/

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;

/// Function for the generation of a sinusoid
struct GenSinusoid : VectorFunctionData<2,2,0> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
	      r[0] = 1.0;
        r[1] = Ariadne::cos(x[0]);
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
    Vector<Float> system_parameters(0);

    /// Build the Hybrid System

    /// Create a HybridAutomaton object
    HybridAutomaton sinusoid;

    /// Create the discrete state
    DiscreteState l1(1);

    /// Create the discrete event
    DiscreteEvent e1(1);

    /// Create the resets
    //VectorAffineFunction reset(Matrix<Float>(2,2,0.0,0.0,0.0,1.0),Vector<Float>(2,0.0,0.0));
    //cout << "reset=" << reset << endl << endl;

    /// Create the guards.
    /// Guards are true when f(x) = Ax + b > 0
    //VectorAffineFunction guard(Matrix<Float>(1,2,1.0,0.0),Vector<Float>(1,-2*pi<Float>()));
    //cout << "guard=" << guard << endl << endl;

    /// Create the dynamics
    VectorUserFunction<GenSinusoid> sinusoid_dynamic(system_parameters);

    cout << "dynamic = " << sinusoid_dynamic << endl << endl;

    /// Build the automaton
    sinusoid.new_mode(l1,sinusoid_dynamic);
    //sinusoid.new_forced_transition(e1,l1,l1,reset,guard);

    /// Finished building the automaton

    cout << "Automaton = " << sinusoid << endl << endl;

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver;

    /// Set the evolution parameters
    evolver.settings().maximum_enclosure_cell = Vector<Float>(2,0.5);
    evolver.settings().hybrid_maximum_step_size[l1] = 0.1*pi<Float>();

    std::cout <<  evolver.settings() << std::endl;
    evolver.verbosity=1;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    std::cout << "Computing evolution starting from location l1, t = 0.0, x =0.0" << std::endl;

    Box initial_box(2, 0.0,0.0,0.0,0.0);
    HybridEnclosureType initial_enclosure(l1,initial_box);

    HybridTime evolution_time(3*pi<Float>(),4);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(sinusoid,initial_enclosure,evolution_time,LOWER_SEMANTICS);
    std::cout << "done." << std::endl << orbit << std::endl << orbit.reach() << std::endl;

    //std::cout << "Orbit="<<orbit<<std::endl;

    Box graphic_box(2, -0.5,3*pi<Float>()+0.5, -2.5,2.5);
    plot("sinusoid_sin_orbit", 0, 1, 2, graphic_box, Colour(0.0,0.5,1.0), orbit, -1);

}
