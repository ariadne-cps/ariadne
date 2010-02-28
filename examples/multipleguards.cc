/***************************************************************************
 *            multipleguards.cc
 *
 *	      Author: Luca Geretti
 *
 * Example for multiple guards: an X-Y clockwise rotation in which a transition
 * occurs as soon as a new quadrant is reached.
 *
 ****************************************************************************/

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;

struct mg_gg : VectorFunctionData<1,3,2> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
        r[0] = min(x[0]-p[0],x[1]-p[1]);
    }
};

struct mg_ll : VectorFunctionData<1,3,2> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
        r[0] = min(p[0]-x[0],p[1]-x[1]);
    }
};

struct mg_gl : VectorFunctionData<1,3,2> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
        r[0] = min(x[0]-p[0],p[1]-x[1]);
    }
};

struct mg_lg : VectorFunctionData<1,3,2> {
    template<class R, class A, class P> static void
    compute(R& r, const A& x, const P& p) {
        r[0] = min(p[0]-x[0],x[1]-p[1]);
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

int main(int argc, char **argv) 
{
  
    float f = 1.0; // Frequency

    /// Set the system parameters
    double w1 = 2*Ariadne::pi<Float>()*f;
    double w2 = 2*Ariadne::pi<Float>()*f;

    /// Guards parameters
    Vector<Float> dp(2);
    
    double A[9]={0.0,w1,0.0, -w2,0.0,0.0, 0.0,0.0,0.0};
    double b[3]={0.0,0.0,1.0};

    float EVOL_TIME = 1.0/f;
    if(argc > 1) {
        EVOL_TIME = atof(argv[1]);
    }
    int VERBOSITY = 1;
    int EVOL_TRANS = 4;
    if(argc > 2) {
        VERBOSITY = atoi(argv[2]);
    }
   
    float MAX_ENCL_WIDTH = 1e-1;
    float MAX_STEP_SIZE = 1e-2;

    /// Build the Hybrid System
  
    /// Create a HybridAutomton object
    HybridAutomaton multipleguards;
  
    /// Create the discrete states
    DiscreteState pospos(1);
    DiscreteState posneg(2);
    DiscreteState negpos(3);
    DiscreteState negneg(4);	

    /// Create the discrete events
    DiscreteEvent pospos2posneg(12);
    DiscreteEvent pospos2negpos(13);
    DiscreteEvent posneg2pospos(21);
    DiscreteEvent posneg2negneg(24);
    DiscreteEvent negpos2pospos(31);
    DiscreteEvent negpos2negneg(34);
    DiscreteEvent negneg2posneg(42);
    DiscreteEvent negneg2negpos(43);
  
    /// Create the dynamics
    VectorAffineFunction dynamic(Matrix<Float>(3,3,A),Vector<Float>(3,b));	
   
    /// Create the guards
    dp(0) = 0.0;
    dp(1) = 0.0;
    VectorUserFunction<mg_gg> pospos_g(dp); 
    VectorUserFunction<mg_gl> posneg_g(dp); 
    VectorUserFunction<mg_lg> negpos_g(dp); 
    VectorUserFunction<mg_ll> negneg_g(dp); 
 
    /// Build the automaton
    multipleguards.new_mode(pospos,dynamic);
    multipleguards.new_mode(posneg,dynamic);
    multipleguards.new_mode(negpos,dynamic);
    multipleguards.new_mode(negneg,dynamic);

    /// Automaton transitions in the case of multiple guards 

    multipleguards.new_forced_transition(pospos2posneg,pospos,posneg,IdentityFunction(3),posneg_g);
    multipleguards.new_forced_transition(pospos2negpos,pospos,negpos,IdentityFunction(3),negpos_g);
    multipleguards.new_forced_transition(negpos2pospos,negpos,pospos,IdentityFunction(3),pospos_g);
    multipleguards.new_forced_transition(negpos2negneg,negpos,negneg,IdentityFunction(3),negneg_g);
    multipleguards.new_forced_transition(posneg2pospos,posneg,pospos,IdentityFunction(3),pospos_g);
    multipleguards.new_forced_transition(posneg2negneg,posneg,negneg,IdentityFunction(3),negneg_g);
    multipleguards.new_forced_transition(negneg2negpos,negneg,negpos,IdentityFunction(3),negpos_g);
    multipleguards.new_forced_transition(negneg2posneg,negneg,posneg,IdentityFunction(3),posneg_g);

/*
    /// Automaton transitions in the case of single guards (transition as soon as the clock reaches period/4)

    VectorAffineFunction quarterperiod_g(Matrix<Float>(1,3,0.0,0.0,1.0),Vector<Float>(1,-0.25/f));
    VectorAffineFunction resettime_r(Matrix<Float>(3,3, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,0.0),Vector<Float>(3,0.0,0.0,0.0));

    multipleguards.new_forced_transition(pospos2posneg,pospos,posneg,resettime_r,quarterperiod_g);
    multipleguards.new_forced_transition(posneg2negneg,posneg,negneg,resettime_r,quarterperiod_g);
    multipleguards.new_forced_transition(negneg2negpos,negneg,negpos,resettime_r,quarterperiod_g);
    multipleguards.new_forced_transition(negpos2pospos,negpos,pospos,resettime_r,quarterperiod_g);
*/

    /// Finished building the automaton

    cout << "Automaton = " << multipleguards << endl << endl;

    /// Compute the system evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver;
    evolver.verbosity = VERBOSITY;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_cell = Vector<Float>(3,MAX_ENCL_WIDTH);
    evolver.parameters().maximum_step_size = MAX_STEP_SIZE;
    evolver.parameters().enable_subdivisions = true;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    Box initial_box(3, 0.0,0.0, 1.0,1.0, 0.0,0.0);
    HybridEnclosureType initial_enclosure(pospos,initial_box);
    Box bounding_box1(2, -1.0,1.0, -1.0,1.0);
    Box bounding_box2(2, 0.0,EVOL_TIME,-1.0,1.0);

    HybridTime evolution_time(EVOL_TIME,EVOL_TRANS); 

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(multipleguards,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.final="<<orbit.final()<<std::endl;

    plot("multipleguards_xy", 0, 1, 3, bounding_box1, Colour(0.0,0.5,1.0), orbit, -1);
    plot("multipleguards_tx", 2, 0, 3, bounding_box2, Colour(0.0,0.5,1.0), orbit, -1);
    plot("multipleguards_ty", 2, 1, 3, bounding_box2, Colour(0.0,0.5,1.0), orbit, -1);

    if (orbit.reach().find(pospos)!=orbit.reach().locations_end())
        plot("multipleguards_xy_pospos", 0, 1, 3, bounding_box1, Colour(0.0,0.5,1.0), orbit.reach()[pospos], -1);
    if (orbit.reach().find(posneg)!=orbit.reach().locations_end())
        plot("multipleguards_xy_posneg", 0, 1, 3, bounding_box1, Colour(0.0,0.5,1.0), orbit.reach()[posneg], -1);
    if (orbit.reach().find(negpos)!=orbit.reach().locations_end())
        plot("multipleguards_xy_negpos", 0, 1, 3, bounding_box1, Colour(0.0,0.5,1.0), orbit.reach()[negpos], -1);
    if (orbit.reach().find(negneg)!=orbit.reach().locations_end())
        plot("multipleguards_xy_negneg", 0, 1, 3, bounding_box1, Colour(0.0,0.5,1.0), orbit.reach()[negneg], -1);
}
