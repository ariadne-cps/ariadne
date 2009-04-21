/***************************************************************************
 *            springs.cc
 ****************************************************************************/

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;

/// Control variables:
/// x1: position of the first mass
/// x2: position of the second mass
/// v1: speed of the first mass
/// v2: speed of the second mass

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
    draw(fig,set); 
    fig.write(filename); 
}

int main() 
{    

    /// Constants
    float m1 = 4.0; // Mass of the first ball
    float m2 = 0.75; // Mass of the second ball
    float k1 = 2.0; // Elastic constant of the first spring
    float k2 = 1.0; // Elastic constant of the second spring
    float p1 = 1.0; // Neutral position for the first spring
    float p2 = 2.0; // Neutral position for the second spring 
    float x1_0 = 0.0; // Initial position for the first spring
    float x2_0 = 3.0; // Initial position for the second spring
    float st = 1.9; // Stickyness
    float EVOL_TIME = 25.0; // Evolution time
    int   EVOL_TRANS = 4; // Evolution transitions
    float MAX_ENCLOSURE_RADIUS = 0.02; // Maximum enclosure radius
    float MAX_STEP_SIZE = 0.05; // Maximum integration step size

    /// Build the Hybrid System
  
    /// Create a HybridAutomaton object
    HybridAutomaton springs;
  
    /// Create the discrete states
    DiscreteState free(1);
    DiscreteState stuck(2);

    /// Create the discrete events
    DiscreteEvent sticking(1);
    DiscreteEvent unsticking(2);
    
    /// Create the dynamics

    /// Free oscillation (x1'=vx1; x2'=vx2; vx1'=k1*(p1-x1)/m1; vx2'=k2*(p2-x2)/m2; t'=1 )
    Matrix<Float> A1 = Matrix<Float>(5,5);
    Vector<Float> b1 = Vector<Float>(5);
    A1[0][2] = 1.0;
    A1[1][3] = 1.0;
    A1[2][0] = -k1/m1;
    b1[2] = k1*p1/m1;
    A1[3][1] = -k2/m2;
    b1[3] = k2*p2/m2;
    b1[4] = 1.0;
    AffineFunction free_d(A1,b1);
    
    /// Stuck oscillation (x1'=vx1; x2'=vx2; vx1'=vx2'=(k1*p1+k2*p2-(k1+k2)*x1)/(m1+m2); t'=1 )
    Matrix<Float> A2 = Matrix<Float>(5,5);
    Vector<Float> b2 = Vector<Float>(5);
    A2[0][2] = 1.0;
    A2[1][3] = 1.0;
    A2[2][0] = -(k1+k2)/(m1+m2);
    b2[2] = (k1*p1+k2*p2)/(m1+m2);
    A2[3][1] = -(k1+k2)/(m1+m2);
    b2[3] = (k1*p1+k2*p2)/(m1+m2);
    b2[4] = 1.0;
    AffineFunction stuck_d(A2,b2);

    /// Create the resets

    /// Stick (x1^=x1; x2^=x2; vx1^=vx2^=(m1*vx1+m2*vx2)/(m1+m2) ) 
    AffineFunction stick_r(Matrix<Float>(5,5, 1.0, 0.0, 0.0,        0.0,        0.0,
                                              0.0, 1.0, 0.0,        0.0,        0.0,
                                              0.0, 0.0, m1/(m1+m2), m2/(m1+m2), 0.0,
                                              0.0, 0.0, m1/(m1+m2), m2/(m1+m2), 0.0,
                                              0.0, 0.0, 0.0,        0.0,        1.0),
                           Vector<Float>(5));
    /// Unstick (do nothing)
    IdentityFunction unstick_r(5);

    /// Create the guards
    
    /// Guard for the transition between free and stuck (x1 >= x2)
    AffineFunction free2stuck_g(Matrix<Float>(1,5,1.0,-1.0,0.0,0.0,0.0),Vector<Float>(1));
    /// Guard for the transition between stuck and just free ((k1-k2)*x1 + k2*p2 - k1*p1 >= st)
    AffineFunction stuck2free_g(Matrix<Float>(1,5,k1-k2,0.0,0.0,0.0,0.0),Vector<Float>(1,k2*p2 -k1*p1 -st));
    //AffineFunction stuck2free_g(Matrix<Float>(1,4,0.0,0.0,-1.0,1.0),Vector<Float>(1,-st)); 


    /// Build the automaton
    
    /// Locations
    springs.new_mode(free,free_d);
    springs.new_mode(stuck,stuck_d);
    /// Invariants
    //springs.new_invariant(free,free2stuck_g);
    /// Events
    springs.new_forced_transition(sticking,free,stuck,stick_r,free2stuck_g);
    springs.new_forced_transition(unsticking,stuck,free,unstick_r,stuck2free_g);

    /// Finished building the automaton

    cout << "Automaton = " << springs << endl << endl;

    /// Compute the system evolution

    /// Create a HybridEvolver object
    StableHybridEvolver evolver;
    evolver.verbosity = 1;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = MAX_ENCLOSURE_RADIUS;
    evolver.parameters().maximum_step_size = MAX_STEP_SIZE;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the system evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    std::cout << "Computing evolution..." << std::endl;

    Box initial_box(5, x1_0,x1_0, x2_0,x2_0, 0.0,0.0, 0.0,0.0, 0.0,0.0);
//    Box initial_box(4, -0.994037,-0.994037, -0.994037,-0.994037, 0.225821,0.225821, 0.225821,0.225821);

    HybridEnclosureType initial_enclosure(free,initial_box);
    
    std::cout << "Initial set=" << initial_enclosure << std::endl;
  
    HybridTime evolution_time(EVOL_TIME,EVOL_TRANS);

    std::cout << "Computing orbit... " << std::flush;
    OrbitType orbit = evolver.orbit(springs,initial_enclosure,evolution_time,UPPER_SEMANTICS);
    std::cout << "done." << std::endl;

    std::cout << "Orbit.initial="<<orbit.initial()<<std::endl;
    std::cout << "Orbit.final="<<orbit.final()<<std::endl;

    Box graphic_box_x1x2(2,min(p1-abs(p1-x1_0),p2-abs(p2-x2_0)),max(p1+abs(p1-x1_0),p2+abs(p2-x2_0)),min(p1-abs(p1-x1_0),p2-abs(p2-x2_0)),max(p1+abs(p1-x1_0),p2+abs(p2-x2_0)));
    Box graphic_box_x1v1(2,p1-abs(p1-x1_0),p1+abs(p1-x1_0),-ceil(sqrt(k1/m1)*abs(p1-x1_0)),ceil(sqrt(k1/m1)*abs(p1-x1_0)));
    Box graphic_box_x2v2(2,p2-abs(p2-x2_0),p2+abs(p2-x2_0),-ceil(sqrt(k2/m2)*abs(p2-x2_0)),ceil(sqrt(k2/m2)*abs(p2-x2_0)));
    Box graphic_box_xv(2, -0.0,3.0, -2.0,2.0);
    Box graphic_box_xt(2, -0.0,EVOL_TIME, -0.0,3.0);
    plot("springs_x1v1_orbit", 0, 2, 5, graphic_box_xv, Colour(0.0,0.5,1.0), orbit, 2);
    plot("springs_x2v2_orbit", 1, 3, 5, graphic_box_xv, Colour(0.0,0.5,1.0), orbit, 2);
    plot("springs_x1x2_orbit", 0, 1, 5, graphic_box_x1x2, Colour(0.0,0.5,1.0), orbit, 2);
    plot("springs_x1t_orbit", 4, 0, 5, graphic_box_xt, Colour(0.0,0.5,1.0), orbit, 2);
    plot("springs_x2t_orbit", 4, 1, 5, graphic_box_xt, Colour(0.0,0.5,1.0), orbit, 2);

}
