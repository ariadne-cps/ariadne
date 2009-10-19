/*******************************************************************************************************
 *            voltagetransition.cc
 *  
 *            by Luca Geretti
 *
 * Provides the verification of a transition between a Vhigh supply voltage and a Vlow supply voltage
 *
 *******************************************************************************************************/

#include <cstdlib>
#include <cstring>
#include <cstdarg>

#include "coconut-demo.h"

/// Function for plotting the orbit and reachability set
template<class SET> void plot(const char* filename, const int& xaxis, const int& yaxis, const int& numVariables, const Box& bbox, const Colour& fc, const SET& set, const int& MAX_GRID_DEPTH) { 
    // Assigns local variables
    Figure fig; 
    array<uint> xy(2,xaxis,yaxis);

    fig.set_projection_map(ProjectionFunction(xy,numVariables)); 
    fig.set_bounding_box(bbox); 
    
    fig.set_x_axis_label("x label");
    fig.set_y_axis_label("y label");
    
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

int print_usage(string cmdname) {
    std::cerr << "Usage: "<<cmdname<<" REACH grid_depth"<<std::endl;
    std::cerr << "       "<<cmdname<<" SYNTH grid_depth parameter min max"<<std::endl;
    return 1;
}

int main(int argc, char** argv) 
{   
    string cmdname(argv[0]);
    cmdname = cmdname.substr( cmdname.find_last_of("/\\") +1 );

    // Check arguments
    if(argc < 3) return print_usage(cmdname);
    
    enum { C_REACH, C_SYNTH } command;
    if(strcmp(argv[1], "REACH") == 0) {
        command = C_REACH;
    } else if(strcmp(argv[1], "SYNTH") == 0) {
        command = C_SYNTH;
    } else {
        return print_usage(cmdname);
    }
    
	MAX_GRID_DEPTH = atoi(argv[2]);
	if(MAX_GRID_DEPTH <= 0) return print_usage(cmdname);
	
	int pvar = 2;
	float pmin = -1.0;  
	float pmax = 1.0; 
	
	if(command == C_SYNTH) {
	    if(argc < 6) return print_usage(cmdname);
	    pvar = atoi(argv[3]);
	    if(pvar < 0) return print_usage(cmdname);
	    pmin = atof(argv[4]);
	    pmax = atof(argv[5]);
	}
	
    std::cout << "command : " << command << std::endl;
    std::cout << "grid depth : " << MAX_GRID_DEPTH << std::endl;
    std::cout << "pvar : " << pvar << std::endl;
    std::cout << "pmin : " << pmin << std::endl;
    std::cout << "pmax : " << pmax << std::endl;
    
    // Build the automaton and initialize evolution parameters
    build_automaton();

    /// Compute the automaton evolution

    /// Create a HybridEvolver object
    HybridEvolver evolver;
    evolver.verbosity = 1;

    /// Set the evolution parameters
    evolver.parameters().maximum_enclosure_radius = MAX_ENCL_RADIUS;
    evolver.parameters().maximum_step_size = MAX_STEP_SIZE;
    std::cout <<  evolver.parameters() << std::endl;

    // Declare the type to be used for the automaton evolution
    typedef HybridEvolver::EnclosureType HybridEnclosureType;
    typedef HybridEvolver::OrbitType OrbitType;
    typedef HybridEvolver::EnclosureListType EnclosureListType;

    /// Create a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.parameters().lock_to_grid_time = LOCK_TOGRID_TIME;
    analyser.parameters().lock_to_grid_steps = LOCK_TOGRID_STEPS;
    analyser.parameters().maximum_grid_depth= MAX_GRID_DEPTH;
    analyser.verbosity = 3;
    automaton.set_grid(grid); 
    std::cout <<  analyser.parameters() << std::endl;

    HybridImageSet initial_set;
    if(command == C_SYNTH) initial_box[pvar] = Interval(pmin,pmax);
    initial_set[initial_state]=initial_box;

    std::cout << "Computing reach set... " << std::endl << std::flush;
    HybridGridTreeSet reach = analyser.chain_reach(automaton,initial_set);	
    std::cout << "done." << std::endl;

    plot("automaton_reach", 1, 2, 4, graphic_box, Colour(0.0,0.5,1.0), reach, -1);

    if(command == C_SYNTH) {
        if (reach[unsafe].empty())
	        cout << std::endl << "The automaton is safe for any provided value of the parameter." << std::endl;
        else if (reach[safe].empty())
        	cout << std::endl << "The automaton is unsafe for any provided value of the parameter." << std::endl;
        else {
            Float lower_unsafe = (reach[unsafe].bounding_box())[pvar].lower();
            Float lower_safe = (reach[safe].bounding_box())[pvar].lower();
            if(lower_safe < lower_unsafe) {
                cout << std::endl << "The maximum value allowed for the parameter is " << min(lower_unsafe,(reach[safe].bounding_box())[pvar].upper()) << "." << std::endl;
            } else {
                cout << std::endl << "The minimum value allowed for the parameter is " << max(lower_safe,(reach[unsafe].bounding_box())[pvar].upper()) << "." << std::endl;
            }              
        }
    }

    return 0;
}
