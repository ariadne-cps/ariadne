/***************************************************************************
 *            3bouncingballs.cc
 *
 *  Copyright  2008  Davide Bresolin
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

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;

/// Function for plotting the orbit and reachability set
template<class SET> void plot(const char* filename, const int& xaxis, const int& yaxis, int& numVariables, const Box& bbox, const Colour& fc, const SET& set, const int& MAX_GRID_DEPTH) { 
    // Assigns local variables
    Figure fig; 
    array<uint> xy(2,xaxis,yaxis);
    fig.set_projection_map(ProjectionFunction(xy,numVariables)); 
    fig.set_bounding_box(bbox); 

    // If the grid must be shown
    if (MAX_GRID_DEPTH >= 0)
    {
	// The rectangle to be drawn
	Box rect = Box(4);
	// Chooses the fill colour
        fig << fill_colour(Colour(1.0,1.0,1.0));

	// Gets the number of times each variable interval would be divided by 2
        int numDivisions = MAX_GRID_DEPTH / numVariables;
	// Gets the step in the x direction, by 1/2^(numDivisions+h), where h is 1 if the step is to be further divided by 2, 0 otherwise
	double step_x = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > xaxis) ? 1 : 0)));
	// Initiates the x position to the bounding box left bound
        double pos_x = bbox[0].lower();
        // Sets the rectangle 2-nd interval to the corresponding bounding box interval (while the >2 intervals are kept at [0,0])
	rect[1] = bbox[1];
        // While between the interval
        while (pos_x < bbox[0].upper())
        {
	    rect[0] = Interval(pos_x,pos_x+step_x); // Sets the rectangle x coordinate
	    pos_x += step_x; // Shifts the x position
	    fig << rect; // Appends the rectangle
        }

	// Repeats for the rectangles in the y direction
	double step_y = 1.0/(1 << (numDivisions + ((MAX_GRID_DEPTH - numDivisions*numVariables > yaxis) ? 1 : 0)));  
        double pos_y = bbox[1].lower();
	rect[0] = bbox[0];
        while (pos_y < bbox[1].upper())
        {
	    rect[1] = Interval(pos_y,pos_y+step_y);
   	    fig << rect;
	    pos_y += step_y;
        }
    }
    // Draws and creates file
    fig.set_fill_colour(fc); 
    draw(fig,set); 
    fig.write(filename); 
}

/* Variables:
*
*	x1: x coordinate of the first ball
*	y1: y coordinate of the first ball
*	x2: x coordinate of the second ball
*	y2: y coordinate of the second ball
*	x3: x coordinate of the third ball
*	y3: y coordinate of the third ball
*	vx1: x-component of the speed of the first ball
*	vy1: y-component of the speed of the first ball
*	vx2: x-component of the speed of the second ball
*	vy2: y-component of the speed of the second ball
*	vx3: x-component of the speed of the third ball
*	vy3: y-component of the speed of the third ball
*/

int main() 
{ 
    /// Sets the system parameters
    int numVariables = 12; // The number of variables involved
    double m1 = 1.0; // Mass of the first ball
    double m2 = 1.0; // Mass of the second ball
    double m3 = 1.0; // Mass of the third ball
    double eps = 0.5; // Coefficient of restitution
    double kx = 0.5; // Fraction of x-speed recovered in a collision
    double ky = 0.5; // Fraction of y-speed recovered in a collision
    double g = 9.81; // Gravity constant
    double table_x = 4; // Extension of the table
    double table_y = 3; // Height of the table

    /// Sets the evolution parameters
    double EVOL_TIME = 4.0;
    int EVOL_TRANS = 15;
    double MAX_ENCLOSURE_RADIUS = 0.02;
    double MAX_STEP_SIZE = 0.01;
    int VERBOSITY = 1;

    /// Sets the analyzer parameters
    double LOCK_TO_GRID_TIME = 3.0;
    int MAX_GRID_DEPTH = 2;

    /// Builds the Hybrid System
  
    /// Creates a MonolithicHybridAutomaton object
    MonolithicHybridAutomaton balls;
  
    /// Creates discrete states
    DiscreteState all_on_pre12collision(1);
    DiscreteState all_on_pre23collision(2);
    DiscreteState all_on_postcollisions(3);
    DiscreteState firstsecond_on_third_off(4);
    DiscreteState first_on_secondthird_off(5);
    DiscreteState all_off(6);   

    /// Creates the discrete events
    DiscreteEvent collide12_e(1);
    DiscreteEvent collide23_e(2);
    DiscreteEvent first_rollover_e(3);
    DiscreteEvent second_rollover_e(4);
    DiscreteEvent third_rollover_e(5);
    DiscreteEvent first_bounce_e(6);
    DiscreteEvent second_bounce_e(7);
    DiscreteEvent third_bounce_e(8);
  
    /// Creates the dynamics

    /// All rolling
    Matrix<Float> A1 = Matrix<Float>(12,12);
    A1[0][6] = 1.0;
    A1[1][7] = 1.0;
    A1[2][8] = 1.0;
    A1[3][9] = 1.0;
    A1[4][10] = 1.0;
    A1[5][11] = 1.0;
    VectorAffineFunction all_rolling_d(A1,Vector<Float>(12));   

    /// First and second rolling, third bouncing
    double b2[12]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-g};    
    VectorAffineFunction firstsecond_rolling_third_bouncing_d(A1,Vector<Float>(12,b2));

    /// First rolling, second and third bouncing
    double b3[12]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-g,0.0,-g};
    VectorAffineFunction first_rolling_secondthird_bouncing_d(A1,Vector<Float>(12,b3));

    /// All bouncing
    double b4[12]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,-g,0.0,-g,0.0,-g};
    VectorAffineFunction all_bouncing_d(A1,Vector<Float>(12,b4));    

    /// Creates the resets

    /// Collision between the first and the second balls
    Matrix<Float> A2 = Matrix<Float>(12,12);
    A2[0][0] = 1.0;
    A2[1][1] = 1.0;
    A2[2][2] = 1.0;
    A2[3][3] = 1.0;
    A2[4][4] = 1.0;
    A2[5][5] = 1.0;
    A2[6][6] = (m1 - eps*m2)/(m1+m2);
    A2[6][8] = m2*(1.0+eps)/(m1+m2);	
    A2[7][7] = 1.0;
    A2[8][6] = m1*(1.0+eps)/(m1+m2);
    A2[8][8] = (m2 - eps*m1)/(m1+m2);	
    A2[9][9] = 1.0;
    A2[10][10] = 1.0;
    A2[11][11] = 1.0;
    VectorAffineFunction collision12_r(A2,Vector<Float>(12));

    /// Collision between the second and the third balls
    Matrix<Float> A3 = Matrix<Float>(12,12);
    A3[0][0] = 1.0;
    A3[1][1] = 1.0;
    A3[2][2] = 1.0;
    A3[3][3] = 1.0;
    A3[4][4] = 1.0;
    A3[5][5] = 1.0;
    A3[6][6] = 1.0;
    A3[7][7] = 1.0;
    A3[8][8] = (m2 - eps*m3)/(m2+m3);
    A3[8][10] = m3*(1.0+eps)/(m2+m3);	
    A3[9][9] = 1.0;
    A3[10][8] = m2*(1.0+eps)/(m2+m3);
    A3[10][10] = (m3 - eps*m2)/(m2+m3);	
    A3[11][11] = 1.0;
    VectorAffineFunction collision23_r(A3,Vector<Float>(12));


    /// First bouncing
    Matrix<Float> A4 = Matrix<Float>(12,12);
    A4[0][0] = 1.0;
    A4[1][1] = 1.0;
    A4[2][2] = 1.0;
    A4[3][3] = 1.0;
    A4[4][4] = 1.0;
    A4[5][5] = 1.0;
    A4[6][6] = kx;
    A4[7][7] = -ky;
    A4[8][8] = 1.0;
    A4[9][9] = 1.0;
    A4[10][10] = 1.0;
    A4[11][11] = 1.0;
    VectorAffineFunction first_bouncing_r(A4,Vector<Float>(12));   

    /// Second bouncing
    Matrix<Float> A5 = Matrix<Float>(12,12);
    A5[0][0] = 1.0;
    A5[1][1] = 1.0;
    A5[2][2] = 1.0;
    A5[3][3] = 1.0;
    A5[4][4] = 1.0;
    A5[5][5] = 1.0;
    A5[6][6] = 1.0;
    A5[7][7] = 1.0;
    A5[8][8] = kx;
    A5[9][9] = -ky;
    A5[10][10] = 1.0;
    A5[11][11] = 1.0;
    VectorAffineFunction second_bouncing_r(A5,Vector<Float>(12));   

    /// Third bouncing
    Matrix<Float> A6 = Matrix<Float>(12,12);
    A6[0][0] = 1.0;
    A6[1][1] = 1.0;
    A6[2][2] = 1.0;
    A6[3][3] = 1.0;
    A6[4][4] = 1.0;
    A6[5][5] = 1.0;
    A6[6][6] = 1.0;
    A6[7][7] = 1.0;
    A6[8][8] = 1.0;
    A6[9][9] = 1.0;
    A6[10][10] = kx;
    A6[11][11] = -ky;
    VectorAffineFunction third_bouncing_r(A6,Vector<Float>(12));   

    /// Drop from table
    IdentityFunction drop_from_table_r(12);

    /// Creates the guards

    /// Collision between the first and the second ball (x1 >= x2)
    VectorAffineFunction collide12_g(Matrix<Float>(1,12, 1.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),Vector<Float>(1,0.0));
    /// Collision between the second and the third ball (x2 >= x3)
    VectorAffineFunction collide23_g(Matrix<Float>(1,12, 0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),Vector<Float>(1,0.0));
    /// First reaching the table extension (x1 >= table_x)
    VectorAffineFunction first_reaching_table_extension_g(Matrix<Float>(1,12, 1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),Vector<Float>(1,-table_x));
    /// Second reaching the table extension (x2 >= table_x)
    VectorAffineFunction second_reaching_table_extension_g(Matrix<Float>(1,12, 0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),Vector<Float>(1,-table_x));
    /// Third reaching the table extension (x3 >= table_x)
    VectorAffineFunction third_reaching_table_extension_g(Matrix<Float>(1,12, 0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),Vector<Float>(1,-table_x));
    /// First touching the floor (y1 <= 0)
    VectorAffineFunction first_touching_floor_g(Matrix<Float>(1,12, 0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),Vector<Float>(1,0.0));
    /// Second touching the floor (y2 <= 0)
    VectorAffineFunction second_touching_floor_g(Matrix<Float>(1,12, 0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),Vector<Float>(1,0.0));
    /// Third touching the floor (y3 <= 0)
    VectorAffineFunction third_touching_floor_g(Matrix<Float>(1,12, 0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0),Vector<Float>(1,0.0));

    /// Builds the automaton
    balls.new_mode(all_on_pre12collision,all_rolling_d);
    balls.new_mode(all_on_pre23collision,all_rolling_d);
    balls.new_mode(all_on_postcollisions,all_rolling_d);
    balls.new_mode(firstsecond_on_third_off,firstsecond_rolling_third_bouncing_d);
    balls.new_mode(first_on_secondthird_off,first_rolling_secondthird_bouncing_d);
    balls.new_mode(all_off,all_bouncing_d);

    /// Creates the transitions

    /// When the first and second balls collide
    balls.new_forced_transition(collide12_e,all_on_pre12collision,all_on_pre23collision,collision12_r,collide12_g);
    /// When the second and third balls collide
    balls.new_forced_transition(collide23_e,all_on_pre23collision,all_on_postcollisions,collision23_r,collide23_g);

    /// When the third ball leaves the extension of the table, while the first and second are on the table
    balls.new_forced_transition(third_rollover_e,all_on_postcollisions,firstsecond_on_third_off,drop_from_table_r,third_reaching_table_extension_g);
    /// When the second ball leaves the extension of the table, the third is already off the table, and the first is still on the table
    balls.new_forced_transition(second_rollover_e,firstsecond_on_third_off,first_on_secondthird_off,drop_from_table_r,second_reaching_table_extension_g);
    /// When the first ball leaves the extension of the table, while the second and third are already off the table
    balls.new_forced_transition(third_rollover_e,first_on_secondthird_off,all_off,drop_from_table_r,first_reaching_table_extension_g);

    /// When the third ball touches the floor, while the first and second are still on the table
    balls.new_forced_transition(third_bounce_e,firstsecond_on_third_off,firstsecond_on_third_off,third_bouncing_r,third_touching_floor_g);
    /// When the third ball touches the floor, while the first is on the table and the second is off the table
    balls.new_forced_transition(third_bounce_e,first_on_secondthird_off,first_on_secondthird_off,third_bouncing_r,third_touching_floor_g);
    /// When the third ball touches the floor, while the first and the second are off the table
    balls.new_forced_transition(third_bounce_e,all_off,all_off,third_bouncing_r,third_touching_floor_g);

    /// When the second ball touches the floor, while the first is on the table and the third is off the table
    balls.new_forced_transition(second_bounce_e,first_on_secondthird_off,first_on_secondthird_off,second_bouncing_r,second_touching_floor_g);
    /// When the second ball touches the floor, while the first and the third are off the table
    balls.new_forced_transition(second_bounce_e,all_off,all_off,second_bouncing_r,second_touching_floor_g);

    /// When the first ball touches the floor and all balls are off the table
    balls.new_forced_transition(first_bounce_e,all_off,all_off,first_bouncing_r,first_touching_floor_g);

    /// Initial parameters

    /// Sets the initial state parameters
//    Box initial_box(12, 2.0,2.0, table_y,table_y, 4.0,4.0,table_y,table_y, 6.0,6.0,table_y,table_y, 4.0,4.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0);
//    HybridEvolver::EnclosureType initial_enclosure(all_on_pre12collision,initial_box);
    Box initial_box(12, 1.0,1.0, table_y,table_y, 2.0,2.0, table_y,table_y, 3.0,3.0,table_y,table_y, 9.0,9.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0);
    HybridEvolver::EnclosureType initial_enclosure(all_on_pre12collision,initial_box);


    /// Shows the automaton

    cout << "Automaton = " << balls << endl << endl;

    /// Plot-related information

    /// Bounding box
    Box bounding_box_pos(2, 0.0, 12.0, 0.0, 4.0);
    Box bounding_box_speed(2, 0.0, 12.0, 0.0, 10.0);

    /// Computes the system evolution

    /// Creates a HybridEvolver object
    HybridEvolver evolver;
    evolver.verbosity = VERBOSITY;

    /// Sets the evolution parameters
    evolver.parameters().maximum_enclosure_radius = MAX_ENCLOSURE_RADIUS;
    evolver.parameters().maximum_step_size = MAX_STEP_SIZE;
    std::cout <<  evolver.parameters() << std::endl;
  
    HybridTime evol_limits(EVOL_TIME,EVOL_TRANS);
 
    std::cout << "Computing orbit... " << std::flush;
    HybridEvolver::OrbitType orbit = evolver.orbit(balls,initial_enclosure,evol_limits,UPPER_SEMANTICS);

    std::cout << std::endl << "Orbit.final.size()="<<orbit.final().size()<<std::endl;

    std::cout << "Plotting result to text file..." << std::flush;
    
//    textplot("3balls_orbit.txt", orbit);
    
    std::cout << " done." << std::endl;

    std::cout << "Plotting result to png files..." << std::flush;
    
    plot("3balls-x1y1_orbit", 0, 1, numVariables, bounding_box_pos, Colour(0.0,0.5,1.0), orbit, -1);
    plot("3balls-x2y2_orbit", 2, 3, numVariables, bounding_box_pos, Colour(0.0,0.5,1.0), orbit, -1);
    plot("3balls-x3y3_orbit", 4, 5, numVariables, bounding_box_pos, Colour(0.0,0.5,1.0), orbit, -1);
    plot("3balls-x1vx1_orbit", 0, 6, numVariables, bounding_box_speed, Colour(0.0,0.5,1.0), orbit, -1);
    plot("3balls-x2vx2_orbit", 2, 8, numVariables, bounding_box_speed, Colour(0.0,0.5,1.0), orbit, -1);
    plot("3balls-x3vx3_orbit", 4, 10, numVariables, bounding_box_speed, Colour(0.0,0.5,1.0), orbit, -1);

    std::cout << " done." << std::endl;

/*
    /// Creates a ReachabilityAnalyser object
    HybridReachabilityAnalyser analyser(evolver);
    analyser.verbosity = 1;
    analyser.parameters().lock_to_grid_time = LOCK_TO_GRID_TIME;
    analyser.parameters().maximum_grid_depth= MAX_GRID_DEPTH;

    HybridImageSet initial_set;
    initial_set[both_on]=initial_box;

    //plot("2balls-initial_set1", xaxis, yaxis, numVariables, bounding_box, Colour(0.0,0.5,1.0), initial_set, MAX_GRID_DEPTH);


    // Computes evolved sets and reach sets using upper semantics.
    // These functions compute over-approximations to the evolved and reachable sets. Subdivision is used
    // as necessary to keep the local errors reasonable. The accumulated global error may be very large.
    std::cout << "Computing upper reach set... " << std::flush;
    HybridGridTreeSet* upper_reach_set_ptr = analyser.upper_reach(balls,initial_set,limits);

    /// Plots the reachability results
    plot("2balls-upper_reach", xaxis, yaxis, numVariables, bounding_box, Colour(0.0,0.5,1.0), *upper_reach_set_ptr, MAX_GRID_DEPTH);
*/
}
