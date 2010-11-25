/***************************************************************************
 *            automaton.cc
 *
 *  Copyright 2009  Davide Bresolin
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

#include "system.h"

void build_automaton() {
    // Define evolution parameters
    MAX_ENCL_RADIUS = 0.2; /// Maximum enclosure radius
    MAX_STEP_SIZE = 0.5e-2; /// Maximum step size
    LOCK_TOGRID_TIME = 100.0;
    LOCK_TOGRID_STEPS = 100;

    // Parameters and constants
    float k = -0.1; /// Steepness of the reference voltage, k
    float Vhigh = 1.2; /// Operating voltage, Vhigh
    float Vlow = 0.8; /// Operating low voltage, Vlow

    /// Constants
    int NH = 24; // The number of High PMOS transistors

    float Vr0 = Vhigh;  // Initial value for the reference voltage
    float Vo0 = 1.164;  // Initial value for the output voltage
    float Vrf = 0.0;    // Final value for the reference voltage (initially 0)
    float P0 = 0.57;    // Default value for the period

    // Gets the values for S 
    Vector<float> S(NH+1);
    S[NH] = 1494.128;
    for (int n=NH-1;n>=1;n--)
		S[n] = S[n+1]/1.18;
    // Gets the values of the output voltages
    Vector<float> Vo(NH+1);
    Vo[0] = 0.0;
    for (int n=NH;n>=1;n--)
    {
    	//Vo[n] = (beta*S[n]*Rl*Vth-1)/(beta*S[n]*Rl) + std::sqrt(std::pow((beta*S[n]*Rl*Vth-1)/(beta*S[n]*Rl),2) + Vhigh*(Vhigh-2*Vth));	    
		Vo[n] = Vo0-0.05*(NH-n);
		if (n < NH-1 && Vo[n+2] > Vlow && Vo[n+1] <= Vlow)
			Vrf = (Vo[n+1] + Vo[n])/2;	
	}

	int DEPTH_ADD[4] = {3,0,0,3};

    /// Create a MonolithicHybridAutomaton object

    /// Create the events
    DiscreteEvent chk(1); // Checks the current voltage
    DiscreteEvent up(2); // Reacts to a positive voltage difference
    DiscreteEvent down(3); // Reacts to a negative voltage difference
    DiscreteEvent to_end(4); // Ends the Vr transition
    DiscreteEvent to_safe(5); // Switches to the safe location
    DiscreteEvent to_unsafe(6); // Switches to the unsafe location

    /// Create the discrete states
    DiscreteLocation states[NH+1][2];
    for (int i=0;i<=NH;i++)
    {
        states[i][0] = DiscreteLocation(10*i+1);
        states[i][1] = DiscreteLocation(10*i+2);		
    }
    DiscreteLocation end(1000);

    /// Create the dynamics 

    // Main dynamics (a'=1,Vr'=k,Vo'=0,P'=0)
    VectorAffineFunction work_d(Matrix<Float>(4,4),Vector<Float>(4,1.0,k,0.0,0.0));
		
    // For each state
    for (int i=0;i<=NH;i++)
    {
    	automaton.new_mode(states[i][0],work_d);
		automaton.new_mode(states[i][1],work_d);
	}
	automaton.new_mode(end,work_d);
	automaton.new_mode(safe,work_d);	
	automaton.new_mode(unsafe,work_d);	
	
    /// Create the invariants

	// Invariant to stop execution (Vr >= Vhigh)
	VectorAffineFunction stop_i(Matrix<Float>(1,4,0.0,-1.0,0.0,0.0),Vector<Float>(1,Vhigh));

	automaton.new_invariant(safe,stop_i);
	automaton.new_invariant(unsafe,stop_i);

    /// Create the resets

	/// Resets the control clock counter and the output voltage to the value at level n (a^=0,Vr^=Vr,Vo^=Vo[n],P^=P)
    Vector<VectorAffineFunction*> vo_r(NH+1);
	for (int n=0;n<=NH;n++)
	    vo_r[n] = new VectorAffineFunction(Matrix<Float>(4,4, 0.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0, 0.0,0.0,0.0,0.0, 0.0,0.0,0.0,1.0), Vector<Float>(4,0.0,0.0,Vo[n],0.0));  
		
    /// Create the guards
   
    /// Guard for the checking of the voltage difference (a >= P)
    VectorAffineFunction chkdiff_g(Matrix<Float>(1,4,1.0,0.0,0.0,-1.0),Vector<Float>(1));
    /// Guard for a positive voltage difference (Vr >= Vo)
    VectorAffineFunction posdiff_g(Matrix<Float>(1,4,0.0,1.0,-1.0,0.0),Vector<Float>(1));
    /// Guard for a negative voltage difference (Vr <= Vo)
    VectorAffineFunction negdiff_g(Matrix<Float>(1,4,0.0,-1.0,1.0,0.0),Vector<Float>(1));
    /// Guard for the end of the reference voltage excursion (Vr <= Vrf)
    VectorAffineFunction end_g(Matrix<Float>(1,4,0.0,-1.0,0.0,0.0),Vector<Float>(1,Vrf));
    /// Guard for a safe result (Vo <= Vlow)
    VectorAffineFunction safe_g(Matrix<Float>(1,4,0.0,0.0,-1.0),Vector<Float>(1,Vlow));
    /// Guard for a unsafe result (Vo >= Vlow)
    VectorAffineFunction unsafe_g(Matrix<Float>(1,4,0.0,0.0,1.0),Vector<Float>(1,-Vlow));

    /// Create the transitions
    automaton.new_transition(end,to_safe,safe,IdentityFunction(4),safe_g,urgent);
    automaton.new_transition(end,to_unsafe,unsafe,IdentityFunction(4),unsafe_g,urgent);
    // NH
    automaton.new_transition(states[NH][0],chk,states[NH][1],IdentityFunction(4),chkdiff_g,urgent);
	automaton.new_transition(states[NH][1],up,states[NH][0],*vo_r[NH],posdiff_g,urgent);
	automaton.new_transition(states[NH][1],down,states[NH-1][0],*vo_r[NH-1],negdiff_g,urgent);
	   
	// From NH-1 downto 1
	for (int n=NH-1;n>0;n--)
	{
	  	automaton.new_transition(states[n][0],chk,states[n][1],IdentityFunction(4),chkdiff_g,urgent);
	   	automaton.new_transition(states[n][1],up,states[n+1][0],*vo_r[n+1],posdiff_g,urgent);
	   	automaton.new_transition(states[n][1],down,states[n-1][0],*vo_r[n-1],negdiff_g,urgent);
       	automaton.new_transition(states[n][0],to_end,end,IdentityFunction(4),end_g,urgent);
	}

	// 0
	automaton.new_transition(states[0][0],chk,states[0][1],IdentityFunction(4),chkdiff_g,urgent);
	automaton.new_transition(states[0][1],up,states[1][0],*vo_r[1],posdiff_g,urgent);
	automaton.new_transition(states[0][1],down,states[0][0],IdentityFunction(4),negdiff_g,urgent);
    automaton.new_transition(states[0][0],to_end,end,IdentityFunction(4),end_g,urgent);
    automaton.new_transition(states[0][0],to_safe,safe,IdentityFunction(4),safe_g,urgent);
    automaton.new_transition(states[0][0],to_unsafe,unsafe,IdentityFunction(4),unsafe_g,urgent);

    /// Finished building the automaton    

    // Definition of the initial set
    initial_box = Box(4, 0.0,0.0, Vr0,Vr0, Vo0,Vo0, P0,P0);
    initial_state = states[24][0];    
    
    // Definition of the grid
    grid = Grid(Vector<Float>(4, std::pow(2,-1-DEPTH_ADD[0]-MAX_GRID_DEPTH/4), Vr0+std::pow(2,-1-DEPTH_ADD[1]-MAX_GRID_DEPTH/4), Vo0+std::pow(2,-1-DEPTH_ADD[2]-MAX_GRID_DEPTH/4), P0+std::pow(2,-1-DEPTH_ADD[3]-MAX_GRID_DEPTH/4)),Vector<Float>(4, std::pow(2,-DEPTH_ADD[0]), std::pow(2,-DEPTH_ADD[1]), std::pow(2,-DEPTH_ADD[2]), std::pow(2,-DEPTH_ADD[3])));

    // Definition of the graphics parameters
    graphic_box = Box(2, Vrf*0.9, Vhigh*1.1, Vrf*0.9, Vhigh*1.1);
    nvar = 4;
    xvar = 1;
    yvar = 2;
    
}
