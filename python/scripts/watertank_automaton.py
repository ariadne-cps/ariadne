#!/usr/bin/python

##############################################################################
#            watertank-automaton.py
#
#  Copyright 2008-9  Davide Bresolin, Pieter Collins
##############################################################################

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

from ariadne import *

def bounding_box_str(self):
    return str(self.bounding_box())

DiscreteState.__repr__=DiscreteState.__str__
TaylorSet.__repr__=bounding_box_str


# Create the system object
watertank=HybridAutomaton();

# Declare some constants. Note that system parameters should be given as variables.
T=4.0;
hmin=5.5;
hmax=8.0;
delta=0.05;
lamb=0.02;
b=0.3;

# Declare the system variables
x=ScalarFunction.variable(2,0);
alpha=ScalarFunction.variable(2,1);

# Declare useful constant functions
zero=ScalarFunction.constant(2,0.0);
recT=ScalarFunction.constant(2,1.0/T);

open=DiscreteState("open");
closed=DiscreteState("closed");
opening=DiscreteState("opening");
closing=DiscreteState("closing");

# Declare the events we use
start_opening=DiscreteEvent("start_opening");
start_closing=DiscreteEvent("start_closing");
finished_opening=DiscreteEvent("finished_opening");
finished_closing=DiscreteEvent("finished_closing");

identity=VectorFunction.identity(2)

# The water level is always given by the same dynamic.
# When alpha is 1, the valve is open and water flows in
# When alpha is 0, the valve is open and no water flows in
# The valve opens and closes at rate 1/T
watertank.new_mode(open,[-lamb*x+b*alpha,zero]);
watertank.new_mode(closed,[-lamb*x+b*alpha,zero]);
watertank.new_mode(opening,[-lamb*x+b*alpha,recT]);
watertank.new_mode(closing,[-lamb*x+b*alpha,-recT]);

# Example syntax for invariants
# watertank.new_invariant(open,x-hmax)

# Specify the condition that the valve starts closing when hmax <= x <= hmax+delta
# using an invariant and guard.
#watertank.new_invariant(open, x-(hmax+delta));
#watertank.new_transition(start_closing,open,closing,identity,x-hmax,False)
#watertank.new_transition(start_closing,opening,closing,identity,x-hmax,False)
watertank.new_transition(start_closing,open,closing,[x,alpha],x-hmax,True)
watertank.new_transition(start_closing,opening,closing,identity,x-hmax,True)

# Specify the condition that the valve starts opening when hmin <= x <= hmin+delta
# using a combined 'invariant and activation'. The event may occur when x<=hmin, and
# must occur while x>=hmin-delta.
watertank.new_transition(start_opening,closed,opening,identity,hmin-x,True)
watertank.new_transition(start_opening,closing,opening,identity,hmin-x,True)

watertank.new_transition(finished_opening,opening,open,identity,alpha-1,True)
watertank.new_transition(finished_closing,closing,closed,identity,-alpha,True)

if __name__=='__main__':


    print dir()
    #evolver=StableHybridEvolver()
    #evolver=ImageSetHybridEvolver()
    evolver=ConstrainedImageSetHybridEvolver()
    evolver=HybridEvolver()

    err=Interval(-1,+1)*0.03125;
    initial_location=opening
    initial_box=Box([1.0+err,0.125+err])

    initial_set=HybridBox(initial_location,initial_box)

    evolution_time=HybridTime(101.0,12)

    orbit=evolver.orbit(watertank,initial_set,evolution_time)
    reach=orbit.reach()
    evolve=orbit.evolve()

    print orbit
    print
    print reach
    print evolve

    print len(reach[open])

    fig=Figure()
    fig.set_bounding_box([{0.0:8.5},{-0.1:+1.1}])
    fig.set_fill_colour(1.0,0.0,0.0)
    for (loc,sets) in reach.items():
        print loc,len(sets)
        for set in sets:
            fig.draw(set)
    fig.set_fill_colour(1.0,1.0,0.0)
    for (loc,sets) in evolve.items():
        for set in sets:
            fig.draw(set.bounding_box()+IntervalVector([{-0.16:0.16},{-0.02:0.02}]))
    fig.write("watertank-automaton")

