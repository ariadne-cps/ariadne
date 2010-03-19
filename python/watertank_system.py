#!/usr/bin/python

##############################################################################
#            watertank-system.py
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

def BooleanConstant(name,value=None):
    if value==None:
        value=name
        name=str(value)
    return DiscretePredicate(value)
def TriboolConstant(name,value=None):
    if value==None:
        value=name
        name=str(value)
    return ContinuousPredicate(value)
def StringConstant(name,value=None):
    if value==None:
        value=name
        name=str(value)
    return StringExpression(value)
def RealConstant(name,value=None):
    if value==None:
        value=name
        name=str(value)
    return RealExpression(value)

if __name__=='__main__':

    all=~EventSet() # Used for guards/resets true for any set
    any=BooleanConstant('any',True) # Used for equations true in any location

    # Create the system object
    watertank=HybridSystem();

    # Declare a discrete variable
    valve=StringVariable("valve"); # Values "closed", "opening", "open", "closing"

    # Declare some constants. Note that system parameters should be given as variables.
    T=RealConstant("T",4.0);
    hmin=RealConstant("hmin",5.5);
    hmax=RealConstant("hmax",8.0);
    delta=RealConstant("delta",0.05);
    lamb=RealConstant("lambda",0.02);
    b=RealConstant("b",0.3);

    # Declare the system variables
    x=RealVariable("x");
    alpha=RealVariable("alpha");

    # Declare the events we use
    start_opening=Event("start_opening");
    start_closing=Event("start_closing");
    finished_opening=Event("finished_opening");
    finished_closing=Event("finished_closing");

    # The water level is always given by the same dynamic.
    # When alpha is 1, the valve is open and water flows in
    # When alpha is 0, the valve is open and no water flows in
    watertank.new_dynamic(any, dot(x) << -lamb*x+b*alpha);

    # Specify the equation for how the valve opens/closes
    watertank.new_dynamic(valve=="opening", dot(alpha) << +1.0/T);
    watertank.new_dynamic(valve=="closing", dot(alpha) << -1.0/T);

    watertank.new_dynamic(valve=="opening", [dot(alpha) << +1.0/T, dot(x) << -lamb*x+b*alpha]);
    # When the valve is open or closed, alpha is constant.
    # Note that since we know alpha=0.0 or alpha=1.0, we should not need to consider alpha as a state variable.
    # This requires some cleverness on the part of the symbolic system analyser.
    # It would be possible to model the system with alpha explicitly a constant, but this would require some
    # cleverness on the part of the modeler.
    watertank.new_dynamic((valve=="closed") | (valve=="open"), dot(alpha) << RealConstant(0.0));

    # Specify the condition that the valve starts opening when hmax <= x <= hmax+delta
    # using an invariant and guard.
    watertank.new_invariant(any, x<=hmax+delta);
    watertank.new_guard((valve=="open") | (valve=="opening"), start_closing, x>=hmax);

    # Specify the condition that the valve starts closing when hmin <= x <= hmin+delta
    # using a combined 'invariant and activation'. The event may occur when x<=hmin, and
    # must occur while x>=hmin-delta.
    watertank.new_guard((valve=="closed") | (valve=="closing"), start_opening, x<=hmin,x>=hmin-delta);

    # Specify the guards for when the valve reaches the desired position
    watertank.new_guard(valve=="opening", finished_opening, alpha>=1.0, alpha<=1.0);
    watertank.new_guard(valve=="closing", finished_closing, alpha<=0.0, alpha>=0.0);

    # Explicitly disallow events when they don't make sense
    watertank.new_guard(~(valve=="opening"), finished_opening, TriboolConstant(False));
    watertank.new_guard(~(valve=="closing"), finished_closing, TriboolConstant(False));
    watertank.new_guard(~((valve=="open") | (valve=="opening")), start_closing, TriboolConstant(False));
    watertank.new_guard(~((valve=="closed") | (valve=="closing")), start_opening, TriboolConstant(False));

    # Specify the discrete resets
    watertank.new_transition(any, finished_opening, next(valve) << StringConstant("open"));
    watertank.new_transition(any, finished_closing, next(valve) << StringConstant("closed"));
    watertank.new_transition(any, start_opening, next(valve) << StringConstant("opening"));
    watertank.new_transition(any, start_closing, next(valve) << StringConstant("closing"));

    # For any event occurring in any location, the value of x and alpha are not updated.
    watertank.new_reset(any, all, next(x) << x);
    watertank.new_reset(any, all, next(alpha) << alpha);


    # Finished building the automaton

    print "Watertank: \n", watertank
