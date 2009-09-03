#!/usr/bin/python

##############################################################################
#            hybrid_system.py
#
#  Copyright 2009  Pieter Collins <Pieter.Collins@cwi.nl>
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

import sys

from formula import *

class Event:
    def __init__(self,name):
        self.name=str(name)
    def __hash__(self):
        return hash(self.name)
    def __eq__(self,other):
        return self.name==other.name
    def __str__(self):
        return "Event('"+self.name+"')"
    def __repr__(self):
        return self.name
    def __invert__(self):
        return EventSet(self).__invert__()
    def __contains__(self,event):
        return (event is self)

class EventSet:
    def __init__(self,events=set()):
        if isinstance(events,Event):
            self.set=set([events])
            self.compl=False
        elif isinstance(events,EventSet):
            self.set=events.set
            self.compl=events.compl
        else:
            self.set=set(events)
            self.compl=False

    def __invert__(self):
        result=EventSet(self)
        result.compl=not result.compl
        return result

    def __contains__(self,event):
        return (event in self.set) ^ self.compl

    def __str__(self):
        if len(self.set)==0:
            if self.compl: return 'all'
        if len(self.set)==1:
            if self.compl: return '~['+repr(list(self.set)[0])+']'
        if self.compl:
            return '~'+str(list(self.set))
        else:
            return str(list(self.set))

    def __repr__(self):
        return str(self)


"""A hybrid system, defined by dynamical rules.

A hybrid system has the following attributes:
    transitions: A list of (Event,State,Variable,Expression) tuples
    equations: A list of (State,Variable,Expression) tuples, giving algebraic equations
    dynamics: A list of (State,Variable,Expression) tuples, giving differential equations
    resets: A list of (Event,State,Variable,Expression) tuples
    guards: A list of (Event,State,Activation,Invariant) tuples,
"""

class HybridSystem:
    def __init__(self):
        self.transitions=[]
        self.equations=[]
        self.dynamics=[]
        self.invariants=[]
        self.guards=[]
        self.resets=[]

    def new_transition(self,events,locations,variable,expression):
        self.transitions.append((events,locations,variable,expression))

    def new_equation(self,locations,variable,expression):
        self.equations.append((locations,variable,expression))

    def new_dynamic(self,locations,variable,expression):
        self.dynamics.append((locations,variable,expression))

    def new_reset(self,events,locations,variable,expression):
        self.resets.append((events,locations,variable,expression))

    def new_guard(self,events,locations,activation,invariant=None):
        #if invariant==None: invariant=activation
        if invariant==None: self.guards.append((events,locations,activation))
        else: self.guards.append((events,locations,activation,invariant))

    def new_invariant(self,locations,predicate):
        self.invariants.append((locations,predicate))

    def order(self,assignments):
        active={}
        variables=set()
        for (variable,expression) in assignments:
            active[variable]=expression
            variables.add(variable)
        # Order relations by variable use
        print active
        print variables
        arguments={}
        for (variable,expression) in active.items():
            arguments[variable]=expression.variables().intersection(variables)
        result=[]
        while arguments:
            ok=False
            for (variable,dependencies) in arguments.items():
                if not dependencies:
                    for var in arguments.keys():
                        arguments[var].discard(variable)
                    del arguments[variable]
                    result.append((variable,active[variable]))
                    ok=True
            if not ok:
                print "Circular dependency in",arguments.keys()
        print "arguments:",arguments
        return result

    def active_variables(self,mode):
        argument_variables=set()
        differential_variables=set()
        algebraic_variables=set()
        for (locations,variable,expression) in self.equations:
            if locations.evaluate(mode):
                assert variable not in algebraic_variables
                algebraic_variables.add(variable)
                argument_variables.update(expression.variables())
        for (locations,variable,expression) in self.dynamics:
            if locations.evaluate(mode):
                assert variable not in algebraic_variables and variable not in differential_variables
                differential_variables.add(variable)
                argument_variables.update(expression.variables())
        #assert argument_variables.issubset(differential_variables.union(algebraic_variables))
        return (algebraic_variables,differential_variables,argument_variables)

    def active_equations(self,mode):
        result=[]
        for (locations,variable,expression) in self.equations:
            if locations.evaluate(mode):
                result.append((variable,expression))
        return self.order(result)

    def active_dynamics(self,mode):
        result=[]
        for (locations,variable,expression) in self.dynamics:
            if locations.evaluate(mode):
                result.append((variable,expression))
        return result

    def active_guards(self,mode):
        result=[]
        for guard in self.guards:
            (events,locations,activation) = guard[0:3]
            if locations.evaluate(mode):
                result.append((events,activation))
        return result

    def active_resets(self,active_event,mode):
        result=[]
        for (events,locations,variable,expression) in self.resets:
            if active_event in events and locations.evaluate(mode):
                result.append((variable,expression))
        return result

    def active_transitions(self,active_event,mode):
        result=[]
        for (events,locations,variable,expression) in self.transitions:
            if active_event in events and locations.evaluate(mode):
                result.append((variable,expression))
        return result

    def __str__(self):
        return "HybridSystem(\n  equations="+str(self.equations) \
                 + ",\n  dynamics="+str(self.dynamics) \
                 + ",\n  invariants="+str(self.invariants) \
                 + ",\n  guards="+str(self.guards) \
                 + ",\n  resets="+str(self.resets)+"\n)"

