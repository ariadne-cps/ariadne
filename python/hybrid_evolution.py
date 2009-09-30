#!/usr/bin/python

##############################################################################
#            hybrid_evolution.py
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
from ariadne import *
from constrained_image_set import *
#from hybrid_system import *

Indeterminate=indeterminate()
Float = float

print dir()

Taylor=ScalarTaylorFunction
Taylor.__repr__=Taylor.__str__
VectorFunction.__repr__=VectorFunction.__str__
ScalarFunction.__repr__=ScalarFunction.__str__
VectorTaylorFunction.__repr__=VectorTaylorFunction.__str__
ScalarTaylorFunction.__repr__=ScalarTaylorFunction.__str__
DiscreteMode.__repr__=DiscreteMode.__str__
DiscreteTransition.__repr__=DiscreteTransition.__str__
ConstrainedImageSet.__repr__=ConstrainedImageSet.__str__


def evolve(system,initial_set,time):
    integrator=TaylorIntegrator(5)
    
    assert(isinstance(system,HybridAutomaton))
    (initial_location,initial_set)=initial_set
    assert(isinstance(initial_location,DiscreteState))
    assert(isinstance(initial_set,Box))
    initial_mode=system.mode(initial_location)

    dynamic=initial_mode.dynamic()
    invariants=initial_mode.invariants().values()
    transitions=sys.transitions(initial_location)
    flow=integrator.flow(dynamic,initial_set,time)
    flow_domain=Box(join(initial_set,Interval(0,time)))
    flow=restrict(flow,flow_domain)
    print flow
    
    event_sequence=[]
    evolved_set=ConstrainedImageSet(Box(join(initial_set,Interval(0,time))),flow,[],[],[],[])
    evolved_sets={tuple(event_sequence):evolved_set}
    for invariant in invariants:
        print flow.range()
        if compose(invariant,flow).range()[0].upper() < 0.0:
            print "invariant",invariant,"is not active"
        elif evaluate(invariant,initial_set)[0].lower() > 0.0:
            print "invariant",invariant,"is initially active"
        else:
            print "invariant",invariant,"may be active"
            evolved_set.maximum_negative_constraints.append(compose(invariant,flow))
    
    for (event,guard) in [ (transition.event(),transition.guard()) for transition in transitions ]:
        if compose(guard,flow).range()[0].upper() < 0.0:
            print "guard",guard,"is not active"
        elif evaluate(guard,initial_set)[0].lower() > 0.0:
            print "guard",guard,"is initially active"
        else:
            print "guard",guard,"may be active"
            new_event_sequence=event_sequence+event
            new_domain=Box(join(evolved_set.domain(),Interval(0,time)))
            new_function=None #Compute new function here
            new_invariants=None #Compute new invariants here
            new_guards=None #Compute new guards here
            new_evolved_set=Box(join(initial_set,join(Interval(0,time),Interval(0,time))))
    
    print
    print evolved_set
    


if __name__=="__main__":
    q1=DiscreteState(1)
    q2=DiscreteState(2)
    e12=DiscreteEvent(12)
    x=ScalarFunction.variable(2,0)
    y=ScalarFunction.variable(2,1)

    sys=HybridAutomaton()
    sys.new_mode(q1,[0.0*x+1.0,x])
    #sys.new_mode(q1,[0.5*x+y,-x-0.5*y])
    sys.new_invariant(q1,y-0.126)
    sys.new_mode(q2,[x,y])
    sys.new_transition(e12,q1,q2,VectorFunction([x+1,y+1]),x-1,False)

    print sys
    
    m1=sys.mode(q1)
    print m1
    print m1.location(), m1.dynamic(),m1.invariants().values()
    t12=sys.transition(e12,q1)
    print t12.event(), t12.source(), t12.target(), t12.reset(), t12.guard(), t12.urgency()
    print
    
    initial_location=q1
    initial_set=Box([{0:0.125},{0:0.125}])

    print sys.modes()
    print sys.transitions()
    print sys.transitions(q1)
    evolve(sys,(initial_location,initial_set),0.25)