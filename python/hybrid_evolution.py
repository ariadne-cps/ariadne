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

# Redifine string function for intervals
def interval_str(self):
    return '{' + str(self.lower()) + ':' + str(self.upper()) + '}'
Interval.__str__=interval_str
del interval_str

Taylor=ScalarTaylorFunction

#Set representation functions
Taylor.__repr__=Taylor.__str__
Interval.__repr__=Interval.__str__
VectorFunction.__repr__=VectorFunction.__str__
ScalarFunction.__repr__=ScalarFunction.__str__
VectorTaylorFunction.__repr__=VectorTaylorFunction.__str__
ScalarTaylorFunction.__repr__=ScalarTaylorFunction.__str__
DiscreteState.__repr__=DiscreteState.__str__
DiscreteEvent.__repr__=DiscreteEvent.__str__
DiscreteMode.__repr__=DiscreteMode.__str__
DiscreteTransition.__repr__=DiscreteTransition.__str__
ConstrainedImageSet.__repr__=ConstrainedImageSet.__str__

# Define extra functionality for TaylorFunctions
def taylor_set(f,j,g):
    """Set the value of the jth variable to g. The new domain is the same as the old domain.
       Either g is a constant, or g is a scalar variable with domain equal to the new domain."""
    n=f.domain().size()
    subdomain=Box( [ f.domain()[i] for i in range(0,j)+range(j+1,n) ] )
    if not isinstance(g,ScalarTaylorFunction):
        g=ScalarTaylorFunction.constant(subdomain,g)
    assert subdomain==Box(g.domain())
    h=VectorTaylorFunction([ ScalarTaylorFunction.variable(subdomain,i) for i in range(0,j) ] + [ g ] + [ ScalarTaylorFunction.variable(subdomain,i) for i in range(j,n-1) ])
    return compose(f,h)

def taylor_substitute(f,j,g):
    """Substitute the value of the jth variable with g. 
       Either g is a constant, or g is a scalar variable with domain equal to the domain of f."""
    domain=Box(f.domain())
    if not isinstance(g,ScalarTaylorFunction):
        g = ScalarTaylorFunction.constant(domain,g)
    assert domain==Box(g.domain())
    h=VectorTaylorFunction.identity(domain)
    h[j]=g
    return compose(f,h)

ScalarTaylorFunction.set=taylor_set
VectorTaylorFunction.set=taylor_set
ScalarTaylorFunction.substitute=taylor_substitute
VectorTaylorFunction.substitute=taylor_substitute
del taylor_set
del taylor_substitute

class Set(set):
    def __str__(set):
        res="{"
        first=True
        for x in set:
            if first: first=False
            else: res+=","
            res+=str(x)
        return res+"}"

class Data:
    def __init__(self):
        pass
    def __repr__(self):
        res="{ "
        for field in dir(self):
            if field[0]!='_': res+=str(field)+":"+str(getattr(self,field))+", "
        return res+"}"
    
def evolve(system,initial_location,initial_set,initial_time,maximum_time):
    integrator=TaylorIntegrator(5)
    
    assert isinstance(system,HybridAutomaton)
    assert isinstance(initial_location,DiscreteState)
    if isinstance(initial_set,IntervalVector):
        initial_set=ConstrainedImageSet(initial_set,VectorFunction.identity(initial_set.size()))
    assert isinstance(initial_set,ConstrainedImageSet)
    if isinstance(initial_time,float):
        initial_time=ScalarTaylorFunction.constant(initial_set.domain,initial_time)
    assert isinstance(initial_time,ScalarTaylorFunction)
    assert initial_time.domain()==initial_set.domain

    initial_bounding_box=Box(initial_set.bounding_box())

    # Extract dynamic and constraints
    initial_mode=system.mode(initial_location)
    dynamic=initial_mode.dynamic()
    invariants=initial_mode.invariants()
    transitions=system.transitions(initial_location)
    guards={}
    activations={}
    resets={}
    for transition in transitions:
        assert isinstance(transition.guard(),ScalarFunction)
        resets[transition.event()]=transition.reset()
        if transition.urgency():
            guards[transition.event()]=transition.guard()
        else:
            activations[transition.event()]=transition.guard()
    print "\ninvariants:",invariants,"\nguards:",guards,"\nactivations:",activations,"\nresets:",resets
    
    # Compute initially active guards
    print
    initial_ranges={}
    no_flow=False
    for (event,constraint) in invariants.items() + guards.items():
        constraint_range=constraint(initial_bounding_box)
        print str(event)+".initial_range():",constraint_range
        initial_ranges[event]=constraint(initial_bounding_box)
        if constraint_range.lower()>0.0:
            no_flow=True
    print
    assert(no_flow==False)
    
    # Compute flow
    flow=integrator.flow(dynamic,initial_set.bounding_box(),maximum_time)
    print "flow:",flow
    flow_domain=Box(join(initial_set.bounding_box(),Interval(0,maximum_time)))
    flow=restrict(flow,flow_domain)
    print "flow:",flow
    flow_range=Box(flow.range())
    print "flow.domain():",flow_domain
    print "flow.range():",flow_range
    print
    #print flow.polynomial()

    # Compute ranges of constraints
    constraint_data={}
    for (event,constraint) in invariants.items()+guards.items()+activations.items():
        data=Data()
        flowed_constraint=compose(constraint,flow)
        constraint_range1=flowed_constraint(flow_domain)
        constraint_range2=constraint(flow(flow_domain))
        constraint_range=intersection(constraint_range1,constraint_range2)
        print event,constraint_range1,constraint_range2,constraint_range
        if constraint_range.upper()>=0.0 or True:
            constraint_derivative=lie_derivative(constraint,dynamic)
            flowed_constraint_derivative=compose(constraint_derivative,flow)
            constraint_derivative_range=flowed_constraint_derivative(flow_domain)
            data.constraint_range=constraint_range
            data.derivative_range=constraint_derivative_range
            data.flowed_constraint=flowed_constraint
            constraint_data[event]=data
    print "constraint_data:",constraint_data
    
    event_sequence=()
    evolved_set=ConstrainedImageSet(flow.domain(),flow,[],[],[],[])
    evolved_sets={tuple(event_sequence):evolved_set}

    # Determine which constraints may limit progress
    progress_constraints={}
    for event in invariants.keys()+guards.keys():
        data=constraint_data[event]
        if data.constraint_range.upper() >= 0.0:
            progress_constraints[event]=data.flowed_constraint
    print "progress_constraints:", [ event for event in progress_constraints ]
    
    # Determine which constraints may activate transitions
    transition_constraints={}
    for event in activations.keys()+guards.keys():
        data=constraint_data[event]
        if data.constraint_range.upper() >= 0.0:
            transition_constraints[event]=data.flowed_constraint
    print "transition_constraints:", [ event for event in transition_constraints ]

    # Compute the reached set under the flow
    reach_set=ConstrainedImageSet(domain=flow_domain,function=flow,maximum_negative_constraints=progress_constraints)
    print "\nreach_set:",reach_set,"\n"

    # Compute the reached set under a single event
    for event in transition_constraints:
        reset=resets[event]
        reset_flow=compose(reset,flow)
        event_progress_constraints=dict(progress_constraints)
        if event in event_progress_constraints:
            # event is urgent
            event_guard_constraint={event:progress_constraints[event]}
            event_progress_constraints.remove(event)
            event_reach_set=ConstrainedImageSet(domain=flow_domain,function=reset_flow,
                                                maximum_negative_constraints=event_progress_constraints,
                                                maximum_equality_constraints=event_guard_constraint)
        else:
            event_activation_constraint={event:transition_constraints[event]}
            event_reach_set=ConstrainedImageSet(domain=flow_domain,function=reset_flow,
                                                maximum_negative_constraints=progress_constraints,
                                                positive_constraints=event_activation_constraint)
        print "\njumped_set:",jumped_set,"\n"

    print "\nDONE\n"

    return { (): reach_set }
    for (event,guard,reset,urgency) in \
            [ (transition.event(),transition.guard(),transition.reset(),transition.urgency()) for transition in transitions ]:
        event_str="g"+str(event)[1:]
        assert isinstance(guard,ScalarFunction)
        guard_derivative=lie_derivative(guard,dynamic)
        guard_range1=compose(guard,flow).range()
        guard_range2=guard(flow.range())
        guard_range=intersection(guard_range1,guard_range2)
        guard_derivative_range=compose(guard_derivative,flow).range()
        print event_str+".range():",guard_range1
        print event_str+".range():",guard_range2
        print "d"+event_str+".range():",guard_derivative_range
        if guard_range.upper() < 0.0:
            print "  guard",event_str,"is not active"
        elif evaluate(guard,initial_set).lower() > 0.0:
            print "  guard",event_str,"is initially active"
        else: 
            print "  guard",event_str,"may be active"
            print "  reset:",reset
            new_event_sequence=event_sequence+(event,)
            new_domain=evolved_set.domain
            new_function=compose(reset,flow)
            new_invariants=None #Compute new invariants here
            new_guards=None #Compute new guards here
            print new_domain,new_function
            new_evolved_set=ConstrainedImageSet(evolved_set)
            new_evolved_set.function=new_function
            new_evolved_set.positive_constraints.append((event,compose(guard,flow)))
            #new_evolved_set.positive_constraints.append((event,compose(guard_derivative,flow)))
            #(join(initial_set,join(Interval(0,time),Interval(0,time))))
            evolved_sets[new_event_sequence]=new_evolved_set
    return evolved_sets


def test_function():
    tx = [ ScalarTaylorFunction.variable([{-2:2},{-1:1},{0:3}],i) for i in range(0,3) ]
    ty = [ ScalarTaylorFunction.variable([{-2:2},{0:3}],i) for i in range(0,2) ]
    tf=tx[0]+2*tx[1]+3*tx[2]
    tg3=tx[0]*tx[1]
    tg2=ty[0]*ty[1]
    print tf
    print tg2
    print tg3

    print tf.set(1,0.375)
    print tf.substitute(1,0.375)
    print tf.set(1,0.125*tg2)
    print tf.substitute(1,0.125*tg3)


if __name__=="__main__":
    #test_function()

    q1=DiscreteState(1)
    q2=DiscreteState(2)
    e12=DiscreteEvent(12)
    z=ScalarFunction.constant(2,0.0)
    o=ScalarFunction.constant(2,1.0)
    x=ScalarFunction.variable(2,0)
    y=ScalarFunction.variable(2,1)

    hsys=HybridAutomaton()
    hsys.new_mode(q1,[o,3*x])
    #hsys.new_mode(q1,[0.0*x+1.0,0.0*y+2.0])
    #sys.new_mode(q1,[0.5*x+y,-x-0.5*y])
    hsys.new_invariant(q1,y-0.126-0.5)
    hsys.new_mode(q2,[x,y])
    #hsys.new_transition(e12,q1,q2,[x+1,y-1],x-0.175,False)
    hsys.new_transition(e12,q1,q2,[x+1,y-1],y-0.5,False)
    #hsys.new_transition(e12,q1,q2,VectorFunction([x+1,y+1]),x-1,False)

    DiscreteState.__repr__=DiscreteState.__str__
    DiscreteEvent.__repr__=DiscreteEvent.__str__
    

    print hsys
    print [ mode.location() for mode in hsys.modes() ]
    print [ (mode.location(), [ (transition.event(),transition.target().location()) for transition in hsys.transitions(mode.location()) ]) for mode in hsys.modes() ]
    print

    m1=hsys.mode(q1)
    print m1
    print m1.location(), m1.dynamic(),m1.invariants()
    t12=hsys.transition(e12,q1)
    print t12.event(), t12.source().location(), t12.target().location(), t12.reset(), t12.guard(), t12.urgency()
    print
    
    initial_location=q1
    initial_set=Box([{0:0.125},{0:0.125}])
    initial_set=Box([{0:0.0625},{0:0.0625}])

    print hsys.modes()
    print hsys.transitions()
    print hsys.transitions(q1)

    initial_time=0.0
    evolution_time=0.5
    
    evolved_sets = evolve(hsys,initial_location,initial_set,initial_time,evolution_time)

    print
    print "evolved_sets=\n"+str(evolved_sets)
    print
    bounding_box=Box([{-1:2},{-2:1}])
    i=0;
    for (events,evolved_set) in evolved_sets.items():
        print events,evolved_set.function,evolved_set.maximum_negative_constraints
        evolved_set.plot("hybrid_evolution-"+str(i)+".png",bounding_box,resolution=6)
        i=i+1