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
#ConstrainedImageSet=TaylorConstrainedFlowSet

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
    """Set the value of the jth variable to g. The new domain is the old domain with the jth variable removed.
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

def lie_derivative(g,f):
    assert(isinstance(g,ScalarFunction))
    assert(isinstance(f,VectorFunction))
    assert(g.argument_size()==f.result_size()==f.argument_size())
    r=ScalarFunction(g.argument_size())
    for i in range(0,g.argument_size()):
        r=r+derivative(g,i)*f[i]
    return r

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

class ConstraintData:
    def __init__(self):
        pass
    def __repr__(self):
        res="{ constraint_range:"+str(self.constraint_range)+", derivative_range:"+str(self.derivative_range)+" }"
        return res

class Composer:
    def __call__(self,f1,f2):
        return compose(f1,f2)


class PrototypeHybridEvolver:
    def __init__(self):
        self.integrator=TaylorIntegrator(5)
        self.maximum_step_size=2.0

    def constraint_nonnegative(self,constraint,image_set):
        return (compose(constraint,image_set.function())(image_set.domain())).upper()>=0

    def orbit(self,system,initial_hybrid_set,maximum_hybrid_time,semantics=None):
        assert isinstance(system,HybridAutomaton)
        initial_location=initial_hybrid_set.location()
        initial_set=initial_hybrid_set.continuous_state_set()
        assert isinstance(initial_location,DiscreteState)
        if isinstance(initial_set,IntervalVector):
            initial_set=ConstrainedImageSet(initial_set,VectorFunction.identity(initial_set.size()))
        assert isinstance(initial_set,ConstrainedImageSet)
        initial_time=ScalarFunction.constant(initial_set.domain().size(),0.0)
        assert isinstance(maximum_hybrid_time,HybridTime)
        maximum_time=maximum_hybrid_time.continuous_time()

        initial_events=()
        working_sets=[(initial_events,initial_location,initial_set,initial_time)]
        reach_sets=[]
        evolve_sets=[]

        while len(working_sets)!=0:
            self.upper_evolution_step(working_sets,reach_sets,evolve_sets,system,maximum_time)
        return (reach_sets,evolve_sets)


    def upper_evolution_step(self,working_sets,reach_sets,evolve_sets,system,maximum_time):

        (starting_events,starting_location,starting_set,starting_time)=working_sets[-1]
        assert isinstance(starting_time,ScalarFunction)
        working_sets.pop()
        # FIXME: Remove this!
        if len(starting_set.domain())>4: return

        n=starting_set.dimension()
        m=starting_set._domain.dimension()
        np=starting_set.domain().size()

        print "UPPER EVOLUTION STEP"

        # Extract dynamic and constraints
        starting_mode=system.mode(starting_location)
        dynamic=starting_mode.dynamic()
        invariants=starting_mode.invariants()
        guards={}
        activations={}
        resets={}
        targets={}
        transitions=system.transitions(starting_location)
        for transition in transitions:
            assert isinstance(transition.guard(),ScalarFunction)
            resets[transition.event()]=transition.reset()
            targets[transition.event()]=transition.target().location()
            if transition.urgency():
                guards[transition.event()]=transition.guard()
            else:
                activations[transition.event()]=transition.guard()
        print "\nstarting_events:",starting_events,
        print "\nstarting_location:",starting_location,
        print "\nstarting_set:",starting_set
        print "\nstarting_time:",starting_time
        print "\ninvariants:",invariants,"\nguards:",guards,"\nactivations:",activations,"\nresets:",resets

        # Compute flow
        maximum_step_size=maximum_time
        flow=self.integrator.flow(dynamic,starting_set.bounding_box(),self.maximum_step_size)
        print "flow:",flow
        flow_domain=Box(flow.domain())
        step_size=flow_domain[-1].upper()
        flow_domain[-1]=Interval(0,step_size)
        flow=restrict(flow,flow_domain)
        print "flow:",flow
        print "flow.domain():",flow.domain()
        print "flow.range():",flow.range()
        flowed_set=ConstrainedImageSet(starting_set)
        flowed_set.apply_flow(flow)
        composed_flow=flowed_set.function()
        print "\ncomposed_flow:",composed_flow
        print "composed_flow.domain():",composed_flow.domain()
        print "composed_flow.range():",composed_flow.range()
        flow=composed_flow
        print
        print flowed_set._domain
        print flowed_set._function
        print "\nflowed_set:",flowed_set
        #print flow.polynomial()

        # Compute the elapsed evolution time and the constraint bounding the final time
        print "\nstarting_time:",starting_time
        starting_time_function=embed(starting_time,1)
        print "starting_time_function:",starting_time_function
        dwell_time_function=ScalarFunction.variable(starting_time_function.argument_size(),m)
        print "dwell_time_function:",dwell_time_function
        evolution_time_function=starting_time_function+dwell_time_function
        print "evolution_time_function:",evolution_time_function
        final_constraint=evolution_time_function-maximum_time
        print "final_constraint:",final_constraint

        # Compute restrictions on continuous evolution
        reached_set=ConstrainedImageSet(flowed_set)
        for (event,constraint) in invariants.items()+guards.items():
            derivative=lie_derivative(constraint,dynamic)
            reached_set.new_invariant(event,constraint,derivative)
        if len(starting_events)>0:
            reached_set.new_raw_invariant(DiscreteEvent(0),final_constraint)
        reached_set._location=starting_location
        print reached_set
        reach_sets.append((starting_events,reached_set,))

        # Compute the set reached at the end of the evolution
        final_set=ConstrainedImageSet(reached_set)
        final_set.new_raw_equation(DiscreteEvent(0),final_constraint)
        final_set._location=starting_location
        print "\nreached_set:",reached_set
        print "\nfinal_set:",final_set
        if not final_set.empty():
            print "nonempty_final_set:",final_set
            evolve_sets.append((starting_events,final_set))


        # Compute the reached set under a single event
        for (event,constraint) in guards.items()+activations.items():
            constraint_range=reached_set.compute_range(constraint)
            print "\nconstraint_range("+str(event)+"):",constraint_range
            if constraint_range.upper()>=0:
                derivative=lie_derivative(constraint,dynamic)
                target=targets[event]
                reset=resets[event]
                jumping_set=ConstrainedImageSet(reached_set)
                if event in guards:
                    jumping_set.new_guard(event,constraint,derivative)
                else:
                    jumping_set.new_activation(event,constraint,derivative)
                if jumping_set.empty():
                    print "empty_jumping_set:",jumping_set
                else:
                    #jumped_set=apply(reset,jumping_set)
                    jumped_set=ConstrainedImageSet(jumping_set)
                    jumped_set.apply_map(reset)
                    jumped_set._location=target
                    print "\njumped_set("+str(event)+"):",jumped_set
                    dwell_time_function=ScalarFunction.variable(flow.argument_size(),flow.argument_size()-1)
                    jumped_time_function=embed(starting_time,1)+dwell_time_function
                    working_sets.append((starting_events+(event,),target,jumped_set,jumped_time_function))
        print "\nDONE STEP\n"

        raw_input("Press any key to continue...")

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

def plot(filename,set,bounding_box=None,resolution=4):
    fig=Figure()
    if bounding_box==None:
        bounding_box=set.bounding_box()
    fig.set_bounding_box(bounding_box)
    set.draw(fig,resolution)
    #for box in self.grid_outer_approximation(resolution):
        #    fig.draw(box)
    fig.write(filename)


def build_system():
    q1=DiscreteState(1)
    q2=DiscreteState(2)
    q3=DiscreteState(3)
    e12=DiscreteEvent(12)
    e23=DiscreteEvent(23)
    z=ScalarFunction.constant(2,0.0)
    o=ScalarFunction.constant(2,1.0)
    x=ScalarFunction.variable(2,0)
    y=ScalarFunction.variable(2,1)

    hsys=HybridAutomaton()
    hsys.new_mode(q1,[o,z])
    hsys.new_mode(q2,[o,z])
    hsys.new_mode(q3,[o,z])
    #hsys.new_invariant(q1,x-y-0.5)
    hsys.new_invariant(q2,x-4*pow(y-1,2)-1.5)
    hsys.new_transition(e12,q1,q2,[x,y+1],x-y-0.5,True)
    hsys.new_transition(e23,q2,q3,[x,y+1],x-4*pow(y-1,2)-1.5,False)

    #hsys.new_mode(q1,[0.0*x+1.0,0.0*y+2.0])
    #sys.new_mode(q1,[0.5*x+y,-x-0.5*y])
    #hsys.new_transition(e12,q1,q2,[x+1,y-1],x-0.175,False)
    #hsys.new_transition(e12,q1,q2,VectorFunction([x+1,y+1]),x-1,False)

    initial_location=q1
    initial_box=Box([{-0.0625:0.0625},{-0.0625:0.0625}])
    initial_box=Box([{-0.125:0.125},{-0.125:0.125}])
    initial_box=Box([{-0.00390625:0.00390625},{-0.125:0.125}])
    initial_box=Box([{-0.00390625:0.00390625},{-0.25:0.25}])
    hsys.initial_set=HybridBox(initial_location,initial_box)
    hsys.evolution_time=HybridTime(2.5,5)
    hsys.evolution_time=HybridTime(2.0,5)
    hsys.bounding_box=Box([{-0.5:3.5},{-0.5:3.5}])

    return hsys


if __name__=="__main__":
    #test_function()

    from watertank_automaton import watertank

    err=Interval(-1,+1)*0.03125;
    initial_location=DiscreteState("opening") #opening
    initial_box=Box([1.0+err,0.0625+err])
    watertank.initial_set=HybridBox(initial_location,initial_box)
    watertank.evolution_time=HybridTime(0.5,12)
    watertank.bounding_box=Box([{0:9},{-0.1:1.1}])

    automaton=build_system()
    resolution=4

    print
    print automaton
    print
    print [ mode.location() for mode in automaton.modes() ]
    print [ (mode.location(), [ (transition.event(),transition.target().location()) for transition in automaton.transitions(mode.location()) ]) for mode in automaton.modes() ]
    print


    print "computing evolution..."
    evolver=PrototypeHybridEvolver()
    (reached_sets,evolved_sets)=evolver.orbit(automaton,automaton.initial_set,automaton.evolution_time)

    #evolver=HybridEvolver()
    #orbit = evolver.orbit(automaton,automaton.initial_set,automaton.evolution_time)
    #print orbit
    #reached_sets=orbit.reach()
    #evolved_sets=orbit.evolve()


    print "\nreached_sets:"
    for (events,reached_set) in reached_sets:
        print str(events+(reached_set._location,))+":\n"+str(reached_set)
    print "\nevolved_sets:"
    for (events,evolved_set) in evolved_sets:
        print str(events)+":\n"+str(evolved_set)
    print


    bounding_box=automaton.bounding_box
    print "reach..."
    i=0;
    for (events,reached_set) in reached_sets:
        i=i+1
        print "   ",events+(reached_set._location,)
        #print events,evolved_set.function,evolved_set.maximum_negative_constraints
        plot("hybrid_evolution-r"+str(i)+".png",reached_set,bounding_box,resolution)
    i=0;
    print "\nevolve..."
    for (events,evolved_set) in evolved_sets:
        i=i+1
        print "   ",events+(evolved_set._location,)
        #print events,evolved_set.function,evolved_set.maximum_negative_constraints
        plot("hybrid_evolution-e"+str(i)+".png",evolved_set,bounding_box,resolution)
