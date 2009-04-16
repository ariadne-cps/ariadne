#!/usr/bin/python

##############################################################################
#            hybrid_evolver.py
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

Indeterminate=indeterminate()

def is_pos(ivl):
    if ivl.lower()>0 : return True
    if ivl.upper()<0 : return False
    return Indeterminate

def is_zero(ivl):
    if ivl.lower()==0 and ivl.upper()==0: return True
    if ivl.lower()>0 or ivl.upper()<0: return False
    return Indeterminate

Interval.is_zero=is_zero
Interval.is_pos=is_pos


class HybridEvolver:

    temporal_order=4;
    spacial_order=2;

    starting_event="starting_event"
    finishing_event="finishing_event"
    blocking_event=-3
    invariant_event=-3


    def possibly_touching(self,guard,set):
        guard_range=guard(set)
        if not isinstance(guard_range,Interval):
            guard_range=guard_range.range()[0]
        if guard_range.lower()>0 or guard_range.upper()<0:
            return False
        else:
            return True


    def is_active(self,guard,set):
        guard_range=guard(set)
        if not isinstance(guard_range,Interval):
            guard_range=guard_range.range()[0]
        return is_pos(guard_range)


    def initially_active(self, guards, starting_set):
        # Determine which events are initially active
        initially_active={}
        initially_active["blocking_event"]=False
        #print guards.items()
        for (event,(guard,urgent)) in guards.items():
            initially_active[event]=self.is_active(guard,starting_set)
            if urgent and definitely(initially_active[event]):
                initially_active["blocking_event"]=True
        return initially_active


    def flow_bounds(self,dynamic, starting_box, suggested_time_step):
        assert(isinstance(dynamic,FunctionInterface))
        time_step=suggested_time_step
        small_neighbourhood=starting_box+(starting_box-midpoint(starting_box))/16
        one_box_neighbourhood=starting_box+(starting_box-midpoint(starting_box))
        assert(subset(starting_box,one_box_neighbourhood))
        assert(subset(small_neighbourhood,small_neighbourhood))
        print starting_box,small_neighbourhood,one_box_neighbourhood
        bound=one_box_neighbourhood+dynamic(one_box_neighbourhood)*(4*Interval(0,suggested_time_step))
        bound=one_box_neighbourhood+dynamic(bound)*(2*Interval(0,suggested_time_step))
        while(not subset(small_neighbourhood+time_step*dynamic(bound),bound)):
            time_step/=2
            assert(time_step>1e-12)
        bound=small_neighbourhood+Interval(0,time_step)*dynamic(bound)
        bound=small_neighbourhood+Interval(0,time_step)*dynamic(bound)
        assert(subset(starting_box+time_step*dynamic(bound),bound))
        print "flow_bounds(dynamic,",starting_box,",",suggested_time_step,")=",(time_step,bound),"mapping to",\
            starting_box+Interval(0,time_step)*dynamic(bound)
        return (time_step,bound)


    def flow_model(self,dynamic,starting_box,suggested_time_step):
        (time_step,flow_bounds)=self.flow_bounds(dynamic,starting_box,suggested_time_step)
        dynamic_model=TaylorFunction(flow_bounds,dynamic)
        return unchecked_flow(dynamic_model,starting_box,Interval(-time_step,+time_step),self.temporal_order)


    def integration_step(self,flow_model,starting_set_model,integration_time_model):
        return compose(flow_model,join(starting_set_model,integration_time_model))


    def reachability_step(self,flow_model,starting_set_model,zero_time_model, finishing_time_model):
        dimension=flow_model.result_size()
        print "reachability_step:"
        print "  ",starting_set_model
        #FIXME: Segfault in printing zero size model
        print "  ",zero_time_model
        print "  ",finishing_time_model
        assert(zero_time_model.domain()==starting_set_model.range())
        composed_zero_time_model=compose(zero_time_model,starting_set_model)
        composed_finishing_time_model=compose(finishing_time_model,starting_set_model)
        step_time_domain=Interval(0,1)
        step_time_model=TaylorVariable.identity(Interval(0,1))
        embedded_step_time_model=embed(step_time_model,finishing_time_model.domain())
        embedded_zero_time_model=embed(composed_zero_time_model,step_time_domain)
        embedded_finishing_time_model=embed(composed_finishing_time_model,step_time_domain)
        # The embedded, composed reachability_time_model is the function giving how far the point with paramter s has to flow
        embedded_reachability_time_model=embedded_zero_time_model \
            +embedded_step_time_model*(embedded_finishing_time_model-embedded_zero_time_model)
        embedded_starting_set_model=embed(starting_set_model,step_time_domain)
        return compose(flow_model,join(starting_set_model,reachability_time_model))


    def compute_touching_time_interval(self,guard,flow_model):
        """Compute a time interval such that any touching of any point of the initial
           set with the guard in [0,time_step] occurs in this interval."""
        guard_flow=compose(guard,flow_model)
        
        # Test for the case of clearly no crossing
        if not self.possibly_touching(guard_flow,flow_model.domain()):
            print "No touching"
            return None
        
        space_domain=flow_model.domain()[0:dimension]
        time_domain=Interval(0,flow_model.domain()[dimension].upper())
        
        # Compute a lower bound for the crossing time
        lower_time=0.0
        upper_time=flow_model.domain()[dimension].upper()
        if self.possibly_touching(guard_flow,join(space_domain,Interval(lower_time))):
            lower_touching_time=lower_time
        else:
            upper_time=(upper_time+lower_time)/2
            for i in range(0,9):
                if self.possibly_touching(guard_flow,join(space_domain,Interval(lower_time,upper_time))):
                    upper_time=(upper_time+lower_time)/2
                else:
                    (lower_time,upper_time)=(upper_time,2*upper_time-lower_time)
                if upper_time==flow_model.domain()[dimension].upper():
                    break
            lower_touching_time=lower_time
        
        # Compute an upper bound for the crossing time
        lower_time=0.0
        upper_time=flow_model.domain()[dimension].upper()
        if self.possibly_touching(guard_flow,join(space_domain,Interval(upper_time))):
            upper_touching_time=upper_time
        else:
            lower_time=(upper_time+lower_time)/2
            for i in range(0,9):
                if self.possibly_touching(guard_flow,join(space_domain,Interval(lower_time,upper_time))):
                    lower_time=(upper_time+lower_time)/2
                else:
                    (lower_time,upper_time)=(2*lower_time-upper_time,lower_time)
            upper_touching_time=upper_time
        
        if(lower_touching_time>=upper_touching_time):
            return None
        else:
            return Interval(lower_touching_time,upper_touching_time)


    def compute_transverse_crossing_time(self,guard,flow_model,touching_time_interval=None):
        """Compute transverse crossing time using the implicit function theorem.
           First perform some approximate steps, and then try to improve error."""
        dimension=guard.argument_size()
        guard_flow=compose(guard,flow_model)
        guard_flow_gradient=derivative(guard_flow,dimension)
        approx_guard_flow=midpoint(guard_flow)
        print "guard_flow =",guard_flow.polynomial()
        print "guard_flow_gradient =",guard_flow_gradient.polynomial(),guard_flow_gradient.range()
        
        space_domain=flow_model.domain()[0:dimension]
        if touching_time_interval:
            time_domain=touching_time_interval
        else:
            time_domain=Interval(0,flow_model.domain()[dimension].upper())
        id=TaylorFunction.identity(space_domain)
        tau=TaylorExpression.constant(space_domain,0.0)
        tau+=time_domain
        
        for i in range(0,7):
            tau_midpoint=midpoint(tau)
            guard_at_time_tau=unchecked_compose(guard_flow,join(id,tau_midpoint))
            gradient_at_time_tau=unchecked_compose(guard_flow_gradient,join(id,tau_midpoint))
            dtau=guard_at_time_tau/gradient_at_time_tau
            tau=tau_midpoint-dtau
        
        tau+=16*tau.error()*Interval(-1,+1)
        tau_midpoint=midpoint(tau)
        guard_at_time_tau=unchecked_compose(guard_flow,join(id,tau_midpoint))
        gradient_at_time_tau=unchecked_compose(guard_flow_gradient,join(id,tau))
        new_tau=tau_midpoint-guard_at_time_tau/gradient_at_time_tau
        assert(refines(new_tau,tau))
        return new_tau


    def compute_crossing_time(self,guard,flow_model):
        print "compute_crossing_time"
        # Compute an approximate crossing time
        guard_flow=compose(guard,flow_model)[0]
        
        dimension=guard.argument_size()
        space_domain=flow_model.domain()[0:dimension]
        step_size=flow_model.domain()[dimension].upper()
        
        domain=join(space_domain,Interval(0,step_size))
        guard_range=guard_flow.evaluate(domain)
        print "guard_range =",guard_range
        if not possibly(guard_range.is_zero()):
            return None
        time_domain=Interval(0)
        initial_guard_range=guard_flow.evaluate(join(flow_space_domain,time_domain))
        time_domain=Interval(flow_time_domain.upper())
        final_guard_range=guard_flow.evaluate(join(flow_space_domain,time_domain))
        print "initial_guard_range =",initial_guard_range
        print "final_guard_range =",final_guard_range
        time_domain=Interval(0,step_size)
        initial_guard_range=guard_flow.evaluate(join(flow_space_domain,time_domain))


    def compute_lower_hitting_time(self,guard,flow_model):
        return self.lower_implicit(compose(guard,flow_model))


    def compute_upper_hitting_time(self,guard,flow_model):
        return self.upper_implicit(compose(guard,flow_model))


    def lower_implicit(self,f):
        """Compute a function h such that f(x,h(x))~0, and that f(x,h(x))*f(x,0)>=0."""
        dimension=f.argument_size()-1
        s1=2.0
        h_domain=join(f.domain()[0:dimension],Interval(0,s1))
        h0=TaylorExpression.constant(h_domain,0.0)
        h=h0
        id=TaylorFunction.identity(h_domain)[0:dimension]
        
        for i in range(0,9):
            h=antiderivative(-unchecked_compose(f,join(id,h)),dimension)+h0
            h.sweep(1e-8)
            print h,unchecked_compose(f,join(id,h))
        h.sweep(1e-8)
        
        h_domain=f.domain()[0:dimension]
        id=TaylorFunction.identity(h_domain)
        s1=TaylorExpression.constant(h_domain,s1)
        h=unchecked_compose(h,join(id,s1))
        return h


    def upper_implicit(self,f):
        """Compute a function h such that f(x,h(x))~0, and that f(x,h(x))*f(x,0)>=0."""
        dimension=f.argument_size()-1
        df=derivative(f,dimension)
        
        s1=0.125
        h_domain=join(f.domain()[0:dimension],Interval(0,s1))
        h0=TaylorExpression.constant(h_domain,0.0)
        h=h0
        id=TaylorFunction.identity(h_domain)[0:dimension]
        s1m=TaylorExpression.constant(h_domain,s1)

        print "\nComputing upper hitting time"

        #Flow forward to hitting time
        print h.error(),h.range(),unchecked_compose(f,join(id,h)).range()
        for i in range(0,20):
            for j in range(0,9):
                h=antiderivative(-unchecked_compose(f,join(id,h)),dimension)+h0
                h.sweep(1e-12)
                #print h,unchecked_compose(f,join(id,h))
            print h.error(),h.range(),unchecked_compose(f,join(id,h)).range()
            h.sweep(1e-12)
            h0=unchecked_compose(h,join(id,s1m))
        print

        # Flow a fixed time step
        h0+=3e-2
        h=h0
        
        #Flow a small step past hitting time
        print h.error(),h.range(),unchecked_compose(f,join(id,h)).range(),unchecked_compose(df,join(id,h)).range()
        for i in range(0,3):
            for j in range(0,9):
                h=antiderivative(unchecked_compose(df,join(id,h)),dimension)+h0
                h.sweep(1e-12)
            print h.error(),h.range(),unchecked_compose(f,join(id,h)).range(),unchecked_compose(df,join(id,h)).range()
            h.sweep(1e-12)
            h0=unchecked_compose(h,join(id,s1m))
        print
        
        #Flow backward to hitting time
        print h.error(),h.range(),unchecked_compose(f,join(id,h)).range()
        for i in range(0,3):
            for j in range(0,9):
                h=antiderivative(-unchecked_compose(f,join(id,h)),dimension)+h0
                h.sweep(1e-10)
                #print h,unchecked_compose(f,join(id,h))
            print h.error(),h.range(),unchecked_compose(f,join(id,h)).range()
            h.sweep(1e-10)
            h0=unchecked_compose(h,join(id,s1m))
        print

        h_domain=f.domain()[0:dimension]
        id=TaylorFunction.identity(h_domain)
        s1m=TaylorExpression.constant(h_domain,s1)
        print h.error()
        h=unchecked_compose(h,join(id,s1m))
        return h


    def crossing_time(self,guard,flow_model):
        touching_time_interval=self.touching_time_interval(guard,flow_model)
        if touching_time_interval==None:
            return None
        crossing_time_model=self.compute_crossing_time(guard,flow_model,touching_time_interval)


    def evolution_events(self, guards, flow_model, zero_time_model):
        print "evolution_events"
        print guards, flow_model
        
        # Output of function
        active_events={}
        
        time_step_model=zero_time_model+time_step
        flow_domain=flow_model.domain()
        crossing_times={ self.starting_event:zero_time_model, self.finishing_event:time_step_model }
        for (event,(guard,urgent)) in guards.items():
            gphi=compose(guard,flow_model)[0]
            gphi_range=gphi(flow_domain)
            print "guard range =",gphi_range
            if(gphi_range.lower()>0):
                # Guard is active over whole time interval
                print "\nguard",guard,"has value",gphi.range()," so is active over entire time interval"
            elif gphi_range.upper()<0:
                # Guard is not active
                pass
            else:
                touching_time=self.crossing_time(guard,flow_model)
                crossing_times[event]=touching_time
        print "crossing_times =",crossing_times
        
        ordered_crossing_times=[]
        for (event,time) in crossing_times.items():
            ordered_crossing_times.append((time,event))
        ordered_crossing_times.sort()
        
        print "ordered_crossing_times =",ordered_crossing_times
        
        if(len(ordered_crossing_times)==2):
            # Just starting and finishing pseudo-events
            return {"finishing_event":time_step_model,"blocking_event":time_step_model}


    def evolution_step(self, dynamic, guards, resets, starting_set, starting_time, time_step):
        # Change time to a model
        if(type(starting_time)==float):
            starting_time=starting_set[0]*0.0+starting_time
        
        zero_time_model=starting_time*0+0
        time_step_model=zero_time_model+time_step
        finishing_time_model=starting_time+time_step
        
        # If an urgent event is definitely initially_active, then there is
        # no continuous evolution.
        initially_active=self.initially_active(guards,starting_set)
        if initially_active["blocking_event"]==True:
            for event in initially_active.keys():
                if possibly(initially_active[event]):
                    active_events[event]=(starting_set,starting_time,None)
        
        # Compute the continuous evolution
        starting_box=starting_set.range()
        [time_step,flow_bounds]=self.flow_bounds(dynamic,starting_box,time_step)
        print "flow_bounds",flow_bounds
        dynamic_model=TaylorFunction(flow_bounds,dynamic)
        flow_time_interval=Interval(0,time_step)
        flow_model=unchecked_flow(dynamic_model,starting_box,Interval(-time_step,+time_step),self.temporal_order)
        flow_domain=join(starting_box,flow_time_interval)
        print "\nflow =",flow_model
        print "flow_domain =",flow_domain
        
        # Compute events
        event_times=self.evolution_events(guards,flow_model,starting_time*0.0)
        
        # Compute the blocking time
        blocking_time_model=time_step_model
        blocking_events={}
        for (event,event_time_model) in event_times.items():
            if guards.has_key(event) and guards[event][1]:
                if blocking_time_model>=event_time_model:
                    blocking_events={event:None}
                    blocking_time_model=event_time_model
                elif blocking_time_model<=event_time_model:
                    pass
                else:
                    blocking_events.insert(event)
                    blocking_time_model=min(blocking_time_model,event_time_model)
        
        #Apply events and continuous reach/evolve steps
        jumped_sets=[]
        print "events =",event_times
        for tup in event_times.items():
            print "tup=",tup
            (event,crossing_time_model)=tup
            
            if event=="blocking_event":
                reachable_set=self.reachability_step(flow_model,initial_set,zero_time_model,crossing_time_model)
            elif event=="finishing_event":
                evolved_set=self.integration_step(flow_model,initial_set,crossing_time_model)
            else:
                print "Event ",event," from time ",crossing_time_model,"to time",other_time
                reset=resets[events]
                assert(other_time==None)
                active_set=self.integration_step(flow_model,initial_set,crossing_time_model)
                jump_set=apply(reset,jump_set)
                jumped_sets.append(jump_set)
        
        print "reachable_set =",reachable_set
        print "evolved_set =",evolved_set
        print "jumped_sets =",jumped_sets
        
        jumped_sets.append(evolved_set)
        return ((reachable_set,),jumped_sets)


class Mode:
    def __init__(self,dynamic,guard,reset=None):
        self.dynamic=dynamic
        self.guard=guard
        self.reset=reset
        assert(isinstance(dynamic,FunctionInterface))
        assert(isinstance(guard,ExpressionInterface))
        assert(dynamic.result_size()==dynamic.argument_size())
        assert(guard.argument_size()==dynamic.argument_size())


def run_example(dynamic,guard,domain,step_size=0.25):
    evolver=HybridEvolver()
    
    initial_set=TaylorFunction.identity(domain)
    space_domain=domain
    flow_model=evolver.flow_model(dynamic,domain,step_size)
    time_domain=flow_model.domain()[dimension]
    print "\ninitial_set =",initial_set
    print "\nflow_model =",flow_model
    final_time=TaylorExpression.constant(space_domain,time_domain.upper())
    print "\nfinal_time =",final_time
    final_set=unchecked_compose(flow_model,join(initial_set,final_time))
    print "\nfinal_set =",final_set
    initial_guard=compose(guard,initial_set)
    print "\ninitial_guard =",initial_guard.range()
    final_guard=compose(guard,final_set)
    print "\nfinal_guard =",final_guard.range()
    
    # Compute the gradient function; only works for affine functions at the moment
    guard_gradient=dynamic[0]*derivative(guard,0)
    for i in range(1,dimension):
        guard_gradient=guard_gradient+dynamic[i]*derivative(guard,i)
    print "\nguard_gradient =",guard_gradient
    
    touching_time_interval=evolver.compute_touching_time_interval(guard,flow_model)
    print "\ntouching_time_interval =",touching_time_interval
    crossing_time_model=evolver.compute_transverse_crossing_time(guard,flow_model)
    print "\ncrossing_time_model =",crossing_time_model
    lower_hitting_time_model = evolver.compute_lower_hitting_time(guard,flow_model)
    print "\nlower_hitting_time_model =",lower_hitting_time_model
    upper_hitting_time_model = evolver.compute_upper_hitting_time(guard,flow_model)
    print "\nupper_hitting_time_model =",upper_hitting_time_model
    print "\ncrossing_time_range =",crossing_time_model.range()
    print "lower_hitting_time_range =",lower_hitting_time_model.range()
    print "upper_hitting_time_range =",upper_hitting_time_model.range()
    print "(crossing_time_model-lower_hitting_time_model).range() =",(crossing_time_model-lower_hitting_time_model).range()
    print "(crossing_time_model-upper_hitting_time_model).range() =",(crossing_time_model-upper_hitting_time_model).range()
    sys.exit()
    lower_hitting_time=evolver.compute_lower_hitting_time(guard,flow_model)
    print "\ncrossing_time =",crossing_time_model
    print "\nlower_hitting_time =",lower_hitting_time
    print initial_set,lower_hitting_time
    print initial_set.domain(),lower_hitting_time.domain()
    print initial_set.domain()[0].lower()-lower_hitting_time.domain()[0].lower()
    print initial_set.domain()==lower_hitting_time.domain()
    print "\ng(lower_hitting_time) =",compose(guard,join(initial_set,lower_hitting_time))[0].range()
    sys.exit()
    
    (reachable_sets,evolved_sets)=evolver.evolution_step(dynamic,{"hit":(guard,True)},["hit",reset],initial_set,initial_time,time_step)
    print "reachable_sets =",reachable_sets
    print "evolved_sets =",evolved_sets




if __name__=="__main__":
    g=9.8 # gravity
    a=0.5 # restitution
    dimension=2
    d=IVector([[-2,3],[-4,5]])
    x=TaylorExpression.variable(d,0)
    y=TaylorExpression.variable(d,1)
    print x+1.2345667e-12*x*x+1.23241312e-13*y*x
    #sys.exit()
    
    ball_dynamic=AffineFunction(Matrix([[0,1],[0,0]]),Vector([0,-g]))
    ball_guard=AffineFunction(Matrix([[-1,0]]),Vector([0]))
    ball_guard=AffineExpression(Vector([-1,0]),0)
    ball_reset=AffineFunction(Matrix([[1,0],[0,-a]]),Vector([0,0]))
    ball_domain=IVector([[0.1,0.11],[-0.5,-0.49]])
    print "ball_dynamic =",ball_dynamic,"\nball_guard =",ball_guard
    run_example(ball_dynamic,ball_guard,ball_domain)
    
    tangency_dynamic=AffineFunction(Matrix([[0,0],[-2,0]]),Vector([1,0]))
    tangency_guard=AffineExpression(Vector([0,1]),0.0)
    tangency_domain=IVector([[-0.25,-0.125],[-3/32.,-1/32.]])
    print "tangency_dynamic =",tangency_dynamic,"\ntangency_guard =",tangency_guard
    run_example(tangency_dynamic,tangency_guard,tangency_domain)

