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
Float = float

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


class HybridEvolverPrototype:

    temporal_order=4;
    spacial_order=2;

    starting_event="starting_event"
    finishing_event="finishing_event"
    blocking_event=-3
    invariant_event=-3

    def __init__(self):
        pass

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
            guard_range=guard_range.range()
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


    def approxmate_newton_implicit_step(self,f,h,order=None):
        n=f.argument_size()-1
        id=TaylorFunction.identity(h.domain())
        df=derivative(f,n)
        df.clobber()
        mh=midpoint(h)
        fmh=unchecked_compose(f,join(id,mh))
        fmh.clobber()
        dfh=unchecked_compose(df,join(id,h))
        if(order): dfh.truncate(order)
        dfh.clobber()
        try:
            rdfh=rec(dfh)
        except:
            rdfh=rec(dfh.value())
        nh=(mh-rdfh*fmh)
        nh.clobber()
        #print nh.error(),nh.range(),unchecked_compose(f,join(id,nh)).range()
        return nh

    def newton_implicit_step(self,f,h,order=None):
        n=f.argument_size()-1
        id=TaylorFunction.identity(h.domain())
        df=derivative(f,n)
        mh=midpoint(h)
        fmh=unchecked_compose(f,join(id,mh))
        dfh=unchecked_compose(df,join(id,h))
        if(order): dfh.truncate(order)
        nh=mh-rec(dfh)*fmh
        #print nh.error(),nh.range(),unchecked_compose(f,join(id,nh)).range()
        return nh


    def krawczyk_implicit_step(self,f,h,order=1):
        n=f.argument_size()-1
        id=TaylorFunction.identity(h.domain())
        df=derivative(f,n)
        mh=midpoint(h)
        fmh=unchecked_compose(f,join(id,mh))
        dfh=unchecked_compose(df,join(id,h))
        mdfh=midpoint(dfh.truncate(order))
        rmdfh=rec(mdfh)
        nh=mh-rmdfh*fmh+(1-rmdfh*dfh)*(h-mh)
        return nh


    def newton_implicit(self,f):
        n=f.argument_size()-1
        df=derivative(f,n)
        h_domain=f.domain()[0:n]
        h0=TaylorExpression.constant(h_domain,0.0)
        id=TaylorFunction.identity(h_domain)

        h=h0
        for i in range(0,9):
            h=self.newton_implicit_step(f,h)
        return h


    def krawczyk_implicit(self,f,order=1):
        n=f.argument_size()-1
        df=derivative(f,n)
        h_domain=f.domain()[0:n]
        h0=TaylorExpression.constant(h_domain,0.0)
        id=TaylorFunction.identity(h_domain)

        h=h0
        for i in range(0,9):
            self.krawczyk_implicit_step(f,h,order)
        return h


    def lower_implicit(self,f):
        """Compute a function h such that f(x,h(x))~0, and that f(x,h(x))*f(x,0)>=0."""
        dimension=f.argument_size()-1
        s1=0.125
        h_domain=join(f.domain()[0:dimension],Interval(0,s1))
        h0=TaylorExpression.constant(h_domain,0.0)
        h=h0
        id=TaylorFunction.identity(h_domain)[0:dimension]
        s1m=TaylorExpression.constant(h_domain,s1)

        #print "\nComputing lower implicit function approximation"
        #print h.error(),h.range(),unchecked_compose(f,join(id,h)).range()
        for i in range(0,20):
            for j in range(0,9):
                h=antiderivative(-unchecked_compose(f,join(id,h)),dimension)+h0
                h.sweep(1e-12)
            #print "  ",h.error(),h.range(),unchecked_compose(f,join(id,h)).range()
            h.sweep(1e-12)
            h-=h.error()
            h.clobber()
            h0=unchecked_compose(h,join(id,s1m))
        #print

        h_domain=f.domain()[0:dimension]
        id=TaylorFunction.identity(h_domain)
        s1=TaylorExpression.constant(h_domain,s1)
        h=unchecked_compose(h,join(id,s1))
        h-=h.error()
        h.clobber()
        #print h.error(),h.range(),unchecked_compose(f,join(id,h)).range()
        return h

    def lower_implicit(self,f):
        """Compute a function h such that f(x,h(x))~0, and that f(x,h(x))*f(x,0)>=0."""
        dimension=f.argument_size()-1
        s=0.75
        h_domain=f.domain()[0:dimension]
        h0=TaylorExpression.constant(h_domain,0.0)
        h=h0
        id=TaylorFunction.identity(h_domain)

        #print "\nComputing lower implicit function approximation"
        #print h.error(),h.range(),unchecked_compose(f,join(id,h)).range()
        fh=unchecked_compose(f,join(id,h))
        for i in range(0,20):
            s*=4; s/=3
            nh=h-s*fh
            nh.clobber()
            nfh=unchecked_compose(f,join(id,nh))
            while nfh.range().upper()>0:
                s*=.75
                nh=h-s*fh
                nh.clobber()
                nfh=unchecked_compose(f,join(id,nh))
                #print s,nfh.range()
            h=nh
            fh=nfh
            #print "  ",s,h.error(),h.range(),fh.range()
        #print

        return h


    def upper_implicit(self,f):
        """Compute a function h such that f(x,h(x))~0, and that f(x,h(x))*f(x,0)>=0."""
        dimension=f.argument_size()-1
        df=derivative(f,dimension)

        s1=0.5
        h_domain=join(f.domain()[0:dimension],Interval(0,s1))
        h0=TaylorExpression.constant(h_domain,0.0)
        h=h0
        id=TaylorFunction.identity(h_domain)[0:dimension]
        s1m=TaylorExpression.constant(h_domain,s1)

        #print "\nComputing upper implicit function approximation"

        #Flow forward to hitting time
        #print h.error(),h.range(),unchecked_compose(f,join(id,h)).range()
        for i in range(0,0):
            for j in range(0,9):
                h=antiderivative(-unchecked_compose(f,join(id,h)),dimension)+h0
                h.sweep(1e-12)
                #print h,unchecked_compose(f,join(id,h))
            print "  ",h.error(),h.range(),unchecked_compose(f,join(id,h)).range()
            h+=3e-2
            h.sweep(1e-12)
            h+=h.error()
            h.clobber()
            h0=unchecked_compose(h,join(id,s1m))
        #print

        # Flow a fixed time step
        #print " ",h.error(),h.range(),unchecked_compose(f,join(id,h)).range(),unchecked_compose(df,join(id,h)).range()
        for i in range(0,0):
            h0+=3e-2
            h=h0
            #print "  ",h.error(),h.range(),unchecked_compose(f,join(id,h)).range(),unchecked_compose(df,join(id,h)).range()
        #print

        #Flow a small step past hitting time
        #print " ",h.error(),h.range(),unchecked_compose(f,join(id,h)).range(),unchecked_compose(df,join(id,h)).range()
        for i in range(0,9):
            for j in range(0,1):
                h=antiderivative(unchecked_compose(df,join(id,h)),dimension)+h0
                h.sweep(1e-12)
            #print "  ",h.error(),h.range(),unchecked_compose(f,join(id,h)).range(),unchecked_compose(df,join(id,h)).range()
            h.sweep(1e-12)
            h+=h.error()
            h.clobber()
            h0=unchecked_compose(h,join(id,s1m))
        #print

        #Flow backward to hitting time
        #print " ",h.error(),h.range(),unchecked_compose(f,join(id,h)).range()
        for i in range(0,18):
            for j in range(0,1):
                h=antiderivative(-unchecked_compose(f,join(id,h)),dimension)+h0
                h.sweep(1e-10)
                #print h,unchecked_compose(f,join(id,h))
            #print "  ",h.error(),h.range(),unchecked_compose(f,join(id,h)).range()
            h.sweep(1e-10)
            h+=h.error()
            h.clobber()
            h0=unchecked_compose(h,join(id,s1m))
        #print

        h_domain=f.domain()[0:dimension]
        id=TaylorFunction.identity(h_domain)
        s1m=TaylorExpression.constant(h_domain,s1)
        h=unchecked_compose(h,join(id,s1m))
        h+=h.error()
        h.clobber()
        #print h.error(),h.range(),unchecked_compose(f,join(id,h)).range()
        return h


    def flow_bounds(self,dynamic, starting_box, suggested_time_step):
        assert(isinstance(dynamic,FunctionInterface))
        time_step=suggested_time_step
        small_neighbourhood=starting_box+(starting_box-midpoint(starting_box))/16
        one_box_neighbourhood=starting_box+(starting_box-midpoint(starting_box))
        assert(subset(starting_box,one_box_neighbourhood))
        assert(subset(small_neighbourhood,small_neighbourhood))
        #print "Computing flow bounds"
        #print starting_box,small_neighbourhood,one_box_neighbourhood
        bound=one_box_neighbourhood+dynamic(one_box_neighbourhood)*(4*Interval(0,suggested_time_step))
        bound=one_box_neighbourhood+dynamic(bound)*(2*Interval(0,suggested_time_step))
        while(not subset(small_neighbourhood+time_step*dynamic(bound),bound)):
            time_step/=2
            assert(time_step>1e-12)
        bound=small_neighbourhood+Interval(0,time_step)*dynamic(bound)
        bound=small_neighbourhood+Interval(0,time_step)*dynamic(bound)
        assert(subset(starting_box+time_step*dynamic(bound),bound))
        #print "flow_bounds(dynamic,",starting_box,",",suggested_time_step,")=",(time_step,bound),"mapping to",\
        #    starting_box+Interval(0,time_step)*dynamic(bound)
        return (time_step,bound)


    def flow_model(self,dynamic,starting_box,suggested_time_step):
        (time_step,flow_bounds)=self.flow_bounds(dynamic,starting_box,suggested_time_step)
        dynamic_model=TaylorFunction(flow_bounds,dynamic)
        return unchecked_flow(dynamic_model,starting_box,Interval(-time_step,+time_step),self.temporal_order)


    def integration_step(self,flow_model,starting_set_model,integration_time_model):
        return compose(flow_model,join(starting_set_model,integration_time_model))


    def reachability_step(self,flow_model,starting_set_model,zero_time_model, finishing_time_model):
        dimension=flow_model.result_size()
        #print "reachability_step:"
        #print "  ",starting_set_model
        #FIXME: Segfault in printing zero size model
        #print "  ",zero_time_model
        #print "  ",finishing_time_model
        #print zero_time_model.domain(),starting_set_model.domain()
        #assert(zero_time_model.domain()==starting_set_model.domain())
        composed_zero_time_model=compose(zero_time_model,starting_set_model)
        composed_finishing_time_model=compose(finishing_time_model,starting_set_model)
        step_time_domain=Interval(0,1)
        step_time_model=TaylorFunction.identity(IVector([Interval(0,1)]))[0]
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
        #print "guard_flow =",guard_flow.polynomial()
        #print "guard_flow_gradient =",guard_flow_gradient.polynomial(),guard_flow_gradient.range()

        space_domain=flow_model.domain()[0:dimension]
        if touching_time_interval:
            time_domain=touching_time_interval
        else:
            time_domain=Interval(0,flow_model.domain()[dimension].upper())
        id=TaylorFunction.identity(space_domain)
        tau=TaylorExpression.constant(space_domain,0.0)
        tau+=time_domain

        self.implicit_step=self.krawczyk_implicit_step;
        self.implicit_step=self.newton_implicit_step;
        self.approximate_implicit_step=self.approxmate_newton_implicit_step;

        for i in range(0,7):
            tau=self.approximate_implicit_step(guard_flow,tau)
        tau=self.implicit_step(guard_flow,tau)

        tau+=16*tau.error()*Interval(-1,+1)
        new_tau=self.implicit_step(guard_flow,tau)
        assert(refines(new_tau,tau))
        return new_tau


    def compute_lower_hitting_time(self,guard,flow_model):
        return self.lower_implicit(compose(guard,flow_model))


    def compute_upper_hitting_time(self,guard,flow_model):
        return self.upper_implicit(compose(guard,flow_model))



    def crossing_time(self,guard,flow_model):
        domain=flow_model.domain()[0:flow_model.result_size()]
        id=TaylorFunction.identity(domain)

        touching_time_interval=self.compute_touching_time_interval(guard,flow_model)
        if touching_time_interval==None:
            return None

        try:
            crossing_time_model=self.compute_transverse_crossing_time(guard,flow_model,touching_time_interval)
        except:
            crossing_time_model=TaylorExpression.constant(domain,touching_time_interval)

        print
        print "  crossing_time_error =",crossing_time_model.error()
        print "  crossing_time_range =",crossing_time_model.range()
        print "  crossing_time_guard_range =",unchecked_compose(compose(guard,flow_model),join(id,crossing_time_model)).range()

        lower_hitting_time_model=self.compute_lower_hitting_time(guard,flow_model)
        print "  lower_hitting_time_range =",lower_hitting_time_model.range()
        print "  lower_hitting_time_guard_range =", unchecked_compose(compose(guard,flow_model),join(id,lower_hitting_time_model)).range()

        upper_hitting_time_model=self.compute_upper_hitting_time(guard,flow_model)
        print "  upper_hitting_time_range =",upper_hitting_time_model.range()
        print "  upper_hitting_time_guard_range =", unchecked_compose(compose(guard,flow_model),join(id,upper_hitting_time_model)).range()
        print

        #print "(upper_hitting_time-lower_hitting_time).range() =",(upper_hitting_time_model-lower_hitting_time_model).range()
        #print "(crossing_time_model-lower_hitting_time_model).range() =",(crossing_time_model-lower_hitting_time_model).range()
        #print "(crossing_time_model-upper_hitting_time_model).range() =",(crossing_time_model-upper_hitting_time_model).range()

        return crossing_time_model


    def evolution_events(self, guards, flow_model, zero_time_model):
        #print "evolution_events"
        print
        print "guards =",guards

        # Output of function
        active_events={}
        dimension=flow_model.result_size()
        time_step=flow_model.domain()[dimension].upper()
        time_step_model=zero_time_model+time_step
        flow_domain=flow_model.domain()
        crossing_times={ self.starting_event:zero_time_model, self.finishing_event:time_step_model }
        for (event,(guard,urgent)) in guards.items():
            crossing_time_model=self.crossing_time(guard,flow_model)
            if crossing_time_model:
                crossing_times[event]=crossing_time_model

        ordered_crossing_times=[]
        for (event,time) in crossing_times.items():
            ordered_crossing_times.append((time,event))
        ordered_crossing_times.sort()

        print "ordered_events = [",
        for (crossing_time,event) in ordered_crossing_times:
            print str(event)+":"+str(crossing_time.range())+",",
        print "]"

        if(len(ordered_crossing_times)==2):
            # Just starting and finishing pseudo-events
            return {"finishing_event":time_step_model,"blocking_event":time_step_model}

        return ordered_crossing_times


    def evolution_step(self, dynamic, guards, resets, starting_set, starting_time, time_step):
        print "\nevolution_step:","t:",starting_time.range(),
        print "c:",starting_set.centre(), "r:",mag(norm(starting_set.range()-starting_set.centre()))


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
        print
        print "flow_bounds",flow_bounds
        dynamic_model=TaylorFunction(flow_bounds,dynamic)
        flow_time_interval=Interval(0,time_step)
        flow_model=unchecked_flow(dynamic_model,starting_box,Interval(-time_step,+time_step),self.temporal_order)
        flow_domain=join(starting_box,flow_time_interval)
        #print "\nflow =",flow_model
        print "flow_domain =",flow_model.domain()
        print "flow_polynomial =",flow_model.polynomial()

        # Compute events
        ordered_event_times=self.evolution_events(guards,flow_model,starting_time*0.0)

        # Compute the blocking time
        blocking_time_model=time_step_model
        blocking_events={}
        for (event_time_model,event) in ordered_event_times:
            if guards.has_key(event) and guards[event][1]:
                if blocking_time_model>=event_time_model:
                    blocking_events={event:None}
                    blocking_time_model=event_time_model
                elif blocking_time_model<=event_time_model:
                    pass
                else:
                    blocking_events.insert(event)
                    blocking_time_model=min(blocking_time_model,event_time_model)
        print
        print "blocking events:",blocking_events
        print "blocking_time_range:",blocking_time_model.range()
        print
        ordered_event_times.append((blocking_time_model,"blocking_event"))

        #Apply events and continuous reach/evolve steps
        jumped_sets=[]
        for (crossing_time_model,event) in ordered_event_times:

            if event=="blocking_event":
                reachable_set=self.reachability_step(flow_model,starting_set,zero_time_model,crossing_time_model)
            elif event=="finishing_event":
                evolved_set=self.integration_step(flow_model,starting_set,crossing_time_model)
            elif event=="starting_event":
                pass
            else:
                assert(guards[event][1]) # Urgent
                print "Event ",event," at time ",crossing_time_model.range()
                reset=resets[event]
                active_set=self.integration_step(flow_model,starting_set,crossing_time_model)
                jump_set=compose(reset,active_set)
                jumped_sets.append((event,jump_set))

        print "reachable_set =",reachable_set.range()
        print "evolved_set =",evolved_set.range()
        jumped_set_ranges=[]
        for (event,set) in jumped_sets:
            jumped_set_ranges.append((event,set.range()))
        print "jumped_sets =",jumped_set_ranges
        sys.exit()

        jumped_sets.append(evolved_set)
        return ((reachable_set,),jumped_sets)

    def orbit(self, system, initial_set, evolution_time, semantics):
        """Returns a tuple of (reach,evolve,intermediate) sets"""
        reach=HybridTaylorSetList([])
        evolve=HybridTaylorSetList([])
        intermediate=HybridTaylorSetList([])

        print "HybridEvolverPrototype.orbit(system,initial_set,evolution_time,semantics) is not implemented"
        return (reach,evolve,intermediate)

class Mode:
    def __init__(self,dynamic,guard,reset=None):
        self.dynamic=dynamic
        self.guard=guard
        self.reset=reset
        assert(isinstance(dynamic,FunctionInterface))
        assert(isinstance(guard,ExpressionInterface))
        assert(dynamic.result_size()==dynamic.argument_size())
        assert(guard.argument_size()==dynamic.argument_size())


def run_example(dynamic,guard,reset,domain,step_size=0.25):
    evolver=HybridEvolver()

    initial_set=TaylorFunction.identity(domain)
    initial_time=TaylorExpression.constant(domain,0.0)

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

    #touching_time_interval=evolver.compute_touching_time_interval(guard,flow_model)
    #print "\ntouching_time_interval =",touching_time_interval
    #crossing_time_model=evolver.compute_transverse_crossing_time(guard,flow_model)
    #print "\ncrossing_time_model =",crossing_time_model
    lower_hitting_time_model = evolver.compute_lower_hitting_time(guard,flow_model)
    #print "\nlower_hitting_time_model =",lower_hitting_time_model
    upper_hitting_time_model = evolver.compute_upper_hitting_time(guard,flow_model)
    #print "\nupper_hitting_time_model =",upper_hitting_time_model
    id=TaylorFunction.identity(space_domain)

    (reachable_sets,evolved_sets)=evolver.evolution_step(dynamic,{"hit":(guard,True)},{"hit":reset},initial_set,initial_time,step_size)
    print "reachable_sets =",reachable_sets
    print "evolved_sets =",evolved_sets

def evolve(automaton,initial,time,semantics):
    print "Embedded Python evolver"
    print automaton,initial,time,semantics
    final=HybridTaylorSetList([])
    reach=HybridTaylorSetList([])
    intermediate=HybridTaylorSetList([])
    print "Returning result"
    return (final,reach,intermediate)


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
    ball_guard=AffineExpression(Vector([-1,0]),Float(0.0))
    ball_reset=AffineFunction(Matrix([[1,0],[0,-a]]),Vector([0,0]))
    ball_domain=IVector([[0.1,0.11],[-0.5,-0.49]])
    print "ball_dynamic =",ball_dynamic,"\nball_guard =",ball_guard
    #run_example(ball_dynamic,ball_guard,ball_domain)

    tangency_dynamic=AffineFunction(Matrix([[0,0],[-2,0]]),Vector([1,0]))
    tangency_guard=AffineExpression(Vector([0,1]),Float(0.0))
    tangency_reset=AffineFunction(Matrix([[1,0],[0,1]]),Vector([0,-2]))
    tangency_domain=IVector([[-0.25,-0.125],[-3/32.,-1/32.]])
    tangency_domain=IVector([[-0.5,-0.375],[-3/32.,-1/32.]])
    print "tangency_dynamic =",tangency_dynamic,"\ntangency_guard =",tangency_guard
    run_example(tangency_dynamic,tangency_guard,tangency_reset,tangency_domain)

