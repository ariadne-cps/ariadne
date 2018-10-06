#!/usr/bin/python

##############################################################################
#            constrained_image_set.py
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

import sys,__builtin__
from ariadne import *
#from hybrid_system import *

Indeterminate=indeterminate()
Float = float


class ComposedFunction:
    def __init__(self,first,second):
        print first,second
        assert first.argument_size() == second.result_size()
        self.first=first
        self.second=second
    def argument_size(self):
        return self.second.argument_size()
    def __call__(self,arguments):
        assert len(arguments) >= self.argument_size()
        arguments=arguments[:self.argument_size()]
        return self.first(self.second(arguments))
    def __str__(self):
        if hasattr(self,"codomain"):
            return str(self.first)+","+str(self.codomain)
        else:
            return str(self.first)
    def __repr__(self):
        return str(self)

class ConstrainedImageSet:
    def __init__(self,domain=None,function=None):
        # Default constructor computes empty set
        if domain==None and function==None:
            self._empty=True

        #Copy constructor
        if function==None and isinstance(domain,ConstrainedImageSet):
            other=domain
            self._empty=other._empty
            self._domain=Box(other._domain)
            self._function=VectorTaylorFunction(other._function)
            self._positive=dict(other._positive)
            self._zero=dict(other._zero)
            self._negative=dict(other._negative)
            self._always_negative=dict(other._always_negative)
            self._always_negative_and_zero=dict(other._always_negative)
            return

        #Construct from a TaylorVectorFunction
        if function==None:
            assert isinstance(domain,VectorTaylorFunction)
            model=domain
            self._domain=Box(model.domain())
            self._function=model
        else:
            #Standard constructor
            if isinstance(domain,IntervalVector):
                domain=Box(domain)
            assert isinstance(domain,Box)
            if isinstance(function,VectorFunction):
                function=VectorTaylorFunction(domain,function)
            assert isinstance(function,VectorTaylorFunction)
            assert function.domain()==domain
            self._domain=domain
            self._function=function
        self._empty=None
        self._positive={}
        self._zero={}
        self._negative={}
        self._always_negative={}
        self._always_negative_and_zero={}
        self._restrictions=[GridTreeSet(Grid(self._domain.size()))]

    def compute_range(self,constraint):
        composed_constraint=compose(constraint,self._function)
        constraint_range=composed_constraint(self._domain)
        return constraint_range

    def empty(self):
        if self._empty == None:
            self._empty = self.grid_outer_approximation(Grid(self.dimension()),6).empty()
        return self._empty

    def dimension(self):
        return self._function.result_size()

    def domain(self):
        return self._domain

    def function(self):
        return self._function

    def bounding_box(self):
        return Box(self._function(self._domain)).bounding_box()

    def disjoint(self,box):
        return self.__disjoint(self._domain,box,box.radius())

    def apply_map(self,map):
        self._function=compose(map,self._function)

    def apply_flow(self,flow):
        time_domain=Interval(flow.domain()[-1])
        self._domain=join(self._domain,time_domain)
        self._function=unchecked_compose(flow,combine(self._function,ScalarTaylorFunction.variable([time_domain],0)))

    def new_activation(self,event,constraint,derivative):
        assert isinstance(constraint,ScalarFunction)
        composed_constraint=ComposedFunction(constraint,self._function)
        constraint_range=composed_constraint(self._domain)
        print str(event)+".range():",constraint_range
        if constraint_range.lower()>=0.0:
            self._positive[event]=composed_constraint*0.0+1.0
        elif constraint_range.upper()<=0.0:
            self._positive[event]=composed_constraint*0.0-1.0
            self._empty=True
        else:
            self._positive[event]=composed_constraint

    def new_guard(self,event,constraint,derivative):
        assert isinstance(constraint,ScalarFunction)
        composed_constraint=ComposedFunction(constraint,self._function)
        constraint_range=composed_constraint(self._domain)
        print str(event)+".range():",constraint_range
        if constraint_range.lower()>0.0:
            self._always_negative_and_zero[event]=composed_constraint*0.0+1.0
            self._empty=True
        elif constraint_range.upper()<0.0:
            pass
            #self._always_negative_and_zero[event]=composed_constraint*0.0-1.0
        else:
            composed_derivative=compose(derivative,self._function)
            derivative_range=composed_derivative(self._domain)
            if derivative_range.lower()>0.0 and False:
                print "Solving for transverse crossing"
                solver=KrawczykSolver(1e-8,4)
                solved_function=solver.implicit(composed_derivative.function(),self._domain[:-1],self._domain[-1])
                print solved_function
                raw_input("")
                self._negative[event]=composed_constraint
            else:
                self._always_negative_and_zero[event]=composed_constraint

    def new_invariant(self,event,constraint,derivative):
        assert isinstance(constraint,ScalarFunction)
        composed_constraint=ComposedFunction(constraint,self._function)
        constraint_range=composed_constraint(self._domain)
        print str(event)+".range():",constraint_range
        if constraint_range.lower()>0.0:
            self._empty=True
        elif constraint_range.upper()<0.0:
            pass
        else:
            composed_derivative=compose(derivative,self._function)
            derivative_range=composed_derivative(self._domain)
            if derivative_range.lower()>0.0:
                # Transverse crossing; invariant holds where negative
                self._negative[event]=composed_constraint
            else:
                self._always_negative[event]=composed_constraint

    def new_equation(self,event,constraint):
        print "new_equation"
        assert isinstance(constraint,ScalarFunction)
        composed_constraint=ComposedFunction(constraint,self._function)
        constraint_range=composed_constraint(self._domain)
        print str(event)+".range():",constraint_range
        if constraint_range.upper()<0.0:
            self._zero[event]=composed_constraint*0.0-1.0
            self._empty=True
        elif constraint_range.lower()>0.0:
            self._zero[event]=composed_constraint*0.0-1.0
            self._empty=True
        else:
            self._zero[event]=composed_constraint

    def new_raw_equation(self,event,constraint):
        print "new_raw_equation"
        assert isinstance(constraint,ScalarFunction)
        assert constraint.argument_size() <= len(self._domain)
        constraint_range=constraint(self._domain)
        print str(event)+".range():",constraint_range
        if constraint_range.upper()<0.0:
            self._zero[event]=constraint*0.0-1.0
            self._empty=True
        elif constraint_range.lower()>0.0:
            self._zero[event]=constraint*0.0+1.0
            self._empty=True
        else:
            self._zero[event]=constraint

    def new_raw_invariant(self,event,constraint):
        print "new_raw_invariant"
        assert isinstance(constraint,ScalarFunction)
        assert constraint.argument_size() <= len(self._domain)
        constraint_range=constraint(self._domain)
        print str(event)+".range():",constraint_range
        if constraint_range.upper()<0.0:
            pass
        elif constraint_range.lower()>0.0:
            self._negative[event]=constraint*0.0+1.0
            self._empty=True
        else:
            self._negative[event]=constraint

    def new_termination(self,event,constraint):
        self.new_raw_invariant(event,constraint)



    def __disjoint(self,domain,box,err):
        if False: #positive_constraints(domain)
            return True
        image=Box(self._function(domain))
        if disjoint(image,box):
            return True;
        if subset(image,box):
            return False
        if image.radius()<err:
            return Indeterminate
        (subdomain1,subdomain2)=split(domain)
        return self.__disjoint(subdomain1,box,err) and self.__disjoint(subdomain2,box,err)


    def grid_outer_approximation(self,grid,depth):
        gts=GridTreeSet(grid)
        for box in self.__outer_approximation(self._domain,depth):
            gts.adjoin_outer_approximation(box,depth+4)
        return gts

    def box_outer_approximation(self,depth):
        return self.__outer_approximation(self._domain,depth)

    def outer_approximation(self,grid,depth):
        return self.grid_outer_approximation(grid,depth)

    def __outer_approximation(self,subdomain,depth):
        eps=__builtin__.pow(2.0,-depth)
        for (id,constraint) in self._zero.items():
            if constraint(subdomain[:constraint.argument_size()]).upper() < 0.0 or constraint(subdomain).lower() > 0.0:
                return []
        for (id,constraint) in self._positive.items():
            if constraint(subdomain[:constraint.argument_size()]).upper() < 0.0:
                return []
        for (id,constraint) in self._negative.items():
            if constraint(subdomain[:constraint.argument_size()]).lower() > 0.0:
                return []
        #for (id,constraint) in self._equations.items():
        #    constraint_range=constraint(subdomain)
        #    if constraint_range.lower() > 0.0 or constraint_range.upper()<0.0:
        #        return []
        for (id,constraint) in self._always_negative_and_zero.items():
            constraint_range=constraint(subdomain)
            if constraint_range.lower() > 0.0 or constraint_range.upper()<0.0:
                return []
            lowdomain=Box(subdomain[:constraint.argument_size()])
            while lowdomain[-1].lower()>=self._domain[-1].lower():
                if constraint(lowdomain).lower() > 0.0:
                    return []
                lowdomain[-1]=lowdomain[-1]-lowdomain[-1].width()
        for (id,constraint) in self._always_negative.items():
            lowdomain=Box(subdomain[:constraint.argument_size()])
            while lowdomain[-1].lower()>=self._domain[-1].lower():
                if constraint(lowdomain).lower() > 0.0:
                    return []
                lowdomain[-1]=lowdomain[-1]-lowdomain[-1].width()
        #for (id,constraint,auxiliary) in self._ranged_negative.items():
        #    assert False
        range=Box(self._function(subdomain))
        if range.radius()<eps:
            return [range]
        else:
            (subdomain1,subdomain2)=split(subdomain)
            return self.__outer_approximation(subdomain1,depth) + \
                self.__outer_approximation(subdomain2,depth)


    def plot(self,filename,bounding_box=None,resolution=4):
        print "PLOT"
        assert(self.function.result_size()==2)
        fig=Figure()
        if bounding_box==None:
            bounding_box=self.bounding_box()
        fig.set_bounding_box(bounding_box)
        approximation=self.grid_outer_approximation(resolution)
        if empty(approximation):
            print self,"is empty",
        print approximation.measure()
        fig.draw(approximation)
        #for box in self.grid_outer_approximation(resolution):
        #    fig.draw(box)
        fig.write(filename)

    def box_draw(self,figure,resolution=4):
        assert(self._function.result_size()==2)
        for box in self.box_outer_approximation(resolution):
            figure.draw(box)

    def grid_draw(self,figure,resolution=4):
        assert(self._function.result_size()==2)
        approximation=self.grid_outer_approximation(Grid(2),resolution)
        print approximation.measure()
        figure.draw(approximation)

    def draw(self,figure,resolution=4):
        self.box_draw(figure,resolution)
        #self.grid_draw(figure,resolution)

    def __str__(self):
        domain=self._domain
        invariants=[ str(id)+":"+str(constraint)+"<=0, "+str(constraint(domain)) \
            for (id,constraint) in self._always_negative.items()]
        increasing_invariants=[ str(id)+":"+str(constraint)+"<=0, "+str(constraint(domain)) \
            for (id,constraint) in self._negative.items()]
        activations=[ str(id)+":"+str(constraint)+"<=0, "+str(constraint(domain)) \
            for (id,constraint) in self._negative.items()]
        res="ConstrainedImageSet(\n"
        res+="( domain="+str(self._domain)
        res+="( function="+str(self._function.function())
        res+=", codomain="+str(self._function(self._domain))
        res+=", invariants="+str(self._always_negative)
        res+=", increasing_invariants="+str(self._negative)
        res+=", activations="+str(self._positive)
        res+=", equations="+str(self._zero)
        res+=" )"
        return res

def lie_derivative(g,f):
    assert(isinstance(g,ScalarFunction))
    assert(isinstance(f,VectorFunction))
    assert(g.argument_size()==f.result_size()==f.argument_size())
    r=ScalarFunction(g.argument_size())
    for i in range(0,g.argument_size()):
        r=r+derivative(g,i)*f[i]
    return r

if __name__=='__main__':

    id=VectorFunction.identity(2)
    x=ScalarFunction.variable(2,0)
    y=ScalarFunction.variable(2,1)
    f=VectorFunction([x,y-x*x])
    g=ScalarFunction(1-(x*x+y*y))
    h=ScalarFunction(x-y+y*y*y)
    d=Box([{-2:+2},{-2:+2}])

    print lie_derivative(g,f)

    fig=Figure()
    ConstrainedImageSet.draw=ConstrainedImageSet.grid_draw
    ConstrainedImageSet.draw=ConstrainedImageSet.box_draw
    resolution=4

    fig.set_bounding_box(f(d))
    set=ConstrainedImageSet(d,f,[],[],[('h',h)],[])
    fig.set_fill_colour(0.0,0.0,1.0)
    set.draw(fig,resolution)
    set=ConstrainedImageSet(d,f,[],[],[],[('h',h)])
    fig.set_fill_colour(0.0,1.0,1.0)
    set.draw(fig,resolution)
    for cell in set.grid_outer_approximation(resolution): fig.draw(cell.box())
    fig.write("constrained_image_set-1")
    fig.clear()


    fig.set_bounding_box(id(d))
    set=ConstrainedImageSet(d,id,[],[],[],[])
    fig.set_fill_colour(0.0,0.0,1.0)
    set.draw(fig,resolution)
    set=ConstrainedImageSet(d,id,[('h',h)],[],[],[])
    fig.set_fill_colour(1.0,0.0,1.0)
    set.draw(fig,resolution)
    set=ConstrainedImageSet(d,id,[],[('h',h)],[],[])
    fig.set_fill_colour(0.0,1.0,1.0)
    set.draw(fig,resolution)
    c=InterpolatedCurve(0.0,[0.0,0.0])
    for i in range(0,64+1):
        s=d[1].lower()+d[1].width()*(i/64.0)
        c.insert(s,[s-s*s*s,s])
    fig.draw(c)
    fig.write("constrained_image_set-2")

