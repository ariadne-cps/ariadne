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

class ConstrainedImageSet:
    def __init__(self,domain,function=None, \
                 positive_constraints=[],equality_constraints=[],
                 maximum_negative_constraints=[],maximum_equality_constraints=[]):
        #Copy constructor
        if function==None:
            assert isinstance(domain,ConstrainedImageSet)
            other=domain
            self.domain=Box(other.domain)
            self.function=VectorTaylorFunction(other.function)
            self.positive_constraints=list(other.positive_constraints)
            self.equality_constraints=list(other.equality_constraints)
            self.maximum_negative_constraints=list(other.maximum_negative_constraints)
            self.maximum_equality_constraints=list(other.maximum_equality_constraints)
            return
        
        #Standard constructor
        if isinstance(domain,IntervalVector):
            domain=Box(domain)
        assert isinstance(domain,Box) 
        self.domain=domain
        
        assert isinstance(function,VectorFunction) or isinstance(function,VectorTaylorFunction)
        assert function.argument_size()==domain.size()
        self.function=function
        
        assert(isinstance(positive_constraints,list))
        for (id,constraint) in positive_constraints:
            assert isinstance(constraint,ScalarFunction) or isinstance(constraint,ScalarTaylorFunction)
            assert constraint.argument_size()==domain.size()
        self.positive_constraints=positive_constraints
        
        assert(isinstance(equality_constraints,list))
        for (id,constraint) in equality_constraints:
            assert(isinstance(constraint,ScalarFunction))
            assert(constraint.argument_size()==domain.size())
        self.equality_constraints=equality_constraints
       
        for (id,constraint) in maximum_negative_constraints:
            assert(isinstance(constraint,ScalarFunction))
            assert(constraint.argument_size()==domain.size())
        self.maximum_negative_constraints=maximum_negative_constraints
        
        assert(isinstance(maximum_equality_constraints,list))
        for (id,constraint) in maximum_equality_constraints:
            assert(isinstance(constraint,ScalarFunction))
            assert(constraint.argument_size()==domain.size())
        self.maximum_equality_constraints=maximum_equality_constraints
       
        self.time_variable=len(domain)-1


    def bounding_box(self):
        return Box(self.function(self.domain)).bounding_box()


    def disjoint(self,box):
        return self.__disjoint(self.domain,box,box.radius())


    def __disjoint(self,domain,box,err):
        if False: #positive_constraints(domain)
            return True
        image=Box(self.function(domain))
        if disjoint(image,box):
            return True;
        if subset(image,box):
            return False
        if image.radius()<err:
            return Indeterminate
        (subdomain1,subdomain2)=split(domain)
        return self.__disjoint(subdomain1,box,err) and self.__disjoint(subdomain2,box,err)


    def grid_outer_approximation(self,depth):
        gts=GridTreeSet(self.function.result_size())
        for box in self.__outer_approximation(self.domain,depth):
            gts.adjoin_outer_approximation(box,depth+4)
        return gts
        
    def outer_approximation(self,depth):
        return self.__outer_approximation(self.domain,depth)

    def __outer_approximation(self,subdomain,depth):
        eps=__builtin__.pow(2.0,-depth)
        for (id,constraint) in self.positive_constraints:
            if constraint(subdomain).upper() < 0.0:
                return []
        for (id,constraint) in self.equality_constraints:
            constraint_range=constraint(subdomain)
            if constraint_range.lower() > 0.0 or constraint_range.upper()<0.0:
                return []
        for (id,constraint) in self.maximum_equality_constraints:
            constraint_range=constraint(subdomain)
            if constraint_range.lower() > 0.0 or constraint_range.upper()<0.0:
                return []
            lowdomain=Box(subdomain)
            while lowdomain[-1].lower()>=self.domain[-1].lower():
                if constraint(lowdomain).lower() > 0.0:
                    return []
                lowdomain[-1]=lowdomain[-1]-lowdomain[-1].width()
        for (id,constraint) in self.maximum_negative_constraints:
            lowdomain=Box(subdomain)
            while lowdomain[-1].lower()>=self.domain[-1].lower():
                if constraint(lowdomain).lower() > 0.0:
                    return []
                lowdomain[-1]=lowdomain[-1]-lowdomain[-1].width()
        range=Box(self.function(subdomain))
        if range.radius()<eps:
            return [range]
        else:
            (subdomain1,subdomain2)=split(subdomain)
            return self.__outer_approximation(subdomain1,depth) + \
                self.__outer_approximation(subdomain2,depth)


    def plot(self,filename,bounding_box=None,resolution=4):
        assert(self.function.result_size()==2)
        fig=Figure()
        if bounding_box==None:
            bounding_box=self.bounding_box()
        fig.set_bounding_box(bounding_box)
        fig.draw(self.grid_outer_approximation(resolution))
        #for box in self.grid_outer_approximation(resolution):
        #    fig.draw(box)
        fig.write(filename)

    def box_draw(self,figure,resolution=4):
        assert(self.function.result_size()==2)
        for box in self.outer_approximation(resolution):
            fig.draw(box)

    def grid_draw(self,figure,resolution=4):
        assert(self.function.result_size()==2)
        fig.draw(self.grid_outer_approximation(resolution))

    def __str__(self):
        res="ConstrainedImageSet["+str(self.function.result_size())+","+str(self.function.argument_size())+"]("
        res+="( domain="+str(self.domain)
        res+=", codomain="+str(self.function(self.domain))
        res+=", maximum_constraints="+str([ (str(id)+"<0",constraint(self.domain)) for (id,constraint) in self.maximum_negative_constraints])
        res+=", positive_constraints="+str([ (str(id)+"<0",constraint(self.domain)) for (id,constraint) in self.positive_constraints])
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

