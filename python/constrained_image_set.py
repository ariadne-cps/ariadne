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

import sys
from ariadne import *
#from hybrid_system import *

Indeterminate=indeterminate()
Float = float

class ConstrainedImageSet:
    def __init__(self,domain,function, \
                 positive_constraints,equality_constraints,
                 maximum_positive_constraints,maximum_equality_constraints):
        assert(isinstance(domain,Box))
        self.domain=domain
        
        assert(isinstance(function,VectorFunction))
        assert(function.argument_size()==domain.size())
        self.function=function
        
        for constraint in positive_constraints:
            assert(isinstance(constraint,ScalarFunction))
            assert(constraint.argument_size()==domain.size())
        self.positive_constraints=positive_constraints
        
        assert(isinstance(equality_constraints,list))
        for constraint in equality_constraints:
            assert(isinstance(constraint,ScalarFunction))
            assert(constraint.argument_size()==domain.size())
        self.equality_constraints=equality_constraints
       
        for constraint in maximum_positive_constraints:
            assert(isinstance(constraint,ScalarFunction))
            assert(constraint.argument_size()==domain.size())
        self.maximum_positive_constraints=maximum_positive_constraints
        
        assert(isinstance(maximum_equality_constraints,list))
        for constraint in maximum_equality_constraints:
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


    def outer_approximation(self,depth):
        return self.__outer_approximation(self.domain,depth)

    def __outer_approximation(self,subdomain,depth):
        eps=__builtins__.pow(2.0,-depth)
        for constraint in self.positive_constraints:
            if constraint(subdomain).upper() < 0.0:
                return []
        for constraint in self.equality_constraints:
            constraint_range=constraint(subdomain)
            if constraint_range.lower() > 0.0 or constraint_range.upper()<0.0:
                return []
        for constraint in self.maximum_equality_constraints:
            constraint_range=constraint(subdomain)
            if constraint_range.lower() > 0.0 or constraint_range.upper()<0.0:
                return []
            lowdomain=Box(subdomain)
            while lowdomain[-1].lower()>=self.domain[-1].lower():
                if constraint(lowdomain).lower() > 0.0:
                    return []
                lowdomain[-1]=lowdomain[-1]-lowdomain[-1].width()
        for constraint in self.maximum_positive_constraints:
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


    def plot(self,filename,resolution=4):
        assert(self.function.result_size()==2)
        fig=Figure()
        fig.set_bounding_box(self.bounding_box())
        for box in self.outer_approximation(resolution):
            fig.draw(box)
        fig.write(filename)

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
    resolution=3

    fig.set_bounding_box(f(d))
    s=ConstrainedImageSet(d,f,[],[],[h],[])
    fig.set_fill_colour(0.0,0.0,1.0)
    for box in s.outer_approximation(resolution): fig.draw(box)
    s=ConstrainedImageSet(d,f,[],[],[],[h])
    fig.set_fill_colour(0.0,1.0,1.0)
    for box in s.outer_approximation(resolution): fig.draw(box)
    fig.write("constrained_image_set-1")

    fig.clear()
    
    fig.set_bounding_box(id(d))
    s=ConstrainedImageSet(d,id,[],[],[],[])
    fig.set_fill_colour(0.0,0.0,1.0)
    for box in s.outer_approximation(resolution): fig.draw(box)
    s=ConstrainedImageSet(d,id,[h],[],[],[])
    fig.set_fill_colour(1.0,0.0,1.0)
    for box in s.outer_approximation(resolution): fig.draw(box)
    s=ConstrainedImageSet(d,id,[],[h],[],[])
    fig.set_fill_colour(0.0,1.0,1.0)
    for box in s.outer_approximation(resolution): fig.draw(box)
    fig.write("constrained_image_set-2")
