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
#from hybrid_system import *

Indeterminate=indeterminate()
Float = float

Polynomial=ScalarPolynomialFunction
Polynomial.__repr__=Polynomial.__str__
Taylor=ScalarTaylorFunction
Taylor.__repr__=Taylor.__str__

def evolve(dynamic):
    pass

class ConstraintSet:
    def __init__(self,domain,function,equality_constraints,positive_constraints):
        self.domain=domain
        self.function=function
        self.equality_constraints=equality_constraints
        self.positive_constraints=positive_constraints

    def bounding_box(self):
        return function(domain)

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

if __name__=="__main__":
    x=ScalarAffineFunction([1,2],3)
    print x

    one=Polynomial.constant(2,1.0)
    x=Polynomial.variable(2,0)
    y=Polynomial.variable(2,1)

    g0=y-one
    g1=x
    g2=y-x/2+x*x*x/4

    f=[one,x]
    print f

    x=RealVariable('x')
    y=RealVariable('y')
    e=x*x+sin(y)
    print e
    s=RealSpace([x,y])
    f=RealFunction(e,s)
    print f

    d=Box([{-2:+2},{-3:+3}])
    t=Taylor(d,f)
    t,t([1,2]),f([1,2])

    cs=ConstraintSet(d,VectorFunction([f]),None,None)
    print cs.disjoint(Box([{-1:+1}]))
