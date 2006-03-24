#!/usr/bin/python

##############################################################################
#            test_zonotope.py
#
#  Copyright 2006  Alberto Casagrande, Pieter Collins 
#  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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


from ariadne.base import *
from ariadne.evaluation import *
from ariadne.geometry import *
from ariadne.linear_algebra import *
from math import *
import sys

r=Rectangle("[9,11]x[5,11]")
r2=Rectangle("[5,6]x[3,4]")


z=Zonotope(r);
z2=Zonotope(r2);

A=Matrix(2,2)
b=Vector(2)

alpha=3.14/4
r_fact=7.0/8

A[0,0]=r_fact*cos(alpha)
A[0,1]=r_fact*sin(alpha)
A[1,0]=r_fact*-sin(alpha)
A[1,1]=r_fact*cos(alpha)

M=AffineMap(A,b)

z2=M(z2)

z3=minkowski_sum(z,z2)

p=Point(2)

p[0]=17
p[1]=15

p2=Point(p)
p3=Point(p)

p2[1]=8
p3[0]=14

p4=Point(p3)
p5=Point(p3)
p6=Point(z3.centre())

p4[1]=p2[1]

p5[1]=p2[1]-3

print z3.contains(p)
print z3.contains(p2)
print z3.contains(p3)
print z3.contains(p4)
print z3.contains(p5)
print z3.interior_contains(p)
print z3.interior_contains(p2)
print z3.interior_contains(p3)
print z3.interior_contains(p4)
print z3.interior_contains(p6)

print z3 
print p
print p2
print p3
print p4
print p5
print p6


