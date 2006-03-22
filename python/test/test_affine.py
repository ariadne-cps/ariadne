#!/usr/bin/python

##############################################################################
#            test_affine.py
#
#  Copyright 2006  Pieter Collins <Pieter.Collins@cwi.nl>
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

A=Matrix(2,2)
b=Vector(2)

alpha=3.14/4
r_fact=7.0/8

A[0,0]=r_fact*cos(alpha)
A[0,1]=r_fact*sin(alpha)
A[1,0]=r_fact*-sin(alpha)
A[1,1]=r_fact*cos(alpha)

M=AffineMap(A,b)

r=Rectangle("[9,11]x[5,11]")
r2=Rectangle("[5,6]x[3,4]")
box=Rectangle("[-14,14]x[-14,14]")

p=Parallelotope(r)
p2=Parallelotope(r2)

print p

eps=EpsPlot("affine.eps",box)
eps.set_fill_colour("blue")
eps.write(p)
eps.write(p2)
eps.set_fill_colour("green")

for i in range(0,100):
	p=M(p)
	p2=M(p2)
	eps.write(p)
	eps.write(p2)

eps.close()
