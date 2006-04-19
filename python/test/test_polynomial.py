#!/usr/bin/python

##############################################################################
#            test_polynomial.py
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

from ariadne import Real
from ariadne.base import *
from ariadne.numeric import *
from ariadne.evaluation import *
from ariadne.geometry import *
from ariadne.linear_algebra import *

a=SizeArray(3)
a[0]=1
a[1]=2
a[2]=3

m=Monomial(Real(2.3),a)
print m
m=Monomial("2.3 x_0^2 * x_2^3  z")
print m

m=Monomial("1.1 x_0  z")
print m
m=Monomial("2.3 x_1^2 z")
print m
m=Monomial("x_0 x_1 * x_2 z")
print m

print

#p=Polynomial("1.1 x_0 + 2.3 * x_1^2  z")
#print p
#q=Polynomial("2.3 * x_1^2 + 1.1 x_0  z")
#print q
p=Polynomial("4.2  - x_0 + x_1 + 2 x_0^2 - 1.0 + 2 x_0 x_1 + x_1^2")
q0=Point("(0,0)")
q1=Point("(1,0)")
q2=Point("(0,1)")
q3=Point("(1,1)")
r3=Rectangle("[0.99,1.01]x[0.99,1.01]")
print p
print p.apply(q0), p.apply(q1), p.apply(q2), p.apply(q3)
print p.apply(r3)

print

pm=PolynomialMap("[4.2 - 1.1 x_0 + 2.3 x_1^2, 1 x_1 ]  z")
print pm
print pm.derivative()
print pm.derivative(q3)
print pm.derivative(r3)

print "Exiting"
