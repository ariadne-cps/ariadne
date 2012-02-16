#!/usr/bin/python

##############################################################################
#            test_calculus.py
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

from ariadne import *
Float=float

d=Float(2.5)
i=Interval(1.5,1.75)

bx=Box([{1:3},{-1:2},{-3:3}])

swp = ThresholdSweeper(1e-8);

c=ScalarTaylorFunction.constant(bx,1.5,swp)
x=ScalarTaylorFunction.coordinate(bx,0,swp)
y=ScalarTaylorFunction.coordinate(bx,1,swp)
v=VectorTaylorFunction.identity(bx,swp)
y=v[1]
t=5+2*x+y

+t; -t; t+t; t-t; t*t; t/t;
+t; -t; t+t; t-t; t*t;
t+d; t-d; t*d; t/d;
d+t; d-t; d*t;
t+i; t-i; t*i; t/i;
i+i; i-t; i*t;

derivative(t,0)
antiderivative(t,0)

f=VectorTaylorFunction([x,c,y])
g=ScalarTaylorFunction(t)
compose(f,f); compose(g,f)

p=RealVectorFunction.identity(3)
q=RealScalarFunction.coordinate(3,1)
compose(p,f); compose(q,f)

