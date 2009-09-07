#!/usr/bin/python

##############################################################################
#            test_function.py
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

c=ScalarPolynomialFunction.constant(3,1.5)
x=ScalarPolynomialFunction.variable(3,0)
y=ScalarPolynomialFunction.variable(3,1)
v=ScalarPolynomialFunction.variables(3)

p=x+y

+p; -p; p+p; p-p; p*p;
p+d; p-d; p*d; p/d;
d+p; d-p; d*p;
p+i; p-i; p*i; p/i;
i+i; i-p; i*p;

derivative(p,0)
antiderivative(p,0)

f=VectorPolynomialFunction([p,p])

