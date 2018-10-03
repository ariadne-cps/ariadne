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

def test_function():

    def dec(x): return Decimal(x)

    dp = DoublePrecision()
    r=Real(Decimal(1.75))

    c=EffectiveScalarMultivariateFunction.constant(3,r)
    x=EffectiveScalarMultivariateFunction.coordinate(3,0)
    y=EffectiveScalarMultivariateFunction.coordinate(3,1)
    id=EffectiveVectorMultivariateFunction.identity(3)

    p=x+y

    +p; -p; p+p; p-p; p*p;
    p+r; p-r; p*r; p/r;
    #r+p; r-p; r*p;

    derivative(p,0)

#    b=ExactBox([{0:dec(0.25)},{dec(0.25):dec(0.50)},{dec(0.50):dec(0.75)}])
    b=FloatDPBoundsVector([1,2,3],dp)
    f=EffectiveVectorMultivariateFunction([c,x,y])
    g=EffectiveScalarMultivariateFunction(p)

    join(g,g); join(f,g); join(g,f); join(f,f)
    compose(g,f); compose(f,f)
    evaluate(f,b);
    evaluate(g,b)

test_function()
