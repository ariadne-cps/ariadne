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

def test_calculus():

    dp=DoublePrecision()

    cy=ValidatedNumber(2)
    cx=FloatDPBounds(2,dp)

    bx=ExactBox([{Dyadic(1):3},{-1:2},{-3:3}])

    swp = ThresholdSweeper(1e-8);

    tc=ScalarTaylorFunction.constant(bx,dec(1.5),swp)
    tx=ScalarTaylorFunction.coordinate(bx,0,swp)
    ty=ScalarTaylorFunction.coordinate(bx,1,swp)
    tid=VectorTaylorFunction.identity(bx,swp)
    ty=tid[1]

    t=5+2*tx+ty

    +t; -t; t+t; t-t; t*t; t/t;
    +t; -t; t+t; t-t; t*t;
    t+cy; t-cy; t*cy; t/cy;
    cy+t; cy-t; cy*t;
    t+cx; t-cx; t*cx; t/cx;
    cx+cx; cx-t; cx*t;

    derivative(t,0)
    antiderivative(t,0)

    f=VectorTaylorFunction([x,c,y])
    g=ScalarTaylorFunction(t)
    compose(f,f); compose(g,f)

    p=RealVectorFunction.identity(3)
    q=RealScalarFunction.coordinate(3,1)
    compose(p,f); compose(q,f)

test_calculus()
