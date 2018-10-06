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

    swp = ThresholdSweeper(dp,1e-8);

    tc=ValidatedScalarMultivariateTaylorFunctionModel.constant(bx,cast_exact(1.5),swp)
    tx=ValidatedScalarMultivariateTaylorFunctionModel.coordinate(bx,0,swp)
    ty=ValidatedScalarMultivariateTaylorFunctionModel.coordinate(bx,1,swp)
    tid=ValidatedVectorMultivariateTaylorFunctionModel.identity(bx,swp)
    ty=tid[1]

    tf=5+2*tx+ty

    +tf; -tf; tf+tf; tf-tf; tf*tf; tf/tf;
    +tf; -tf; tf+tf; tf-tf; tf*tf;
    tf+cy; tf-cy; tf*cy; tf/cy;
    cy+tf; cy-tf; cy*tf;
    tf+cx; tf-cx; tf*cx; tf/cx;
    cx+cx; cx-tf; cx*tf;

    derivative(tf,0)
    antiderivative(tf,0)

    vtf=ValidatedVectorMultivariateTaylorFunctionModel([tx,tc,ty])
    vtf=tid
    stf=ValidatedScalarMultivariateTaylorFunctionModel(tf)
    compose(vtf,vtf)
    compose(stf,vtf)

    vf=EffectiveVectorMultivariateFunction.identity(3)
    sf=EffectiveScalarMultivariateFunction.coordinate(3,1)
    compose(sf,vtf); compose(vf,vtf)
