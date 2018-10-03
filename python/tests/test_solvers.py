#!/usr/bin/python

##############################################################################
#            test_solvers.py
#
#  Copyright 2012  Pieter Collins <Pieter.Collins@cwi.nl>
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


fo=EffectiveScalarMultivariateFunction(2)
fx=EffectiveScalarMultivariateFunction.coordinate(2,0)
fy=EffectiveScalarMultivariateFunction.coordinate(2,1)

fn=ValidatedVectorMultivariateFunction([fx+fy,fy])
b=ExactBox([{-1:+1},{-1:+1}])

solver=IntervalNewtonSolver(1e-8,12)
solver.solve(fn,b)

solver=KrawczykSolver(1e-8,12)
solver.solve(fn,b)

d=ExactBox([{-1:+1},{-1:+1}])
h=Dyadic(0.25)
h=FloatDP(h,DoublePrecision())

vf=EffectiveVectorMultivariateFunction([fo,fx])
vf=ValidatedVectorMultivariateFunction(vf)
integrator=TaylorPicardIntegrator(1e-8)
integrator.flow_step(vf,d,h)

vf=ValidatedVectorMultivariateFunction([fo,fx])
integrator=TaylorSeriesIntegrator(1e-8)
integrator.flow_step(vf,d,h)

