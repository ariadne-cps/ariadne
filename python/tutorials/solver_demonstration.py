#!/usr/bin/python3

##############################################################################
#            solver_demonstration.py
#
#  Copyright  2009-21  Pieter Collins
##############################################################################

# This file is part of Ariadne.

# Ariadne is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Ariadne is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Ariadne. If not, see <https://www.gnu.org/licenses/>.

# Import all classes in the ariadne module
from pyariadne import *


def algebraic_solver_demonstration():
    #! [Algebraic Solver demonstration]

    # Compute the solution h to the vector equation f(x,h(x))=0
    slv=IntervalNewtonSolver(1e-8,6)
    dom=BoxDomainType([{-1:+1},{-1:+1},{-1:+1}])
    x=EffectiveScalarMultivariateFunction.coordinate(3,0)
    y0=EffectiveScalarMultivariateFunction.coordinate(3,1)
    y1=EffectiveScalarMultivariateFunction.coordinate(3,2)
    f=join(x+4*y0+y1,y0+y1)
    print("f:",f)
    hf=slv.implicit(f,BoxDomainType([{-1:+1}]),BoxDomainType([{-1:+1},{-1:+1}]))
    print("implicit(f):",hf)

    # Compute the solution h to the scalar equation g(x,h(x))=0
    # with f(x,y)=4+x-y^2, so y=sqrt(4+x)
    dom=BoxDomainType([{-1:+1},{-1:+1}])
    x=EffectiveScalarMultivariateFunction.coordinate(2,0)
    y=EffectiveScalarMultivariateFunction.coordinate(2,1)
    g=x-4*y+y*y
    print("g:",g)
    hg=slv.implicit(g, BoxDomainType([{-1:+1}]),IntervalDomainType({-1:+1}))
    print("implicit(g):",hg)
    #! [Algebraic Solver demonstration]


def differential_solver_demonstration():
    #! [Differential Solver demonstration]

    # Compute the flow of the Taylor function f starting in the domain dom for time interval [-h,+h]
    integrator=GradedTaylorSeriesIntegrator(1e-8)

    dom=BoxDomainType([{-1:+1}])
    bbx=BoxDomainType([{-4:+4}])
    h=1/two
    f=ValidatedVectorMultivariateFunction.identity(1)

    print(f,dom,h)
    phis=integrator.flow(f,dom,h)
    assert(len(phis)==1)
    phi=phis[0]
    print("phi:",phi,type(phi))

    # Compute two time steps of the flow of the Taylor function f starting in domain D for the interval [h,2h]
    phi0=phi
    print("phi.domain():",phi.domain())
    print("h:",h)
    phi0h=partial_evaluate(phi,1,FloatDPBounds(h,dp))
    dom1=phi0h.codomain()
    phi=integrator.flow(f,dom1,h)[0]
    print("phi:",phi)
    swp=GradedSweeperDP(dp,6);
    tr=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate([{0:2*h}],0,swp)-h
    tr=ValidatedScalarMultivariateFunctionModelDP(tr)
    phi1=compose(phi,combine(phi0h,tr))
    print("phi0:",phi0)
    print("phi1:",phi1)
    #! [Differential Solver demonstration]


if __name__=='__main__':
    algebraic_solver_demonstration()
    print()
    differential_solver_demonstration()
