#!/usr/bin/python

##############################################################################
#            test_polyhedron.py
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
from ariadne.linear_algebra import *
from ariadne.combinatoric import *
from ariadne.geometry import *
from ariadne.system import *
from ariadne.evaluation import *
from ariadne.output import *
from ariadne.models import *
import sys

import ariadne.models
print "Available models:",dir(ariadne.models)

mu=Real(0.25)

vdp=VanDerPolEquation(mu)
interval_newton=IntervalNewtonSolver(Real(0.00001),64)
lohner=LohnerIntegrator(0.125,0.5,0.1) #(maximum_step_size,lock_to_grid_time,maximum_set_size)

subdivisions=128
grid_extent=Rectangle("[--4,4]x[-2,2]") # grid bounding box
finite_grid=FiniteGrid(grid_extent,128)
grid=finite_grid.grid()
initial_set=Rectangle("[1.00,1.002]x[05,0.501]")
initial_set=Parallelotope(initial_set)
initial_set=Zonotope(initial_set)

h=Rational(0.125)
t=Rational(0.25)
print "Initial set: ", initial_set
print "Integrate initial parallelotope for one time step"
intermediate_set=lohner.integration_step(vdp,initial_set,h)
print "Integrate initial parallelotope for time 1"
intermediate_set=lohner.integrate(vdp,initial_set,t)
print intermediate_set

print initial_set

#bounding_set.adjoin(over_approximation(grid_extent,grid))
time=Rational(5)
print "Computing chain-reachable set..."
reach_set=lohner.reach(vdp,initial_set,time)
print reach_set

print "Exporting to postscript output...",
epsbb=Rectangle("[-4.1,4.1]x[-2.1,2.1]") # eps bounding box
eps=EpsPlot("van_der_pol_oscillator-1.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(reach_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()
print " done."
