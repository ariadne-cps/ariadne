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

mu=Real(0.25)

vdp=VanDerPolEquation(mu)
interval_newton=IntervalNewtonSolver(Real(0.00001),64)

subdivisions=128
grid_extent=Rectangle("[--4,4]x[-2,2]") # grid bounding box
finite_grid=FiniteGrid(grid_extent,128)
grid=finite_grid.grid()
initial_set=Rectangle("[1.00,1.000002]x[0,0.000001]") # initial state
initial_set=Parallelotope(initial_state)
initial_cell=over_approximation(initial_state,g)

cb=Rectangle(gbb) # cutoff box
initial_set=GridMaskSet(fg)
initial_set.adjoin(over_approximation(fixed_point,g))
bounding_set=GridMaskSet(fg)
bounding_set.adjoin(over_approximation(gbb,g))
time=Rational(5)
print "Computing chain-reachable set..."
reach=reach(vdp,initial_set,time)
print reach
#print "Found",cr.size(),"cells in grid with",cr.capacity(),"cells."

ptscr=PartitionTreeSet(gmscr)
print "Reduced to",ptscr.size(),"cells in partition tree with",ptscr.capacity(),"cells."

print "Exporting to postscript output...",
epsbb=Rectangle("[-4.1,4.1]x[-2.1,2.1]") # eps bounding box
eps=EpsPlot("van_der_pol_oscillator-1.eps",epsbb)
eps.set_line_style(False)
eps.set_fill_colour("red")
eps.write(difference(gmscr.neighbourhood(),gmscr))
eps.set_fill_colour("blue")
eps.write(gmscr.adjoining())
eps.set_fill_colour("green")
eps.write(gmscr)
eps.set_fill_colour("blue")
eps.write(fixed_point_cell)

eps=EpsPlot("henon_attractor-2.eps",epsbb)
eps.set_pen_colour("black")
eps.set_fill_colour("white")
eps.write(cb)
eps.set_line_style(0)
eps.set_fill_colour("green")
eps.write(gmscr)
eps.set_line_style(1)
eps.set_fill_style(0)
eps.write(ptscr.partition_tree())
eps.set_fill_style(1)
eps.set_fill_colour("blue")
eps.write(fixed_point_cell)
eps.close()
print " done."
