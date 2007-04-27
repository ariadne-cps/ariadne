#!/usr/bin/python

##############################################################################
#            henon_attractor.py
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

from ariadne import *
import sys

print "\nAvailable classes and functions\n",dir(),"\n\n"

a=Float(1.5)
b=Float(0.875)
henon_map=HenonMap(a,b)
print henon_map


grid_extent=Rectangle("[-10,6]x[-8,8]") 
number_of_subdivisions=128
finite_grid=FiniteGrid(grid_extent,number_of_subdivisions)
grid=finite_grid.grid()
initial_guess=IntervalPoint("([0,2],[0,2])") # initial state

interval_newton=IntervalNewtonSolver(0.00001,64)
fixed_point=interval_newton.fixed_point(henon_map,initial_guess)
fixed_point_cell=outer_approximation(fixed_point,grid)
print "Found fixed point in",fixed_point

initial_set=GridMaskSet(finite_grid)
initial_set.adjoin(outer_approximation(fixed_point,grid))
bounding_set=GridMaskSet(finite_grid)
bounding_set.adjoin(over_approximation(grid_extent,grid))

apply=Applicator();

print "Computing chain-reachable set..."

chain_reach_set = apply.chainreach(henon_map,initial_set,bounding_set)
print "Found", chain_reach_set.size(), "cells in grid with", chain_reach_set.capacity(), "cells."

chain_reach_tree_set = PartitionTreeSet(chain_reach_set)
print "Reduced to", chain_reach_tree_set.size()," cells" \
    "in partition tree with", chain_reach_tree_set.capacity(),"cells."

print "Exporting to postscript output...",
epsbb=Rectangle("[-4.1,4.1]x[-4.1,4.1]") # eps bounding box
eps=EpsPlot()
eps.open("henon_attractor-1.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(chain_reach_set)
eps.set_fill_colour("blue")
eps.write(fixed_point_cell)
eps.close()

eps.open("henon_attractor-2.eps",epsbb)
eps.set_pen_colour("black")
eps.set_fill_colour("green")
eps.set_line_style(0)
eps.write(chain_reach_tree_set)
eps.set_line_style(1)
eps.set_fill_style(0)
eps.write(chain_reach_tree_set.partition_tree())
eps.set_fill_style(1)
eps.set_fill_colour("blue")
eps.write(fixed_point_cell)
eps.close()
print " done."
