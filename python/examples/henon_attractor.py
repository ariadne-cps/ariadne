#!/usr/bin/python

##############################################################################
#            henon_attractor.py
#
#  Copyright 2006  Pieter Collins
#
# Maintainer: Pieter Collins <Pieter.Collins@cwi.nl>
#
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

#Set parameters for Henon attractor
a=Float(1.5)
b=Float(-0.375)

#Construct Python-based Henon map.
print "Constructing and evaluating Henon map defined in Python"
def henon_function(x):
    return [ a-x[0]*x[0]-b*x[1], x[0] ]
henon_function.result_size=2
henon_function.argument_size=2
henon_function=AriadneFunction(henon_function)
print "henon_function =",henon_function
python_henon_map=Map(henon_function)
pt=IntervalPoint([1.2,1.2])
print python_henon_map(pt)

#Construct built-in Henon map
print "Constructing built-in Henon map"
builtin_henon_map=HenonMap(Point([a,b]))
print builtin_henon_map

#Select which version of the Henon map to use
#henon_map=python_henon_map
henon_map=builtin_henon_map

# Find a fixed point
initial_guess=IntervalPoint( [Interval(0,2),Interval(0,2)] ) # initial state
interval_newton=IntervalNewtonSolver(0.00001,64)
fixed_point=interval_newton.fixed_point(henon_map,initial_guess)
print "Found fixed point in",fixed_point

# Construct initial set and bounding set
initial_set=RectangularSet(fixed_point)
bounding_set=RectangularSet([[-10,6],[-8,8]])

# Construct evolver
parameters=EvolutionParameters()
#parameters.set_grid_length(0.125)
parameters.set_grid_length(0.0625)
parameters.set_bounding_domain_size(4.0)
evolver=MapEvolver(parameters)

# Compute chain-reachable set
print "Computing chain-reachable set..."
chain_reach_set = evolver.chain_reach(henon_map,ConstraintSet(initial_set))
print "Found", chain_reach_set.size(), "cells in grid with", chain_reach_set.capacity(), "cells."

# Reduce chain-reachable set to partition tree set
print "Computing tree representation of chain-reachable set..."
chain_reach_tree_set = PartitionTreeSet(chain_reach_set)
print "Reduced to", chain_reach_tree_set.size(),"cells " \
    "in partition tree with", chain_reach_tree_set.capacity(),"cells."

# Export to postscript output
print "Exporting to postscript output...",
epsbb=Box([[-4.1,4.1],[-4.1,4.1]]) # eps bounding box
eps=EpsPlot()
eps.open("henon_attractor-mask_set.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(chain_reach_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

eps.open("henon_attractor-tree_set.eps",epsbb)
eps.set_line_colour("black")
eps.set_fill_colour("green")
eps.set_line_style(0)
eps.write(chain_reach_tree_set)
eps.set_line_style(1)
eps.set_fill_style(0)
eps.write(chain_reach_tree_set.partition_tree())
eps.set_fill_style(1)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()
print " done."
