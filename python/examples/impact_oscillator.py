#!/usr/bin/python

##############################################################################
#            impact_oscillator.py
#
#  Copyright 2008  Pieter Collins
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

#Set parameters for impact oscillator
z=1.0
T=1.0
l=0.75
d=0.125

parameters=Point([z,T,l,d])
p=parameters

#Construct Python-based impact osciallator
print "Constructing impact osciallator"
def dynamic_function(x):
    twopi=6.2831853071795862;
    return [ x[1], -p[0]*x[1]-x[0]+cos(twopi*x[2]/p[1]), +1.0 ];
dynamic_function.result_size=3
dynamic_function.argument_size=3
dynamic_function=AriadneFunction(dynamic_function)

def guard_function(x):
    return [ p[3] - r[1] ]
guard_function.result_size=1
guard_function.argument_size=3
guard_function=AriadneFunction(guard_function)

def reset_function(x):
    return [ x[0], -p[2]*x[1], x[2] ]
reset_function.result_size=3
reset_function.argument_size=3
reset_function=AriadneFunction(reset_function)

impact_oscillator = ImpactSystem(dynamic_function, guard_function, reset_function)


# Choose initial point
initial_point=Point([0.5,0.0,0.0])

# Construct initial set and bounding set
initial_set=ImageSet(initial_point)
bounding_set=ConstraintSet(Box([[0,3],[-2,2],[0,32]]))

# Construct evolver
parameters=EvolutionParameters()
#parameters.set_grid_length(0.125)
parameters.set_grid_length(0.0625)
parameters.set_bounding_domain_size(4.0)
evolver=ImpactSystemEvolver(parameters)

# Compute chain-reachable set
print "Computing chain-reachable set..."
chain_reach_set = evolver.chain_reach(impact_system,ConstraintSet(initial_set))
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
eps.open("impact_oscillator-mask_set.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(chain_reach_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

eps.open("impact_oscillator-tree_set.eps",epsbb)
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
