#!/usr/bin/python

##############################################################################
#            lorenz_system.py
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

#Set standard parameters for the Lorenz system
beta=Float(8./3.)
rho=Float(10)
sigma=Float(28)

#Give the evolution time for finite-time evolution.
time=Rational(2)

#Construct built-in Lorenz system
print "Constructing built-in Lorenz system"
lorenz_system=LorenzSystem(Point([beta,rho,sigma]))
print lorenz_system

# Start at the origin
initial_point=Point([0,0,0]);

# Construct concrete initial sets
initial_box=Box([[0.9,1.1],[-0.1,0.1],[-0.1,0.1]])
initial_zonotope=Zonotope(initial_box)

# Construct initial set and bounding set
initial_set=ImageSet(initial_point)

# Define evolution parameters
parameters=EvolutionParameters()
#parameters.set_grid_length(0.125)
parameters.set_grid_length(0.0625)
parameters.set_maximum_enclosure_radius(0.25)
parameters.set_bounding_domain_size(8.0)

evolver=default_evolver(lorenz_system,initial_zonotope,parameters)
approximator=default_approximator(Box(0),Zonotope(0))

analyser=VectorFieldReachabilityAnalyser(parameters,evolver,approximator)

# Compute finite-time reachable and evolved set
print "Computing evolve set..."
evolve_set = evolver.evolve(lorenz_system,initial_zonotope,time)
print "Computing reach set..."
reach_set = evolver.reach(lorenz_system,initial_zonotope,time)
print "evolve_set has",evolve_set.size(),"enclosures"
print "reach_set has",reach_set.size(),"enclosures"


# Compute chain-reachable set
print "Computing chain-reachable set..."
#chain_reach_set = analyser.chain_reach(lorenz_system,initial_set)
chain_reach_set=initial_set
#print "Found", chain_reach_set.size(), "cells in grid with", chain_reach_set.capacity(), "cells."

# Reduce chain-reachable set to partition tree set
print "Computing tree representation of chain-reachable set..."
#tree_chain_reach_set = PartitionTreeSet(chain_reach_set)
#print "Reduced to", chain_reach_tree_set.size(),"cells " \
#    "in partition tree with", chain_reach_tree_set.capacity(),"cells."

# Export to postscript output
print "Exporting to postscript output...",
epsbb=Box([[-8.1,8.1],[-8.1,8.1],[-8.1,8.1]]) # eps bounding box
eps=EpsPlot()
eps.open("lorenz_attractor-reach_set.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(reach_set)
eps.set_fill_colour("yellow")
eps.write(evolve_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

eps.open("lorenz_attractor-tree_set.eps",epsbb)
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
