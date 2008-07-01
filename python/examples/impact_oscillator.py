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

# Set the evolution time
time=Rational(15)

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
initial_box=Box(initial_point)
initial_zonotope=Zonotope(initial_box)
initial_set=ImageSet(initial_point)
bounding_set=ConstraintSet(Box([[0,3],[-2,2],[0,32]]))

# Construct evolver
parameters=EvolutionParameters()
#parameters.set_grid_length(0.125)
parameters.set_grid_length(0.0625)
parameters.set_bounding_domain_size(4.0)

# Construct the evolver. The first two arguments are dummys which are
# used to infer the type of the evolver to be constructed; the values
# are irrelevant.
evolver=default_evolver(impact_oscillator,initial_zonotope,parameters)

# Compute reachable set
print "Computing reachable set..."
reach_set = evolver.reach(impact_oscillator,initial_zonotope,time)
print "Computing evolved set..."
reach_set = evolver.evolve(impact_system,initial_zonotope,time)

print "Found", reach_set.size(), "enclosures."

# Export to postscript output
print "Exporting to postscript output...",
epsbb=Box([[-4.1,4.1],[-4.1,4.1]]) # eps bounding box
eps=EpsPlot()
eps.open("impact_oscillator-reach_set.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(reach_set)
eps.set_fill_colour("yellow")
eps.write(evolve_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

