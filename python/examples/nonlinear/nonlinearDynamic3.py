#!/usr/bin/python

##############################################################################
#            nonlinear_integration.py
#
#  Copyright 2008  Davide Bresolin
#
# Maintainer: Davide Bresolin <davide.bresolin@univr.it>
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
#from ariadne.numeric import *

import sys

#Prova di utilizzo di funzioni trigonometriche
#import math 


# Construct Python-based nonlinear function
print "Constructing and evaluating nonlinear integration system defined in Python"
def fun(x):
     secondaVariabile = sin(x[1])*cos(x[0])
     terzaVariabile = 3*cos(x[2]) + 2*sin(x[1])
     quartaVariabile = 2*cos(x[0])
     return [ 1.5*x[0]-x[0]*x[1]+0.375*x[1], secondaVariabile, terzaVariabile, quartaVariabile ]
     
fun.result_size=4
fun.argument_size=4
nl_fun=AriadneFunction(fun)
print "nl_function =",nl_fun

# Compare the behaviour of fun e nl_fun
for i in range(11):
  x = Float(i*0.1)
  print "fun(",x,") = ",fun([x, x, x, x])
  #print "nl_fun(",x,") = ",nl_fun(Vector([x, x, x, x]))

#sys.exit(0)

# Define the corresponding VectorField
vector_field=VectorField(nl_fun)

# Construct initial set and bounding set
initial_set=RectangularSet([[0.0,0.00001],[0.0,0.000001],[0.0,0.00001],[1.0,1.00001]])
bounding_set=RectangularSet([[0,2],[0,2],[0,2],[0,2]])

# Construct integrator (use Kuhn integrator for nonlinear systems)
integrator=KuhnIntegrator(3,4)

# Construct evolver
parameters=EvolutionParameters()
parameters.set_grid_length(0.1)
parameters.set_bounding_domain_size(2.0)
parameters.set_maximum_step_size(0.125)

evolver=VectorFieldEvolver(parameters,integrator)

# Compute chain-reachable set
print "Computing chain-reachable set..."
chain_reach_set = evolver.upper_reach(vector_field,initial_set,Rational(1.0))
print "Found", chain_reach_set.size(), "cells in grid with", chain_reach_set.capacity(), "cells."

# Reduce chain-reachable set to partition tree set
#print "Computing tree representation of chain-reachable set..."
#chain_reach_tree_set = PartitionTreeSet(chain_reach_set)
#print "Reduced to", chain_reach_tree_set.size(),"cells " \
#    "in partition tree with", chain_reach_tree_set.capacity(),"cells."

# Export to postscript output
print "Exporting to postscript output...",
epsbb=Box([[-4.1,4.1],[-4.1,4.1],[-4.1,4.1],[-4.1,4.1]]) # eps bounding box
eps=EpsPlot()
eps.open("nonlinear_system-mask_set.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(chain_reach_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

# eps.open("nonlinear_system-tree_set.eps",epsbb)
# eps.set_pen_colour("black")
# eps.set_fill_colour("green")
# eps.set_line_style(0)
# eps.write(chain_reach_tree_set)
# eps.set_line_style(1)
# eps.set_fill_style(0)
# eps.write(chain_reach_tree_set.partition_tree())
# eps.set_fill_style(1)
# eps.set_fill_colour("blue")
# eps.write(initial_set)
# eps.close()
print " done."
