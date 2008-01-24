#!/usr/bin/python

##############################################################################
#            tank_hybrid_automaton.py
#
#  Copyright 2006  Alberto Casagrande <casagrande@dimi.uniud.it>
#            2007  Pieter Collins <Pieter.Collins@cwi.nl>
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

#Use the set-based hybrid system framework
HybridEvolver=SetBasedHybridEvolver

#Create hybrid automaton object
automaton=HybridAutomaton("Affine hybrid automaton")


#Create first location
#Discrete state
q1=DiscreteState(1);
#Dynamic
dynamic=AffineVectorField(Matrix([[-10,1],[0,-1]]),Vector([0,18]));
#Invariant
space=RectangularSet([[0,3],[0,40]]);

mode1=automaton.new_mode(q1,dynamic,space);

#Create second location
#Discrete state
q2=DiscreteState(2);
#Dinamico
dynamic=AffineVectorField(Matrix([[0,0],[0,-1]]),Vector([0,18]));
#Invariant
space=RectangularSet([[2.9,3],[0,40]]);

mode2=automaton.new_mode(q2,dynamic,space);

#Create transition
e=DiscreteEvent(3);
jump_set = RectangularSet([[2.9,3],[30,40]]);
reset = AffineMap(Matrix([[1,0],[0,1]]),Vector([0,0]));


transition=automaton.new_transition(e,q1,q2,reset,jump_set)

#Initial condition in continuous space
initial_rectangle=RectangularSet([[0,0.001],[0,36]]);

print "\nCreating initial hybrid set"
initial_set=HybridSet()
initial_set.new_location(q1,initial_rectangle)
initial_set.new_location(q2,EmptySet(2))
print initial_set

print "\nCreating bounding hybrid set"
bounding_box=Rectangle([[-0.1,3.4],[-0.1,42]]);
bounding_set=HybridSet()
bounding_set.new_location(q1,bounding_box)
bounding_set.new_location(q2,bounding_box)
print bounding_set

print "\nCreating evolution parameters"
params=EvolutionParameters()
params.set_lock_to_grid_time(0.5)
params.set_grid_length(0.5)
print params

print "\nCreating hybrid evolver"
applicator=StandardApplicator()
integrator=AffineIntegrator()
hybrid_evolver=HybridEvolver(params,applicator,integrator);

print "\nComputing chain reachable set"
res=hybrid_evolver.chainreach(automaton,initial_set)

print "\nExporting to postscript output...",
#
epsbb=Box([[0,5],[0,50]]) # eps bounding box
eps=EpsPlot()
eps.open("affine_hybrid_automaton-1.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("cyan")
eps.write(jump_set)
eps.set_fill_colour("yellow")
eps.write(res[q1])
eps.set_fill_colour("blue")
eps.write(initial_set[q1])
eps.close()

eps2=EpsPlot()
eps2.open("affine_hybrid_automaton-2.eps",epsbb)
eps2.set_line_style(True)
eps2.set_fill_colour("cyan")
eps2.write(jump_set)
eps2.set_fill_colour("yellow")
eps2.write(res[q2])
eps2.close()

print " done."

