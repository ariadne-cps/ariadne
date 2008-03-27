#!/usr/bin/python

##############################################################################
#            affine_hybrid_automaton.py
#
#  Copyright 2005-7  Alberto Casagrande, Pieter Collins
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

print dir()


space=RectangularSet([[-7.5,7.5],[-7.5,7.5]])
identity=AffineMap(Matrix([[1,0],[0,1]]),Vector([0,0]))

invariant1=RectangularSet([[-7.5,7.5],[-7.5,2.5]])
invariant2=RectangularSet([[-7.5,7.5],[-7.5,7.5]])
dynamic1=AffineVectorField(Matrix([[-2,-1],[1,-2]]),Vector([-1,0]))
dynamic2=AffineVectorField(Matrix([[-2,-1],[1,-2]]),Vector([1,0]))
activation11 = RectangularSet([[-7.5,7.5],[-3,-2]])
activation21 = RectangularSet([[-7.5,7.5],[-7,-6]])
activation12 = RectangularSet([[-7.5,7.5],[1,2]])
reset11 = AffineMap(Matrix([[1,0],[0,-1]]),Vector([0,0]))
reset21 = AffineMap(Matrix([[1,0],[0,-1]]),Vector([0,0]))
reset12 = AffineMap(Matrix([[1,0],[0,-1]]),Vector([0,0]))

automaton=HybridAutomaton("Affine hybrid automaton")
mode1_id=DiscreteState(2)
mode2_id=DiscreteState(3)
#cinvariant1=ConstraintSet(invariant1)
vfdynamic1=VectorField(dynamic1)
#mode1=automaton.new_mode(mode1_id,vfdynamic1,cinvariant1)
mode1=automaton.new_mode(mode1_id,dynamic1,invariant1)
mode2=automaton.new_mode(mode2_id,dynamic2,invariant2)
event1_id=DiscreteEvent(5)
event2_id=DiscreteEvent(7)
#transition=automaton.new_transition(event1_id,mode1_id,mode1_id,reset11,activation11)
#transition=automaton.new_transition(event2_id,mode1_id,mode2_id,reset12,activation12)
#transition=automaton.new_transition(event1_id,mode2_id,mode1_id,reset21,activation21)
print automaton
#print automaton.invariant()

initial_rectangle1=RectangularSet([[-6.96875,-6.9375],[-6.96875,-6.9375]]);
initial_rectangle2=RectangularSet([[6.9375,6.96875],[6.9375,6.96875]])


print "Creating initial hybrid set"
initial_set=HybridSet()
initial_set.new_location(mode1_id,initial_rectangle1)
initial_set.new_location(mode2_id,initial_rectangle2)
print "initial_set.locations() =",initial_set.locations()

print "Creating hybrid grid"
grid=Grid(Vector([0.25,0.25]))
hgrid=HybridGrid()
hgrid.new_location(mode1_id,grid)
hgrid.new_location(mode2_id,grid)

parameters=EvolutionParameters()
#parameters.set_grid_length(0.125)
parameters.set_grid_length(0.05)
parameters.set_hybrid_grid(hgrid)
parameters.set_lock_to_grid_time(10);
parameters.set_maximum_step_size(0.125)
#parameters.set_maximum_enclosure_radius(0.25);
parameters.set_maximum_enclosure_radius(2.5);
parameters.set_verbosity(2);
parameters.set_bounding_domain_size(8);
print parameters
applicator=KuhnApplicator(3)
integrator=AffineIntegrator();
hybrid_evolver=SetBasedHybridEvolver(parameters,applicator,integrator);

time=5.0

set_geometry_verbosity(0)

print "Computing lower reach set"
print initial_set
lower_reach_set=hybrid_evolver.lower_reach(automaton,initial_set,time)

#sys.exit(0)

lower_evolve_set=hybrid_evolver.lower_evolve(automaton,initial_set,time)

print "Exporting to postscript output...",
epsbb=RectangularSet([[-8.1,8.1],[-8.1,8.1]]) # eps bounding box
eps=EpsPlot()
eps.open("affine_hybrid_automaton-lower_reach.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(lower_reach_set[mode1_id])
eps.write(lower_reach_set[mode2_id])
eps.set_fill_colour("red")
eps.write(lower_evolve_set[mode1_id])
eps.write(lower_evolve_set[mode2_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.write(initial_set[mode2_id])
eps.close()
print

print "Exporting to txt output..."
txt=TextFile()
txt.open("affine_hybrid_automaton-lower_reach.txt")
txt.write(lower_reach_set[mode1_id])
txt.write(lower_reach_set[mode2_id])
txt.close()


#sys.exit(0)

print "Computing upper reach set"
upper_reach_set=hybrid_evolver.upper_reach(automaton,initial_set,time)

print "Exporting to postscript output...",
eps.open("affine_hybrid_automaton-upper_reach.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("red")
eps.write(upper_reach_set[mode2_id])
eps.set_fill_colour("green")
eps.write(upper_reach_set[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.write(initial_set[mode2_id])
eps.close()
print

print "Exporting to txt output..."
txt=TextFile()
txt.open("affine_hybrid_automaton-upper_reach.txt")
txt.write(upper_reach_set[mode1_id])
txt.write(upper_reach_set[mode2_id])
txt.close()


print "Computing chainreach set"
chainreach_set=hybrid_evolver.chain_reach(automaton,initial_set)
#chainreach_set=initial_set

print "Exporting to postscript output...",
eps.open("affine_hybrid_automaton-chain_reach.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("cyan")
eps.write(activation11)
eps.set_fill_colour("cyan")
eps.write(activation21)
eps.set_fill_colour("magenta")
eps.write(activation12)
eps.set_fill_colour("red")
eps.write(chainreach_set[mode2_id])
eps.set_fill_colour("green")
eps.write(chainreach_set[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.write(initial_set[mode2_id])
eps.close()

print "Exporting to txt output..."
txt=TextFile()
txt.open("affine_hybrid_automaton-chain_reach.txt")
txt.write(chainreach_set[mode1_id])
txt.write(chainreach_set[mode2_id])
txt.close()


print " done."
