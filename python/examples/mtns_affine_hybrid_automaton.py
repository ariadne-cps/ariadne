#!/usr/bin/python

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

print dir()

space=RectangularSet([[-2.0,2.0],[-2.0,2.0]])
identity=AffineMap(Matrix([[1,0],[0,1]]),Vector([0,0]))

invariant=RectangularSet([[-1.0,1.0],[-1.0,1.0]])
dynamic=AffineVectorField(Matrix([[-0.5,-1.0],[ 1.0,-0.5]]),Vector([0.00,0.00]))
activation12 = RectangularSet([[-0.2,-0.1],[-0.2,-0.1]])
activation21 = RectangularSet([[0.1,0.2],[0.1,0.2]])
reset = AffineMap(Matrix([[-7,0],[0,-7]]),Vector([0,0]))

automaton=HybridAutomaton("Affine hybrid automaton")
mode1_id=DiscreteState(0)
mode2_id=DiscreteState(1)
mode1=automaton.new_mode(mode1_id,dynamic,invariant)
mode2=automaton.new_mode(mode2_id,dynamic,invariant)
event_id=DiscreteEvent(2)
transition=automaton.new_transition(event_id,mode1_id,mode2_id,reset,activation12)
transition=automaton.new_transition(event_id,mode2_id,mode1_id,reset,activation21)
print automaton
#print automaton.invariant()

initial_rectangle=RectangularSet([[-0.8,-0.7],[0.7,0.8]])
bounding_box=RectangularSet([[-2,2],[-2,2]])
epsbb=RectangularSet([[-2.1,2.1],[-2.1,2.1]]) # eps bounding box

v = Vector(2);
v[0] = 1.0/8;
v[1] = 1.0/8;
grid=Grid(v)

print "Creating initial hybrid set"
initial_set=HybridSet()
initial_set.new_location(mode1_id,initial_rectangle)
initial_set.new_location(mode2_id,EmptySet(2))
print "initial_set.locations() =",initial_set.locations()

print "Creating bounding hybrid set",
bounding_set=HybridSet();
bounding_set.new_location(mode1_id,bounding_box)
bounding_set.new_location(mode2_id,bounding_box)
print "bounding_set.locations() =",bounding_set.locations()
print

par = EvolutionParameters()
par.set_maximum_step_size(1.0/32);
par.set_minimum_step_size(1.0/128);
par.set_lock_to_grid_time(1);
par.set_grid(grid)
par.set_hybrid_bounding_domain(bounding_set)
par.set_verbosity(2)

apply=StandardApplicator()
integrator=AffineIntegrator();
hybrid_evolver=SetBasedHybridEvolver(par,apply,integrator);


print "Computing chainreach set"
chainreach_set=hybrid_evolver.chainreach(automaton,initial_set)

print "Exporting to postscript output...",
eps=EpsPlot()
eps.open("mtns_affine_hybrid_automaton-2.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("cyan")
eps.write(activation21)
eps.set_fill_colour("magenta")
eps.write(activation12)
eps.set_fill_colour("red")
eps.write(chainreach_set[mode2_id])
eps.set_fill_colour("yellow")
eps.write(chainreach_set[mode1_id])
eps.set_fill_colour("green")
eps.write(continuous_chainreach_set[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.close()
print " done."
