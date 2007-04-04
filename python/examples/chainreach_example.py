#!/usr/bin/python
#
# continuous_chainreach example
#
# Written by Davide Bresolin on April, 2007
#

import sys
from ariadne import *

# Definition of the automaton
automaton=HybridAutomaton("Stupid Automaton")

# Only one discrete mode
dyn=AffineVectorField(Matrix("[0, 0; 0, 0]"), Vector("[1,1]"))
inv=Rectangle("[0,1]x[0,1]")
m0=automaton.new_mode(0, dyn, PolyhedralSet(inv))

print automaton

#bounding box used for Postscript output
bounding_box=Rectangle("[-1,4]x[-1,4]")

# Defintion of the underlying grid
grid=Grid(Vector("[0.25,0.25]"))
block=LatticeBlock("[-2,14]x[-2,14]")
fgrid=FiniteGrid(grid,block)

# Initial set 
init=Rectangle("[0,0.1]x[0,0.1]")
initial_set=HybridGridMaskSet()
initial_set.new_location(m0.id(),fgrid)
initial_set[m0.id()].adjoin_over_approximation(init)

print "initial_set.locations() =",initial_set.locations()
print "initial_set[m0.id].size(),capacity()=",initial_set[m0.id()].size(),initial_set[m0.id()].capacity()

# Definition of the Hybrid Evolver
maximum_step_size=0.125;
lock_to_grid_time=0.5;
maximum_set_radius=0.25;

apply=Applicator()
integrator=AffineIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
#integrator=LohnerIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
hybrid_evolver=HybridEvolver(apply,integrator);

# FIRST example: continuous chainreach with bounding_box [0,5]x[0,5]
bound=Rectangle("[0,3]x[0,3]")
bounding_set=HybridGridMaskSet()
bounding_set.new_location(m0.id(),fgrid)
bounding_set[m0.id()].adjoin_over_approximation(bound)


print "Computing continuous chainreachable set with bounding box [0,5]x[0,5]..."
reach_set=hybrid_evolver.continuous_chainreach(automaton,initial_set,bounding_set)

print "Done."

# Eps output
eps=EpsPlot()
eps.open("chainreach_example1.eps",bounding_box,0,1)

# Write the bound
eps.set_fill_colour("cyan")
eps.write(bound)

# Write the invariant
eps.set_fill_colour("yellow")
eps.write(inv)

# Write the reached set
eps.set_fill_colour("green")
eps.write(reach_set[m0.id()])

# Write the initial set
eps.set_fill_colour("blue")
eps.write(initial_set[m0.id()])

eps.write(fgrid)

eps.close()

print "Computing chainreachable set with bounding box [0,5]x[0,5]..."
reach_set=hybrid_evolver.chainreach(automaton,initial_set,bounding_set)

print "Done."

# Eps output
eps.open("chainreach_example2.eps",bounding_box,0,1)

# Write the bounding box
#eps.set_fill_colour("white")
#eps.write(bounding_box)

# Write the bound
eps.set_fill_colour("cyan")
eps.write(bound)

# Write the invariant
eps.set_fill_colour("yellow")
eps.write(inv)

# Write the reached set
eps.set_fill_colour("green")
eps.write(reach_set[m0.id()])

# Write the initial set
eps.set_fill_colour("blue")
eps.write(initial_set[m0.id()])
eps.write(init)

eps.write(fgrid)

eps.close()

print "Computing chainreachable set with bounding box",inv,"using Integrator..."
block=LatticeBlock("[-2,6]x[-2,6]")
fgrid=FiniteGrid(grid,block)
ginit=GridMaskSet(fgrid)
ginit.adjoin_over_approximation(init)
ginv=GridMaskSet(fgrid)
ginv.adjoin_over_approximation(inv)

set_evaluation_verbosity(4)

reach=integrator.chainreach(dyn,ginit,ginv)

# Eps output
eps.open("chainreach_example3.eps",bounding_box,0,1)

# Write the invariant
eps.set_fill_colour("yellow")
eps.write(ginv)
eps.write(inv)

# Write the reached set
eps.set_fill_colour("green")
eps.write(reach)

# Write the initial set
eps.set_fill_colour("blue")
eps.write(ginit)
eps.write(init)

#eps.write(fgrid)

eps.close()

print
