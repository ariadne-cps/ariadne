#!/usr/bin/python
#
# continuous_chainreach example
#
# Written by Davide Bresolin on April, 2007
#

import sys
from ariadne import *

set_evaluation_verbosity(2)

# Definition of the automaton
automaton=HybridAutomaton("Stupid Automaton")

# Only one discrete mode
m0_id=DiscreteState(0)
dyn=AffineVectorField(Matrix([[0, 0],[0, 0]]), Vector([1,1]))
inv=RectangularSet([[0,1],[0,1]])
inv=PolyhedralSet(Matrix([[-1,0], [0,-1]]),Vector([1,1]));
m0=automaton.new_mode(m0_id, dyn, inv)

print automaton

#bounding box used for Postscript output
bounding_box=RectangularSet([[-2,4],[-2,4]])

# Defintion of the underlying grid
grid=Grid(Vector([0.25,0.25]))

# Initial set 
init=RectangularSet([[0,0.1],[0,0.1]])
initial_set=HybridSet()
initial_set.new_location(m0_id,init)

print "initial_set.locations() =",initial_set.locations()

# Bounding set
bound=RectangularSet([[-1.5,3],[-1.5,3]])
bounding_set=HybridSet()
bounding_set.new_location(m0_id,bound)


# Definition of the Hybrid Evolver
par=EvolutionParameters();
par.set_maximum_step_size(0.125);
par.set_lock_to_grid_time(10);
par.set_grid(grid);
par.set_hybrid_bounding_domain(bounding_set);

apply=StandardApplicator()
integrator=AffineIntegrator();
#integrator=LohnerIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
hybrid_evolver=SetBasedHybridEvolver(par, apply,integrator);

# FIRST example: continuous chainreach with bounding_box [0,5],[0,5]


print "Computing chainreachable set with bounding box [0,5],[0,5]..."
reach_set=hybrid_evolver.chainreach(automaton,initial_set)

print "Done."

# Eps output
eps=EpsPlot()
eps.open("chainreach_example1.eps",bounding_box,0,1)

# Write the bound
eps.set_fill_colour("cyan")
eps.write(bound)

# Write the invariant
eps.set_fill_colour("yellow")
#eps.write(closed_intersection(inv,bounding_box))

# Write the reached set
eps.set_fill_colour("green")
eps.write(reach_set[m0_id])

# Write the initial set
eps.set_fill_colour("blue")
eps.write(initial_set[m0_id])

#eps.write(fgrid)

eps.close()

print "Computing lower reachable set for 5 seconds..."
reach_set=hybrid_evolver.lower_reach(automaton,initial_set,5)

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
#eps.write(closed_intersection(inv,Polyhedron(bound)))

# Write the reached set
eps.set_fill_colour("green")
eps.write(reach_set[m0_id])

# Write the initial set
eps.set_fill_colour("blue")
eps.write(initial_set[m0_id])
eps.write(init)

#eps.write(fgrid)

eps.close()

print "Computing chainreachable set with bounding box [0,5],[0,5] using Integrator..."

vfevolver = VectorFieldEvolver(par,integrator);

#redirect_log("chainreach_example.log")
reach=vfevolver.chainreach(dyn,init,Box([[0,5],[0,5]]))


# Eps output
eps.open("chainreach_example3.eps",bounding_box,0,1)
# Write the invariant
eps.set_fill_colour("yellow")
eps.write(inv)
#eps.write(closed_intersection(inv,bound))

# Write the reached set
eps.set_fill_colour("green")
eps.write(reach)

# Write the initial set
eps.set_fill_colour("blue")
#eps.write(ginit)
eps.write(init)

#eps.write(fgrid)
eps.close()

print "Done."

