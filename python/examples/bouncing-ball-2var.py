#!/usr/bin/python
#
# Bouncing ball example
#
# Written by Davide Bresolin on March, 2007
#

from ariadne import *
import sys

# Definition of the automaton
automaton=HybridAutomaton("Bouncing ball")

# Location fly:
mode_id=DiscreteState(0)
dyn=AffineVectorField(Matrix([[0, 1],[ 0, 0]]), Vector([0,-1.0]))
inv=RectangularSet([[0,10],[-10,10]])
# inv=PolyhedralSet(Polyhedron(Matrix([[-1,0]]),Vector([0])))
l1=automaton.new_mode(mode_id, dyn, inv)

# Transition l1 -> l1
event_id=DiscreteEvent(1)
act=RectangularSet([[0,0.001],[-10,0]])
res=AffineMap(Matrix([[1,0],[0,-0.5]]), Vector([0,0]))
e12=automaton.new_transition(event_id,mode_id,mode_id,res,act)

print automaton


# Initial set 
init=Rectangle([[1.999,2.0],[0,0.001]])
initial_set=HybridSet()
initial_set.new_location(mode_id,init)

print "initial_set.locations() =",initial_set.locations()

# Bounding set 
bounding_box=RectangularSet([[-1,2.5],[-2.5,2.5]])
#bounding_set=HybridSet()
#bounding_set.new_location(mode_id,bounding_box)
#print "bounding_set.locations() =",bounding_set.locations()

# Definition of the Hybrid Evolver
par = EvolutionParameters()
par.set_maximum_step_size(0.5)
par.set_lock_to_grid_time(4)
par.set_grid_length(1.0/32)
par.set_bounding_domain_size(10.0)
par.set_verbosity(2)


print par

apply=StandardApplicator()
integrator=AffineIntegrator();
#integrator=LohnerIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
#hybrid_evolver=HybridEvolver(apply,integrator);
hybrid_evolver=SetBasedHybridEvolver(par,apply,integrator);

time=4

print "initial set=",initial_set

print "Computing upper reachable set for",time,"seconds..."
chain_reach_set=hybrid_evolver.upper_reach(automaton,initial_set,time)
evolve_set=hybrid_evolver.upper_evolve(automaton,initial_set,time)


print "Exporting to eps..."
eps=EpsPlot()
eps.open("bouncing-ball-upper.eps",bounding_box,0,1)

# Print the bounding_box
eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(bounding_box)

# Write the activation area
eps.set_fill_colour("yellow")
eps.write(act)

# Write the reached set
eps.set_line_style(False)
eps.set_fill_colour("green")
eps.write(chain_reach_set[mode_id])
eps.set_fill_colour("yellow")
eps.write(evolve_set[mode_id])

# Write the initial set
eps.set_line_style(False)
eps.set_fill_colour("blue")
eps.write(initial_set[mode_id])
#eps.write(act12)
#eps.write(act23)
#eps.write(act30)
#eps.write(discrete_step[mode_id])
#eps.write(discrete_step[l2.id()])
#eps.write(discrete_step[l3.id()])

eps.close()


print "initial set=",initial_set

print "Computing lower reachable set for",time,"seconds..."
lower_reach_set=hybrid_evolver.lower_reach(automaton,initial_set,time)
lower_evolve_set=hybrid_evolver.lower_evolve(automaton,initial_set,time)

print "Exporting to eps..."
eps=EpsPlot()
eps.open("bouncing-ball-lower.eps",bounding_box,0,1)

# Print the bounding_box
eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(bounding_box)

# Write the activation area
eps.set_fill_colour("yellow")
eps.write(act)

# Write the reached set
eps.set_line_style(False)
eps.set_fill_colour("green")
eps.write(lower_reach_set[mode_id])
eps.set_fill_colour("yellow")
eps.write(lower_evolve_set[mode_id])

# Write the initial set
eps.set_line_style(False)
eps.set_fill_colour("blue")
eps.write(initial_set[mode_id])
#eps.write(act12)
#eps.write(act23)
#eps.write(act30)
#eps.write(discrete_step[mode_id])
#eps.write(discrete_step[l2.id()])
#eps.write(discrete_step[l3.id()])

eps.close()

print "Computing chain reachable set..."
chainreach_set=hybrid_evolver.chain_reach(automaton,initial_set)

print "Exporting to eps..."

# Eps output
eps.open("bouncing-ball-chain.eps",bounding_box,0,1)

# Print the bounding_box
eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(bounding_box)

# Write the activation area
eps.set_fill_colour("yellow")
eps.write(act)

# Write the reached set
eps.set_line_style(False)
eps.set_fill_colour("green")
eps.write(chainreach_set[mode_id])

# Write the initial set
eps.set_line_style(False)
eps.set_fill_colour("blue")
eps.write(initial_set[mode_id])
#eps.write(act12)
#eps.write(act23)
#eps.write(act30)
#eps.write(discrete_step[mode_id])
#eps.write(discrete_step[l2.id()])
#eps.write(discrete_step[l3.id()])

eps.close()




   

  	

