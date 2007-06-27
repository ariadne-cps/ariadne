#!/bin/python
#
# Bouncing ball example
#
# Written by Davide Bresolin on March, 2007
#

from ariadne import *

# Definition of the automaton
automaton=HybridAutomaton("Bouncing ball")

# Location fly:
dyn=AffineVectorField(Matrix("[0, 1; 0, 0]"), Vector("[0,-1.0]"))
inv=RectangularSet("[0,10]x[-10,10]")
# inv=PolyhedralSet(Polyhedron(Matrix("[-1,0]"),Vector("[0]")))
l1=automaton.new_mode(0, dyn, inv)

# Transition l1 -> l1
act=RectangularSet("[0,0.001]x[-10,0]")
res=AffineMap(Matrix("[1,0;0,-0.5]"), Vector("[0,0]"))
e12=automaton.new_transition(0,l1.id(),l1.id(),res,act)

print automaton

# Defintion of the underlying grid
v=Vector(2)
v[0]=1.0/32
v[1]=1.0/32
grid=Grid(v)
block=LatticeBlock("[-80,200]x[-200,200]")
fgrid=FiniteGrid(grid,block)

# Initial set 
init=Rectangle("[1.999,2.0]x[0,0.001]")
initial_set=HybridGridMaskSet()
initial_set.new_location(l1.id(),fgrid)
initial_set[l1.id()].adjoin_over_approximation(init)

print "initial_set.locations() =",initial_set.locations()
#print "initial_set[l0.id].size(),capacity()=",initial_set[l0.id()].size(),initial_set[l0.id()].capacity()

# Bounding set 
bounding_box=Rectangle("[-1,2.5]x[-2.5,2.5]")
bounding_set=HybridGridMaskSet()
bounding_set.new_location(l1.id(),fgrid)
bounding_set[l1.id()].adjoin_over_approximation(bounding_box)
print "bounding_set.locations() =",bounding_set.locations()
#print "bounding_set[mode1_id].size(),capacity()=",bounding_set[m1.id()].size(),bounding_set[m1.id()].capacity()
#print "bounding_set[mode2_id].size(),capacity()=",bounding_set[m2.id()].size(),bounding_set[m2.id()].capacity()

# Definition of the Hybrid Evolver
maximum_step_size=1.0/8
lock_to_grid_time=0.5
maximum_set_radius=0.25

print maximum_step_size,lock_to_grid_time,maximum_set_radius;

apply=Applicator()
integrator=AffineIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
#integrator=LohnerIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
hybrid_evolver=HybridEvolver(apply,integrator);

set_hybrid_evolver_verbosity(4)
#set_applicator_verbosity(7)
#set_integrator_verbosity(7)
#set_geometry_verbosity(7)

#print "Computing continuous chainreachable set..."
#continuous_chainreach_set=hybrid_evolver.continuous_chainreach(automaton,initial_set,bounding_set)
continuous_chainreach_set=initial_set

print "Computing chainreachable set..."
chainreach_set=hybrid_evolver.chainreach(automaton,initial_set,bounding_set)

#print "Computing discrete transition..."
#discrete_init=HybridGridCellListSet(reach_set)
#discrete_step=hybrid_evolver.discrete_step(automaton,discrete_init)

print "Done."

print "continuous_chainreach_set =",continuous_chainreach_set[l1.id()]
print "chainreach_set=",chainreach_set[l1.id()]

# Eps output
eps=EpsPlot()
eps.open("bouncing-ball.eps",bounding_box,0,1)

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
eps.write(chainreach_set[l1.id()])
eps.set_fill_colour("yellow")
eps.write(continuous_chainreach_set[l1.id()])

# Write the initial set
eps.set_line_style(False)
eps.set_fill_colour("blue")
eps.write(initial_set[l1.id()])
#eps.write(act12)
#eps.write(act23)
#eps.write(act30)
#eps.write(discrete_step[l1.id()])
#eps.write(discrete_step[l2.id()])
#eps.write(discrete_step[l3.id()])

eps.close()

print

# Txt output
txt=TxtPlot()
txt.open("bouncing-ball.txt")
txt.write(chainreach_set[l1.id()])
txt.close()


   

  	

