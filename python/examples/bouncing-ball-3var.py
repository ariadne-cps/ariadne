#!/bin/python
#
# Bouncing ball example
#
# Written by Davide Bresolin on March, 2007
#

from ariadne import *

set_hybrid_evolver_verbosity(4)
#set_integrator_verbosity(6)

max_time="5"

# Definition of the automaton
automaton=HybridAutomaton("Bouncing ball")

# Location fly:
dyn=AffineVectorField(Matrix("[0, 1, 0; 0, 0, 0; 0, 0, 0]"), Vector("[0,-9.8,1]"))
inv=RectangularSet("[0,15]x[-20,20]x[0,100]")
# inv=PolyhedralSet(Polyhedron(Matrix("[-1,0]"),Vector("[0]")))
l1=automaton.new_mode(0, dyn, inv)

# Transition l1 -> l1
#act=Polyhedron(Matrix("[1,0; 0,1]"),Vector("[0,0]"))
act=RectangularSet("[0,0.01]x[-20,20]x[0,100]")
res=AffineMap(Matrix("[1,0,0;0,-1,0;0,0,1]"), Vector("[0,0,0]"))
e12=automaton.new_transition(0,l1.id(),l1.id(),res,act)

print automaton

# Defintion of the underlying grid
grid=Grid(Vector("[0.25,0.25,0.25]"))
bounding_box=Rectangle("[0,25]x[-25,25]x[0,"+max_time+"]")
bounding_box=bounding_box.neighbourhood(Float(0.5));
extent=bounding_box.neighbourhood(Float(1));
fgrid=FiniteGrid(grid,extent)

# Initial set 
init=Rectangle("[10,10.1]x[0,0.1]x[0,0.1]")
initial_set=HybridGridMaskSet()
initial_set.new_location(l1.id(),fgrid)
initial_set[l1.id()].adjoin_over_approximation(init)

print "initial_set.locations() =",initial_set.locations()
#print "initial_set[l0.id].size(),capacity()=",initial_set[l0.id()].size(),initial_set[l0.id()].capacity()

# Bounding set 
bounding_set=HybridGridMaskSet()
bounding_set.new_location(l1.id(),fgrid)
bounding_set[l1.id()].adjoin_over_approximation(bounding_box)
print "bounding_set.locations() =",bounding_set.locations()
#print "bounding_set[mode1_id].size(),capacity()=",bounding_set[m1.id()].size(),bounding_set[m1.id()].capacity()
#print "bounding_set[mode2_id].size(),capacity()=",bounding_set[m2.id()].size(),bounding_set[m2.id()].capacity()

# Definition of the Hybrid Evolver
maximum_step_size=0.25;
lock_to_grid_time=1.0;
maximum_set_radius=0.5;

apply=Applicator()
integrator=AffineIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
#integrator=LohnerIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
hybrid_evolver=HybridEvolver(apply,integrator);

print "Computing chainreachable set..."
#reach_set=hybrid_evolver.continuous_chainreach(automaton,initial_set,bounding_set)
reach_set=hybrid_evolver.chainreach(automaton,initial_set,bounding_set)


#print "Computing discrete transition..."
#discrete_init=HybridGridCellListSet(reach_set)
#discrete_step=hybrid_evolver.discrete_step(automaton,discrete_init)

print "Done."

print "reach_set =",reach_set
print "reach_set.size() =",reach_set.size()
print "reach_set.capacity() =",reach_set.capacity()

# Eps output
eps=EpsPlot()
eps.open("bouncing-ball.eps",bounding_box,0,2)

# Defintion of the reference grid
#eps.set_line_style(False)
#eps.set_fill_colour("yellow")
#eps.write(Rectangle("[0,5]x[0,5]x[0,5]"))

# Write the safe area
#eps.set_fill_colour("yellow")
#eps.write(Rectangle("[0,20]x[1,12]"))

eps.set_fill_colour("green")
eps.write(reach_set[l1.id()])

# Write the initial set
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


  	

