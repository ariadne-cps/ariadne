#!/usr/bin/python
#
# Bouncing ball example
#
# Written by Davide Bresolin on March, 2007
#

from ariadne import *

#set_hybrid_evolver_verbosity(6)
#set_integrator_verbosity(6)

max_time="5"
g="-1"

# Definition of the automaton
automaton=HybridAutomaton("Bouncing ball")

# Location fly:
mode_id=0
dyn=AffineVectorField(Matrix("[0, 1, 0; 0, 0, 0; 0, 0, 0]"), Vector("[0,"+g+",1]"))
inv=RectangularSet("[0,2]x[-2,1]x[0,"+max_time+"]")
# inv=PolyhedralSet(Polyhedron(Matrix("[-1,0]"),Vector("[0]")))
l1=automaton.new_mode(mode_id, dyn, inv)

# Transition l1 -> l1
#act=Polyhedron(Matrix("[1,0; 0,1]"),Vector("[0,0]"))
event_id=1
act=RectangularSet("[-0.1,0.0]x[-2,0]x[0,"+max_time+"]")
res=AffineMap(Matrix("[1,0,0;0,-0.5,0;0,0,1]"), Vector("[0,0,0]"))
e12=automaton.new_transition(event_id,l1.id(),l1.id(),res,act)

print automaton

# Defintion of the underlying grid
grid=Grid(Vector("[0.01,0.01,0.01]"))
block=LatticeBlock("[-20,220]x[-220,120]x[-20,520]")
fgrid=FiniteGrid(grid,block)

# Initial set 
init=Rectangle("[1.999,2]x[0,0.001]x[0,0.001]")
#init=Rectangle("[0,0.001]x[0.999,1.0]x[0,0.001]")
initial_set=HybridGridMaskSet()
initial_set.new_location(l1.id(),fgrid)
initial_set[l1.id()].adjoin_over_approximation(init)

print "initial_set.locations() =",initial_set.locations()
#print "initial_set[l0.id].size(),capacity()=",initial_set[l0.id()].size(),initial_set[l0.id()].capacity()

# Bounding set
bounding_box=Rectangle("[-0.1,2.1]x[-2.1,1.1]x[-0.1,"+max_time+".1]") 
bounding_set=HybridGridMaskSet()
bounding_set.new_location(l1.id(),fgrid)
bounding_set[l1.id()].adjoin_over_approximation(bounding_box)
print "bounding_set.locations() =",bounding_set.locations()
#print "bounding_set[mode1_id].size(),capacity()=",bounding_set[m1.id()].size(),bounding_set[m1.id()].capacity()
#print "bounding_set[mode2_id].size(),capacity()=",bounding_set[m2.id()].size(),bounding_set[m2.id()].capacity()

# Definition of the Hybrid Evolver
maximum_step_size=0.1;
lock_to_grid_time=2.0;
maximum_set_radius=0.5;

apply=Applicator()
integrator=AffineIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
#integrator=LohnerIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
hybrid_evolver=HybridEvolver(apply,integrator);

print "Computing reach set..."
#reach_set=hybrid_evolver.continuous_chainreach(automaton,initial_set,bounding_set)
reach_set=hybrid_evolver.chainreach(automaton,initial_set,bounding_set)

print "Done."

print "reach_set[l1.id].size() =",reach_set[l1.id()].size()
print "reach_set[l1.id].capacity() =",reach_set[l1.id()].capacity()

#print "Computing discrete transition..."
#trans_set=hybrid_evolver.discrete_step(automaton,reach_set)
#print "Done."

#print "trans_set[l1.id].size() =",trans_set[l1.id()].size()
#print "trans_set[l1.id].capacity() =",trans_set[l1.id()].capacity()

#print "Computing second continuous reach set..."
#reach_set2=hybrid_evolver.continuous_chainreach(automaton,trans_set,bounding_set)
#print "Done."

#print "reach_set2[l1.id].size() =",reach_set2[l1.id()].size()
#print "reach_set2[l1.id].capacity() =",reach_set2[l1.id()].capacity()

# Txt output
print "Creating txt output..."
txt=TxtPlot()
txt.open("bouncing-ball-3var.txt")

# Write the reached set
txt.write(reach_set[l1.id()])
txt.close()

print "Done."

# Eps output
print "Creating eps output..."
eps=EpsPlot()
eps.open("bouncing-ball-3var.eps", bounding_box,2,0)

# Write the reached set
eps.write(reach_set[l1.id()])
eps.close()

print "Done."





  	

