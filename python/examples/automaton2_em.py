#!/usr/bin/python

from ariadne import *
import sys

#print dir()
invariant1=RectangularSet([[0,8],[0,30]])
dynamic1=AffineVectorField(Matrix([[0,0],[0,-1]]),Vector([1,25]))

#locazione 2
invariant2=RectangularSet([[6,20],[0,30]])
dynamic2=AffineVectorField(Matrix([[0,0],[0,-2]]),Vector([1,5]))

#creazione automa ibrido
automaton=HybridAutomaton("Affine hybrid automaton")

mode1_id=DiscreteState(1)
mode1=automaton.new_mode(mode1_id,dynamic1,invariant1)

mode2_id=DiscreteState(2)
mode2=automaton.new_mode(mode2_id,dynamic2,invariant2)

activation12 = RectangularSet([[7,20],[0,30]])
reset12 = AffineMap(Matrix([[1,0],[0,1]]),Vector([0,0]))

event1_id=DiscreteEvent(3)
transition=automaton.new_transition(event1_id,mode1_id,mode2_id,reset12,activation12)

print automaton

#definizione condizioni iniziali
initial_rectangle1=Rectangle([[7,7.01],[22,22.01]]);
initial_rectangle2=Rectangle([[12,12.01],[10,10.01]]);

#definizione spazio di lavoro e griglia.
bounding_box=Box([[-1,30],[-1,30]])

grid=Grid(Vector([0.15,0.15]))
#block=LatticeBlock([[-25,500],[-25,500]])
fgrid=FiniteGrid(grid,bounding_box)
print fgrid


print "Creating initial hybrid set"
initial_set=HybridSet()

initial_set.new_location(mode1_id,initial_rectangle1)
initial_set.new_location(mode2_id,initial_rectangle2)
#initial_set.new_location(mode2_id,EmptySet(2))
print "initial_set.locations() =",initial_set.locations()


# set evolution parameters
par=EvolutionParameters()
par.set_maximum_enclosure_radius(1.25);
par.set_maximum_step_size(1.25);
par.set_lock_to_grid_time(1.5);
par.set_grid_length(0.5)
par.set_argument_grid_length(0.25)
par.set_result_grid_length(0.25)
par.set_verbosity(2)

print par

apply=StandardApplicator()
#integrator=AffineIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
integrator=LohnerIntegrator();
#hybrid_evolver=HybridEvolver(apply,integrator);
hybrid_evolver=SetBasedHybridEvolver(par,apply,integrator);


#print "Computing continuous chainreach set"
#set_applicator_verbosity(4)
set_integrator_verbosity(0)
set_hybrid_evolver_verbosity(3)
#continuous_chainreach_set=hybrid_evolver.continuous_chainreach(automaton,initial_set,bounding_set)
print

#print "Computing single discrete step..."
#discrete_step_set=hybrid_evolver.discrete_step(automaton,initial_set)
print

print "Computing lower evolve set for 1 seconds..."
final_set=hybrid_evolver.lower_evolve(automaton,initial_set,1)
print "Computing lower reach set for 1 seconds..."
chainreach_set=hybrid_evolver.lower_reach(automaton,initial_set,1)

print "Sizes of the computed regions:"
print "Mode 1:",chainreach_set[mode1_id].size()
print "Mode 2:",chainreach_set[mode2_id].size()

print "Exporting to postscript output...",
epsbb=Box([[-5,30],[-5,30]]) # eps bounding box
eps=EpsPlot()
eps.open("automaton2_em_lower.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("red")
eps.write(invariant1)
eps.set_fill_colour("yellow")
eps.write(invariant2)
eps.set_fill_colour("green")
eps.write(chainreach_set[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.set_fill_colour("cyan")
eps.write(chainreach_set[mode2_id])
eps.set_fill_colour("black")
eps.write(final_set[mode1_id])
eps.write(final_set[mode2_id])


eps.close()

print " done."


print "Computing upper reach set for 1 seconds..."
chainreach_set=hybrid_evolver.upper_reach(automaton,initial_set,1)

print "Sizes of the computed regions:"
print "Mode 1:",chainreach_set[mode1_id].size()
print "Mode 2:",chainreach_set[mode2_id].size()

print "Exporting to postscript output...",
epsbb=Box([[-5,30],[-5,30]]) # eps bounding box
eps=EpsPlot()
eps.open("automaton2_em_upper.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("red")
eps.write(invariant1)
eps.set_fill_colour("yellow")
eps.write(invariant2)
eps.set_fill_colour("green")
eps.write(chainreach_set[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.set_fill_colour("cyan")
eps.write(chainreach_set[mode2_id])
eps.set_fill_colour("black")
#eps.write(initial_set[mode2_id])


eps.close()

print " done."
