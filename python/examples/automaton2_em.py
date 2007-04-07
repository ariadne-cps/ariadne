#!/projects/ariadne/bin/python2.4

from ariadne import *

#print dir()

space=Rectangle("[-5,30]x[-5,30]")

#locazione 1
invariant1=RectangularSet("[0,8]x[0,30]")
dynamic1=AffineVectorField(Matrix("[0,0;0,-1]"),Vector("[1,25]"))

#locazione 2
invariant2=RectangularSet("[8,20]x[0,30]")
dynamic2=AffineVectorField(Matrix("[0,0;0,-2]"),Vector("[1,5]"))

#creazione automa ibrido
automaton=HybridAutomaton("Affine hybrid automaton")

mode1_id=1
mode1=automaton.new_mode(mode1_id,dynamic1,invariant1)

mode2_id=2
mode2=automaton.new_mode(mode2_id,dynamic2,invariant2)

activation12 = RectangularSet("[8,20]x[0,30]")
reset12 = AffineMap(Matrix("[1,0;0,1]"),Vector("[0,0]"))

event1_id=3
transition=automaton.new_transition(event1_id,mode1_id,mode2_id,reset12,activation12)

print automaton

#definizione condizioni iniziali
initial_rectangle1=Rectangle("[2,2.01]x[5,5.01]");
initial_rectangle2=Rectangle("[12,12.01]x[10,10.01]");

#definizione spazio di lavoro e griglia.
bounding_box=Rectangle("[-1,30]x[-1,30]")

grid=Grid(Vector("[0.15,0.15]"))
block=LatticeBlock("[-25,500]x[-25,500]")
fgrid=FiniteGrid(grid,block)
print fgrid


print "Creating initial hybrid set"
initial_set=HybridGridMaskSet()

initial_set.new_location(mode1_id,fgrid)
initial_set[mode1_id].adjoin_over_approximation(initial_rectangle1)
print "initial_set.locations() =",initial_set.locations()
print "initial_set[mode1_id].size(),capacity()=",initial_set[mode1_id].size(),initial_set[mode1_id].capacity()

initial_set.new_location(mode2_id,fgrid)
initial_set[mode1_id].adjoin_over_approximation(initial_rectangle2)
print "initial_set.locations() =",initial_set.locations()
print "initial_set[mode1_id].size(),capacity()=",initial_set[mode2_id].size(),initial_set[mode2_id].capacity()


print "Creating initial hybrid cell list set"
initial_cell_list_set=HybridGridCellListSet(initial_set)

print "initial_cell_list_set.locations() =",initial_cell_list_set.locations()

print "Creating bounding hybrid set",
bounding_set=HybridGridMaskSet();
bounding_set.new_location(mode1_id,fgrid)
bounding_set.new_location(mode2_id,fgrid)
bounding_set[mode1_id].adjoin_over_approximation(bounding_box);
bounding_set[mode2_id].adjoin_over_approximation(bounding_box);
print "bounding_set.locations() =",bounding_set.locations()
print "bounding_set[mode1_id].size(),capacity()=",bounding_set[mode1_id].size(),bounding_set[mode1_id].capacity()
print

print "bounding_set.locations() =",bounding_set.locations()
print "bounding_set[mode1_id].size(),capacity()=",bounding_set[mode2_id].size(),bounding_set[mode2_id].capacity()
print


maximum_step_size=0.05;
lock_to_grid_time=1.5;
maximum_set_radius=0.5;

apply=Applicator()
#integrator=AffineIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
integrator=LohnerIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
hybrid_evolver=HybridEvolver(apply,integrator);


print "Computing continuous chainreach set"
#set_applicator_verbosity(4)
set_integrator_verbosity(6)
set_hybrid_evolver_verbosity(4)
continuous_chainreach_set=hybrid_evolver.continuous_chainreach(automaton,initial_set,bounding_set)
print

print "Computing single discrete step..."
discrete_step_set=hybrid_evolver.discrete_step(automaton,initial_cell_list_set)
print

print "discrete_step_set",discrete_step_set
print

print "Computing chainreach set..."
chainreach_set=hybrid_evolver.chainreach(automaton,initial_set,bounding_set)

print "Exporting to postscript output...",
epsbb=space # eps bounding box
eps=EpsPlot()
eps.open("automaton2_em.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("red")
eps.write(invariant1)
eps.set_fill_colour("yellow")
eps.write(invariant2)
eps.set_fill_colour("green")
eps.write(continuous_chainreach_set[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.set_fill_colour("brown")
eps.write(continuous_chainreach_set[mode2_id])
eps.set_fill_colour("black")
eps.write(initial_set[mode2_id])


eps.close()

print " done."
