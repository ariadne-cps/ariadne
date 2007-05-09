#!/projects/ariadne/bin/python2.4

from ariadne import *

#print dir()

#Creazione automa ibrido
automaton=HybridAutomaton("Affine hybrid automaton")


#Creazione prima locazione
#sistema dinamico
dynamic=AffineVectorField(Matrix("[-10,1;0,-1]"),Vector("[0,18]"));
#invariante
space=PolyhedralSet(Rectangle("[0,3]x[0,40]"));

mode1_id=1;
mode2_id=2;
mode1=automaton.new_mode(mode1_id,dynamic,space);

#Creazione seconda locazione
#sistema dinamico
dynamic=AffineVectorField(Matrix("[0,0; 0,-1]"),Vector("[0,18]"));
#invariante
space=PolyhedralSet(Rectangle("[2.9,3]x[0,40]"));

mode2=automaton.new_mode(mode2_id,dynamic,space);

#Creazione transizione
jump_set = PolyhedralSet(Rectangle("[2.9,3]x[30,40]"));
reset = AffineMap(Matrix("[1,0;0,1]"),Vector("[0,0]"));

event_id=1;

transition=automaton.new_transition(event_id,mode1_id,mode2_id,reset,jump_set)

#Condizioni iniziali nello spazio continuo
initial_rectangle=Rectangle("[0,0.001]x[0,36]");

#definizione dello spazio continuo da analizzare
grid=Grid(Vector("[0.02,0.1]"));
block=LatticeBlock("[-10,200]x[-10,450]");

#definizione della griglia nel Lattice
fgrid=FiniteGrid(grid,block)

print "Creating initial hybrid set",
initial_set=HybridGridMaskSet()
initial_set.new_location(mode1_id,fgrid)
initial_set.new_location(mode2_id,fgrid)
initial_set[mode1_id].adjoin_over_approximation(initial_rectangle)
print "initial set"
print initial_set[mode1_id].size(),initial_set[mode1_id].capacity()

bounding_box=Rectangle("[-0.1,3.4]x[-0.1,42]");

bounding_box_set=HybridGridMaskSet()
bounding_box_set.new_location(mode1_id,fgrid)
bounding_box_set.new_location(mode2_id,fgrid)
bounding_box_set[mode1_id].adjoin_over_approximation(bounding_box)
bounding_box_set[mode2_id].adjoin_over_approximation(bounding_box)

print "bounding box set"
print bounding_box_set[mode1_id].size(),bounding_box_set[mode1_id].capacity()

apply=Applicator()
integrator=AffineIntegrator(0.125,0.5,0.25);
hybrid_evolver=HybridEvolver(apply,integrator);

#chain reach
res=hybrid_evolver.chainreach(automaton,initial_set,bounding_box_set)

print "Exporting to postscript output...",
#
epsbb=Rectangle("[0,5]x[0,50]") # eps bounding box
eps=EpsPlot()
eps.open("affine_hybrid_automaton-1.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("cyan")
eps.write(jump_set)
eps.set_fill_colour("yellow")
eps.write(res[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.close()

eps2=EpsPlot()
eps2.open("affine_hybrid_automaton-2.eps",epsbb)
eps2.set_line_style(True)
eps2.set_fill_colour("cyan")
eps2.write(jump_set)
eps2.set_fill_colour("yellow")
eps2.write(res[mode2_id])
eps2.close()

print " done."

