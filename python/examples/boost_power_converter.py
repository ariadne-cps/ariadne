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
vector=extract_vector
matrix=extract_matrix

print dir()

E=1.0
C=1.0
L=1.0
R=1.0
T=1.0
S=0.75

# variables I, Q, t

# mode1; switch open, I>=0, t<=T*S
invariant1=PolyhedralSet(matrix([[-1,0,0],[0,0,1]]),vector([0,T*S]))
dynamic1=AffineVectorField(matrix([[0,-1/(C*L),0],[1,-1/(R*C),0],[0,0,0]]),vector([E/L,0,1]))

# mode2; switch open, I=0 Q>=EC t<=T*S
invariant2=PolyhedralSet(matrix([[0,-1,0],[0,0,1]]),vector([-E*C,T*S]))
dynamic2=AffineVectorField(matrix([[0,0,0],[0,-1/(R*C),0],[0,0,0]]),vector([0,0,1]))

# mode3; switch closed, Q>=0 t<=1
invariant3=PolyhedralSet(matrix([[0,0,1]]),vector([1]))
dynamic3=AffineVectorField(matrix([[0,0,0],[0,-1/(R*C),0],[0,0,0]]),vector([E/L,0,1]))

identity=AffineMap(Matrix("[1,0,0;0,1,0;0,0,1]"),Vector("[0,0,0]"))
clock_reset=AffineMap(Matrix("[1,0,0;0,1,0;0,0,0]"),Vector("[0,0,0]"))

# switch if I<=0
activation12 = PolyhedralSet(matrix([[0,0,0]]),vector([0]))
reset12 = identity
# switch if t>=T*S
activation13 = PolyhedralSet(matrix([[0,0,-1]]),vector([-T*S]))
activation23 = activation13
reset13 = identity
reset23 = identity
# switch if t>=T
activation31 = PolyhedralSet(matrix([[0,0,-1]]),vector([-T]))
reset31 = clock_reset

automaton=HybridAutomaton("Affine hybrid automaton")
mode1_id=1
mode2_id=2
mode3_id=3
mode1=automaton.new_mode(mode1_id,dynamic1,invariant1)
mode2=automaton.new_mode(mode2_id,dynamic2,invariant2)
mode3=automaton.new_mode(mode3_id,dynamic3,invariant3)
event1_id=1
event2_id=2
event3_id=3
automaton.new_transition(event1_id,mode1_id,mode2_id,reset12,activation12)
automaton.new_transition(event2_id,mode1_id,mode3_id,reset13,activation13)
automaton.new_transition(event2_id,mode2_id,mode3_id,reset23,activation23)
automaton.new_transition(event3_id,mode3_id,mode1_id,reset31,activation31)
print automaton

initial_point=Point("(0.96875,0.28125,0.0)")
initial_rectangle=Rectangle(initial_point,initial_point)
bounding_box=Rectangle("[0,2]x[0,2]x[0,"+str(T)+"]")

print initial_rectangle
print bounding_box

grid=Grid(Vector("[0.0625,0.0625,0.03125]"))
block=LatticeBlock("[-1,33]x[-1,33]x[-1,33]")
fgrid=FiniteGrid(grid,block)
print fgrid

print "Creating hybrid set locations"
initial_set=HybridGridMaskSet()
initial_set.new_location(mode1_id,fgrid)
initial_set.new_location(mode2_id,fgrid)
initial_set.new_location(mode3_id,fgrid)
print initial_set.locations()
print "Adjoining initial rectangle"
initial_set[mode1_id].adjoin_over_approximation(initial_rectangle)
print "initial_set.locations() =",initial_set.locations()
print "initial_set[mode1_id].size(),capacity()=",initial_set[mode1_id].size(),initial_set[mode1_id].capacity()
print "initial_set[mode2_id].size(),capacity()=",initial_set[mode2_id].size(),initial_set[mode2_id].capacity()

print "Creating initial hybrid cell list set"
initial_cell_list_set=HybridGridCellListSet(initial_set)
#initial_cell_list_set.new_location(mode1_id,grid)
#initial_cell_list_set.new_location(mode2_id,grid)
#initial_cell_list_set[mode2_id].adjoin_over_approximation(initial_rectangle2)
print "initial_cell_list_set.locations() =",initial_cell_list_set.locations()

print "Creating bounding hybrid set",
bounding_set=HybridGridMaskSet();
bounding_set.new_location(mode1_id,fgrid)
bounding_set.new_location(mode2_id,fgrid)
bounding_set.new_location(mode3_id,fgrid)
bounding_set[mode1_id].adjoin_over_approximation(bounding_box);
bounding_set[mode2_id].adjoin_over_approximation(bounding_box);
bounding_set[mode3_id].adjoin_over_approximation(bounding_box);
print "bounding_set.locations() =",bounding_set.locations()
print "bounding_set[mode1_id].size(),capacity()=",bounding_set[mode1_id].size(),bounding_set[mode1_id].capacity()
print "bounding_set[mode2_id].size(),capacity()=",bounding_set[mode2_id].size(),bounding_set[mode2_id].capacity()
print


maximum_step_size=0.25;
lock_to_grid_time=0.25;
maximum_set_radius=0.25;

apply=Applicator()
integrator=AffineIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
hybrid_evolver=HybridEvolver(apply,integrator);

print "Computing continuous chainreach set"
continuous_chainreach_set=hybrid_evolver.continuous_chainreach(automaton,initial_set,bounding_set)
print "Exporting to postscript output...",
epsbb=Rectangle("[-0.1,2.1]x[-0.1,2.1]x[-0.1,1.1]") # eps bounding box
epspmap0=PlanarProjectionMap(3,2,0) # eps projection map
epspmap1=PlanarProjectionMap(3,2,1) # eps projection map
eps=EpsPlot()
eps.open("boost_power_converter-1.eps",epsbb,epspmap0)
eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(bounding_box)
eps.set_fill_colour("green")
eps.write(continuous_chainreach_set[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.close()

eps.open("boost_power_converter-2.eps",epsbb,epspmap1)
eps.set_fill_colour("white")
eps.write(bounding_box)
eps.set_fill_colour("green")
eps.write(continuous_chainreach_set[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.close()
print

