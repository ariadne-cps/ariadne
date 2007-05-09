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

print dir()

space=PolyhedralSet(Rectangle("[-7.5,7.5]x[-7.5,7.5]"))
identity=AffineMap(Matrix("[1,0;0,1]"),Vector("[0,0]"))

invariant1=PolyhedralSet(Rectangle("[-7.5,7.5]x[-7.5,2.5]"))
invariant2=PolyhedralSet(Rectangle("[-7.5,7.5]x[-7.5,7.5]"))
dynamic1=AffineVectorField(Matrix("[-2,-1; 1,-2]"),Vector("[-1,0]"))
dynamic2=AffineVectorField(Matrix("[-2,-1; 1,-2]"),Vector("[1,0]"))
activation11 = PolyhedralSet(Rectangle("[-7.5,7.5]x[-3,-2]"))
activation21 = PolyhedralSet(Rectangle("[-7.5,7.5]x[-7,-6]"))
activation12 = PolyhedralSet(Rectangle("[-7.5,7.5]x[1,2]"))
reset11 = AffineMap(Matrix("[1,0;0,-1]"),Vector("[0,0]"))
reset21 = AffineMap(Matrix("[1,0;0,-1]"),Vector("[0,0]"))
reset12 = AffineMap(Matrix("[1,0;0,-1]"),Vector("[0,0]"))

automaton=HybridAutomaton("Affine hybrid automaton")
mode1_id=2
mode2_id=3
mode1=automaton.new_mode(mode1_id,dynamic1,invariant1)
mode2=automaton.new_mode(mode2_id,dynamic2,invariant2)
event1_id=5
event2_id=7
transition=automaton.new_transition(event1_id,mode1_id,mode1_id,reset11,activation11)
transition=automaton.new_transition(event2_id,mode1_id,mode2_id,reset12,activation12)
transition=automaton.new_transition(event1_id,mode2_id,mode1_id,reset21,activation21)
print automaton
#print automaton.invariant()

initial_rectangle1=Rectangle("[-6.96875,-6.9375]x[-6.96875,-6.9375]");
initial_rectangle2=Rectangle("[6.9375,6.96875]x[6.9375,6.96875]")
bounding_box=Rectangle("[-8,8]x[-8,8]")


v = Vector(2);
v[0] = 1.0/32;
v[1] = 1.0/32;
grid=Grid(v)
block=LatticeBlock("[-300,300]x[-300,300]")
fgrid=FiniteGrid(grid,block)
print fgrid

print "Creating initial hybrid set"
initial_set=HybridGridMaskSet()
initial_set.new_location(mode1_id,fgrid)
initial_set.new_location(mode2_id,fgrid)
initial_set[mode1_id].adjoin_over_approximation(initial_rectangle1)
initial_set[mode2_id].adjoin_over_approximation(initial_rectangle2)
print "initial_set.locations() =",initial_set.locations()
print "initial_set[mode1_id].size(),capacity()=",initial_set[mode1_id].size(),initial_set[mode1_id].capacity()
print "initial_set[mode2_id].size(),capacity()=",initial_set[mode2_id].size(),initial_set[mode2_id].capacity()

print "Creating bounding hybrid set",
bounding_set=HybridGridMaskSet();
bounding_set.new_location(mode1_id,fgrid)
bounding_set.new_location(mode2_id,fgrid)
bounding_set[mode1_id].adjoin_over_approximation(bounding_box);
bounding_set[mode2_id].adjoin_over_approximation(bounding_box);
print "bounding_set.locations() =",bounding_set.locations()
print "bounding_set[mode1_id].size(),capacity()=",bounding_set[mode1_id].size(),bounding_set[mode1_id].capacity()
print "bounding_set[mode2_id].size(),capacity()=",bounding_set[mode2_id].size(),bounding_set[mode2_id].capacity()
print

maximum_step_size=1.0/32;
lock_to_grid_time=0.25;
maximum_set_radius=0.5;

apply=Applicator()
integrator=AffineIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
#integrator=LohnerIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
hybrid_evolver=HybridEvolver(apply,integrator);

print "Computing continuous chainreach set"
continuous_chainreach_set=hybrid_evolver.continuous_chainreach(automaton,initial_set,bounding_set)
print "Exporting to postscript output...",
epsbb=Rectangle("[-8.1,8.1]x[-8.1,8.1]") # eps bounding box
eps=EpsPlot()
eps.open("affine_hybrid_automaton-1.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("red")
eps.write(continuous_chainreach_set[mode2_id])
eps.set_fill_colour("green")
eps.write(continuous_chainreach_set[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.write(initial_set[mode2_id])
eps.close()
print

print "Computing single discrete step"
discrete_step_set=hybrid_evolver.discrete_step(automaton,initial_set)

print "Exporting to postscript output...",
eps.open("affine_hybrid_automaton-2.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("red")
eps.write(discrete_step_set[mode2_id])
eps.set_fill_colour("green")
eps.write(discrete_step_set[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.write(initial_set[mode2_id])
eps.close()
print

print "Computing chainreach set"
chainreach_set=hybrid_evolver.chainreach(automaton,initial_set,bounding_set)

print "Exporting to postscript output...",
eps.open("affine_hybrid_automaton-3.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("cyan")
eps.write(activation11)
eps.write(activation21)
eps.set_fill_colour("magenta")
eps.write(activation12)
eps.set_fill_colour("red")
eps.write(chainreach_set[mode2_id])
eps.set_fill_colour("yellow")
eps.write(chainreach_set[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.write(initial_set[mode2_id])
eps.close()
print " done."
