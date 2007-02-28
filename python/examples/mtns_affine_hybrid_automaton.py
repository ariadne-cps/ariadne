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

space=PolyhedralSet(Rectangle("[-2.0,2.0]x[-2.0,2.0]"))
identity=AffineMap(Matrix("[1,0;0,1]"),Vector("[0,0]"))

invariant=PolyhedralSet(Rectangle("[-1.0,1.0]x[-1.0,1.0]"))
dynamic=AffineVectorField(Matrix("[-0.5,-1.0; 1.0,-0.5]"),Vector("[0.00,0.00]"))
activation12 = PolyhedralSet(Rectangle("[-0.2,-0.1]x[-0.2,-0.1]"))
activation21 = PolyhedralSet(Rectangle("[0.1,0.2]x[0.1,0.2]"))
reset = AffineMap(Matrix("[-7,0;0,-7]"),Vector("[0,0]"))

automaton=HybridAutomaton("Affine hybrid automaton")
mode1_id=0
mode2_id=1
mode1=automaton.new_mode(mode1_id,dynamic,invariant)
mode2=automaton.new_mode(mode2_id,dynamic,invariant)
event_id=0
transition=automaton.new_transition(event_id,mode1_id,mode2_id,reset,activation12)
transition=automaton.new_transition(event_id,mode2_id,mode1_id,reset,activation21)
print automaton
#print automaton.invariant()

initial_rectangle=Rectangle("[-0.8,-0.7]x[0.7,0.8]")
bounding_box=Rectangle("[-2,2]x[-2,2]")
epsbb=Rectangle("[-2.1,2.1]x[-2.1,2.1]") # eps bounding box


#grid=RegularGrid(Vector("[0.03125,0.03125]"))
#block=LatticeBlock("[-32,32]x[-32,32]")
grid=RegularGrid(Vector("[0.0625,0.0625]"))
block=LatticeBlock("[-32,32]x[-32,32]")
fgrid=FiniteGrid(grid,block)
print fgrid

print "Creating initial hybrid set"
initial_set=HybridGridMaskSet()
initial_set.new_location(mode1_id,fgrid)
initial_set.new_location(mode2_id,fgrid)
initial_set[mode1_id].adjoin_over_approximation(initial_rectangle)
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


maximum_step_size=0.125;
lock_to_grid_time=0.25;
maximum_set_radius=0.5;

apply=Applicator()
integrator=AffineIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius);
hybrid_evolver=HybridEvolver(apply,integrator);

print "Computing continuous chainreach set"
continuous_chainreach_set=hybrid_evolver.continuous_chainreach(automaton,initial_set,bounding_set)
print "Exporting to postscript output...",
eps=EpsPlot("mtns_affine_hybrid_automaton-1.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(continuous_chainreach_set[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.close()
print

print "Computing chainreach set"
chainreach_set=hybrid_evolver.chainreach(automaton,initial_set,bounding_set)

print "Exporting to postscript output...",
eps=EpsPlot("mtns_affine_hybrid_automaton-2.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("cyan")
eps.write(activation21)
eps.set_fill_colour("magenta")
eps.write(activation12)
eps.set_fill_colour("red")
eps.write(chainreach_set[mode2_id])
eps.set_fill_colour("yellow")
eps.write(chainreach_set[mode1_id])
eps.set_fill_colour("green")
eps.write(continuous_chainreach_set[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.close()
print " done."
