#!/usr/bin/python

##############################################################################
#            affine_hybrid_automaton.py
#
#  Copyright 2005-7  Alberto Casagrande, Pieter Collins
#
# Maintainer: Pieter Collins <Pieter.Collins@cwi.nl>
#
##############################################################################

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
import sys 

print dir()


space=RectangularSet([[-7.5,7.5],[-7.5,7.5]])
identity=AffineMap(Matrix([[1,0],[0,1]]),Vector([0,0]))

invariant1=PolyhedralSet(Matrix([[-1,0],[0,-1]]),Vector([0,0]))
#invariant1=RectangularSet([[0,2],[0,2]])
dynamic1=AffineVectorField(Matrix([[0,0],[0,0]]),Vector([1,1]))

automaton=HybridAutomaton("Affine hybrid automaton")
mode1_id=DiscreteState(1)
mode1=automaton.new_mode(mode1_id,dynamic1,invariant1)
print automaton

initial_rectangle1=RectangularSet([[0.0001,0.001],[0.0001,0.001]]);
bounding_box=Box([[-8,8],[-8,8]])


print "Creating initial hybrid set"
initial_set=HybridSet()
initial_set.new_location(mode1_id,initial_rectangle1)
print "initial_set.locations() =",initial_set.locations()

print "Creating hybrid grid"
grid=Grid(Vector([0.25,0.25]))
hgrid=HybridGrid()
hgrid.new_location(mode1_id,grid)

print "Creating hybrid bounding domain"
bounding_set=HybridSet()
bounding_set.new_location(mode1_id,space)
print "bounding_set.locations() =",bounding_set.locations()


parameters=EvolutionParameters()
#parameters.set_grid_length(0.125)
parameters.set_grid_length(0.05)
parameters.set_hybrid_grid(hgrid)
parameters.set_lock_to_grid_time(4.0);
parameters.set_maximum_step_size(0.125)
#parameters.set_maximum_enclosure_radius(0.25);
parameters.set_maximum_enclosure_radius(2.5);
parameters.set_verbosity(2);
parameters.set_bounding_domain_size(0.5);
parameters.set_hybrid_bounding_domain(bounding_set);
print parameters

applicator=KuhnApplicator(3)
integrator=AffineIntegrator();
hybrid_evolver=SetBasedHybridEvolver(parameters,applicator,integrator);

time=1

print "Computing lower reach set"
print initial_set
lower_reach_set=hybrid_evolver.upper_reach(automaton,initial_set,time)

print "Exporting to postscript output...",
epsbb=RectangularSet([[-8.1,8.1],[-8.1,8.1]]) # eps bounding box
eps=EpsPlot()
eps.open("unbounded_automaton-upper_reach.eps",epsbb)
eps.set_fill_colour("yellow")
eps.write(invariant1)
eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(lower_reach_set[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.close()
print

#sys.exit(0)

print "Computing chainreach set"
chainreach_set=hybrid_evolver.chain_reach(automaton,initial_set)
#chainreach_set=initial_set

print "Exporting to postscript output...",
eps.open("unbounded_hybrid_automaton-chain_reach.eps",epsbb)
eps.set_fill_colour("yellow")
eps.write(invariant1)
eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(chainreach_set[mode1_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode1_id])
eps.close()
print " done."
