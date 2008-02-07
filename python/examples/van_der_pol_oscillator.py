#!/usr/bin/python

##############################################################################
#            van_der_pol_oscillator.py
#
#  Copyright 2006  Pieter Collins 
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


mu=0.25

vdp=VanDerPolEquation(Point([mu]))

parameters=EvolutionParameters()
parameters.set_lock_to_grid_time(1.0)
parameters.set_grid_length(0.25)
parameters.set_maximum_step_size(0.125)
parameters.set_maximum_enclosure_radius(2.5)
integrator=LohnerIntegrator()
evolver=VectorFieldEvolver(parameters,integrator)
print "Available integration methods:",dir(evolver),"\n"

initial_set=RectangularSet([[0.99,1.01],[0.49,0.51]])
print "initial_set =", initial_set,"\n"

step_size=Rational(0.125)
flow_steps=20
reach_steps=30

flow_time=flow_steps*step_size
reach_time=reach_steps*step_size

print "Computing lower reachable set to time",reach_time
lower_reach_set=evolver.lower_reach(vdp,initial_set,reach_time)
print lower_reach_set,"\n\n"

print "Computing reachable set to time",reach_time
reach_set=evolver.upper_reach(vdp,initial_set,reach_time)
print reach_set,"\n\n"

print "Integrate initial set for time",flow_time
evolve_set=evolver.upper_evolve(vdp,initial_set,flow_time)
print evolve_set

print "Computing reachable set from time",flow_time,"to time",flow_time,"+",reach_time
reach_set2=evolver.upper_reach(vdp,evolve_set,reach_time)
print reach_set,"\n\n"



print "Exporting to postscript output...",
epsbb=Box([[-4.1,4.1],[-2.1,2.1]]) # eps bounding box
#epsbb=Box([[0.4,1.4],[-0.4,0.6]]) # eps bounding box
eps=EpsPlot()
eps.open("van_der_pol_oscillator-1.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(reach_set2)
eps.set_fill_colour("cyan")
eps.write(reach_set)
eps.set_fill_colour("yellow")
eps.write(evolve_set)
eps.set_fill_colour("magenta")
eps.write(lower_reach_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()
print " done."
