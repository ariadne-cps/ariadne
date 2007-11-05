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


mu=Float(0.25)

vdp=VanDerPolEquation(mu)

parameters=EvolutionParameters()
parameters.set_maximum_step_size(Rational(0.125))
parameters.set_lock_to_grid_time(Rational(1.0))
parameters.set_maximum_basic_set_radius(Float(0.0625))

integrator=C1LohnerIntegrator()
evolver=VectorFieldEvolver(parameters,integrator)

print "Available integration methods:",dir(evolver),"\n"

subdivisions=128
initial_set=RectangularSet([[0.99,1.01],[0.49,0.51]])

print "initial_set =", initial_set,"\n"

step_size=Rational(0.125)
flow_steps=2
reach_steps=3

flow_time=flow_steps*step_size
reach_time=reach_steps*step_size

print "Integrate initial set for time",flow_time
flowed_set=evolver.integrate(vdp,initial_set,flow_time)
print flowed_set
print "  ",flowed_set,"\n"

print "Computing reachable set from time",flow_time,"to time",flow_time,"+",reach_time
reach_set=evolver.reach(vdp,flowed_set,reach_time)
print reach_set,"\n\n"


intermediate_sets=ZonotopeListSet(initial_set)
current_set=initial_set
for i in range(0,flow_steps+reach_steps):
  current_set=evolver.integrate(vdp,current_set,step_size)
  intermediate_sets.adjoin(ZonotopeListSet(current_set))

print "initial set:",initial_set
print "flowed set:",flowed_set
print
print "reach set:",reach_set

print "Exporting to postscript output...",
epsbb=Rectangle("[-4.1,4.1]x[-2.1,2.1]") # eps bounding box
epsbb=Rectangle("[0.4,1.4]x[-0.4,0.6]") # eps bounding box
eps=EpsPlot()
eps.open("van_der_pol_oscillator-1.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(reach_set)
eps.set_fill_colour("white")
eps.write(intermediate_sets)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()
print " done."
