#!/usr/bin/python
#
# integration example
#
# Written by Davide Bresolin on May, 2008
#

from ariadne import *
from math import *


# y' = -y
dyn=AffineVectorField(Matrix([[0,0],[0,-1]]),Vector([1,0]))

grid=Grid(Vector([0.01,0.01]))
box=Box([[-0.02,1.02],[-0.02,1.02]])
fgrid=FiniteGrid(grid,box)

print "Creating reference set"
reference_set=GridMaskSet(fgrid)
for i in range(101):
	x = 0.01*i+0.005
	y = exp(-x)
	reference_set.adjoin_outer_approximation(RectangularSet([[x,x],[y,y]]))
	
print "reference_set.size(),capacity()=",reference_set.size(),reference_set.capacity()

print "Creating inital set"
initial_set=RectangularSet([[0.0,0.0],[1.0,1.0]])

time = 1.0
# Evolution parameters
par=EvolutionParameters()
par.set_maximum_step_size(0.25);
par.set_lock_to_grid_time(time);
par.set_grid_length(0.01);


integrator=AffineIntegrator()
evolver=VectorFieldEvolver(par,integrator)

print "Computing upper reach set with Affine Integrator..."
upper_reach_set=evolver.upper_reach(dyn,initial_set,Rational(time))

print "Exporting to postscript output...",
eps=EpsPlot()
eps.open("e-y-affine.eps",box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(upper_reach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done.\n"

integrator=LohnerIntegrator()
evolver=VectorFieldEvolver(par,integrator)

print "Computing upper_reach sets with Lohner Integrator..."
upper_reach_set=evolver.upper_reach(dyn,initial_set,Rational(time))

print "Exporting to postscript output...",
eps.open("e-y-lohner.eps",box,0,1)
eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(upper_reach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."

integrator=KuhnIntegrator(3,3)
evolver=VectorFieldEvolver(par,integrator)

print "Computing upper_reach sets with Kuhn Integrator..."
upper_reach_set=evolver.upper_reach(dyn,initial_set,Rational(time))

eps.open("e-y-kuhn.eps",box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(upper_reach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."







  	

