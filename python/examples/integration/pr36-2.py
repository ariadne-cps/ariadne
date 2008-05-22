#!/usr/bin/python
#
# integration example
#
# Written by Davide Bresolin on May, 2007
#

from ariadne import *
from math import *


# y' = z & z' = 3z + 2y
dyn=AffineVectorField(Matrix([[0,0,0],[0,0,1],[0,2,3]]),Vector([1,0,0]))

grid=Grid(Vector([0.1,0.1,0.1]))
box=Box([[-0.2,1.2],[-1.2,2.2],[-0.2,8.2]])
fgrid=FiniteGrid(grid,box)

print "Creating reference set"
reference_set=GridMaskSet(fgrid)
for i in range(11):
	x = 0.1*i+0.001
	y = exp(2*x)-2*exp(x)
	reference_set.adjoin_outer_approximation(RectangularSet([[x,x],[y,y],[0,0]]))
	


print "Creating inital set"
initial_set=RectangularSet([[0.0,0.0],[-1.0,-1.0],[0.0,0.0]])

time = 1.0
# Evolution parameters
par=EvolutionParameters()
par.set_maximum_step_size(0.25);
par.set_lock_to_grid_time(time);
par.set_grid_length(0.1);


integrator=AffineIntegrator()
evolver=VectorFieldEvolver(par,integrator)

print "Computing upper reach set with Affine Integrator..."
upper_reach_set=evolver.upper_reach(dyn,initial_set,Rational(time))

print "Exporting to postscript output...",
eps=EpsPlot()
eps.open("pr36-2-affine.eps",box,0,1)

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
eps.open("pr36-2-lohner.eps",box,0,1)
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

eps.open("pr36-2-kuhn.eps",box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(upper_reach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."



  	

