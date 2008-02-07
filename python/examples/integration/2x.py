#!/usr/bin/python
#
# Integration example
#
# Written by Davide Bresolin on May, 2007
#

from ariadne import *
from math import *

# y' = 2x
dyn=AffineVectorField(Matrix([[0,0],[2,0]]),Vector([1,0]))

grid=Grid(Vector([0.01,0.01]))
block=Box([[-0.01,1.01],[-0.01,1.01]])
fgrid=FiniteGrid(grid,block)
print fgrid

print "Creating reference set"
reference_set=GridMaskSet(fgrid)
for i in range(101):
	x = 0.01*i+0.005
	y = x*x
	pstr = "["+str(x)+","+str(x)+"]x["+str(y)+","+str(y)+"]"
	# print pstr
	point = Box([[x,x],[y,y]])
	reference_set.adjoin_outer_approximation(point)
	
print "reference_set.size(),capacity()=",reference_set.size(),reference_set.capacity()

print "Creating inital set"
initial_set=RectangularSet([[0.0,0.01],[0.0,0.01]])

# Evolution parameters
par=EvolutionParameters()
par.set_maximum_step_size(0.005);
par.set_lock_to_grid_time(2);
par.set_grid_length(0.01);

integrator=AffineIntegrator()
evolver=VectorFieldEvolver(par,integrator)

print "Computing upper reach set with Affine Integrator..."
chainreach_set=evolver.upper_reach(dyn,initial_set,Rational(1))

print "Computing lower reach set with Affine Integrator..."
lower_reach_set=evolver.lower_reach(dyn,initial_set,Rational(1))

print "Exporting to postscript output...",
bounding_box=Box([[-0.01,1.01],[-0.01,1.01]])
eps=EpsPlot()
eps.open("2x-affine-upper.eps",bounding_box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

eps.open("2x-affine-lower.eps",bounding_box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("cyan")
eps.write(lower_reach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."

integrator=LohnerIntegrator()
evolver=VectorFieldEvolver(par,integrator)

print "Computing chainreach sets with Lohner Integrator..."
chainreach_set=evolver.upper_reach(dyn,initial_set,Rational(1))

print "Exporting to postscript output...",
eps.open("2x-lohner.eps",bounding_box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(bounding_box)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."

integrator=KuhnIntegrator(3,3)
evolver=VectorFieldEvolver(par,integrator)

print "Computing chainreach sets with Kuhn Integrator..."
chainreach_set=evolver.upper_reach(dyn,initial_set,Rational(1))

eps.open("2x-kuhn.eps",bounding_box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(bounding_box)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."





  	

