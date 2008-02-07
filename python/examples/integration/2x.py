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
initial_set=RectangularSet([[0.0,0.0],[0.0,0.0]])

# Evolution parameters
par=EvolutionParameters()
par.set_maximum_step_size(0.005);
par.set_lock_to_grid_time(2);
par.set_grid_length(0.01);

integrator=AffineIntegrator()
evolver=VectorFieldEvolver(par,integrator)

print "Computing reach set with Affine Integrator..."
chainreach_set=evolver.upper_reach(dyn,initial_set,Rational(1))

print "Exporting to postscript output...",
bounding_box=Box([[-0.01,1.01],[-0.01,1.01]])
eps=EpsPlot()
eps.open("2x-affine.eps",bounding_box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."

integrator=LohnerIntegrator()
evolver=VectorFieldEvolver(par,integrator)

print "Computing chainreach sets with Lohner Integrator..."
chainreach_set=evolver.chainreach(dyn,initial_set)

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

integrator=EulerIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius)

print "Computing chainreach sets with Euler Integrator..."
chainreach_set=integrator.chainreach(dyn,initial_set,bounding_set)
print "chainreach_set.size(),capacity()=",chainreach_set.size(),chainreach_set.capacity()

print "Exporting to postscript output...",
eps.open("2x-euler.eps",bounding_box,0,1)

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



  	

