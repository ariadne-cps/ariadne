#!/bin/python
#
# Integration example
#
# Written by Davide Bresolin on May, 2007
#

from ariadne import *
from math import *

vector=extract_vector
matrix=extract_matrix

# y' = 2x
dyn=AffineVectorField(matrix([[0,0],[2,0]]),vector([1,0]))

grid=Grid(Vector("[0.01,0.01]"))
block=LatticeBlock("[-2,102]x[-2,102]")
fgrid=FiniteGrid(grid,block)
print fgrid

print "Creating reference set"
reference_set=GridMaskSet(fgrid)
for i in range(101):
	x = 0.01*i+0.005
	y = x*x
	pstr = "["+str(x)+","+str(x)+"]x["+str(y)+","+str(y)+"]"
	# print pstr
	point = Rectangle(pstr)
	reference_set.adjoin_outer_approximation(point)
	
print "reference_set.size(),capacity()=",reference_set.size(),reference_set.capacity()

print "Creating inital set"
initial_rectangle=Rectangle("[0.0,0.0]x[0.0,0.0]")
initial_set=GridMaskSet(fgrid)
print "Adjoining initial rectangle"
initial_set.adjoin_outer_approximation(initial_rectangle)
print "initial_set.size(),capacity()=",initial_set.size(),initial_set.capacity()

bounding_box=Rectangle("[-0.01,1.01]x[-0.01,1.01]")

print "Creating bounding set"
bounding_set=GridMaskSet(fgrid)
print "Adjoining bounding rectangle"
bounding_set.adjoin_over_approximation(bounding_box)
print "bounding_set.size(),capacity()=",bounding_set.size(),bounding_set.capacity()

maximum_step_size=0.005;
lock_to_grid_time=2;
maximum_set_radius=0.1;

integrator=AffineIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius)

print "Computing chainreach sets with Affine Integrator..."
chainreach_set=integrator.chainreach(dyn,initial_set,bounding_set)
print "chainreach_set.size(),capacity()=",chainreach_set.size(),chainreach_set.capacity()

print "Exporting to postscript output...",
eps=EpsPlot()
eps.open("2x-affine.eps",bounding_box,0,1)

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

integrator=LohnerIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius)

print "Computing chainreach sets with Lohner Integrator..."
chainreach_set=integrator.chainreach(dyn,initial_set,bounding_set)
print "chainreach_set.size(),capacity()=",chainreach_set.size(),chainreach_set.capacity()

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



  	

