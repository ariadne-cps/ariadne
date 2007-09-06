#!/usr/bin/python
#
# esempio di integrazione
#
# Written by Davide Bresolin on May, 2007
#

from ariadne import *
from math import *

vector=extract_vector
matrix=extract_matrix

# dinamica x' = -x -y & y' = x - y
dyn=AffineVectorField(matrix([[-1,-1],[1,-1]]),vector([0,0]))

grid=Grid(Vector("[0.01,0.01]"))
block=LatticeBlock("[-102,102]x[-102,102]")
fgrid=FiniteGrid(grid,block)
print fgrid

print "Creating inital set"
initial_set=GridMaskSet(fgrid)
print "Adjoining initial rectangles"
initial_set.adjoin_outer_approximation(Rectangle("[0.0,0.0]x[1.0,1.0]"))
initial_set.adjoin_outer_approximation(Rectangle("[1.0,1.0]x[1.0,1.0]"))
initial_set.adjoin_outer_approximation(Rectangle("[1.0,1.0]x[0.0,0.0]"))
print "initial_set.size(),capacity()=",initial_set.size(),initial_set.capacity()

bounding_box=Rectangle("[-1.01,1.01]x[-1.01,1.01]")

print "Creating bounding set"
bounding_set=GridMaskSet(fgrid)
print "Adjoining bounding rectangle"
bounding_set.adjoin_over_approximation(bounding_box)
print "bounding_set.size(),capacity()=",bounding_set.size(),bounding_set.capacity()

# crea l'integratore
maximum_step_size=0.005;
lock_to_grid_time=2;
maximum_set_radius=0.1;

integrator=AffineIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius)

print "Computing chainreach sets with Affine Integrator..."
chainreach_set=integrator.chainreach(dyn,initial_set,bounding_set)
print "chainreach_set.size(),capacity()=",chainreach_set.size(),chainreach_set.capacity()

print "Exporting to text output...",
txt=TxtPlot()
txt.open("spirals-affine.txt")
txt.write(chainreach_set)
txt.close()
print "Done."
  
print "Exporting to postscript output...",
eps=EpsPlot()
eps.open("spirals-affine.eps",bounding_box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(bounding_box)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."

integrator=LohnerIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius)

print "Computing chainreach sets with Lohner Integrator..."
#chainreach_set=integrator.chainreach(dyn,initial_set,bounding_set)
print "chainreach_set.size(),capacity()=",chainreach_set.size(),chainreach_set.capacity()

print "Exporting to postscript output...",
eps.open("spirals-lohner.eps",bounding_box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(bounding_box)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."

lock_to_grid_time=0.2; 

integrator=EulerIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius)

print "Computing chainreach sets with Euler Integrator..."
chainreach_set=integrator.chainreach(dyn,initial_set,bounding_set)
print "chainreach_set.size(),capacity()=",chainreach_set.size(),chainreach_set.capacity()

print "Exporting to postscript output...",
eps.open("spirals-euler.eps",bounding_box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(bounding_box)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."




  	

