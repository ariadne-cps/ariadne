#!/usr/bin/python
#
# esempio di integrazione
#
# Written by Davide Bresolin on May, 2007
#

from ariadne import *
from math import *

# dinamica x' = -x -y & y' = x - y
dyn=AffineVectorField(Matrix([[-1,-1],[1,-1]]),Vector([0,0]))

grid=Grid(Vector([0.01,0.01]))
box=Box([[-1.02,1.02],[-1.02,1.02]])
fgrid=FiniteGrid(grid,box)


print "Creating inital set"
initial_set=GridMaskSet(fgrid)
print "Adjoining initial rectangles"
initial_set.adjoin_outer_approximation(RectangularSet([[0.0,0.0],[1.0,1.0]]))
initial_set.adjoin_outer_approximation(RectangularSet([[1.0,1.0],[1.0,1.0]]))
initial_set.adjoin_outer_approximation(RectangularSet([[1.0,1.0],[0.0,0.0]]))
print "initial_set.size(),capacity()=",initial_set.size(),initial_set.capacity()


# Evolution parameters
par=EvolutionParameters()
par.set_maximum_step_size(0.25);
par.set_lock_to_grid_time(2);
par.set_grid_length(0.01);


integrator=AffineIntegrator()
evolver=VectorFieldEvolver(par,integrator)


print "Computing chainreach sets with Affine Integrator..."
chainreach_set=evolver.chainreach(dyn,initial_set,box)
print "chainreach_set.size(),capacity()=",chainreach_set.size(),chainreach_set.capacity()

print "Exporting to text output...",
txt=TextFile()
txt.open("spirals-affine.txt")
txt.write(chainreach_set)
txt.close()
print "Done."
  
print "Exporting to postscript output...",
eps=EpsPlot()
eps.open("spirals-affine.eps",box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(box)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."

integrator=LohnerIntegrator()
evolver=VectorFieldEvolver(par,integrator)

print "Computing chainreach sets with Lohner Integrator..."
chainreach_set=evolver.chainreach(dyn,initial_set,box)
print "chainreach_set.size(),capacity()=",chainreach_set.size(),chainreach_set.capacity()

print "Exporting to postscript output...",
eps.open("spirals-lohner.eps",box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(box)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."

integrator=KuhnIntegrator(3,3)
evolver=VectorFieldEvolver(par,integrator)

print "Computing chainreach sets with Kuhn Integrator..."
chainreach_set=evolver.chainreach(dyn,initial_set,box)
print "chainreach_set.size(),capacity()=",chainreach_set.size(),chainreach_set.capacity()

print "Exporting to postscript output...",
eps.open("spirals-kuhn.eps",box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(box)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

print "Done."




  	

