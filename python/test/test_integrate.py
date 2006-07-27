#!/usr/bin/python

##############################################################################
#            test_integrate.py
#
#  Copyright 2006  Pieter Collins <Pieter.Collins@cwi.nl>
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

from ariadne import Real
from ariadne.base import *
from ariadne.numeric import *
from ariadne.linear_algebra import *
from ariadne.geometry import *
from ariadne.system import *
from ariadne.evaluation import *
import sys

def plot(fn,bb,set):
  eps=EpsPlot(fn+".eps",bb,0,1)
  eps.set_fill_colour("green")
  for bs in set:
    eps.write(bs)
  eps.set_fill_colour("blue")
  eps.write(set[0])
  eps.set_fill_colour("red")
  eps.write(set[-1])
  eps.close()

grid=RegularGrid(2,Real(0.03125))
bb=Rectangle("[-2,2]x[-2,2]")
bbb=Rectangle("[0.8,1.05]x[0.93,1.05]x[-2,2]")

A=Matrix("[-0.5,-0.75;+1.25,-0.5]")
b=Vector("[0,0]")
avf=AffineVectorField(A,b)

lohner=C1LohnerIntegrator(0.125,0.5,0.0625)
  
print "Testing integration of affine map on parallelotope"
initial_rectangle=Rectangle("[0.96,1.04]x[0.46,0.54]")
initial_parallelotope=Parallelotope(initial_rectangle)
initial_set=ParallelotopeListSet(initial_parallelotope)
lohner.integration_step(avf,initial_set[-1],Real(0.125))
for i in range(0,16):
  initial_set.push_back(lohner.integrate(avf,initial_set[-1],Real(0.125)))
initial_set.push_back(lohner.integrate(avf,initial_set[0],Real(2.0)))
plot("integrate1",bb,initial_set)
print "Done\n"

print "Testing integration of affine map on parallelotope list set"
initial_set=ParallelotopeListSet(initial_parallelotope.subdivide())
final_set=lohner.integrate(avf,initial_set,Real(2.0))
plot("integrate2",bb,[initial_set,final_set])
print "Done\n"

print "Testing reach of affine map on parallelotope list set"
initial_set=ParallelotopeListSet(initial_parallelotope.subdivide())
reach_set=lohner.reach(avf,initial_set,Real(2.0))
plot("integrate3",bb,reach_set)
print "Done\n"

print "Testing integration of affine map on grid mask set"
grid=RegularGrid(2,Real(0.02))
bounding_box=Rectangle("[-2,2]x[-2,2]")
bounds=over_approximation(bounding_box,grid).lattice_set()
gcl=over_approximation(initial_parallelotope,grid)
finitegrid=FiniteGrid(grid,bounds)
#print "Cell list =",gcl
initial_set=GridMaskSet(finitegrid)
initial_set.adjoin(gcl)
reach_set=GridMaskSet(finitegrid)
reach_set.adjoin(initial_set)
#print "init_set =", init_set
bounds_set=GridMaskSet(finitegrid)
bounds_set.adjoin(GridRectangle(grid,bounds))
#print "bounding set =",bounds_set
for i in range(0,2):
  initial_set=lohner.integrate(avf,initial_set,bounds_set,Real(2.0))
  reach_set.adjoin(initial_set)
plot("integrate4",bb,reach_set)
print "Done\n"

print "Testing reach of affine map on grid mask set"
grid=RegularGrid(2,Real(0.0625))
bounding_box=Rectangle("[-2,2]x[-2,2]")
bounds=over_approximation(bounding_box,grid).lattice_set()
finitegrid=FiniteGrid(grid,bounds)
gcl=over_approximation(initial_parallelotope,grid)
#print "Cell list =",gcl
initial_set=GridMaskSet(finitegrid)
initial_set.adjoin(gcl)
reach_set=GridMaskSet(finitegrid)
reach_set.adjoin(initial_set)
#print "init_set =", init_set
bounds_set=GridMaskSet(finitegrid)
bounds_set.adjoin(GridRectangle(grid,bounds))
#print "bounding set =",bounds_set
reach_set=lohner.reach(avf,initial_set,bounds_set,Real(1.5))
final_set=lohner.integrate(avf,initial_set,bounds_set,Real(1.5))
plot("integrate5",bb,[initial_set,reach_set,final_set])
print "Done\n"


lohner=C1LohnerIntegrator(0.05,0.1,0.1)
print "Testing integration of lorenz system on rectangle"
h=Real(1./32)
ls=LorenzSystem(Real(8./3.),Real(28.0),Real(10.0))
r0=Rectangle("[1.0,1.01]x[1.0,1.01]x[1.0,1.01]")
r=[r0]
for i in range(0,4):
  r.append(lohner.integration_step(ls,r[-1],h))
plot("integrate6",bbb,r)
print "Done\n"

print "Testing integration of lorenz system on parallelotope"
p0=Parallelotope(r0)
p=[p0]
for i in range(0,1):
  p.append(lohner.integration_step(ls,p[-1],h))
plot("integrate7",bbb,p)
print "Done\n"

print "Testing integration of lorenz system on list set"
initial=ParallelotopeListSet(p0)
final=lohner.integrate(ls,initial,Real(0.14))
#print final
eps=EpsPlot("integrate8.eps",bbb,0,1)
eps.set_fill_colour("blue")
eps.write(initial)
eps.set_fill_colour("green")
eps.write(final)
eps.close
print ls.derivative(r0.centre())
print "Done\n"

print "Exiting"
sys.exit()

print "Testing integration of lorenz system on grid mask set"
grid=InfiniteGrid(3,Real(0.1));
lr=LatticeRectangle("[-100,100]x[-100,100]x[-100,100]")
initial=GridMaskSet(grid,lr)
initial.adjoin(over_approximation(r0,grid))
bounding_set=GridMaskSet(grid,lr)
bounding_set.adjoin(GridRectangle(grid,lr))

final=lohner.integrate(ls,initial,bounding_set,Real(1.0))
eps=EpsPlot("integrate8.eps",bbb,0,1)
eps.set_fill_colour("green")
eps.write(final)
eps.close
print "Done\n"
