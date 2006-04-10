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
from ariadne.evaluation import *
from ariadne.geometry import *
from ariadne.linear_algebra import *
import sys

a=IndexArray(2)
a[0]=1
a[1]=2

lcls=LatticeCellListSet(2)
lc=LatticeCell(a)
lr=LatticeRectangle("[0,2]x[0,3]")
lms=LatticeMaskSet(lr)
lms.adjoin(lc)
print lms
lms.adjoin(lr)
print lms
lcls.adjoin(lms)
print lcls


def plot(fn,bb,set):
  eps=EpsPlot(fn+".eps",bb)
  eps.set_fill_colour("green")
  for bs in set:
    eps.write(bs)
  eps.set_fill_colour("blue")
 # eps.write(set[0])
  eps.set_fill_colour("red")
#  eps.write(set[-1])
  eps.close()

bb=Rectangle("[-2,2]x[-2,2]")

A=Matrix("[-0.25,-1;+1,-0.25]")
b=Vector("[0,0]")
avf=AffineVectorField(A,b)

init_rect=Rectangle("[0.96,1.04]x[0.46,0.54]")
init_paral=Parallelotope(init_rect)

init_set=ParallelotopeListSet(init_paral)
for i in range(0,16):
  print
  init_set.push_back(integrate(avf,init_set[len(init_set)-1],Real(0.1),Real(0.1)))
plot("integrate1",bb,init_set)

print "\n\n\n\n"

init_set=ParallelotopeListSet(init_paral.subdivide())
reach_set=ParallelotopeListSet(init_set)
for i in range(0,11):
  final_set=integrate(avf,init_set,Real(0.2*i),Real(0.1))
  reach_set.adjoin(final_set)
plot("integrate2",bb,reach_set)

print "\n\n\n\n"

grid=InfiniteGrid(2,Real(0.02))
print "grid =", grid
bounding_box=Rectangle("[-2,2]x[-2,2]")
bounds=over_approximation(bounding_box,grid).position()
gcl=over_approximation(init_paral,grid)
print "Cell list =",gcl
init_set=GridMaskSet(grid,bounds)
init_set.adjoin(gcl)
reach_set=GridMaskSet(grid,bounds)
reach_set.adjoin(init_set)
#print "init_set =", init_set
bounds_set=GridMaskSet(grid,bounds)
bounds_set.adjoin(GridRectangle(grid,bounds))
#print "bounding set =",bounds_set
for i in range(0,2):
  init_set=integrate(avf,init_set,bounds_set,Real(2.0),Real(0.2))
  reach_set.adjoin(init_set)
plot("integrate3",bb,reach_set)
#print reach_set

print "Done"
sys.exit()


h=Real(1./32)
ls=LorenzSystem(Real(8./3.),Real(28.0),Real(10.0))
r0=Rectangle("[1.0,1.1]x[1.0,1.1]x[1.0,1.1]")
r=[r0]
for i in range(0,4):
  r.append(integration_step(ls,r[-1],h))
plot("integrate4",bb,r)


print "\n\n\n\n"
p0=Parallelotope(r0)
p=[p0]
for i in range(0,1):
  p.append(integration_step(ls,p[-1],h))
  
print p[-1]
plot("integrate5",bb,p)

print "\n\n\n\n"
print "Computing integration of list set\n\n"
initial=ParallelotopeListSet(p0)
final=integrate(ls,initial,Real(0.3),Real(0.001))
print final
eps=EpsPlot("integrate6.eps",bb)
eps.set_fill_colour("blue")
eps.write(initial)
eps.set_fill_colour("green")
eps.write(final)
eps.close

sys.exit()

print "\n\n\n\n"
print "Computing integration of denotable set\n\n"
grid=InfiniteGrid(3,Real(0.05));
lr=LatticeRectangle("[-400,400]x[-400,400]x[-400,400]")
initial=GridMaskSet(grid,lr)
initial.adjoin(over_approximation(r0,grid))
bounding_set=GridMaskSet(grid,lr)
bounding_set.adjoin(GridRectangle(grid,lr))

final=integrate(ls,initial,bounding_set,Real(1.0),Real(0.1))
eps=EpsPlot("integrate7.eps",bb)
eps.set_fill_colour("green")
eps.write(final)
eps.close
