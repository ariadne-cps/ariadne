#!/usr/bin/python

##############################################################################
#            test_polyhedron.py
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

from ariadne.base import *
from ariadne.evaluation import *
from ariadne.geometry import *
from ariadne.linear_algebra import *
import sys

bb=Rectangle("[-3,3]x[-3,3]")
print bb,"\n"

g=FiniteGrid(bb,16);
print g
print

r1=Rectangle("[-0.2,0.3]x[0.1,0.9]")
print over_approximation(r1,g)
r2=Rectangle("[0.1,0.5]x[0.2,0.4]")
print over_approximation(r2,g)

rls=RectangleListSet(2)
rls.push_back(r1)
rls.push_back(r2)
print rls

c=Point(2)
A=Matrix(2,2)
c[0]=0.5
c[1]=0.1
A[0,0]=1
A[0,1]=0.5
A[1,0]=0.5
A[1,1]=0.6

p=Parallelotope(c,A)

print "\n\n", p, "\n\n"
sp=p.subdivide()
print "Subdivided: ",sp
ap=over_approximation(p,g)
print sp
print ap
asp=over_approximation(sp,g)
print asp
afp=GridMaskSet(over_approximation(sp[0],g))
print afp, asp
cafp=difference(asp,afp)

print "Stepping through RectangleListSet using index"
for i in range(0,len(rls)):
  print rls[i]
print "Iterating through RectangleListSet"
for rect in rls:
  print rect

gms=GridMaskSet(asp)
print "Stepping through GridMaskSet using index"
for i in range(0,len(gms)):
  print gms[i]
print "Iterating through GridMaskSet"
i=0
for cell in gms:
  print cell
  i+=1
print i, len(gms)

gcls=GridCellListSet(ap)
print "Stepping through GridCellListSet using index"
for i in range(0,len(gcls)):
  print gcls[i]
print "Iterating through GridCellListSet"
for cell in gcls:
  print cell

grls=GridRectangleListSet(rls)
print "Stepping through GridRectangleListSet using index"
for i in range(0,len(grls)):
  print grls[i]
print "Iterating through GridRectangleListSet"
for rect in grls:
  print rect

  
eps=EpsPlot("gr1.eps",bb)
eps.set_fill_colour("blue")
eps.write(over_approximation(p.bounding_box(),g))
eps.set_fill_colour("red")
eps.write(asp)
eps.set_fill_colour("green")
eps.write(sp)
eps.close()

eps=EpsPlot("gr2.eps",bb)
eps.set_fill_colour("red")
eps.write(afp)
eps.set_fill_colour("green")
eps.write(cafp)
eps.set_fill_colour("white")
eps.write(sp)
eps.close()

eps=EpsPlot("gr3.eps",bb)
eps.set_fill_colour("white")
eps.write(gms.bounding_box())
eps.set_fill_colour("magenta")
eps.write(gms.adjoining().adjoining())
eps.set_fill_colour("red")
eps.write(gms.neighbourhood())
eps.set_fill_colour("blue")
eps.write(gms.adjoining())
eps.set_fill_colour("green")
eps.write(gms)
eps.close()

pts=PartitionTreeSet(gms.neighbourhood())
grls=GridRectangleListSet(pts)
gms=GridMaskSet(grls)
gms.clear()
gms.adjoin(grls[0*(len(grls)-6)])
eps=EpsPlot("gr4.eps",bb)
eps.set_fill_colour("magenta")
eps.write(grls)
eps.set_fill_colour("green")
eps.write(gms)
eps.close()
