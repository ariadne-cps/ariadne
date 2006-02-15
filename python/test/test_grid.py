#!/usr/bin/python

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

r=Rectangle("[-0.2,0.3]x[0.1,0.9]")
print over_approximation(r,g)

c=Point(2)
A=Matrix(2,2)
c[0]=0.5
c[1]=0.1
A[0,0]=1
A[0,1]=0.5
A[1,0]=0.5
A[1,1]=0.6

p=Parallelopiped(c,A)

print "\n\n", p, "\n\n"
sp=p.subdivide()
ap=over_approximation(p,g)
print sp
print ap
asp=over_approximation(sp,g)
print asp

for i in range(0,len(asp)):
  print asp[i]

eps=EpsPlot("tg.eps",bb)
eps.set_fill_colour("blue")
eps.write(over_approximation(p.bounding_box(),g))
eps.set_fill_colour("red")
eps.write(asp)
eps.set_fill_colour("green")
eps.write(sp)

afp=GridMaskSet(over_approximation(sp[0],g))
print afp, asp
cafp=difference(asp,afp)
eps.set_fill_colour("red")
eps.write(afp)
eps.set_fill_colour("green")
eps.write(cafp)
eps.set_fill_colour("white")
eps.write(sp)

eps.close()
