#!/usr/bin/python

from ariadne import *
from evaluation import *
from geometry import *
from linear_algebra import *
import sys

h=HenonMap(Dyadic(1.5),Dyadic(0.275))

x=Point(2)
#print x,h(x),h(h(x)),"\n\n"

r=Rectangle("[1.4,1.6]x[-0.9,1.1]")
#print r,h(r),h(h(r)),"\n\n"
print h(r)

p=Parallelopiped(r)
sp=p.subdivide()
hp=h(p)
hsp=apply(h,sp)
#eps=EpsPlot("henon.eps",Rectangle("[-2,3]x[-2,3]"))
eps=EpsPlot("henon.eps",h(r))
eps.set_fill_colour("blue")
eps.write(hp)
eps.set_fill_colour("green")
eps.write(hsp)

print subset(hsp[0],hp)
print subset(hsp[1],hp)
print subset(hsp[2],hp)
print subset(hsp[3],hp)
print subset(hsp[3].subdivide()[0],hp)
eps.set_fill_colour("red")
eps.write(hsp[3].subdivide()[0])

eps.close()

rb=Rectangle("[-1,1]x[-1,1]")
print rb

g=FiniteGrid(rb,12);
print g
print over_approximation(g,p)

c=Point(2)
A=Matrix(2,2)
c[0]=1.0
c[1]=1.0
A[0,0]=2.0
A[1,1]=1.0
A[0,1]=1.1
A[1,0]=1.5
p1=Parallelopiped(c,A)
p2=Parallelopiped(c,inverse(A))
pls=ParallelopipedListSet(2)
pls.push_back(p1)
pls.push_back(p2)
