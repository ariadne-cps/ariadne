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

A=Matrix(2,2)
b=Vector(2)
A[0,0]=-0.25;
A[0,1]=-1;
A[1,0]=+1;
A[1,1]=-0.25;
avf=AffineVectorField(A,b)

init_rect=Rectangle("[0.95,1.05]x[0.45,0.55]")
init_paral=Parallelotope(init_rect)

ar=[init_rect]
for i in range(0,16):
  ar.append(integrate(avf,ar[-1],Real(0.1),Real(0.1)))
plot("integrate1",bb,ar)

print "\n\n\n\n"
print avf.name()

ap=[init_paral]
for i in range(0,128):
  print
  ap.append(integrate(avf,ap[-1],Real(0.1),Real(0.1)))
#for p in ap:
#  print p
plot("integrate2",bb,ap)

sys.exit()

h=Real(1./64)
ls=LorenzSystem(Real(8./3.),Real(28.0),Real(10.0))
r0=Rectangle("[1.0,1.1]x[1.0,1.1]x[1.0,1.1]")
r=[r0]
for i in range(0,n):
  r.append(integration_step(ls,r[-1],h))


print "\n\n\n\n"
p=[Parallelotope(r0)]
#p.append(reach_step(ls,p[-1],h/4))
for i in range(0,n):
  p.append(integration_step(ls,p[-1],h))
#p.append(reach_step(ls,p[-1],Real(h/4)))
  
print p[-1]

bb
eps=EpsPlot("integrate2.eps",bb)
eps.set_fill_colour("green")
for i in range(0,n):
  eps.write(p[i])
eps.set_fill_colour("blue")
eps.write(p[0])
eps.set_fill_colour("red")
eps.write(p[-1])
eps.close()
