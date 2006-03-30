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

from ariadne.base import *
from ariadne.evaluation import *
from ariadne.geometry import *
from ariadne.linear_algebra import *
import sys

n=20
h=Dyadic(1./256)
ls=LorenzSystem(Dyadic(8./3.),Dyadic(28.0),Dyadic(10.0))
r0=Rectangle("[1.0,1.1]x[1.0,1.1]x[1.0,1.1]")
r=[r0]
for i in range(0,n):
  r.append(integrate(ls,r[-1],h))

b=Vector(2)
A=Matrix(2,2)
A[0,1]=-1;
A[1,0]=+1;
avf=AffineVectorField(A,b)

ar=[Rectangle("[0.95,1.05]x[-0.05,0.05]")]
for i in range(0,n):
  print
  ar.append(integrate(avf,ar[-1],Dyadic(0.2)))

p=[Parallelotope(r0)]
for i in range(0,n):
  p.append(integrate(ls,p[-1],h))

bb=Rectangle("[-4,4]x[-2,6]")
eps=EpsPlot("integrate.eps",bb)
eps.set_fill_colour("green")
for i in range(0,n):
  eps.write(p[i])
eps.set_fill_colour("blue")
eps.write(p[0])
eps.close()
