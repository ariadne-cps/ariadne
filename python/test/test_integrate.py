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

n=6

ls=LorenzSystem(Dyadic(8./3.),Dyadic(28.0),Dyadic(10.0))

r=[Rectangle("[1.0,1.1]x[1.0,1.1]x[1.0,1.1]")]
for i in range(0,1):
  r.append(integrate(ls,r[-1],Dyadic(0.1)))

b=Vector(2)
A=Matrix(2,2)
A[0,1]=-1;
A[1,0]=+1;
avf=AffineVectorField(A,b)

ar=[Rectangle("[0.95,1.05]x[-0.05,0.05]")]
for i in range(0,n):
  ar.append(integrate(avf,ar[-1],Dyadic(0.2)))

bb=Rectangle("[-2,2]x[-2,2]")
eps=EpsPlot("integrate.eps",bb)
eps.set_fill_colour("green")
for i in range(0,n):
  eps.write(ar[i])
eps.set_fill_colour("blue")
eps.write(ar[0])
eps.close()
