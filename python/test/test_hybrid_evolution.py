#!/usr/bin/python

##############################################################################
#            test_hybrid_evolution.py
#
#  Copyright 2006  Alberto Casagrande, Pieter Collins 
#  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

from ariadne.hybrid_system import *
from ariadne.hybrid_system_evolution import *

import sys

eps=EpsPlot("hybrid_evolution.eps",Rectangle("[-1,1]x[-1,1]"),0,1,"Provadatw","b")
eps.set_fill_colour("yellow")

b=Vector(2)
A=Matrix(2,2)

A[0,0]=-0.25
A[0,1]=-1
A[1,0]=+1
A[1,1]=-0.25
dyn=AffineVectorField(A,b)
inv=Zonotope(Rectangle("[-1,1]x[-1,1]"))
n1=DiscreteNode(dyn,inv)

n2=DiscreteNode(dyn,inv)

A[0,0]=-7
A[0,1]=0
A[1,0]=0
A[1,1]=-7
res=AffineMap(A,b)
act=Zonotope(Rectangle("[-0.2,0]x[-0.2,0]"))

eps.write(act)

e1=DiscreteTransition(n1,n2,res,act)

act=Zonotope(Rectangle("[0,0.2]x[0,0.2]"))

eps.write(act)

e2=DiscreteTransition(n2,n1,res,act)

h=HybridAutomaton("H")

h.add(e1)
h.add(e2)

c_init=Zonotope(Rectangle("[-0.8,-0.7]x[0.7,0.8]"))
reach=bounded_time_reachability(h, n1, c_init, 0.25, 0.125, 2**8, 2**8, 0, 'yes') 

eps.set_fill_colour("blue")
for bs in reach[0].reached_regions():
  eps.write(bs.continuous_set())

eps.set_fill_colour("red")
for bs in reach[1].reached_regions():
  eps.write(bs.continuous_set())
eps.close()
