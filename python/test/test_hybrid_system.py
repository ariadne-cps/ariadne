#!/usr/bin/python

##############################################################################
#            test_hybrid_system.py
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


from ariadne.base import *
from ariadne.numeric import *
from ariadne.linear_algebra import *
from ariadne.geometry import *
from ariadne.system import *
from ariadne.evaluation import *

from ariadne.hybrid_system import *

import sys

b=Vector(2)
A=Matrix(2,2)

A[0,1]=-1;
A[1,0]=+1;
dyn=AffineVectorField(A,b)
inv=Parallelotope(Rectangle("[-1,1]x[-1,1]"))
n1=DiscreteNode(dyn,inv)

A[1,1]=-3;
A[1,0]=+1;
dyn=AffineVectorField(A,b)
inv=Parallelotope(Rectangle("[-1,1]x[-1,1]"))
n2=DiscreteNode(dyn,inv)

A=Matrix(2,2)

A[0,0]=1
A[1,1]=1
res=AffineMap(A,b)
act=Rectangle("[0.5,1]x[0.5,1]")

e1=DiscreteTransition(n1,n2,res,act)

act=Rectangle("[-1,-0.5]x[-1,-0.5]")
e2=DiscreteTransition(n2,n1,res,act)

h=HybridAutomaton("H")

h.add(e1)
h.add(e2)

print h


