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
from ariadne.geometry import *
from ariadne.linear_algebra import *
from ariadne.evaluation import *

sl=PointList()
s=Point(2)
s[0]=0.125
s[1]=0.125
sl.append(s)
s[0]=2.375
sl.append(s)
s[1]=2.375
sl.append(s)
r=Simplex(sl)
s[0]=1.5
sl.append(s)
q=Polyhedron(sl)

A=Matrix(2,2);
b=Vector(2);

A[0,0]=-1.0;
A[1,1]=-1.0;

M=AffineMap(A,b)

q2=M(q);

c=convex_hull(q,q2);

eps=EpsPlot("polyhedron.eps",Rectangle("[-3,3]x[-3,3]"))
eps.set_fill_colour("yellow")

eps.write(c)

eps.set_fill_colour("blue")

eps.write(q)
eps.write(q2)

eps.close()

print sl
print r
print q
print q.vertices()
