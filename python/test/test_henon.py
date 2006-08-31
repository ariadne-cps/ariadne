#!/usr/bin/python

##############################################################################
#            test_henon.py
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
from ariadne.output import *
import sys

h=HenonMap(Real(1.5),Real(0.275))

x=Point(2)
#print x,h(x),h(h(x)),"\n\n"

r=Rectangle("[1.4,1.6]x[-0.9,1.1]")
#print r,h(r),h(h(r)),"\n\n"
print h(r)

p=Parallelotope(r)
sp=p.subdivide()
hp=h(p)
hsp=apply(h,sp)

print subset(hsp[0],hp)
print subset(hsp[1],hp)
print subset(hsp[2],hp)
print subset(hsp[3],hp)
print subset(hsp[3].subdivide()[0],hp)

eps=EpsPlot("henon.eps",h(r))
eps.set_fill_colour("blue")
eps.write(hp)
eps.set_fill_colour("green")
eps.write(hsp)
eps.close()
