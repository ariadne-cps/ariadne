#!/usr/bin/python

##############################################################################
#            test_zonotope.py
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
from ariadne.evaluation import *
from ariadne.geometry import *
from ariadne.linear_algebra import *
from math import *
import sys

r=Rectangle("[9,11]x[5,11]x[0,0]")
r2=Rectangle("[5,6]x[3,4]x[-0.5,0.5]")

z=Zonotope(r);
z2=Zonotope(r2);

z3=minkowski_sum(z,z2)

p=Point(3)

p[0]=17
p[1]=15
p[2]=-0.5

p2=Point(p)
p3=Point(p)

p2[1]=8
p3[0]=14

p4=Point(p3)
p5=Point(p3)
p6=Point(3)

p6[0]=15.5
p6[1]=11.5
p6[2]=0

p4[1]=p2[1]

p5[1]=p2[1]-3

print "Z1 = ",z 
print "Z2 = ",z2 
print "Z3 = Z1 + Z2 = ",z3 

if (z3.contains(p)):
	print "Z contains the point ",p
else:
	print "Z does not contain the point ",p

if (z3.contains(p2)):
	print "Z contains the point ",p2
else:
	print "Z does not contain the point ",p2
	
if (z3.contains(p3)):
	print "Z contains the point ",p3
else:
	print "Z does not contain the point ",p3
	
if (z3.contains(p4)):
	print "Z contains the point ",p4
else:
	print "Z does not contain the point ",p4

if (z3.contains(p5)):
	print "Z contains the point ",p5
else:
	print "Z does not contain the point ",p5
	

if (z3.interior_contains(p)):
	print "Z's interior contains the point ",p
else:
	print "Z's interior does not contain the point ",p

if (z3.interior_contains(p2)):
	print "Z's interior contains the point ",p2
else:
	print "Z's interior does not contain the point ",p2
	
if (z3.interior_contains(p3)):
	print "Z's interior contains the point ",p3
else:
	print "Z's interior does not contain the point ",p3
	
if (z3.interior_contains(p4)):
	print "Z's interior contains the point ",p4
else:
	print "Z's interior does not contain the point ",p4

if (z3.interior_contains(p5)):
	print "Z's interior contains the point ",p5
else:
	print "Z's interior does not contain the point ",p5
	
if (z3.interior_contains(p6)):
	print "Z's interior contains the point ",p6
else:
	print "Z's interior does not contain the point ",p6

