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
from ariadne.evaluation import *
from ariadne.geometry import *
from ariadne.linear_algebra import *
from ariadne.output import *
import sys

bb=Rectangle("[-3,3]x[-3,3]")
print bb,"\n"

ps=PartitionScheme(bb)
print ps

#btl=[0,0,0,1,0,0,1,0,1,1,0,1,1,0,1,1,1]
#bal=[1,1,1,0,1,0,1,0,1]
btl=[0,0,1,0,1,1,0,1,1]
bal=[1,1,1,1,1]
btl=[0,1,0,0,1,0,1,1,0,1,1]
bal=[0,1,1,0,1,1]
#bal=[0,1,0,1,0,1]

bt=BooleanArray(len(btl))
for i in range(0,len(btl)):
  bt[i]=btl[i]
bt=BinaryTree(bt)
ba=BooleanArray(len(bal))
for i in range(0,len(bal)):
  ba[i]=bal[i]
print bt, ba

pt=PartitionTree(ps,bt)
print pt
print "Iterating through PartitionTree"
for cell in pt:
  print cell

pts=PartitionTreeSet(pt,ba)
print pts
print "Iterating through PartitionTreeSet"
for cell in pts:
  print cell

rls=RectangleListSet(pts)
print "RectangleListSet"
for rect in rls:
  print rect

print

grls=GridRectangleListSet(pts)
print "GridRectangleListSet"
for grect in grls:
  print grect

print

gms=GridMaskSet(grls)
print "GridMaskSet"
print gms


npts=PartitionTreeSet(gms)
print pts
print npts
 

eps=EpsPlot("pt.eps",bb,0,1,"a","b")
eps.set_fill_colour("green")
eps.write(bb)
#eps.write(pts)
#eps.set_fill_style(0)
#eps.write(pts.partition_tree())
#eps.set_fill_colour("blue")
#eps.write(grls)
#eps.write(gms)
#eps.set_fill_colour("red")
#eps.write(npts)
eps.close()
