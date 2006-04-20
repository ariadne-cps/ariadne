#!/usr/bin/python

##############################################################################
#            test_lattice.py
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

ary=IndexArray(3)
ary[0]=-1
ary[1]=2
ary[2]=3
lc1=LatticeCell(ary)
assert str(lc1)=="[-1,0]x[2,3]x[3,4]"
ary[0]=-5
ary[2]=4
lc2=LatticeCell(ary)
assert lc1!=lc2
assert lc1>lc2
assert not (lc1<lc2)
print lc1,lc2, lc1<lc2, lc2>lc1, lc2<lc1, lc1>lc2
lc3=LatticeCell(lc1)
print lc1, lc3
assert lc3==lc1
print
lcls=LatticeCellListSet(lc1.dimension())
lcls.adjoin(lc1)
lcls.adjoin(lc2)
print lcls[0],lcls[1]
ary[0]=0
lc3=LatticeCell(ary)
ary[1]=6
lc4=LatticeCell(ary)
ary[0]=0
ary[1]=0
ary[2]=0
lc5=LatticeCell(ary)
lcls.adjoin(lc3)
lcls.adjoin(lc1)
lcls.adjoin(lc4)
lcls.adjoin(lc2)
lcls.adjoin(lc1)
print lcls
lcls.unique_sort()
print lcls
sys.exit()

imglcls=LatticeCellListSet(2)
imgary=IndexArray(2)
imgary[0]=4
imgary[1]=2
imglcls.adjoin(LatticeCell(imgary))
imgary[0]=5
imglcls.adjoin(LatticeCell(imgary))
lm=LatticeMap(3,2)
lm.adjoin_to_image(lc1,imglcls)
assert str(lc3)=="[0,1]x[2,3]x[4,5]"
assert str(imglcls)=="[[4,5]x[2,3],[5,6]x[2,3]]"
lm.adjoin_to_image(lc3,imglcls)
imgary[1]=3
lm.adjoin_to_image(lc3,LatticeCell(imgary))
print lm
print

print lc1, lm(lc1)
print lc2, lm(lc2)
print lc3, lm(lc3)
print lc4, lm(lc4)
print
lcls=LatticeCellListSet(3)
lcls.adjoin(lc1)
lcls.adjoin(lc2)
lcls.adjoin(lc3)
print lcls, lm(lcls)
print
print lm
print "Passed"
