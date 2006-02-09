#!/usr/bin/python

from ariadne import *
from evaluation import *
from geometry import *
from linear_algebra import *

ba=BooleanArray(7)
ba[2]=1
ba[4]=1
ba[5]=1
ba[6]=1
bma=BooleanArray(4)
bma[0]=1
bma[1]=1
bma[3]=1

bt=BinaryTree(ba)

r1=Rectangle(2)
r1[0]=Interval(1.1,2.1)
r1[1]=Interval(0.5,1.5)
r2=Rectangle(2)
r2[0]=Interval(0,4)
r2[1]=Interval(0,1)
r=Rectangle(2)
r[0]=Interval(0,1)
r[1]=Interval(0,1)

r=Rectangle("")
print r
r=Rectangle("[0,1]x[1,2]x[2,3]z")
print r
r=Rectangle("[0,1]x[1,2] x [2,3]")
print r
r2=Rectangle("[0.5,7]x[-1,5] x [5,7]")
r=Rectangle(" [0,1] x[1,2], x [2,3]")
print r

print r1,r2,regular_intersection(r1,r2)
print disjoint(r1,r2), interiors_intersect(r1,r2), inner_subset(r1,r2)

ps=PartitionScheme(r2)
pt=PartitionTree(ps,bt)
pts=PartitionTreeSet(ps,bt,bma)

p = Polyhedron(3)

print ps
print pt
print pts
