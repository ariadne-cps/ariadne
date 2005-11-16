#!/usr/bin/python

from ariadne import *
from geometry import *

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

r=Rectangle(2)
r[0]=Interval(1.1,2.1)
r[1]=Interval(0.5,1.5)
r[0]=Interval(0,1)
r[1]=Interval(0,1)

ps=PartitionScheme(r)
pt=PartitionTree(ps,bt)
pts=PartitionTreeSet(ps,bt,bma)

print ps
print pt
print pts
