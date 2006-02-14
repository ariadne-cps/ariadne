#!/usr/bin/python
from ariadne import *
from geometry import *

sl=StateList()
s=State(2)
s[0]=Dyadic(0.125)
s[1]=Dyadic(0.125)
sl.append(s)
s[0]=Dyadic(2.375)
sl.append(s)
s[1]=Dyadic(2.375)
sl.append(s)
r=Simplex(sl)
s[0]=Dyadic(1.5)
sl.append(s)
q=Polyhedron(sl)

print sl
print r
print q
print q.vertices()
