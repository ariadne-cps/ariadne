#!/usr/bin/python

##############################################################################
#            test_simplex.py
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
from ariadne.numeric import * 
from ariadne.linear_algebra import *
from ariadne.geometry import EpsPlot, Point, Rectangle, Parallelotope, disjoint

import sys

A=Matrix(3,2)
b=Vector(3)
c=Vector(2)
A[0,0]=2.
A[0,1]=1.
A[1,0]=1.
A[1,1]=1.
A[2,0]=1.
A[2,1]=3.
b[0]=1.
b[1]=1.
b[2]=1.
c[0]=4.
c[1]=5.

lp=LinearProgram(A,b,c)
print lp
lp.solve()
print lp.optimizing_point()
print lp.optimal_value()


T=Matrix(4,4)
T[0,0]=0.5
T[0,1]=0.25
T[0,3]=4.0
T[1,0]=1.0
T[1,1]=3.0
T[1,2]=-1.0
T[1,3]=20.0
T[2,0]=1.0
T[2,1]=1.0
T[2,3]=10.0
T[2,0]=1.0
T[3,0]=-2.0
T[3,1]=-4.0
T[3,2]=1.0
T[3,3]=-30.0

lp=LinearProgram(T)
print lp
lp.solve()
print lp.optimizing_point()
print lp.optimal_value()

c=Point(2)
c[0]=-1.29931640625
c[1]=1.5390625
A=Matrix(2,2)
A[0,0]=-0.16796875
A[0,1]=-0.0478515625
A[1,0]=0.0546875
A[1,1]=0.0

r=Rectangle("[-1.359375,-1.25]x[1.421875,1.53125]")
a=r.lower_corner()
b=r.upper_corner()
p=Parallelotope(c,A)
print r,p

n=2
o=Vector(n)
for i in range(0,n):
  o[i]=1;
Ao=A*o

#Construct matrix for testing intersection of rectangle and point
#  Rectangle  a<=x<=b
#  Parallelotope  x==c+Ae,  -1<=e<=1
#
# Translate x'=x-a,  e'=e+1
#   0<=x'<=b-a
#   0<=e'<=2
#   x'+a==c+A(e'-1) ->  x'-Ae' == c-a-A1
#  
#  Introduce slack variables for first two inequalities
#  Introduce auxiliary variables for last equality, changing sign of RHS if necessary
#
#  Need to minimise sum of auxiliary variables -> add sum of last rows to get
#  value function.

T=QMatrix(3*n+1,2*n+1)
a=r.lower_corner()
b=r.upper_corner()
rh=(c-a)-Ao
for i in range(0,n):
  T[i,i]=1
  T[i,2*n]=Rational(b[i])-Rational(a[i])
  T[n+i,n+i]=1
  T[n+i,2*n]=2
  if(rh[i]>=0):
    T[2*n+i,i]=1
    for j in range(0,n):
      T[2*n+i,n+j]=-A[i,j]
    T[2*n+i,2*n]=rh[i]
  else:
    T[2*n+i,i]=-1
    for j in range(0,n):
      T[2*n+i,n+j]=A[i,j]
    T[2*n+i,2*n]=-rh[i]
for i in range(0,n):
  for j in range(0,2*n):
    T[3*n,j]-=T[2*n+i,j]
  T[3*n,2*n]-=T[2*n+i,2*n]
  
print T
TT=QMatrix(T)

lp=QLinearProgram(T)
print lp
lp.solve()
u=lp.optimizing_point()

v=Vector(2)
v[0]=float(eval(str(u[0])+".0"))
v[1]=float(eval(str(u[1])+".0"))

print a,u
print lp.optimal_value()

print
print
print

print "Rectangle",r
print p

print a,b,c,A,"  ",Ao
print TT
print lp.tableau()
print a,v,a+v
print r.contains(a+v), p.contains(a+v)

print
print T
print disjoint(r,p)

bb=Rectangle("[-2,-1]x[1,2]")
eps=EpsPlot("lp.eps",bb)
eps.set_fill_colour("green")
eps.write(r)
eps.set_fill_colour("blue")
eps.write(p)
eps.close()
