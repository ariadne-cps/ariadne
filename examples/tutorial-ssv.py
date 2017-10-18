#!/usr/bin/python

from ariadne import *
print dir()

a=Real(Rational(3,2))
x=RealVariable("x")
y=RealVariable("y")
z=RealVariable("z")

B=RealVariablesBox([x|{2:3},y|{3:4},4<=(z<=5)])

#print make_set(B)
print make_box([x,y,z],B)
print make_set([z,x,y],B)

c1=ContinuousPredicate(x*y<=a)
print c1
#C=RealExpressionConstraintSet([ContinuousPredicate(x*y<=a)])
#C=RealExpressionConstraintSet([x*y<=1, x+z>=2, (x**2 + y**2) <= 5])
#print C

