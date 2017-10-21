#!/usr/bin/python

from ariadne import *
print dir()

print

a=Real(Rational(8,5))
b=Rational(1,3)
x=RealVariable("x")
y=RealVariable("y")

f=make_function([x,y],[a-x*x-b*y,x])
print "f:",f

e=Dyadic(1,4) # 1/2^4
bx0=make_box([x,y],[x|{1-e:1+e},(-e<=y)<=e])
print "bx0:",bx0

ss=make_set([x,y],[x*x+y*y<=2])
print "ss:",ss

z0=ValidatedAffineConstrainedImageSet(bx0)
print "z0:",z0

z1=image(z0,f)
bbx1=z1.bounding_box()
print "z1:",z1,bbx1
z2=image(z1,f)
print "z2:",z2


bx1=ExactBox([{FloatDPValue(0.375):FloatDPValue(0.875)},{FloatDPValue(0.875):FloatDPValue(1.125)}])
print "separated(z1,bx1):",z1.separated(bx1)
print "inside(z1,bx1):",z1.inside(bx1)
print "separated(z2,bx1):",z2.separated(bx1)
print "inside(z2,bx1):",z2.inside(bx1)

fig=Figure()
fig.set_bounding_box([{-0.5:1.5},{-0.5:1.5}])
fig.set_fill_colour(1,1,1)
fig.draw(bx1)
fig.set_fill_colour(0.5,0,0)
fig.draw(image(bbx1,f))
fig.set_fill_colour(1,0,0)
fig.draw(z0.bounding_box())
fig.draw(z1.bounding_box())
fig.draw(z2.bounding_box())
fig.set_fill_colour(1,1,0)
fig.draw(z0)
fig.draw(z1)
fig.draw(z2)
fig.write("tutorial-ssv.png")
