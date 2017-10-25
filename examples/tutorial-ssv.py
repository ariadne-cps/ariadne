#!/usr/bin/python

from ariadne import *
print dir()
print "\n"

# Define system constants and variables
a=Real(Rational(8,5))
b=Rational(1,3)
x=RealVariable("x")
y=RealVariable("y")

# Define system transition function
f=make_function([x,y],[a-x*x-b*y,x])
print "f:",f

# Define initial set as a box
e=Dyadic(1,6) # 1/2^4
bx0=make_box([x,y],[x|{1-e:1+e},(-1-e<=y)<=-1+e])
print "bx0:",bx0

# Define safe set using constraints
ss=make_set([x,y],[x*x+y*y<=2])
print "ss:",ss

# Compute first few iterates of initial set
z0=ValidatedAffineConstrainedImageSet(bx0)
print "z0:",z0
z1=image(z0,f)
bbx1=z1.bounding_box()
print "z1:",z1,bbx1
z2=image(z1,f)
print "z2:",z2

# Test inclusion and disjointness of sets
bx1=ExactBox([{FloatDPValue(0.375):FloatDPValue(0.875)},{FloatDPValue(0.875):FloatDPValue(1.125)}])
print "separated(z1,bx1):",z1.separated(bx1)
print "inside(z1,bx1):",z1.inside(bx1)
print "separated(z2,bx1):",z2.separated(bx1)
print "inside(z2,bx1):",z2.inside(bx1)

# Draw sets
fig=Figure()
fig.set_bounding_box([{-0.5:1.5},{-0.5:1.5}])
fig.set_fill_opacity(0.5)
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
fig.set_fill_colour(0.5,1,1)
fig.draw(intersection(ss,[{-2:3},{-2:3}]))
fig.write("tutorial-ssv.png")


# Compute over-approximation to reachable set
grid=Grid(2)
dpth=4
curr_rch=GridTreeSet(grid)
wrk_rch=outer_approximation(z0,grid,dpth)
encl=[[z0]]
rch=[curr_rch]
k=0
while k<20 and not wrk_rch.is_empty():
    fnd_rch=GridTreeSet(grid)
    fnd_encl=[]
    for c in wrk_rch:
        z=ValidatedAffineConstrainedImageSet(c.box())
        fz=image(z,f)
        fnd_encl.append(fz)
        fnd_rch.adjoin_outer_approximation(fz,dpth)
    wrk_rch=difference(fnd_rch,curr_rch)
    curr_rch=union(curr_rch,wrk_rch)
    encl.append(fnd_encl)
    rch.append(curr_rch)

safe=GridTreeSet(grid)
safe.adjoin_inner_approximation(ss,4,dpth)
safe=inner_approximation(ss,grid,4,dpth)

print "Safe?",subset(curr_rch,safe)

fig=Figure()
for k in range(1,len(rch)):
    fig.set_bounding_box([{-1.5:2.5},{-1.5:2.5}])
    fig.set_fill_opacity(0.5)
    fig.set_fill_colour(0.5,0,0)
    fig.draw(rch[k-1])
    fig.set_fill_colour(1,0 ,0)
    fig.draw(difference(rch[k],rch[k-1]))
    enclk=encl[k]
    m=len(enclk)
    for e in encl[k]: fig.draw(e)
    strk=str(k)
    if k<10:  strk="0"+strk
    fig.write("tutorial-ssv-r"+strk+".png")
    fig.clear()


