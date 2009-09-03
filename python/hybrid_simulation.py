#!/usr/bin/python

##############################################################################
#            hybrid_simulation.py
#
#  Copyright 2009  Pieter Collins <Pieter.Collins@cwi.nl>
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

import sys

from hybrid_system import *

try:
    import matplotlib.pyplot as pyplot
except:
    print "No matplotlib library available"


def sorted_str(dct):
    keys=dct.keys()
    keys.sort()
    key=keys[0]
    res="{"+str(key)+":"+str(dct[key])
    for key in keys[1:]:
        res+=","+str(key)+":"+str(dct[key])
    res+="}"
    return res


def rk4step(relations,dynamics,valuation,step_size):
    x0val=dict(valuation)
    for (var,expr) in relations:
        x0val[var]=expr.evaluate(x0val)

    k1val={}; x1val={}
    for (var,expr) in dynamics:
        k1val[var]=expr.evaluate(x0val)
        x1val[var]=x0val[var]+step_size*k1val[var]
    for (var,expr) in relations:
        x1val[var]=expr.evaluate(x1val)



    k2val={}; x2val={}
    for (var,expr) in dynamics:
        k2val[var]=expr.evaluate(x1val)
        x2val[var]=x0val[var]+(0.5*step_size)*k2val[var]
    for (var,expr) in relations:
        x2val[var]=expr.evaluate(x2val)

    dx=0.0; df=0.0
    for var in k2val:
        dx=max(dx,abs(x0val[var]-x1val[var]))
        df=max(df,abs(k1val[var]-k2val[var]))
    L=df/dx
    raw_input("Estimated Lipschitz constant "+str(df)+"/"+str(dx)+" = "+str(df/dx)),

    if step_size*L>0.5:
        return rk4step(relations,dynamics,valuation,step_size/2)

    k3val={}; x3val={}
    for (var,expr) in dynamics:
        k3val[var]=expr.evaluate(x2val)
        x3val[var]=x0val[var]+(0.5*step_size)*k3val[var]
    for (var,expr) in relations:
        x3val[var]=expr.evaluate(x3val)

    k4val={};
    for (var,expr) in dynamics:
        k4val[var]=expr.evaluate(x3val)

    dxval={}
    for (var,expr) in dynamics:
        dxval[var]=(step_size/6.0)*(k1val[var]+2*k2val[var]+2*k3val[var]+k4val[var])

    return dxval


def simulate(system,initial_valuation,final_time):

    result={}

    time=0.0
    steps=0
    step_size=1./8

    val=dict(initial_valuation)
    print val

    print system

    relations=system.active_relations(val)
    dynamics=system.active_dynamics(val)
    guards=system.active_guards(val)
    while time<final_time:
        xval=dict(val)
        for (var,expr) in relations:
            xval[var]=expr.evaluate(xval)
        result[(time,steps)]=xval
        print xval
        active_event=None
        for (event,activation) in guards:
            if activation.evaluate(xval):
                active_event=event
                break
        if active_event:
            transitions=system.active_transitions(active_event,xval)
            resets=system.active_resets(active_event,xval)
            uval={}
            for (var,expr) in transitions:
                uval[var]=expr.evaluate(xval)
            for (var,expr) in resets:
                uval[var]=expr.evaluate(xval)
            steps+=1
            raw_input('active event'+str(active_event)+' at '+str(val)+' to '+str(uval))
            val=uval
            relations=system.active_relations(val)
            dynamics=system.active_dynamics(val)
            guards=system.active_guards(val)
        else:
            dxval=rk4step(relations,dynamics,val,step_size)
            for var in dxval:
                val[var]+=dxval[var]
            time+=step_size

    # Compute extended evaluation at final time
    xval=dict(val)
    for (var,expr) in relations:
        val[var]=expr.evaluate(val)
    result[(time,steps)]=xval

    return result


def plot(orbit,vars):
    for var in vars:
        pyplot.xlabel('t')
        pyplot.ylabel(var.name)
        # Plot a list of (x, y) pairs (tuples or a numpy array would
        # also be OK):
        times=[]
        trace=[]
        for (time,steps) in sorted(orbit.keys()):
            valuation=orbit[(time,steps)]
            if valuation.has_key(var):
                times.append(time)
                trace.append(valuation[var])
        pyplot.plot(times,trace)
    pyplot.show()
    #raw_input('Please press return to continue...\n')


if __name__=='__main__':
    e=Event("TurnOff")
    q=StringVariable('q')
    x=RealVariable('x')
    y=RealVariable('y')
    z=RealVariable('z')

    g=(x>1)
    print g,type(g)

    h=HybridSystem()
    h.new_relation(Constant(True),z,x*x+y*y)
    h.new_dynamic(q=="On",x,+1.25*x-y)
    h.new_dynamic(q=="Off",x,-1.25*x-y)
    h.new_dynamic(Constant(True),y,x-0.25*y+0.25)
    h.new_guard(e,q=="On",x>1.0)
    h.new_reset(e,Constant(True),x,x-1.0)
    h.new_reset(e,Constant(True),y,0.5*y)
    h.new_transition(e,Constant(True),q,Constant("Off"))
    print h

    #print h.variables({q:"On"})
    print h.active_relations({q:"On",})

    init={q:"On", x:0.2, y:0.3}
    orb=simulate(h,init,12.0)
    print
    print sorted_str(orb)

    plot(orb,[x,y])

