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


def rk4step(equations,dynamics,valuation,step_size):
    x0val=dict(valuation)
    for (var,expr) in equations:
        x0val[var]=expr.evaluate(x0val)

    k1val={}; x1val={}
    for (var,expr) in dynamics:
        k1val[var]=expr.evaluate(x0val)
        x1val[var]=x0val[var]+step_size*k1val[var]
    for (var,expr) in equations:
        x1val[var]=expr.evaluate(x1val)



    k2val={}; x2val={}
    for (var,expr) in dynamics:
        k2val[var]=expr.evaluate(x1val)
        x2val[var]=x0val[var]+(0.5*step_size)*k2val[var]
    for (var,expr) in equations:
        x2val[var]=expr.evaluate(x2val)

    dx=0.0; df=0.0
    for var in k2val:
        dx=max(dx,abs(x0val[var]-x1val[var]))
        df=max(df,abs(k1val[var]-k2val[var]))
    L=df/dx
    #raw_input("Estimated Lipschitz constant "+str(df)+"/"+str(dx)+" = "+str(df/dx)),

    if step_size*L>0.5:
        return rk4step(equations,dynamics,valuation,step_size/2)

    k3val={}; x3val={}
    for (var,expr) in dynamics:
        k3val[var]=expr.evaluate(x2val)
        x3val[var]=x0val[var]+(0.5*step_size)*k3val[var]
    for (var,expr) in equations:
        x3val[var]=expr.evaluate(x3val)

    k4val={};
    for (var,expr) in dynamics:
        k4val[var]=expr.evaluate(x3val)

    dxval={}
    for (var,expr) in dynamics:
        dxval[var]=(step_size/6.0)*(k1val[var]+2*k2val[var]+2*k3val[var]+k4val[var])

    dxval[RealVariable("time")]=step_size

    nval=dict(valuation)
    for var in dxval:
        nval[var]+=dxval[var]

    return nval



def compute_crossing(guard,equations,dynamics,initial_valuation,final_valuation,time_step):
    # Use bisection
    lower_time_step=0.0
    upper_time_step=time_step
    for i in range(0,6):
        new_time_step=(lower_time_step+upper_time_step)/2
        new_valuation=rk4step(equations,dynamics,initial_valuation,new_time_step)
        new_guard_value=guard.evaluate(new_valuation)
        if new_guard_value:
            upper_time_step=new_time_step
        else:
            lower_time_step=new_time_step
    return new_valuation


def simulate(system,initial_valuation,final_time,final_steps=6):

    result=[]

    time=0.0
    steps=0
    step_size=1./32

    val=dict(initial_valuation)
    print val

    time_var=RealVariable("time")
    steps_var=IntegerVariable("steps")

    val[time_var]=0.0
    val[steps_var]=0

    print system

    equations=system.active_equations(val)
    dynamics=system.active_dynamics(val)
    guards=system.active_guards(val)
    while val[time_var]<final_time and val[steps_var]<final_steps:
        xval=dict(val)
        for (var,expr) in equations:
            xval[var]=expr.evaluate(xval)
        result.append(xval)
        print xval
        active_event=None
        for (event,activation) in guards:
            if activation.evaluate(xval):
                active_event=event
                active_guard=activation
                break
        if active_event:
            if len(result)>1:
                pval=result[-2]
                if pval[steps_var]==xval[steps_var]:
                    last_time_step=xval[time_var]-pval[time_var]
                    xval=compute_crossing(active_guard,equations,dynamics,pval,xval,last_time_step)
                    result.pop()
                    result.append(xval)

            transitions=system.active_transitions(active_event,xval)
            resets=system.active_resets(active_event,xval)
            uval={}
            for (var,expr) in transitions:
                uval[var]=expr.evaluate(xval)
            for (var,expr) in resets:
                uval[var]=expr.evaluate(xval)
            uval[time_var]=xval[time_var]
            uval[steps_var]=xval[steps_var]+1
            steps+=1
            #raw_input('active event'+str(active_event)+' at '+str(val)+' to '+str(uval))
            val=uval
            equations=system.active_equations(val)
            dynamics=system.active_dynamics(val)
            guards=system.active_guards(val)
        else:
            nval=rk4step(equations,dynamics,val,step_size)
            val=nval

    # Compute extended evaluation at final time
    xval=dict(val)
    for (var,expr) in equations:
        val[var]=expr.evaluate(val)
    result.append(xval)

    return result


def plot_orbit(orbit,vars):
    for var in vars:
        pyplot.xlabel('t')
        pyplot.ylabel(var.name)
        # Plot a list of (x, y) pairs (tuples or a numpy array would
        # also be OK):
        times=[]
        trace=[]
        for valuation in orbit:
            time=valuation[RealVariable("time")]
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
    h.new_equation(Constant(True),z,x*x+y*y)
    h.new_dynamic(q=="On",x,+1.25*x-y)
    h.new_dynamic(q=="Off",x,-1.25*x-y)
    h.new_dynamic(Constant(True),y,x-0.25*y+0.25)
    h.new_guard(e,q=="On",x>1.0)
    h.new_reset(e,Constant(True),x,x-1.0)
    h.new_reset(e,Constant(True),y,0.5*y)
    h.new_transition(e,Constant(True),q,Constant("Off"))
    print h

    #print h.variables({q:"On"})
    print h.active_equations({q:"On",})

    init={q:"On", x:0.2, y:0.3}
    orb=simulate(h,init,12.0)
    print
    print orb

    plot_orbit(orb,[x,y])

