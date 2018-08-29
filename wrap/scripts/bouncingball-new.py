#!/usr/bin/python

##############################################################################
#            bouncingball.py
#
#  Copyright 2008-9  Davide Bresolin, Pieter Collins
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


from hybrid_system import *
from hybrid_simulation import *

if __name__=='__main__':

    all=~EventSet()
    any=BooleanConstant('any',True)

    # Set the system parameters
    a=RealConstant("a",0.5)
    g=RealConstant("g",9.8)

    bouncingball=HybridSystem()

    x=RealVariable("x");
    v=RealVariable("v");

    bounce=Event("bounce");
    print bounce,~bounce

    bouncingball.new_dynamic(any,x,v);
    bouncingball.new_dynamic(any,v,-g);

    bouncingball.new_invariant(any,x>=0);

    bouncingball.new_guard(bounce,any,(x<=0) & (v<0));
    bouncingball.new_reset(bounce,any,v,-a*v);
    bouncingball.new_reset(bounce,any,x,Constant(0.0));

    bouncingball.new_reset(~all,any,x,x);
    bouncingball.new_reset(~bounce,any,v,v);

    print "BouncingBall = \n",bouncingball;


    initial_state={x:1.0, v:0.0}
    print bouncingball.active_equations({})

    orb=simulate(bouncingball,initial_state,3.0)
    print
    print orb

    plot_orbit(orb,[x])
