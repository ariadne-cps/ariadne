#!/usr/bin/python

##############################################################################
#            solver.py
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
from ariadne import *

Indeterminate=indeterminate()
Float = float

Interval.__repr__=Interval.__str__
Point.__repr__=Point.__str__
TaylorExpression.__repr__=TaylorExpression.__str__


def box_radius(box):
    r=0.0
    for i in range(0,len(box)):
        r=__builtins__.max(r,box[i].radius())
    return r
Box.radius=box_radius
del box_radius

def box_centre(box):
    assert len(box)==2
    return Box([box[0].midpoint(),box[1].midpoint()])
Box.centre=box_centre
del box_centre


def my_join(*args):
    if len(args)==1:
        return args[0]
    else:
        r=join(args[0],args[1])
        for i in range(2,len(args)):
            r=join(r,args[i]);
        return r

    if isinstance(args[0],Interval) or isinstance(args[0],IntervalVector):
        result=[]
        for arg in args:
            assert isinstance(arg,Interval) or isinstance(arg,IntervalVector)
            if isinstance(arg,Interval):
                result.append(arg)
            else:
                for i in range(0,len(arg)):
                    result.append(arg[i])
        return Box(result)
    elif isinstance(args[0],TaylorExpression) or isinstance(args[0],TaylorFunction):
        result=[]
        for arg in args:
            assert isinstance(arg,TaylorExpression) or isinstance(arg,TaylorFunction)
            if isinstance(arg,TaylorExpression):
                result.append(arg)
            else:
                for elmt in arg:
                    result.append(elmt)
        return TaylorFunction(result)
    else:
        print "Unknown type",type(args[0]),"of",args[0]

def lie_derivative(f,g):
    assert(f.argument_size()==g.argument_size())
    j=f.argument_size()-f.result_size()
    result = derivative(g,j) * f[j]
    j+=1
    while j!=f.argument_size():
        result = derivative(g,j) * f[j]
        j!=1

def gradient_vector(f):
    r=[]
    for j in range(0,f.argument_size()):
        r.append(derivative(f,j))
    return PolynomialFunction(r)



def gradient(f,D):
    return f.gradient(D)

range=__builtins__.range


def implicit_solve(function, bounds, variable):
    assert(function.argument_size()>variable)
    n=Box(bounds).dimension()
    f=function
    b=bounds
    j=variable
    print "f:",f,"b:",b,"j:",j
    range=bounds[variable]
    domain=Box(join(bounds[0:j],bounds[j+1:n]))
    print "  bounds=",bounds," range=",range," domain=",domain
    h=TaylorExpression.constant(domain,range)
    print "h:",h
    id=[]
    for ii in __builtins__.range(0,n-1):
        id.append(TaylorExpression.variable(domain,ii))
    id=TaylorFunction(id)
    #id=TaylorExpression.variables(domain)
    print "i:",id
    argument=TaylorFunction(my_join(id[0:j],h,id[j:n-1]))
    print "argument:",argument

    dlst=[]
    for ii in __builtins__.range(0,n):
        dlst.append(IntervalDifferential.constant(1,1,bounds[ii]))
    dlst[j]=IntervalDifferential.variable(1,1,bounds[j],0)
    print dlst
    d=IntervalDifferentialVector(len(dlst),1,1);
    for i in __builtins__.range(0,n):
        d[i]=dlst[i]
    print "d:",d
    steps=4
    while h.error()>1e-4:
        steps=steps-1
        assert(steps>0)
        print d,f(d)
        d1f=f(d)[MultiIndex((1,0))]
        print "diff(f):",f(d),"\nd1f:",d1f
        J=1/d1f.midpoint()
        sf=1-d1f/d1f.midpoint()
        print "sf:",sf
        print "he:",h.error()
        new_error=sf*Interval(-h.error(),+h.error())
        print "ne:",new_error
        m=midpoint(h)
        new_h=m-J*compose(f,join(m,id[0]))+new_error
        h=intersection(h,new_h)
        print "h:",h
    return h

def krawczyk_step(function,derivative,solution,parameters,graph_bound):
    f=function
    id=parameters
    h=solution
    m=midpoint(h)
    print "    m.range():",m.range(),
    d1f=derivative(Box(graph_bound))
    print " d1f:",d1f,
    J=1/d1f.midpoint()
    sf=1-d1f/d1f.midpoint()
    print " he:",h.error(),
    print " sf:",sf.upper(),
    ne=sf*Interval(-h.error(),+h.error())
    print " ne:",ne.upper(),
    dh=J*compose(f,join(m,id[0]))
    print " dh.range():",dh.range()
    new_h=m-dh+ne
    return new_h

def implicit_solve_zeroth_explicit(function, derivative, range, domain):
    assert(derivative(Box([range,domain[0]])).lower()>0)
    f=function
    print "f:",f
    print "g:",derivative
    id=TaylorExpression.variables(domain)
    print "i:",id
    h=TaylorExpression.constant(domain,range)
    print "h:",h
    steps=0
    validated=False
    while h.error()>1e-4:
        steps=steps+1
        assert(steps<64)
        new_h=krawczyk_step(function,derivative,h,id,Box([range,domain[0]]))
        if refines(new_h,h):
            validated=True
        try:
            h=intersection(h,new_h)
        except:
            print "Warning:",h,"and",new_h,"are disjoint"
            h=new_h
            range=h.range()
        print "h:",h
    if not validated:
        h.set_error(h.error()*1.5)
    return h


def implicit_solve_zeroth_centre(function, derivative, range, domain):
    assert(type(domain)==Box)
    assert(domain.dimension()==1)
    print "implicit(",function,",",derivative,",",range,",",domain,")"
    graph_space=Box([range,domain[0]])
    initial_centre_derivative=function(Box([range.lower(),domain[0].midpoint()]))
    final_centre_derivative=function(Box([range.upper(),domain[0].midpoint()]))
    assert(initial_centre_derivative.upper()<0.0)
    assert(final_centre_derivative.lower()>0.0)
    print "    ",graph_space,initial_centre_derivative,final_centre_derivative
    working_range=Interval(range)
    f=function
    print "  f:",f
    print "  g:",derivative
    id=TaylorExpression.variables(domain)
    print "  i:",id
    h=TaylorExpression.constant(domain,range)
    print "  h:",h
    steps=16
    validated=False
    while h.error()>1e-4:
        steps=steps-1
        assert(steps>0)
        new_h=krawczyk_step(function,derivative,h,id,Box([range,domain[0]]))
        if refines(new_h,h):
            validated=True
        else:
            print "    Validation not proved:",(h-new_h).clobber().range(),new_h.error(),h.error()
        try:
            h=intersection(h,new_h)
        except:
            print "Warning:",h,"and",new_h,"are disjoint"
            h=new_h
            print "  ",range,new_h.value()
            assert range.contains(new_h.value())
            working_range=h.range()
        if not validated:
            h.set_error(1.5*h.error())
            print "  h:",h
    print graph_space,initial_centre_derivative,final_centre_derivative,h
    assert(validated)
    return h


def draw_circle():
    (x0,x1)=Polynomial.variables(2)
    g=x0*x0+x1*x1-1
    bbox=Box([[0,1],[1,2]])
    bbox=Box([{0.5:2.5},{-0.5:1.5}])
    var=0
    r=implicit_solve(g,bbox,var)
    fig=Figure()
    fig.draw(bbox)
    fig.set_fill_colour(0,1,0)
    fig.set_bounding_box(bbox)
    set.plot(fig,bbox)
    fig.write("circle")


if __name__=="__main__":
    draw_circle()

