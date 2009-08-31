#!/usr/bin/python

##############################################################################
#            hybrid_evolver.py
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

def is_pos(ivl):
    if ivl.lower()>0 : return True
    if ivl.upper()<0 : return False
    return Indeterminate

def is_zero(ivl):
    if ivl.lower()==0 and ivl.upper()==0: return True
    if ivl.lower()>0 or ivl.upper()<0: return False
    return Indeterminate

Interval.is_zero=is_zero
Interval.is_pos=is_pos

class BooleanExpression
    pass

class And(BooleanExpression):
    def __init__(self,first,second):
        self.first=first
        self.second=second
    def evaluate(self,x):
        return self.first.evaluate(x) and self.second.evaluate(x)
    
class Or(BooleanExpression):
    def __init__(self,first,second):
        self.first=first
        self.second=second
    def evaluate(self,x):
        return self.first.evaluate(x) or self.second.evaluate(x)

class Not(BooleanExpression):
    def __init__(self,first):
        self.first=first
        self.second=None
    def evaluate(self,x):
        return -self.first.evaluate(x)
    
def evolve(dynamic, 


Polynomial=PolynomialExpression
if __name__=="__main__":
    one=Polynomial.constant(2,1.0)
    x=Polynomial.variable(2,0)
    y=Polynomial.variable(2,1)

    g0=y-one
    g1=x
    g2=y-x/2+x*x*x/4

    f=[one,x]

    



    