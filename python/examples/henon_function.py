#!/usr/bin/python

##############################################################################
#            henon_function.py
#
#  Copyright 2007  Pieter Collins <Pieter.Collins@cwi.nl>
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

from ariadne import *
import sys

henon_function=InterpretedFunction()
henon_function.read("henon_function.dat")
print henon_function

v=IntervalVector("[0.125,0.25,1.5,0.375]")
print henon_function(v)
print henon_function.jacobian(v)
print

param=IntervalPoint("(1.5,0.375)")
henon_map=FunctionMap(henon_function,param)
print henon_map
print

pt=IntervalPoint("(0.125,0.25)")
print pt
print henon_map(pt)
print henon_map.image(pt)
print henon_map.jacobian(pt)
print
print

square_function=InterpretedFunction("function square output Real y; input Real x; algorithm y:=x^2; end square;")
print square_function
print
v=IntervalVector("[3]")
print v
print square_function(v)
print square_function.jacobian(v)
