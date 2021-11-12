#!/usr/bin/python3

##############################################################################
#            algebra_demonstration.py
#
#  Copyright  2009-21  Pieter Collins
##############################################################################

# This file is part of Ariadne.

# Ariadne is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Ariadne is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Ariadne. If not, see <https://www.gnu.org/licenses/>.

# Import all classes in the ariadne module
from pyariadne import *


def linear_algebra_demonstration():
    #! [Linear Algebra demonstration]

    # Create an interval vector
    b=FloatDPBoundsVector([1,{2:3},{x_(3.875):x_(4.125)}],dp)
    print("b:",b)

    Aq=RationalMatrix([[1,2,4],[3,1,2],[0,0,1]])
    print("Aq:",Aq)

    # Create an interval matrix
    A=FloatDPBoundsMatrix([[1,2,4],[3,x_(1.5),2],[0,0,1]],dp)
    print("A:",A)

    # Solve the linear equation Ax=b
    x=solve(A,b)
    print("A\\b:",x)
    #! [Linear Algebra demonstration]




if __name__=='__main__':
    linear_algebra_demonstration()
