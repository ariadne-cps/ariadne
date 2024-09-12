#!/usr/bin/python3

##############################################################################
#            test_verification.py
#
#  Copyright  2019-24  Luca Geretti
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

from pyariadne import *

def u_control():

    deltat = dec_(0.1)
    v = 3

    pi_ = 3.141592653589793115997963468544185161590576171875

    pi_fraction = 4
    grid_lengths = FloatDPVector([exact(pi_/pi_fraction)],double_precision)
    grid_origin = FloatDPVector([0],double_precision)

    x = RealVariable("x")
    y = RealVariable("y")
    theta = RealVariable("theta")
    u = RealVariable("u")

    dynamics = IteratedMap({next(x):x+deltat*v*cos(theta),next(y):y+deltat*v*sin(theta),next(theta):theta+u,next(u):u})

    control_grid = Grid(grid_origin,grid_lengths)
    control_domain = RealBox([{exact(pi_-pi_/pi_fraction):exact(pi_+pi_/pi_fraction)}])

    return dynamics, control_grid, control_domain

def test_reach_avoid():

    dynamics, control_grid, control_domain = u_control()

    print("dynamics=",dynamics)
    print("control_grid=",control_grid)
    print("control_domain=",control_domain)

if __name__  == '__main__':
    test_reach_avoid()
