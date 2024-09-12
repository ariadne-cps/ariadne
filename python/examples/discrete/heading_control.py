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

def set_workspace():

    pi_ = 3.141592653589793115997963468544185161590576171875

    grid_lengths = FloatDPVector([exact(pi_/4)],double_precision)
    grid_origin = FloatDPVector([0],double_precision)

    control_grid = Grid(grid_origin,grid_lengths)
    control_domain = RealBox([{exact(-pi_-pi_/4):exact(pi_+pi_/4)}])

    print("control_grid=",control_grid)
    print("control_domain=",control_domain)

    state_grid_lengths = FloatDPVector([exact(0.5),exact(0.5),exact(pi_/4)],double_precision)
    state_grid_origin = FloatDPVector([0,0,0],double_precision)
    state_grid = Grid(state_grid_origin,state_grid_lengths)

    theta_domain = RealInterval({0:exact(2*pi_)})
    state_domain = RealBox([{0:5},{0:5},theta_domain])

    x = RealVariable("x")
    y = RealVariable("y")
    theta = RealVariable("theta")
    u = RealVariable("u")

    control_variables = {u}

    deltat = RealConstant("deltat",dec_(0.1))
    v = RealConstant("v",3)

    func = make_function([x,y,theta,u],[x+deltat*v*cos(theta),y+deltat*v*sin(theta),theta+u,1])

    depth = 0
    eps = exact(1e-10)
    probability_threshold = 1e-4

    dynamics = func.zeros(func.argument_size()-len(control_variables),func.argument_size())
    for idx in range(0,dynamics.result_size()):
        dynamics[idx] = func[idx]

    print("dynamics=",dynamics)

    ra = ReachAvoid("heading",dynamics,state_grid,state_domain,control_grid,control_domain,depth,eps,probability_threshold)

    ra.add_obstacle([{1:exact(3.5)},{exact(4.5):5},theta_domain])
    ra.add_obstacle([{0:1},{2:3},theta_domain])
    ra.add_obstacle([{exact(2.5):5},{2:3},theta_domain])
    ra.add_obstacle([{0:5},{0:exact(0.5)},theta_domain])

    ra.add_goal([{4:5},{exact(4.5):5},theta_domain])

    return ra

def test_reach_avoid():

    ra = set_workspace()

    print("state size:", ra.state_size())
    print("control size:", ra.control_size())

    print("obstacles size:", ra.obstacles_size())
    print("goals size:", ra.goals_size())
    print("unverified size:", ra.unverified_size())

    ra.compute_reachability_graph()
    print("num_transitions after computing reachability graph:", ra.num_transitions())

    ra.refine_to_safety_graph()
    print("num_transitions after refining to safety graph:", ra.num_transitions())
    print("safe goal-reachable abstract states:", ra.num_sources())

    ra.update_unverified()
    print("unverified abstract states:", ra.unverified_size(), "(", ra.unverified_percentage(), "% left)")


if __name__  == '__main__':
    test_reach_avoid()
