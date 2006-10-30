##############################################################################
#            system.py
#
#  Copyright 2006  Pieter Collins, Alberto Casagrande
#  Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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

from ariadne.base import *
from ariadne.geometry import *
from ariadne.system import *

class DiscreteMode:
  _number = 0
  def __init__(self, dynamic, invariant):
    self._dynamic=dynamic
    self._invariant=invariant
    self._name = "n"+str(DiscreteMode._number)
    DiscreteMode._number = DiscreteMode._number + 1
  def dynamic(self):
    return self._dynamic
  def invariant(self):
    return self._invariant
  def name(self):
    return self._name
  def set_name(self, name):
    self._name=name
  def __repr__(self):
    return "%s" % self._name
  def __str__(self):
    return "DiscreteMode \'%s\' = [ Dynamics = [ %s ] Invariant = [ %s ] ]" % (self._name, self._dynamic, self._invariant)
  def is_a_discrete_mode(self):
    return 1
                
class DiscreteTransition:
  _number = 0
  def __init__(self, source, destination, reset, activation):
    try:
      source.is_a_discrete_mode()
    except AttributeError:
      raise "DiscreteTransition::DiscreteTransition() Wrong source type"
    try:
      source.is_a_discrete_mode()
    except AttributeError:
      raise "DiscreteTransition::DiscreteTransition() Wrong destination type"
    self._source=source
    self._destination=destination
    self._reset=reset
    self._activation=activation
    self._name= "e" + str(DiscreteTransition._number)
    DiscreteTransition._number = DiscreteTransition._number + 1 
  def source(self):
    return self._source
  def destination(self):
    return self._destination
  def reset(self):
    return self._reset
  def activation(self):
    return self._activation
  def name(self):
    return self._name
  def set_name(self, name):
    self._name=name
  def __repr__(self):
    return "%s" % self._name
  def __str__(self):
    return "DiscreteTransition \'%s\' = [ From/To = < %s, %s >  Reset = [ %s ] Activation = [ %s ] ]" % (self._name, self._source.name(), self._destination.name(), self._reset, self._activation)
  def is_a_discrete_transition(self):
    return 1

class HybridAutomaton:
  def __init__(self,name):
    self._name=name
    self._discrete_modes=[]
    self._discrete_transitions=[]
  def name(self):
    return self._name
  def discrete_modes(self):
    return self._discrete_modes
  def discrete_transitions(self):
    return self._discrete_transitions
  def add(self,obj):
    try:
      obj.is_a_discrete_mode()
      self.add_discrete_mode(obj)
    except AttributeError:
      try:
        obj.is_a_discrete_transition()
        self.add_discrete_transition(obj)
      except AttributeError:
        raise 'HybridAutomaton::add() Wrong obj type'
  def add_discrete_mode(self,mode):
    try:
      mode.is_a_discrete_mode()
    except AttributeError:
      raise "HybridAutomaton::add_discretemode() Wrong mode type"
    if not(mode in self._discrete_modes):
      self._discrete_modes.append(mode)
  def add_discrete_transition(self,edge):
    try:
      edge.is_a_discrete_transition()
    except AttributeError:
      raise "HybridAutomaton::add_edge() Wrong edge type"
    if not(edge in self._discrete_transitions):
      self._discrete_transitions.append(edge)
      self.add_discrete_mode(edge.source())
      self.add_discrete_mode(edge.destination())
  def get_discrete_transitions_leaving(self,mode):
    try:
      mode.is_a_discrete_mode()
    except AttributeError:
      raise "HybridAutomaton::get_discrete_transitions_leaving() Wrong mode type"
    discrete_transitions=[]
    for e in self._discrete_transitions:
      if (e.source()==mode):
        discrete_transitions.append(e)
    return discrete_transitions
  def get_discrete_transitions_reaching(self,mode):
    try:
      mode.is_a_discrete_mode()
    except AttributeError:
      raise "HybridAutomaton::get_discrete_transitions_reaching() Wrong mode type"
    discrete_transitions=[]
    for e in self._discrete_transitions:
      if (e.destination()==mode):
        discrete_transitions.append(e)
    return discrete_transitions 
  def __repr__(self):
    discrete_modes=""
    discrete_transitions=""
    for i in self._discrete_modes:
      discrete_modes = discrete_modes + i.__str__() + "\n"
    for i in self._discrete_transitions:
      discrete_transitions = discrete_transitions + i.__str__() + "\n"
    return "HybridAutomaton \'%s\' = [ \n%s%s]" % (self._name, discrete_modes, discrete_transitions)

class HybridRegion:
  def __init__(self,discrete_mode, set):
    self._discrete_mode=discrete_mode
    self._set=set
  def discrete_mode(self):
    return self._discrete_mode
  def continuous_set(self):
    return self._set
  def continuous_adjoin(self,c_set):
    self._set.adjoin(c_set)
  def __repf__(self):
    return "HybridRegion = [ \n%s%s]" % (self._discrete_mode.__str__(), self._set.__str__())

class HybridReachSet:
  def __init__(self):
    self._discrete_list=[]
  def add_a_hybrid_region(self,d_mode, c_set):
    for i in self._discrete_list:
      if (i.discrete_mode()==d_mode):
        i.continuous_adjoin(c_set)
        return
    region=HybridRegion(d_mode, c_set)
    self._discrete_list.append(region)
  def reached_region_on(self,d_mode):
    for i in self._discrete_list:
      if (i.discrete_mode()==d_mode):
        return i.continuous_set()
  def reached_regions(self):
    return self._discrete_list


class HybridGridRegion:
  def __init__(self,discrete_mode, grid_dimension):
    self._discrete_mode=discrete_mode
    invariant=discrete_mode.invariant()
    grid=IrregularGrid(invariant.bounding_box(),grid_dimension)
    bounds=over_approximation(invariant.bounding_box(),grid).lattice_set()
    self._finite_grid=FiniteGrid(grid,bounds)
    self._bounds_set=GridMaskSet(self._finite_grid)
    (self._bounds_set).adjoin(GridRectangle(grid,bounds))
    self._set=GridMaskSet(self._finite_grid)
    self._grid=grid
  def discrete_mode(self):
    return self._discrete_mode
  def continuous_set(self):
    return self._set
  def join_over_approximation(self,c_set):
#    self._set.adjoin(over_approximation(c_set,self._finite_grid))
    self._set=join_over_approximation(self._set,c_set)
  def join_under_approximation(self,c_set):
#    self._set.adjoin(under_approximation(c_set,self._finite_grid))
    self._set=join_under_approximation(self._set,c_set)
  def grid_over_approximation(self,c_set):
    output=GridMaskSet(self._finite_grid)
    output.adjoin(over_approximation(c_set,self._finite_grid))
    return output
  def grid(self):
    return self._grid
  def finite_grid(self):
    return self._finite_grid;
  def bounds_set(self):
    return self._bounds_set
  def __repf__(self):
    return "HybridRegion = [ \n%s%s]" % (self._discrete_mode.__str__(), self._set.__str__())

class HybridGridReachSet:
  def __init__(self, H, grid_dimension):
    self._discrete_list=[]
    for dmode in H.discrete_modes():
      (self._discrete_list).append(HybridGridRegion(dmode,grid_dimension))
  def grid_over_approximation(self,d_mode, c_set):
    for i in self._discrete_list:
      if (i.discrete_mode()==d_mode):
        output=i.grid_over_approximation(c_set)
        return output
  def join_over_approximation(self,d_mode, c_set):
    for i in self._discrete_list:
      if (i.discrete_mode()==d_mode):
        i.join_over_approximation(c_set)
        return
  def join_under_approximation(self,d_mode, c_set):
    for i in self._discrete_list:
      if (i.discrete_mode()==d_mode):
        i.join_under_approximation(c_set)
        return
  def reached_region_on(self,d_mode):
    for i in self._discrete_list:
      if (i.discrete_mode()==d_mode):
        return i
  def reached_regions(self):
    return self._discrete_list
