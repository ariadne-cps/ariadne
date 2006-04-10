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


class DiscreteNode:
  _number = 0
  def __init__(self, dynamic, invariant):
    self._dynamic=dynamic
    self._invariant=invariant
    self._name = "n"+str(DiscreteNode._number)
    DiscreteNode._number = DiscreteNode._number + 1
  def dynamic(self):
    return self._dynamic
  def invariant(self):
    return self._invariant
  def name(self):
    return self._name
  def __repr__(self):
    return "%s" % self._name
  def __str__(self):
    return "DiscreteNode \'%s\' = [ Dynamics = [ %s ] Invariant = [ %s ] ]" % (self._name, self._dynamic, self._invariant)
  def is_a_discrete_node(self):
    return 1
		
class DiscreteTransition:
  _number = 0
  def __init__(self, source, destination, reset, activation):
    try:
      source.is_a_discrete_node()
    except AttributeError:
      raise "DiscreteTransition::DiscreteTransition() Wrong source type"
    try:
      source.is_a_discrete_node()
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
  def __repr__(self):
    return "%s" % self._name
  def __str__(self):
    return "DiscreteTransition \'%s\' = [ From/To = < %s, %s >  Reset = [ %s ] Activation = [ %s ] ]" % (self._name, self._source.name(), self._destination.name(), self._reset, self._activation)
  def is_a_discrete_transition(self):
    return 1

class HybridAutomaton:
  def __init__(self,name):
    self._name=name
    self._discrete_nodes=[]
    self._discrete_transitions=[]
  def name(self):
    return self._name
  def discrete_nodes(self):
    return self._discrete_nodes
  def discrete_transitions(self):
    return self._discrete_transitions
  def add(self,obj):
    try:
      obj.is_a_discrete_node()
      self.add_discrete_node(obj)
    except AttributeError:
      try:
        obj.is_a_discrete_transition()
        self.add_discrete_transition(obj)
      except AttributeError:
        raise 'HybridAutomaton::add() Wrong obj type'
  def add_discrete_node(self,node):
    try:
      node.is_a_discrete_node()
    except AttributeError:
      raise "HybridAutomaton::add_discretenode() Wrong node type"
    if not(node in self._discrete_nodes):
      self._discrete_nodes.append(node)
  def add_discrete_transition(self,edge):
    try:
      edge.is_a_discrete_transition()
    except AttributeError:
      raise "HybridAutomaton::add_edge() Wrong edge type"
    if not(edge in self._discrete_transitions):
      self._discrete_transitions.append(edge)
      self.add_discrete_node(edge.source())
      self.add_discrete_node(edge.destination())
  def get_discrete_transitions_leaving(self,node):
    try:
      node.is_a_discrete_node()
    except AttributeError:
      raise "HybridAutomaton::get_discrete_transitions_leaving() Wrong node type"
    discrete_transitions=[]
    for e in self._discrete_transitions:
      if (e.source()==node):
        discrete_transitions.append(e)
    return discrete_transitions
  def get_discrete_transitions_reaching(self,node):
    try:
      node.is_a_discrete_node()
    except AttributeError:
      raise "HybridAutomaton::get_discrete_transitions_reaching() Wrong node type"
    discrete_transitions=[]
    for e in self._discrete_transitions:
      if (e.destination()==node):
        discrete_transitions.append(e)
    return discrete_transitions 
  def __repr__(self):
    discrete_nodes=""
    discrete_transitions=""
    for i in self._discrete_nodes:
      discrete_nodes = discrete_nodes + i.__str__() + "\n"
    for i in self._discrete_transitions:
      discrete_transitions = discrete_transitions + i.__str__() + "\n"
    return "HybridAutomaton \"%s\" = [ \n%s%s]" % (self.name, discrete_nodes, discrete_transitions)

class HybridRegion:
  def __init__(self,discrete_node, set):
    self._discrete_node=discrete_node
    self._set=[set]
  def discrete_node(self):
    return self._discrete_node
  def continuous_set(self):
    return self._set
  def __repf__(self):
    return "HybridRegion = [ \n%s%s]" % (self._discrete_node.__str__(), self._set.__str__())

class HybridReachSet:
  def __init__(self):
    self._discrete_list=[]
  def add_a_hybrid_region(self,d_node, c_set):
    for i in self._discrete_list:
      if (i.discrete_node()==d_node):
        list=i.continuous_set()
	list.append(c_set)
	return
    region=HybridRegion(d_node, c_set)
    self._discrete_list.append(region)
  def reached_region(self,d_node):
    for i in self._discrete_list:
      if (i.discrete_node()==d_node):
        return i.continuous_set()
  def reached_region(self):
    reach_set=[]
    for i in self._discrete_list:
      reach_set.extend(i.continuous_set())
    return reach_set
