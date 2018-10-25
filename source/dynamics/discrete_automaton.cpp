/***************************************************************************
 *            discrete_automaton.cpp
 *
 *  Copyright  2018  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "discrete_automaton.tpl.hpp"

#include "numeric/integer.hpp"

namespace Ariadne {

template class FiniteAutonomousSystem<Char,Char>;
template class FiniteAutonomousSystem<Pair<Char,Char>,Pair<Char,Char>>;

template class FiniteDeterministicAutomaton<Char,Char,Char>;
template class FiniteNondeterministicAutomaton<Char,Char,Char>;

template Map<Pair<Char,Char>,Set<Pair<Char,Char>>> kripke(FiniteNondeterministicAutomaton<Char,Char,Char>);

template Map<Pair<Char,Char>,RegularExpression<Char>> kleene_reduce(Map<Pair<Char,Char>,RegularExpression<Char>> transitions, Set<Char> initial);
template RegularExpression<Pair<Char,Char>> kleene_behaviour(FiniteNondeterministicAutomaton<Char,Char,Char>);
template RegularExpression<Pair<Char,Char>> kleene_behaviour(FiniteAutonomousSystem<Pair<Char,Char>,Pair<Char,Char>>);

template RegularExpression<Pair<Char,Char>> brzozowski_behaviour(FiniteNondeterministicAutomaton<Char,Char,Char>);

template FiniteAutonomousSystem<Pair<Char,Char>,Pair<Char,Char>>
    parallel_composition(FiniteNondeterministicAutomaton<Char,Char,Char>, FiniteNondeterministicAutomaton<Char,Char,Char>);


template class FiniteNondeterministicAutomaton<Int16,Int16,Int32>;
template class FiniteAutonomousSystem<Pair<Int16,Int16>,Pair<Int32,Int32>>;

template Map<Pair<Int32,Int32>,Set<Pair<Int16,Int16>>> kripke(FiniteNondeterministicAutomaton<Int16,Int16,Int32>);

template Map<Pair<Int32,Int32>,RegularExpression<Int16>> kleene_reduce(Map<Pair<Int32,Int32>,RegularExpression<Int16>> transitions, Set<Int32> initial);
template RegularExpression<Pair<Int16,Int16>> kleene_behaviour(FiniteNondeterministicAutomaton<Int16,Int16,Int32>);
template RegularExpression<Pair<Int16,Int16>> kleene_behaviour(FiniteAutonomousSystem<Pair<Int16,Int16>,Pair<Int32,Int32>>);

template RegularExpression<Pair<Int16,Int16>> brzozowski_behaviour(FiniteNondeterministicAutomaton<Int16,Int16,Int32>);

template FiniteAutonomousSystem<Pair<Int16,Int16>,Pair<Int32,Int32>>
    parallel_composition(FiniteNondeterministicAutomaton<Int16,Int16,Int32>, FiniteNondeterministicAutomaton<Int16,Int16,Int32>);



} // namespace Ariadne

