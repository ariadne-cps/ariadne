/***************************************************************************
 *            discrete_automaton.hpp
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

/*! \file discrete_automaton.hpp
 *  \brief Discrete finite automata.
*/

#ifndef ARIADNE_DISCRETE_AUTOMATON_HPP
#define ARIADNE_DISCRETE_AUTOMATON_HPP

#include <vector>
#include <variant>
#include <memory>
#include <iosfwd>

#include "utility/container.hpp"
#include "utility/stlio.hpp"

namespace Ariadne {

template<class T, class F> Set<ResultOf<F(T)>> apply_f(F const& f, Set<T> const& s) {
    Set<ResultOf<F(T)>> r; for(auto t : s) { r.insert(f(t)); } return r; }

template<class T1, class T2> Set<Pair<T1,T2>> cartesian_product(Pair<Set<T1>,Set<T2>> const& sets);
template<class T1, class T2> Set<Pair<T1,T2>> cartesian_product(Set<T1> const& set1, Set<T2> const& set2);

template<class T> class RegularExpression;

template<class Y, class X> class FiniteAutonomousSystem {
  private: public:
    Map<X,Set<X>> _transitions;
    Map<X,Y> _outputs;
    Set<X> _initial;
  public:
    FiniteAutonomousSystem(Map<X,Set<X>> transitions, Map<X,Y> outputs, Set<X> initial)
        : _transitions(transitions), _outputs(outputs), _initial(initial) { }
    Set<X> const state_set() const { return this->_outputs.keys(); }
    Set<Y> const output_set() const { return this->_outputs.values(); }
    Map<X,Set<X>> transition_map() const& { return _transitions; }
    Map<X,Y> output_map() const& { return _outputs; }
    Set<X> reach_states() const;
    Set<Y> reach_outputs() const;
    List<Word<Y>> outputs(SizeType steps) const;
    List<Word<X>> trajectories(SizeType steps) const;
    friend OutputStream& operator<<(OutputStream& os, FiniteAutonomousSystem<Y,X> const& sys) {
        return sys._write(os); }
  private:
    OutputStream& _write(OutputStream& os) const;
};

template<class U, class Y, class X> class FiniteNondeterministicAutomaton {
  private: public:
    Map<Pair<X,U>,Set<X>> _transitions;
    Map<X,Y> _outputs;
    Set<X> _initial;
  public:
    FiniteNondeterministicAutomaton(Map<Pair<X,U>,Set<X>> transitions, Map<X,Y> outputs, Set<X> initial)
        : _transitions(transitions), _outputs(outputs), _initial(initial) { }
    Map<X,Y> output_map() const& { return _outputs; }
    Map<Pair<X,U>,Set<X>> const& transition_map() const { return this->_transitions; }
    Set<X> const state_set() const { return this->_outputs.keys(); }
    Set<U> const input_set() const { return apply_f([](Pair<X,U>const& pr){return pr.second;},this->_transitions.keys()); }
    Set<Y> const output_set() const { return this->_outputs.values(); }

    SizeType number_of_transitions() const;
    List<Word<Y>> outputs(Word<U> input) const;
    List<Word<X>> trajectories(Word<U> input) const;
    friend OutputStream& operator<<(OutputStream& os, FiniteNondeterministicAutomaton<U,Y,X> const& aut) {
        return aut._write(os); }
  private:
    OutputStream& _write(OutputStream& os) const;
};

template<class U, class Y, class X> class FiniteDeterministicAutomaton {
    Map<Pair<X,U>,X> _transitions;
    Map<X,Y> _outputs;
    X _initial;
  public:
    FiniteDeterministicAutomaton(Map<Pair<X,U>,X> transitions, Map<X,Y> outputs, X initial)
        : _transitions(transitions), _outputs(outputs), _initial(initial) { }
    Word<X> trajectory(Word<U> input) const;
    Word<Y> output(Word<U> input) const;
};




template<class Y, class X> auto
FiniteAutonomousSystem<Y,X>::trajectories(SizeType steps) const -> List<Word<X>> {
    List<Word<X>> result;
    for(auto x0 : this->_initial) { result.append(Word<X>({x0})); }

    List<Word<X>> tmp;
    for(SizeType step = 0u; step!=steps; ++step) {
        tmp.swap(result);
        result.clear();
        for(auto trj : tmp) {
            Set<X> nxt = this->_transitions[trj.back()];
            for (auto x : nxt) {
                trj.push_back(x);
                result.append(trj);
                trj.pop_back();
            }
        }
    }

    return result;
}

template<class Y, class X> auto
FiniteAutonomousSystem<Y,X>::outputs(SizeType steps) const -> List<Word<Y>> {
    List<Word<Y>> result;
    List<Word<X>> trajs = this->trajectories(steps);
    for(auto traj : trajs) {
        Word<Y> outpt;
        outpt.reserve(traj.size());
        for(auto state : traj) {
            outpt.append(this->_outputs[state]);
        }
        result.append(outpt);
    }
    return result;
}



template<class U, class Y, class X> auto
FiniteNondeterministicAutomaton<U,Y,X>::trajectories(Word<U> input) const -> List<Word<X>> {
    List<Word<X>> result;
    for(auto x0 : this->_initial) { result.append(Word<X>({x0})); }

    List<Word<X>> tmp;
    for(auto u : input) {
        tmp.swap(result);
        result.clear();
        for(auto trj : tmp) {
            Set<X> nxt = this->_transitions[make_pair(trj.back(),u)];
            for (auto x : nxt) {
                trj.push_back(x);
                result.append(trj);
                trj.pop_back();
            }
        }
    }

    return result;
}

template<class U, class Y, class X> auto
FiniteNondeterministicAutomaton<U,Y,X>::outputs(Word<U> input) const -> List<Word<Y>> {
    List<Word<Y>> result;
    List<Word<X>> trajs = this->trajectories(input);
    for(auto traj : trajs) {
        Word<Y> outpt;
        outpt.reserve(traj.size());
        for(auto state : traj) {
            outpt.push_back(this->_outputs[state]);
        }
        result.append(outpt);
    }
    return result;
}


template<class U, class Y, class X> auto
FiniteDeterministicAutomaton<U,Y,X>::trajectory(Word<U> input) const -> Word<X> {
    Word<X> result;
    result.reserve(input.size()+1u);
    X state=this->_initial;
    result.append(state);
    for(SizeType i=0; i!=input.size(); ++i) {
        state=_transitions[make_pair(state,input[i])];
        result.append(state);
    }
    return result;
}

template<class U, class Y, class X> auto
FiniteDeterministicAutomaton<U,Y,X>::output(Word<U> input) const -> Word<Y> {
    Word<Y> result;
    result.reserve(input.size()+1u);
    X state=this->_initial;
    result.append(this->_outputs[state]);
    for(SizeType i=0; i!=input.size(); ++i) {
        state=_transitions[make_pair(state,input[i])];
        result.append(this->_outputs[state]);
    }
    return result;
}

template<class X, class Y> Map<Pair<X,X>,RegularExpression<Y>> kleene_reduce(Map<Pair<X,X>,RegularExpression<Y>> transitions, Set<X> initial);
template<class X, class Y> Map<X,Map<X,RegularExpression<Y>>> brzozowski_reduce(Map<X,Map<X,RegularExpression<Y>>> transitions, Set<X> initial);

template<class U, class Y, class X> Map<Pair<X,X>,Set<Pair<Y,U>>> kripke(FiniteNondeterministicAutomaton<U,Y,X> automaton);
template<class Y, class X> Map<Pair<X,X>,Set<Y>> kripke(FiniteAutonomousSystem<Y,X> system);

template<class U, class Y, class X> RegularExpression<Pair<Y,U>> kleene_behaviour(FiniteNondeterministicAutomaton<U,Y,X> automaton);
template<class Y, class X> RegularExpression<Y> kleene_behaviour(FiniteAutonomousSystem<Y,X> system);

template<class U, class Y, class X> RegularExpression<Pair<Y,U>> brzozowski_behaviour(FiniteNondeterministicAutomaton<U,Y,X> automaton);
template<class Y, class X> RegularExpression<Y> brzozowski_behaviour(FiniteAutonomousSystem<Y,X> system);

template<class Y1, class Y2, class X1, class X2> auto
parallel_composition(FiniteNondeterministicAutomaton<Y2,Y1,X1>, FiniteNondeterministicAutomaton<Y1,Y2,X2>)
    -> FiniteAutonomousSystem<Pair<Y1,Y2>,Pair<X1,X2>>;


} // namespace Ariadne

#endif /* ARIADNE_BINARY_WORD_HPP */
