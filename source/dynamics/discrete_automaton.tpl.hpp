/***************************************************************************
 *            discrete_automaton.tpl.hpp
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

#include "utility/metaprogramming.hpp"
#include "utility/regular_expression.hpp"

#include "discrete_automaton.hpp"

namespace Ariadne {

template<class T> using IsIntegral = std::is_integral<T>;

namespace {

enum SystemVariableKind { INPUT, STATE, OUTPUT };

template<SystemVariableKind K, class T> class SystemVariable {
    T _t;
  public:
    SystemVariable<K,T>(T const& t) : _t(t) { }
    friend OutputStream& operator<<(OutputStream& os, SystemVariable<K,T> const& var) {
        switch (K) { case INPUT: os << 'u'; break; case STATE: os << 'x'; break; case OUTPUT: os << 'y'; break; default: break; }
        return os << var._t; }
};

} // namespace

template<class T> class Class {
  protected:
    T _v;
  public:
    Class(T v) : _v(v) { }
    operator T const& () const { return _v; }
    friend Bool operator==(Class<T> const& c1, Class<T> const& c2) { return c1._v == c2._v; }
    friend Bool operator< (Class<T> const& c1, Class<T> const& c2) { return c1._v <  c2._v; }
    friend OutputStream& operator<<(OutputStream& os, Class<T> const& c) { return os << c._v; }
};

template<class T> struct Input : public Class<T> {
    Input<T>(T t) : Class<T>(t) { }
    friend OutputStream& operator<<(OutputStream& os, Input<T> const& u) { return os << "u" << static_cast<Class<T>const&>(u); }
};

template<class T> struct State : Class<T> {
    State<T>(T t) : Class<T>(t) { }
};
template<class T> inline OutputStream& operator<<(OutputStream& os, State<T> const& x) { return os << "x" << static_cast<Class<T>const&>(x); }
inline OutputStream& operator<<(OutputStream& os, State<Pair<Char,Char>> const& x) {
    Pair<Char,Char>const& xc=x; return os << "x" << xc.first << xc.second; }

template<class T> struct Output : public Class<T> {
    Output<T>(T t) : Class<T>(t) { }
};
template<class T> inline OutputStream& operator<<(OutputStream& os, Output<T> const& y) { return os << "y" << static_cast<Class<T>const&>(y); }
inline OutputStream& operator<<(OutputStream& os, Output<Pair<Char,Char>> const& y) {
    Pair<Char,Char>const& yc=y; return os << "y" << yc.first << yc.second; }


template<class T1, class T2> Set<Pair<T1,T2>> cartesian_product(Set<T1> const& set1, Set<T2> const& set2) {
    Set<Pair<T1,T2>> result;
    for (auto elmt1 : set1) {
        for (auto elmt2 : set2) {
            result.insert(Pair<T1,T2>(elmt1,elmt2));
        }
    }
    return result;
}

template<class T1, class T2> Set<Pair<T1,T2>> cartesian_product(Pair<Set<T1>,Set<T2>> const& sets) {
    return cartesian_product(sets.first,sets.second);
}


template<class U, class Y, class X> SizeType
FiniteNondeterministicAutomaton<U,Y,X>::number_of_transitions() const {
    SizeType r=0u;
    for (auto transition : this->_transitions) {
        r+=transition.second.size();
    }
    return r;
}

template<class U, class Y, class X> OutputStream&
FiniteNondeterministicAutomaton<U,Y,X>::_write(OutputStream& os) const {
    os << "FiniteNondeterministicAutomaton(\n";
    os << "  transitions=" << Map<Pair<State<X>,Input<U>>,Set<State<X>>>(this->_transitions) << ",\n";
    os << "  output=" << Map<State<X>,Output<Y>>(this->_outputs) << ",\n";
    os << "  initial=" << Set<State<X>>(this->_initial) <<"\n";
    os << ")\n";
    return os;
}


template<class Y, class X> OutputStream&
FiniteAutonomousSystem<Y,X>::_write(OutputStream& os) const {
    os << "FiniteAutonomousSystem(\n";
    os << "  transitions=" << Map<State<X>,Set<State<X>>>(this->_transitions) << ",\n";
    os << "  output=" << Map<State<X>,Output<Y>>(this->_outputs) << ",\n";
    os << "  initial=" << Set<State<X>>(this->_initial) <<"\n";
    os << ")\n";
    return os;
}

template<class Y, class X> Set<X>
FiniteAutonomousSystem<Y,X>::reach_states() const {
    Set<X> reach;
    Set<X> new_reach = this->_initial;
    Set<X> found;

    while (not new_reach.empty()) {
        for (auto source : new_reach) {
            found.adjoin(this->_transitions[source]);
        }
        reach.adjoin(new_reach);
        new_reach = difference(found,reach);
        found.clear();
    }

    return reach;
}

template<class Y, class X> Set<Y>
FiniteAutonomousSystem<Y,X>::reach_outputs() const {
    Set<Y> reach_outputs;
    Set<X> rch = reach_states();
    for (auto state : rch) {
        reach_outputs.insert(_outputs[state]);
    }
    return reach_outputs;
}



template<class Y1, class Y2, class X1, class X2> auto
parallel_composition(FiniteNondeterministicAutomaton<Y2,Y1,X1> aut1, FiniteNondeterministicAutomaton<Y1,Y2,X2> aut2)
    -> FiniteAutonomousSystem<Pair<Y1,Y2>,Pair<X1,X2>>
{
    Map<Pair<X1,X2>,Set<Pair<X1,X2>>> transitions;
    Map<Pair<X1,X2>,Pair<Y1,Y2>> outputs;
    Set<Pair<X1,X2>> initial;

    typedef Y2 U1;
    typedef Y1 U2;

    Map<Pair<X1,U1>,Set<X1>> const& transitions1 = aut1._transitions;
    Map<Pair<X2,U2>,Set<X2>> const& transitions2 = aut2._transitions;
    Map<X1,Y1> const& outputs1 = aut1._outputs;
    Map<X2,Y2> const& outputs2 = aut2._outputs;

    Set<X1> states1=aut1.state_set();
    Set<X2> states2=aut2.state_set();

    for(auto x1 : states1) {
        Y1 const y1 = outputs1[x1];
        for(auto x2 : states2) {
            Y2 const y2 = outputs2[x2];

            U1 const& u1=y2;
            U2 const& u2=y1;

            Set<X1> t1 = transitions1[Pair<X1,U1>(x1,u1)];
            Set<X2> t2 = transitions2[Pair<X2,U2>(x2,u2)];

            transitions[Pair<X1,X2>(x1,x2)]=cartesian_product(t1,t2);
            outputs[Pair<X1,X2>{x1,x2}]=Pair<Y1,Y2>{y1,y2};
        }
    }

    initial = cartesian_product(aut1._initial,aut2._initial);

    return FiniteAutonomousSystem<Pair<Y1,Y2>,Pair<X1,X2>>(transitions,outputs,initial);
}


template<class U, class Y, class X> Map<Pair<X,X>,Set<Pair<Y,U>>> kripke(FiniteNondeterministicAutomaton<U,Y,X> automaton) {
    Map<Pair<X,X>, Set<Pair<Y,U>>> result;

    Map<Pair<X,U>,Set<X>> const& transitions = automaton._transitions;
    Map<X,Y> const& outputs = automaton._outputs;

    for(auto t : transitions) {
        X const& source=t.first.first;
        Y output = outputs[source];
        U const& input=t.first.second;
        Pair<Y,U> output_input(output,input);
        Set<X> const& targets = t.second;
        for(auto target : targets) {
            Pair<X,X> source_target(source,target);
            result[source_target].insert(output_input);
        }
    }
    return result;
}

template<class Y, class X> Map<Pair<X,X>,Set<Y>> kripke(FiniteAutonomousSystem<Y,X> system) {
    Map<Pair<X,X>, Set<Y>> result;

    Map<X,Set<X>> const& transitions = system._transitions;
    Map<X,Y> const& outputs = system._outputs;

    for(auto t : transitions) {
        X const& source=t.first;
        Set<X> const& targets = t.second;
        Y output = outputs[source];
        for(auto target : targets) {
            Pair<X,X> source_target(source,target);
            result[source_target].insert(output);
       }
    }
    return result;
}




template<class X, class Y>
Map<Pair<X,X>,RegularExpression<Y>> kleene_reduce(Map<Pair<X,X>,RegularExpression<Y>> transitions, Set<X> initial) {
    typedef RegularExpression<Y> R;
    X u = transitions.begin()->first.first;
    for (auto t : transitions) {
        u = t.first.first;
        if (not initial.contains(u)) { break; }
    }

    Map<X,R> sources, targets;

    R self;
    Map<Pair<X,X>,R> others;

    // Categorise transitions
    for (auto t : transitions) {
        X const& source = t.first.first;
        X const& target = t.first.second;
        R const& label = t.second;
        if (source == u) {
            if (target==u) { self=label; }
            else { targets[target]=label; }
        } else {
            if (target == u) { sources[source]=label; }
            else { others.insert(t); }
        }
    }

    if (self.template holds_alternative<Repetition<Y>>()) {
    } else if (self.template holds_alternative<Empty<Y>>()) {
        self=Simple<Y>({});
    } else {
        self=self++;
    }

    for (auto from : sources) {
        X const& source = from.first;
        for (auto to : targets) {
            X const& target = to.first;
            R& other = others[make_pair(source,target)];
            other = simplify(static_cast<R>(other | (from.second, self, to.second)));
        }
    }

    return others;
}

template<class X, class Y> RegularExpression<Y> kleene_regular_expression(Map<Pair<X,X>,Set<Y>> const& kripke, Set<X> const& initial) {
    Map<Pair<X,X>,RegularExpression<Y>> result;
    for (auto transition : kripke) {
        Pair<X,X> const& edge = transition.first;
        Set<Y> labels = transition.second;
        if (labels.size()==0) { result.insert(edge,Empty<Y>()); }
        else if (labels.size()==1) { result.insert(edge,Simple<Y>(*labels.begin())); }
        else { result.insert(edge,Alternation<Y>(List<Y>(labels))); }
    }
    std::cerr<<"sz="<<result.size()<<"\n";
    while(result.size()>1) {
        result=kleene_reduce(result,initial);
        std::cerr<<"sz="<<result.size()<<"\n";
    }
    return result.begin()->second;
}




template<class X, class Y>
Map<X,Map<X,RegularExpression<Y>>> brzozowski_reduce(Map<X,Map<X,RegularExpression<Y>>> transitions, Set<X> initial) {
    typedef RegularExpression<Y> R;
    X const& source = transitions.back().first;
    Map<X,RegularExpression<Y>> const& edges = transitions.back().second;

    RegularExpression<Y> loop = (edges.has_key(source)) ? RegularExpression<Y>(edges[source]++) : RegularExpression<Y>(Empty<Y>());

    for(auto other_source : transitions.keys()) {
        if (transitions[other_source].has_key(source)) {
            for (auto next : edges.keys()) {
                if (edges.has_key(source)) {
                    transitions[other_source][next] = transitions[other_source][next] | (transitions[other_source][source],loop,transitions[source][next]);
                } else {
                    transitions[other_source][next] = transitions[other_source][next] | (transitions[other_source][source],transitions[source][next]);
                }
            }
            transitions[other_source].erase(transitions[other_source].find(source));
        }
    }

    transitions.erase(--transitions.end());

    return transitions;
}

template<class X, class Y> RegularExpression<Y> brzozowski_regular_expression(Map<Pair<X,X>,Set<Y>> const& kripke, Set<X> const& initial) {
    Map<X,Map<X,RegularExpression<Y>>> result;
    for (auto transition : kripke) {
        X const& source = transition.first.first;
        X const& target = transition.first.second;
        Set<Y> labels = transition.second;
        RegularExpression<Y> expression;
        if (labels.size()==0) { expression=Empty<Y>(); }
        else if (labels.size()==1) { expression=Simple<Y>(labels.front()); }
        else { expression=Alternation<Y>(List<Y>(labels)); }
        result[source][target]=expression;
    }
    std::cerr<<"sz="<<result.size()<<"\n";
    for(auto transition : result) { std::cerr << transition.first << ":" << transition.second.keys() << " "; } std::cerr<<"\n";

    assert(initial.size()==1);
    while(result.size()>1) {
        result=brzozowski_reduce(result,initial);
        std::cerr<<"sz="<<result.size()<<"\n";
    }
    assert(result.begin()->second.size()==1);
    return result.begin()->second.begin()->second;
}



template<class U, class Y, class X> RegularExpression<Pair<Y,U>> kleene_behaviour(FiniteNondeterministicAutomaton<U,Y,X> automaton) {
    Set<X> i = automaton._initial;
    Map<Pair<X,X>,Set<Pair<Y,U>>> a = kripke(automaton);
    return kleene_regular_expression(a,i);
}

template<class Y, class X> RegularExpression<Y> kleene_behaviour(FiniteAutonomousSystem<Y,X> system) {
    Set<X> i = system._initial;
    Map<Pair<X,X>,Set<Y>> a = kripke(system);
    return kleene_regular_expression(a,i);
}

template<class U, class Y, class X> RegularExpression<Pair<Y,U>> brzozowski_behaviour(FiniteNondeterministicAutomaton<U,Y,X> automaton) {
    Set<X> i = automaton._initial;
    Map<Pair<X,X>,Set<Pair<Y,U>>> a = kripke(automaton);
    return brzozowski_regular_expression(a,i);
}

template<class Y, class X> RegularExpression<Y> brzozowski_behaviour(FiniteAutonomousSystem<Y,X> system) {
    Set<X> i = system._initial;
    Map<Pair<X,X>,Set<Y>> a = kripke(system);
    return brzozowski_regular_expression(a,i);
}

} // namespace Ariadne

