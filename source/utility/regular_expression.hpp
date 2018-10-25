/***************************************************************************
 *            regular_expression.hpp
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

/*! \file regular_expression.hpp
 *  \brief Regular expressions and finite automata.
*/

#ifndef ARIADNE_REGULAR_EXPRESSION_HPP
#define ARIADNE_REGULAR_EXPRESSION_HPP

#include <vector>
#include <variant>
#include <memory>
#include <iosfwd>

#include "utility/container.hpp"
#include "utility/stlio.hpp"

namespace Ariadne {

using Char = char;
template<class... TS> using Variant = std::variant<TS...>;
template<class T> using SharedPointer = std::shared_ptr<T>;

template<class T> class Word : public List<T> {
  public:
    using List<T>::List;
    T const& head() const { return this->front(); }
    Word<T> const tail() const { return Word<T>(this->begin()+1,this->end()); }
    friend OutputStream& operator<<(OutputStream& os, Word<T> const& wrd) {
        for(auto chr : wrd) { os << chr; } return os; }
};

template<class T> class Empty;
template<class T> class Simple;
template<class T> class Catenation;
template<class T> class Alternation;
template<class T> class Repetition;


template<class T> class RegularExpression {
    typedef Variant<Empty<T>,Simple<T>,Catenation<T>,Alternation<T>,Repetition<T>> DataType;
    SharedPointer<DataType> _ptr;
  public:
    inline RegularExpression();
    inline RegularExpression(DataType var);
    inline RegularExpression(Empty<T> var);
    inline RegularExpression(Simple<T> var);
    inline RegularExpression(Alternation<T> var);
    inline RegularExpression(Catenation<T> var);
    inline RegularExpression(Repetition<T> var);
    inline RegularExpression(T t);
    inline RegularExpression(InitializerList<T> cat);
    inline RegularExpression(InitializerList<RegularExpression<T>> cat);
    operator DataType const& () const;
    template<class R> Bool holds_alternative() const;
    template<class R> R const& extract_alternative() const;
    template<class TT> friend RegularExpression<TT> simplify(RegularExpression<TT> re);
    friend OutputStream& operator<<(OutputStream& os, RegularExpression<T> const& re) { return re._write(os); }
    friend Alternation<T> operator|(RegularExpression<T> const&, RegularExpression<T> const&);
    friend Alternation<T> operator|(Alternation<T>, RegularExpression<T> const&);
    friend Catenation<T> operator,(RegularExpression<T> const&, RegularExpression<T> const&);
    friend Catenation<T> operator,(Catenation<T>, RegularExpression<T> const&);
    Repetition<T> operator++(int) const;
  private:
    OutputStream& _write(OutputStream& os) const;
};

template<class T> class InfiniteRegularExpression;

template<class T> class InfiniteRepetition {
    InfiniteRepetition(RegularExpression<T>);
};

template<class T> class InfiniteAlternation {
    InfiniteAlternation(List<InfiniteRegularExpression<T>>);
};

template<class T> class InfiniteRegularExpression {
    InfiniteRegularExpression(InfiniteRepetition<T>);
    InfiniteRegularExpression(InfiniteAlternation<T>);
};



template<class T> class Empty {
  public:
    Empty() { }
    Bool models(Word<T> const& wrd) const { return false; }
    friend OutputStream& operator<<(OutputStream& os, Empty<T> const& re) {
        os << "?"; return os; }
};

template<class T> class Simple {
    List<T> _elements;
  public:
    Simple(List<T> elmts) : _elements(elmts) { }
    Bool models(Word<T> const& wrd) const { return this->_elements == wrd; }
    friend OutputStream& operator<<(OutputStream& os, Simple<T> const& re) {
        for(auto t : re._elements) { os << t; } return os; }
};

template<class T> class Alternation {
    List<RegularExpression<T>> _possibilities;
  public:
    Alternation(List<RegularExpression<T>> possibilities) : _possibilities(possibilities) { }
//    Alternation(Alternation<T> const&) = default;
//    Alternation(Alternation<T>&&) = default;
    List<RegularExpression<T>> const& possibilities() const { return _possibilities; }
    Void adjoin(RegularExpression<T> re) { _possibilities.append(re); }
    Bool models(Word<T> const& wrd) const {
        for(auto re : _possibilities) { if (re.models(wrd)) { return true; } } return false; }
    friend Alternation<T> operator|(RegularExpression<T> const& r1, RegularExpression<T> const& r2) {
        return Alternation<T>({r1,r2}); }
    friend Alternation<T> operator|(Alternation<T> r1, RegularExpression<T> const& r2) {
        r1._possibilities.push_back(r2); return r1; }
    friend OutputStream& operator<<(OutputStream& os, Alternation<T> const& re) {
        return write_sequence(os,re._possibilities.begin(),re._possibilities.end(),"|"); }
};

template<class T> class Catenation {
    List<RegularExpression<T>> _elements;
  public:
    Catenation(List<RegularExpression<T>> elmts) : _elements(elmts) { }
    List<RegularExpression<T>> const& elements() const { return _elements; }
    Void append(RegularExpression<T> re) { _elements.append(re); }
    Bool models(Word<T> const& wrd) const;
    friend Catenation<T> operator,(RegularExpression<T> const& r1, RegularExpression<T> const& r2) {
        return Catenation<T>({r1,r2}); }
    friend Catenation<T> operator,(Catenation<T> r1, RegularExpression<T> const& r2) {
        r1._elements.push_back(r2); return r1; }
    friend OutputStream& operator<<(OutputStream& os, Catenation<T> const& re) {
        for (auto elmt : re.elements()) {
            if (elmt.template holds_alternative<Alternation<T>>()) { os << "(" << elmt << ")"; }
            else { os << elmt; }
        } return os; }
};

template<class T> class Repetition {
    RegularExpression<T> _single;
  public:
    Repetition(RegularExpression<T> single) : _single(single) { }
    Bool models(Word<T> const& wrd) const;
    friend OutputStream& operator<<(OutputStream& os, Repetition<T> const& re) {
        return os << "(" << re._single << ")*"; }
};

template<class T> inline Repetition<T> RegularExpression<T>::operator++(int) const { return Repetition<T>(*this); }


template<class T> inline RegularExpression<T>::RegularExpression() : RegularExpression(Empty<T>()) { }
template<class T> inline RegularExpression<T>::RegularExpression(DataType var) : _ptr(new DataType(std::move(var))) { }
template<class T> inline RegularExpression<T>::RegularExpression(Empty<T> var) : _ptr(new DataType(std::move(var))) { }
template<class T> inline RegularExpression<T>::RegularExpression(Simple<T> var) : _ptr(new DataType(std::move(var))) { }
template<class T> inline RegularExpression<T>::RegularExpression(Alternation<T> var) : _ptr(new DataType(std::move(var))) { }
template<class T> inline RegularExpression<T>::RegularExpression(Catenation<T> var) : _ptr(new DataType(std::move(var))) { }
template<class T> inline RegularExpression<T>::RegularExpression(Repetition<T> var) : _ptr(new DataType(std::move(var))) { }
template<class T> inline RegularExpression<T>::RegularExpression(T t) : RegularExpression(Simple<T>(t)) { }
template<class T> inline RegularExpression<T>::RegularExpression(InitializerList<T> cat) : RegularExpression(Simple<T>(cat)) { }
template<class T> inline RegularExpression<T>::RegularExpression(InitializerList<RegularExpression<T>> cat) : RegularExpression(Catenation<T>(cat)) { }

template<class T> inline RegularExpression<T>::operator DataType const& () const { return *_ptr; }
template<class T> template<class R> inline Bool RegularExpression<T>::holds_alternative() const { return std::holds_alternative<R>(*_ptr); }
template<class T> template<class R> inline R const& RegularExpression<T>::extract_alternative() const { return std::get<R>(*_ptr); }

template<class T> inline OutputStream& RegularExpression<T>::_write(OutputStream& os) const {
    std::visit([&os](auto const& re){os<<re;},*this->_ptr); return os; }


} // namespace Ariadne

#endif /* ARIADNE_BINARY_WORD_HPP */
