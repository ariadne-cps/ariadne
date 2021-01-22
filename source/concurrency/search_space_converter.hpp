/***************************************************************************
 *            concurrency/search_space_converter.hpp
 *
 *  Copyright  2011-20  Luca Geretti
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

/*! \file concurrency/search_space_converter.hpp
 *  \brief Classes for converting data from/to the interger search space.
 */

#ifndef ARIADNE_SEARCH_SPACE_CONVERTER_HPP
#define ARIADNE_SEARCH_SPACE_CONVERTER_HPP

namespace Ariadne {

//! \brief Interface for conversion from/into the integer search space
template<class T>
struct SearchSpaceConverterInterface {
    virtual int to_int(T const& value) const = 0;
    virtual T to_value(int i) const = 0;

    virtual SearchSpaceConverterInterface* clone() const = 0;
    virtual ~SearchSpaceConverterInterface() = default;
};

template<class T> struct Log10SearchSpaceConverter;
template<class T> struct Log2SearchSpaceConverter;
template<class T> struct LinearSearchSpaceConverter;

template<> struct Log10SearchSpaceConverter<ExactDouble> : SearchSpaceConverterInterface<ExactDouble> {
    int to_int(ExactDouble const& value) const override {
        if (value == inf) return std::numeric_limits<int>::max();
        else if (value == -inf) return std::numeric_limits<int>::min();
        else return std::round(log(value.get_d())/log(10.0)); }
    ExactDouble to_value(int i) const override { return ExactDouble(exp(log(10.0)*i)); }
    SearchSpaceConverterInterface* clone() const override { return new Log10SearchSpaceConverter(*this); }
};

template<> struct Log2SearchSpaceConverter<ExactDouble> : SearchSpaceConverterInterface<ExactDouble> {
    int to_int(ExactDouble const& value) const override {
        if (value == inf) return std::numeric_limits<int>::max();
        else if (value == -inf) return std::numeric_limits<int>::min();
        return std::round(log(value.get_d())/log(2.0)); }
    ExactDouble to_value(int i) const override { return ExactDouble(exp(log(2.0)*i)); }
    SearchSpaceConverterInterface* clone() const override { return new Log2SearchSpaceConverter(*this); }
};

template<> struct LinearSearchSpaceConverter<ExactDouble> : SearchSpaceConverterInterface<ExactDouble> {
    int to_int(ExactDouble const& value) const override {
        if (value == inf) return std::numeric_limits<int>::max();
        else if (value == -inf) return std::numeric_limits<int>::min();
        return round(value).get<int>(); }
    ExactDouble to_value(int i) const override { return ExactDouble(i); }
    SearchSpaceConverterInterface* clone() const override { return new LinearSearchSpaceConverter(*this); }
};

template<> struct LinearSearchSpaceConverter<int> : SearchSpaceConverterInterface<int> {
    int to_int(int const& value) const override { return value; }
    int to_value(int i) const override { return i; }
    SearchSpaceConverterInterface* clone() const override { return new LinearSearchSpaceConverter(*this); }
};

} // namespace Ariadne

#endif // ARIADNE_SEARCH_SPACE_CONVERTER_HPP
