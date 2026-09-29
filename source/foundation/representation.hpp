/***************************************************************************
 *            foundation/representation.hpp
 *
 *  Copyright  2026  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 */

#ifndef ARIADNE_FOUNDATIONS_REPRESENTATION_HPP
#define ARIADNE_FOUNDATIONS_REPRESENTATION_HPP

#include "utility/typedefs.hpp"

namespace Ariadne {

template<class T> struct Representation {
    T const* pointer;
    explicit Representation(T const& value) : pointer(&value) { }
    T const& reference() const { return *pointer; }
};

template<class T> inline Representation<T> representation(T const& value) {
    return Representation<T>(value);
}

template<class T> concept HasRepresentation = requires(OutputStream& os, T const& value) {
    value._repr(os);
};

template<class T> inline OutputStream& operator<<(OutputStream& os, Representation<T> const& object) {
    if constexpr (HasRepresentation<T>) {
        object.reference()._repr(os);
        return os;
    } else {
        return os << object.reference();
    }
}

} // namespace Ariadne

#endif /* ARIADNE_FOUNDATIONS_REPRESENTATION_HPP */
