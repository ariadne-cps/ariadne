/***************************************************************************
 *            operators.h
 *
 *  Copyright 2008-9  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


/*! \file operators.h
 *  \brief Numerical operator classes
 */

#ifndef ARIADNE_OPERATORS_H
#define ARIADNE_OPERATORS_H

#include <cstdarg>
#include <iostream>
#include <string>


namespace Ariadne {




struct GtrZero {}; struct LessZero {};

struct Gtr {
    tribool operator()(const Float& x1, const Float& x2) const {
        return (x1==x2) ? indeterminate : tribool(x1>x2); }
    tribool operator()(const Interval& x1, const Interval& x2) const {
        if(x1.lower()>x2.upper()) { return true; } else if(x1.upper()<x2.lower()) { return false; } else { return indeterminate; } }
};

struct Less {
    tribool operator()(const Float& x1, const Float& x2) const {
        return (x1==x2) ? indeterminate : tribool(x1<x2); }
    tribool operator()(const Interval& x1, const Interval& x2) const {
        if(x1.lower()>x2.upper()) { return false; } else if(x1.upper()<x2.lower()) { return true; } else { return indeterminate; } }
};

struct Equal {
    template<class T1, class T2> bool operator()(const T1& a1, const T2& a2) const { return a1 == a2; }
};

struct And { template<class T> T operator()(const T& a1, const T& a2) const { return a1 && a2; } };
struct Or { template<class T> T operator()(const T& a1, const T& a2) const { return a1 || a2; } };
struct Not { template<class T> T operator()(const T& a) const { return !a; } };

struct Add { template<class T> T operator()(const T& a1, const T& a2) const { return a1+a2; } };
struct Sub { template<class T> T operator()(const T& a1, const T& a2) const { return a1-a2; } };
struct Mul { template<class T> T operator()(const T& a1, const T& a2) const { return a1*a2; } };
struct Div { template<class T> T operator()(const T& a1, const T& a2) const { return a1/a2; } };

//struct Pow { template<class T, class N> T operator()(const T& a, const N& n) const { return Ariadne::pow(a,n); } };
struct Pow { Pow(int n) : n(n) { } template<class T> T operator()(const T& a) const { return Ariadne::pow(a,n); } int n; };

struct Neg { template<class T> T operator()(const T& a) const { return Ariadne::neg(a); } };
struct Rec { template<class T> T operator()(const T& a) const { return Ariadne::rec(a); } };
struct Sqr { template<class T> T operator()(const T& a) const { return Ariadne::sqr(a); } };
struct Sqrt { template<class T> T operator()(const T& a) const { return Ariadne::sqrt(a); } };

struct Exp { template<class T> T operator()(const T& a) const { return Ariadne::exp(a); } };
struct Log { template<class T> T operator()(const T& a) const { return Ariadne::log(a); } };
struct Sin { template<class T> T operator()(const T& a) const { return Ariadne::sin(a); } };
struct Cos { template<class T> T operator()(const T& a) const { return Ariadne::cos(a); } };
struct Tan { template<class T> T operator()(const T& a) const { return Ariadne::tan(a); } };


inline std::ostream& operator<<(std::ostream& os, const Less& v) { return os << "<="; }
inline std::ostream& operator<<(std::ostream& os, const Gtr& v) { return os << ">="; }
inline std::ostream& operator<<(std::ostream& os, const Equal& v) { return os << "=="; }

inline std::ostream& operator<<(std::ostream& os, const And& v) { return os << "&&"; }
inline std::ostream& operator<<(std::ostream& os, const Or& v) { return os << "||"; }
inline std::ostream& operator<<(std::ostream& os, const Not& v) { return os << "!"; }

inline std::ostream& operator<<(std::ostream& os, const Add& v) { return os << "+"; }
inline std::ostream& operator<<(std::ostream& os, const Sub& v) { return os << "-"; }
inline std::ostream& operator<<(std::ostream& os, const Mul& v) { return os << "*"; }
inline std::ostream& operator<<(std::ostream& os, const Div& v) { return os << "/"; }

inline std::ostream& operator<<(std::ostream& os, const Neg& op) { return os << "neg"; }
inline std::ostream& operator<<(std::ostream& os, const Rec& op) { return os << "rec"; }
inline std::ostream& operator<<(std::ostream& os, const Pow& op) { return os << "pow"; }
inline std::ostream& operator<<(std::ostream& os, const Sqr& op) { return os << "sqr"; }
inline std::ostream& operator<<(std::ostream& os, const Sqrt& op) { return os << "sqrt"; }

inline std::ostream& operator<<(std::ostream& os, const Exp& op) { return os << "exp"; }
inline std::ostream& operator<<(std::ostream& os, const Log& op) { return os << "log"; }
inline std::ostream& operator<<(std::ostream& os, const Sin& op) { return os << "sin"; }
inline std::ostream& operator<<(std::ostream& os, const Cos& op) { return os << "cos"; }
inline std::ostream& operator<<(std::ostream& os, const Tan& op) { return os << "tan"; }




} // namespace Ariadne

#endif // ARIADNE_OPERATORS_H
