/***************************************************************************
 *            numeric/symbolic_real.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file numeric/symbolic_real.h
 *  \brief 
 */



#ifndef ARIADNE_SYMBOLIC_REAL
#define ARIADNE_SYMBOLIC_REAL

#include "../function/symbolic_formula.h"

namespace Ariadne {

class SymbolicReal : NumberObject<SymbolicReal>
{
    SymbolicFormula _s;
  private:
    explicit SymbolicReal(SymbolicFormula s) : _s(s) { }
  public:
    static SymbolicReal create(SymbolicFormula s) { return SymbolicReal(s); }
  public:
    SymbolicReal(int c) : _s(ValidFloat64(c)) { }
    SymbolicReal(double c) : _s(ValidFloat64(c)) { }
    SymbolicReal(ValidFloat64 c) : _s(c) { }
    Void data() const { }
    SymbolicFormula formula() const { return _s; }
};
inline SymbolicReal operator+(SymbolicReal f1, SymbolicReal f2) {
    return SymbolicReal::create(add(f1.formula(),f2.formula())); }
inline SymbolicReal operator-(SymbolicReal f1, SymbolicReal f2) {
    return SymbolicReal::create(sub(f1.formula(),f2.formula())); }
inline SymbolicReal operator*(SymbolicReal f1, SymbolicReal f2) {
    return SymbolicReal::create(mul(f1.formula(),f2.formula())); }
inline SymbolicReal operator/(SymbolicReal f1, SymbolicReal f2) {
    return SymbolicReal::create(div(f1.formula(),f2.formula())); }
inline OutputStream& operator<<(OutputStream& os, SymbolicReal f) {
    return os << f.formula(); }

} // namespace Ariadne

#endif
