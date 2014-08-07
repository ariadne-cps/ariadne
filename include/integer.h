/***************************************************************************
 *            integer.h
 *
 *  Copyright 2008-10  Pieter Collins
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

/*! \file integer.h
 *  \brief Exact integer class and fixed-precision integer functions.
 */
#ifndef ARIADNE_INTEGER_H
#define ARIADNE_INTEGER_H

#ifdef HAVE_GMPXX_H
#include <gmpxx.h>
#endif // HAVE_GMPXX_H

#include <iostream>
#include <string>
#include <cstdlib>
#include <stdint.h>

// Simplifying typedefs for unsigned types
// These may be inclused in other headers,
// but repeating a typedef is not an error
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

namespace Ariadne {

#ifdef DOXYGEN
//! \ingroup NumericModule
//! \brief Integers of arbitrary size with exact arithmetic.
//! (Only available if the Gnu Multiple Precision library (GMP) is installed.)
//! \details
//! Unlike C++ and the Python 2, integer division is performed exactly and returns a rational.
//! The operations \c quot(Integer,Integer) and \c rem(Integer,Integer) can be used to perform integer division.
class Integer { };
#endif // DOXYGEN


#ifdef HAVE_GMPXX_H
class Integer : public mpz_class {
  public:
    using mpz_class::mpz_class;
};
#else
class Integer {
  public:
    Integer() : _value() { }
    Integer(const int& n) : _value(n) { }
    Integer(const std::string& s) : _value(std::atoi(s.c_str())) { }
    int get_i() const  { return _value; }
  private:
    int _value;
};
inline Integer operator+(const Integer& z) {
    return Integer(+z.get_i()); }
inline Integer operator-(const Integer& z) {
    return Integer(-z.get_i()); }
inline Integer operator+(const Integer& z1, const Integer& z2) {
    return Integer(z1.get_i()+z2.get_i()); }
inline Integer operator-(const Integer& z1, const Integer& z2) {
    return Integer(z1.get_i()-z2.get_i()); }
inline Integer operator*(const Integer& z1, const Integer& z2) {
    return Integer(z1.get_i()*z2.get_i()); }
inline bool operator==(const Integer& z1, const Integer& z2) {
    return z1.get_i()==z2.get_i(); }
inline bool operator!=(const Integer& z1, const Integer& z2) {
    return z1.get_i()!=z2.get_i(); }
inline bool operator<=(const Integer& z1, const Integer& z2) {
    return z1.get_i()<=z2.get_i(); }
inline bool operator>=(const Integer& z1, const Integer& z2) {
    return z1.get_i()>=z2.get_i(); }
inline bool operator< (const Integer& z1, const Integer& z2) {
    return z1.get_i()< z2.get_i(); }
inline bool operator> (const Integer& z1, const Integer& z2) {
    return z1.get_i()> z2.get_i(); }
inline std::ostream& operator<<(std::ostream& os, const Integer& z) {
    return os << z.get_i(); }
#endif // HAVE_GMPXX_H

uint32_t fac(uint8_t n);
uint16_t fac(uint16_t n);
uint32_t fac(uint32_t n);
uint64_t fac(uint64_t n);
uint32_t bin(uint8_t n, uint8_t k);
uint16_t bin(uint16_t n, uint16_t k);
uint32_t bin(uint32_t n, uint32_t k);
uint64_t bin(uint64_t n, uint64_t k);

using std::min;
using std::max;

} // namespace Ariadne

#endif
