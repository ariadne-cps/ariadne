/***************************************************************************
 *            rational.code.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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

#include <iostream>
#include <sstream>

#include "numeric/rational.h"

namespace Ariadne { 


double
Rational::get_d() const 
{
  double r; 
  set_(r,*this,round_approx); 
  return r;
}
  

Rational::Rational(const std::string& str) 
{
  mpq_init(_value); 
  std::stringstream ss(str); 
  ss>>*this; 
}


std::ostream& 
operator<<(std::ostream& os, const Rational& q) 
{
  return os << q._value;
}

// FIXME: Allow decimal input without leading zero e.g.  ".25" or "-.25"
std::istream& 
operator>>(std::istream& is, Rational& q) 
{
  Integer intz;
  char sep;
  bool neg_;

  // Test if input is neg_ative
  is >> sep;
  neg_ = (sep=='-');
  is.putback(sep);

  // Input leading part
  is >> intz;

  sep=is.get();
  if(sep=='.') {
    // Decimal input
    Integer numz;
    Integer denz;
    std::stringstream numstr;
    std::stringstream denstr;
    if(neg_) {
      numstr << '-';
    }
    denstr << '1';
    char digit=is.get();
    while(std::isdigit(digit)) {
      numstr << digit;
      denstr << '0';
      digit=is.get();
    }
    is.putback(digit);
    numstr >> numz;
    denstr >> denz;
    q=Rational(intz+Rational(numz,denz));
  } else if (sep=='/') {
    // Fraction input
    Integer numz=intz;
    Integer denz;
    is >> denz;
    q=Rational(numz,denz);
  } else {
    // Integer input 
    is.putback(sep);
    q=intz;
  }
  q.canonicalize();
  return is;
}



} // namespace Ariadne
