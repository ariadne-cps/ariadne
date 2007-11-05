/***************************************************************************
 *            rational.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

std::ostream& 
Numeric::operator<<(std::ostream& os, const Rational& q) 
{
  //std::cerr<<"ostream& operator<<(ostream& os, const Rational& q)"<<std::flush;
  //std::cerr<<": q="<<std::flush;
  //std::cerr<<q.get_base()<<std::endl;
  os << q.get_base();
  return os;
}

// FIXME: Allow decimal input without leading zero e.g.  ".25" or "-.25"
std::istream& 
Numeric::operator>>(std::istream& is, Rational& q) 
{
  mpz_class intz;
  char sep;
  bool neg;

  // Test if input is negative
  is >> sep;
  neg = (sep=='-');
  is.putback(sep);

  // Input leading part
  is >> intz;

  sep=is.get();
  if(sep=='.') {
    // Decimal input
    mpz_class numz;
    mpz_class denz;
    std::stringstream numstr;
    std::stringstream denstr;
    if(neg) {
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
    q=intz+mpq_class(numz,denz);
  } else if (sep=='/') {
    // Fraction input
    mpz_class numz=intz;
    mpz_class denz;
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

}
