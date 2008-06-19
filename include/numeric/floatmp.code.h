/***************************************************************************
 *            floatmp.code.h
 *
 *  Copyright  2006-7 Pieter Collins
 * 
 * Contains code from gmpfrxx by Jon Wilkening 
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

#include "numeric/floatmp.h"
#include "numeric/interval.h"


std::istream& 
operator>>(std::istream& is, mpfr_ptr x)
{
  mpf_t f;
  mpf_init2(f,mpfr_get_prec(x)+64u);
  is >> f;
  mpfr_set_f(x,f,GMP_RNDN);
  return is;
}


std::ostream& 
operator<<(std::ostream& os, mpfr_srcptr x)
{
  mpf_t f;
  mpf_init2(f,mpfr_get_prec(x)+1);
  mpfr_get_f(f,x,GMP_RNDN);
  return os << f;
}


/*

std::istream& 
operator>>(std::istream& is, mpfr_ptr x) 
{
  std::stringstream ss;
  char c;
  is >> c;
  while( ((c>='0' && c<='9') || c=='e' || c=='.' || c=='-') && !is.eof() ) {
    ss << c;
    is.get(c);
  }
  is.putback(c);
  std::string str; 
  ss >> str;
  const char* cstr=str.c_str();
  mpfr_set_str(x,cstr,10,GMP_RNDN);
  return is;
}



std::ostream& 
operator<<(std::ostream& os, mpfr_srcptr x)
{

  char *s, *c, *t;
  mp_exp_t  e;

  const int base=10;
  const int prec=0;

  // for debugging:
  // mpfr_out_str(stdout, 10, 0, x, RND); printf("\n");
  
  if (mpfr_nan_p(x)) {
    os << "nan";
    return os;
  }
  
  if (mpfr_inf_p(x)) {
    if (MPFR_SIGN(x) > 0)
      os << "inf";
    else
      os << "-inf";
    return os;
  }

  if (mpfr_zero_p(x)) {
    if (MPFR_SIGN(x) > 0)
      os << "0.";
    else
      os << "-0.";
    return os;
  }
  
  s = mpfr_get_str (NULL, &e, base, prec, x, GMP_RNDN);
  
  c = s; // Pointer to current value

  // for a=3.1416 we have s = "31416" and e = 1 
  
  if (*c == '-')
    os.put(*c++);
  
  // outputs mantissa
  os.put(*c++); e--; // leading digit 
  os.put('.');
  
  // Find the last nonzero digit
  t=c;
  while(*t != '\0') { ++t; }
  --t;
  while(*t == '0') { --t; }
  ++t; *t='\0';

  // print rest of mantissa
  while (*c != '\0') { os.put(*c++); }  

  mpfr_free_str(s);
  
  // outputs exponent 
  if (e) {
    os << (base <= 10 ? 'e' : '@') << (long) e;
  }
  
  return os;
}

*/



namespace Ariadne {
  


    FloatMP::Float(const std::string& str) 
    {
      std::stringstream ss(str); ss>>*this; 
    }

    std::ostream& 
    operator<<(std::ostream& os, const FloatMP& x)
    {
      return operator<<(os,x._value); 
    }

    std::istream& 
    operator>>(std::istream& is, FloatMP& x) {
      return is >> x._value;
    }

      
  } 
}
