/***************************************************************************
 *            sequence_io.h
 *
 *  Copyright  2005-7  Alberto Casagrande
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
 
/*! \file write_sequence.h
 *  \brief Input-output utilities 
 */

#ifndef ARIADNE_ARIADNE_SEQUENCE_IO_H
#define ARIADNE_ARIADNE_SEQUENCE_IO_H

#include <iostream>

namespace Ariadne {


template<class InputIterator>
std::ostream&
write_sequence(std::ostream& os, InputIterator first, InputIterator last, 
               char opening='[', char closing=']', char separator=',') 
{
  os << opening;
  while(first!=last) {
    os << (*first);
    ++first;
    if(first!=last) {
      os << separator;
    }
  }
  os << closing;
  return os;
}


template<class InputIterator>
std::ostream&
write_pointer_sequence(std::ostream& os, InputIterator first, InputIterator last, 
                       char opening='[', char closing=']', char separator=',') 
{
  os << opening;
  while(first!=last) {
    os << (**first);
    ++first;
    if(first!=last) {
      os << separator;
    }
  }
  os << closing;
  return os;
}


template<class InputIterator>
std::ostream&
write_map_sequence(std::ostream& os, InputIterator first, InputIterator last, 
                   char opening='{', char closing='}', char separator=',', char descriptor=':')
{
  os << opening;
  while(first!=last) {
    os << first->first << descriptor << first->second;
    ++first;
    if(first != last) {
      os << separator;
    }
  }
  os << closing;
  return os;
}

template<class InputIterator>
std::ostream&
write_ariadne_map_sequence(std::ostream& os, InputIterator first, InputIterator last, 
                           char opening='{', char closing='}', char separator=',', char descriptor=':') 
{
  os << opening;
  while(first!=last) {
    os << first->key() << descriptor << first->data();
    ++first;
    if(first != last) {
      os << separator;
    }
  }
  os << closing;
  return os;
}


template<class Container>
std::istream&
read_sequence(std::istream& is, Container& v, 
              char opening='[', char closing=']', char separator=',') 
{
  typedef typename Container::value_type T;

  T x;
  char c;
  
  v.clear();
  std::streampos pos = is.tellg();
  
  try {
    is >> c;
    if(c != opening) {
      throw std::ios_base::failure("Ariadne::Base::read_vector: Input must begin with "+opening);
    }
    
    /* Handle case of empty list */
    is >> c;
    if(c != closing) {
      is.putback(c);
      c=separator;
    }
    
    while(c != closing) {
      if(is.eof()) {
        throw std::ios_base::failure("Ariadne::Base::read_vector: End-of-file reached");
      }
      if(c!=separator) {
        throw std::ios_base::failure("Ariadne::Base::read_vector: Items in list must be separated by "+separator);
      }
      is >> x;
      if(is.fail()) {
        throw std::ios_base::failure("Ariadne::Base::read_vector: Error inputting value in list");
      }
      v.push_back(x);
      is >> c;
    }
  }
  catch(...) {
    // is.seekg(pos);
    throw; 
  }
  
  return is;
}


} // namespace Ariadne


#endif /* ARIADNE_SEQUENCE_IO_H */
