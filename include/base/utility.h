/***************************************************************************
 *            utility.h
 *
 *  2 May 2005
 *  Copyright  2005  Pieter Collins, Alberto Casagrande
 *  Email: Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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
 
/*! \file utility.h
 *  \brief Input-output utilities 
 */

#ifndef _ARIADNE_UTILITY_H
#define _ARIADNE_UTILITY_H

#include <iostream>
#include <stdexcept>

#include <vector>
#include <valarray>
#include <set>

#include "base/array.h"

namespace Ariadne {
  namespace Base { 
    template<typename InputIterator>
    inline
    std::ostream&
    write_sequence(std::ostream& os, InputIterator first, InputIterator last, char opening='[', char closing=']', char separator=',') 
    {
      os << opening;
      if(first!=last) {
        os << (*first);
        ++first;
        while(first != last) {
          os.flush();
          os << separator << (*first);
          ++first;
        }
      }
      os << closing;
      return os;
    }
    
    template<typename T>
    inline
    std::istream&
    read_vector(std::istream& is, std::vector<T>& v, char opening='[', char closing=']', char separator=',') 
    {
      T x;
      char c;
      
      v.clear();
      std::streampos pos = is.tellg();
      
      try {
        is >> c;
        if(c != opening) {
          throw std::invalid_argument("input must begin with "+opening);
        }
        
        /* Handle case of empty list */
        is >> c;
        if(c != closing) {
          is.putback(c);
          c=separator;
        }
        
        while(c != closing) {
          if(is.eof()) {
            throw std::invalid_argument("End-of-file reached");
          }
          if(c!=separator) {
            throw std::invalid_argument("Items in list must be separated by "+separator);
          }
          is >> x;
          if(is.fail()) {
            throw std::invalid_argument("Error inputting value in list");
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
    
    template<typename Iter> inline
    std::ostream&
    operator<<(std::ostream& os, const range<Iter>& a) {
      write_sequence(os,a.begin(),a.end());
      return os;
    }
    
    template<typename T> inline
    std::ostream&
    operator<<(std::ostream& os, const array<T>& a) {
      write_sequence(os,a.begin(),a.end());
      return os;
    }
    
    template<typename T> inline
    std::ostream&
    operator<<(std::ostream& os, const array_vector<T>& a) {
      os << "[ ";
      for(size_t i=0; i!=a.size(); ++i) {
        if(i!=0) {
          os << ", ";
        }
        write_sequence(os,a[i].begin(),a[i].end());
        os.flush();
      }
      os << " ]";
      return os;
    }
  } // namespace Base
} // namespace Ariadne


/* FIXME: This is a hack to allow io of STL classes.
   But really we should not modify namespace std.
   Unfortunately, we need to include the code in 
   any namespace using operator<<.
*/
namespace std {
  template <typename T> 
  inline
  std::ostream& 
  operator<< (std::ostream &os, const std::vector<T>& v) 
  {
    return Ariadne::Base::write_sequence(os,v.begin(),v.end());
  }
  
  template <typename T> 
  inline 
  std::ostream& 
  operator<<(std::ostream &os, const std::set<T>& s) 
  {
    return Ariadne::Base::write_sequence(os,s.begin(), s.end(), '{', '}');
  }
  
  template <typename T> 
  inline
  ostream& 
  operator<< (ostream &os, const std::valarray<T>& v) {
    return Ariadne::Base::operator<<(os,v);
  }
  
  template <typename T> 
  inline
  istream& 
  operator>> (istream &is, vector<T>& v) {
    return Ariadne::Base::read_vector(is,v);
  }
} // namespace std


#endif /* _ARIADNE_UTILITY_H */


