/***************************************************************************
 *            stlio.h
 *
 *  Copyright  2005-6  Pieter Collins, Alberto Casagrande
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
 
/*! \file stlio.h
 *  \brief Input-output utilities 
 */

#ifndef _ARIADNE_STLIO_H
#define _ARIADNE_STLIO_H

#include <iostream>
#include <stdexcept>

#include <vector>
#include <list>
#include <deque>
#include <valarray>
#include <set>
#include <map>

#include <cassert>

#include "../base/array.h"

namespace Ariadne {
  namespace Utility { 
    template<typename InputIterator>
    inline
    std::ostream&
    write_sequence(std::ostream& os, InputIterator first, InputIterator last, char opening='[', char closing=']', char separator=',') 
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
    
    template<typename InputIterator>
    inline
    std::ostream&
    write_map_sequence(std::ostream& os, InputIterator first, InputIterator last, char opening='{', char closing='}', char separator=',') 
    {
      os << opening;
      while(first!=last) {
        os << first->first << ":" << first->second;
        ++first;
        if(first != last) {
          os << separator;
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
  }
}

namespace Ariadne {
  namespace Base { 

    template<typename Iter> inline
    std::ostream&
    operator<<(std::ostream& os, const range<Iter>& a) {
      Utility::write_sequence(os,a.begin(),a.end());
      return os;
    }

    
    template<typename T> inline
    std::ostream&
    operator<<(std::ostream& os, const array<T>& a) {
      Utility::write_sequence(os,a.begin(),a.end());
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
  template <typename S,typename T> 
  inline 
  std::ostream& 
  operator<<(std::ostream &os, const std::pair<S,T>& s) 
  {
    return os << '(' << s.first << ',' << s.second << ')';
  }
  
  template <typename T> 
  inline
  std::ostream& 
  operator<< (std::ostream &os, const std::vector<T>& v) 
  {
    return Ariadne::Utility::write_sequence(os,v.begin(),v.end());
  }
  
  template <typename T> 
  inline
  std::ostream& 
  operator<< (std::ostream &os, const std::list<T>& l) 
  {
    return Ariadne::Utility::write_sequence(os,l.begin(),l.end());
  }
  
  template <typename T> 
  inline
  std::ostream& 
  operator<< (std::ostream &os, const std::deque<T>& d) 
  {
    return Ariadne::Utility::write_sequence(os,d.begin(),d.end());
  }
  
  template <typename T> 
  inline
  ostream& 
  operator<< (std::ostream &os, const std::valarray<T>& v) {
    return Ariadne::Utility::write_sequence(os,&(v[0]),&(v[v.size()-1]));
  }
  
  template <typename T, typename C> 
  inline 
  std::ostream& 
  operator<<(std::ostream &os, const std::set<T,C>& s) 
  {
    return Ariadne::Utility::write_sequence(os,s.begin(), s.end(), '{', '}');
  }
  
  template <typename K, typename T, typename C> 
  inline 
  std::ostream& 
  operator<<(std::ostream &os, const std::map<K,T,C>& m) 
  {
    return Ariadne::Utility::write_map_sequence(os,m.begin(), m.end(), '{', '}');
  }
  
  template <typename T> 
  inline
  istream& 
  operator>> (std::istream &is, std::vector<T>& v) {
    return Ariadne::Utility::read_vector(is,v);
  }
} // namespace std


#endif /* _ARIADNE_STLIO_H */
