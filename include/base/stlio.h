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

#include "../base/array.h"

namespace Ariadne {
  namespace Base { 
    template<class InputIterator>
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
    
    template<class InputIterator>
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
    
    template<class T>
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

    template<class T>
    inline
    std::istream&
    read_array(std::istream& is, Base::array<T>& a, char opening='[', char closing=']', char separator=',') 
    {
      T x;
      char c;
      
      a.reallocate(0);
      std::streampos pos = is.tellg();
      
      try {
        is >> c;
        if(c != opening) {
          throw std::ios_base::failure("Ariadne::Base::read_array: Input must begin with "+opening);
        }
        
        /* Handle case of empty list */
        is >> c;
        if(c != closing) {
          is.putback(c);
          c=separator;
        }
        
        while(c != closing) {
          if(is.eof()) {
            throw std::ios_base::failure("Ariadne::Base::read_array: End-of-file reached");
          }
          if(c!=separator) {
            throw std::ios_base::failure("Ariadne::Base::read_array: Items in list must be separated by "+separator);
          }
          is >> x;
          if(is.fail()) {
            throw std::ios_base::failure("Ariadne::Base::read_array: Error inputting value in list");
          }
          a.resize(a.size()+1);
          a[a.size()-1]=x;
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

    template<class Iter> inline
    std::ostream&
    operator<<(std::ostream& os, const range<Iter>& a) {
      Base::write_sequence(os,a.begin(),a.end());
      return os;
    }
    
    template<class T> inline
    std::ostream&
    operator<<(std::ostream& os, const array<T>& a) {
      return Base::write_sequence(os,a.begin(),a.end());
    }
    
    template<class T> inline
    std::ostream&
    operator<<(std::ostream& os, const array_vector<T>& a) {
      os << "[ ";
      for(size_t i=0; i!=a.size(); ++i) {
        if(i!=0) {
          os << ", ";
        }
        Base::write_sequence(os,a[i].begin(),a[i].end());
        os.flush();
      }
      os << " ]";
      return os;
    }

    template<class T> inline
    std::istream&
    operator>>(std::istream& is, array<T>& a) {
      return Base::read_array(is,a);
    }
    
  } // namespace Base
} // namespace Ariadne


/* FIXME: This is a hack to allow io of STL classes.
   But really we should not modify namespace std.
   Unfortunately, we need to include the code in 
   any namespace using operator<<.
*/
namespace std {
  template<class S, class T> 
  inline 
  std::ostream& 
  operator<<(std::ostream &os, const std::pair<S,T>& s) 
  {
    return os << '(' << s.first << ',' << s.second << ')';
  }
  
  template<class T> 
  inline
  std::ostream& 
  operator<< (std::ostream &os, const std::vector<T>& v) 
  {
    return Ariadne::Base::write_sequence(os,v.begin(),v.end());
  }
  
  template<class T> 
  inline
  std::ostream& 
  operator<< (std::ostream &os, const std::list<T>& l) 
  {
    return Ariadne::Base::write_sequence(os,l.begin(),l.end());
  }
  
  template<class T> 
  inline
  std::ostream& 
  operator<< (std::ostream &os, const std::deque<T>& d) 
  {
    return Ariadne::Base::write_sequence(os,d.begin(),d.end());
  }
  
  template<class T> 
  inline
  ostream& 
  operator<< (std::ostream &os, const std::valarray<T>& v) {
    return Ariadne::Base::write_sequence(os,&(v[0]),&(v[v.size()-1]));
  }
  
  template<class T, class C> 
  inline 
  std::ostream& 
  operator<<(std::ostream &os, const std::set<T,C>& s) 
  {
    return Ariadne::Base::write_sequence(os,s.begin(), s.end(), '{', '}');
  }
  
  template<class K, class T, class C> 
  inline 
  std::ostream& 
  operator<<(std::ostream &os, const std::map<K,T,C>& m) 
  {
    return Ariadne::Base::write_map_sequence(os,m.begin(), m.end(), '{', '}');
  }
  
  template<class T> 
  inline
  istream& 
  operator>> (std::istream &is, std::vector<T>& v) {
    return Ariadne::Base::read_vector(is,v);
  }
} // namespace std


#endif /* _ARIADNE_STLIO_H */
