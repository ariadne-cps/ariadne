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

#ifndef _UTILITY_H
#define _UTILITY_H

#include <iostream>
#include <stdexcept>

#include <vector>
#include <valarray>
#include <set>

namespace Ariadne {
    namespace IO {
	template <typename T> 
	std::ostream& 
	operator<< (std::ostream &os, const std::vector<T>& v) 
	{
	    os << "[";
	    if(v.size() > 0) {
		os << v[0];
		for (size_t i=1; i<v.size(); ++i) {
		    os << ", " << v[i];
		}
	    }
	    os << "]" ;
	    
	    return os;
	}
	
	template <typename T> 
	std::ostream& 
	operator<< (std::ostream &os, const std::valarray<T>& va)
	{
	    os << "(";
	    if(va.size() > 0) {
		os << va[0];
		for (size_t i=1; i<va.size(); ++i) {
		    os << ", " << va[i];
		}
	    }
	    os << ")" ;
	    
	    return os;
	}
	
	template <typename T> 
	std::ostream& 
	operator<< (std::ostream &os, const std::set<T>& s) 
	{
	    os << "[";
	    if(s.size() > 0) {
		typename std::set<T>::iterator iter=s.begin();
		os << *iter++;
		while(iter!=s.end()) {
		    os << ", " << *iter++;
		}
	    }
	    os << "]" ;
	    
	    return os;
	}
	
	template <typename T> 
	std::istream& 
	operator>> (std::istream &is, std::vector<T> &v)
	{
	    T x;
	    char c;
	    
	    v.clear();
	    std::streampos pos = is.tellg();
	    
	    try {
		is >> c;
		if(c != '[') {
		    throw std::invalid_argument("std::vector input must begin with '['");
		}
		
		/* Handle case of empty list */
		is >> c;
		if(c != ']') {
		    is.putback(c);
		    c=',';
		}
		
		while(c != ']') {
		    if(is.eof()) {
			throw std::invalid_argument("End-of-file reached");
		    }
		    if(c!=',') {
			throw std::invalid_argument("Items in list must be separated by ','");
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
	
    } // namespace IO
} // namespace Ariadne


/* FIXME: This is a hack to allow io of STL classes.
          But really we should not modify namespace std.
	  Unfortunately, we need to include the code in 
	  any namespace using operator<<.
*/
  namespace std {
  template <typename T> 
  inline
  ostream& 
  operator<< (ostream &os, const vector<T>& v) {
    return Ariadne::IO::operator<<(os,v);
  }

  template <typename T> 
  inline
  istream& 
  operator>> (istream &is, vector<T>& v) {
    return Ariadne::IO::operator>>(is,v);
  }

  template <typename T> 
  inline
  ostream& 
  operator<< (ostream &os, const valarray<T>& v) {
    return Ariadne::IO::operator<<(os,v);
  }

}


#endif /* _UTILITY_H */
