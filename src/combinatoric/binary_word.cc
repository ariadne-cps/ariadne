/***************************************************************************
 *            binary_word.cc
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
#include <string>
#include <sstream>
#include <cassert>

#include "combinatoric/binary_word.h"

#include "base/stlio.h"

namespace Ariadne {
  namespace Combinatoric {

    BinaryWord::BinaryWord(const std::string& str)
    {
      std::stringstream ss(str);
      ss >> *this;
    }

    std::istream& operator>>(std::istream& is, BinaryWord& b)
    {
      char c;
      std::vector<bool> v;
      is >> c;
      if(c=='[') {
        is.putback(c);
        is >> v;
      } else {
        while(is && (c=='0' || c=='1')) {
          v.push_back( c=='0' ? 0 : 1);
          is.get(c);
        }
        is.putback(c);
      }
      b=BinaryWord(v);
      return is;
    }
    
    std::ostream& operator<<(std::ostream& os, const BinaryWord& b) 
    {    
      if(b.empty()) {
        os << "e";
      }
      for(size_type i=0; i!=b.size(); ++i) {
        if(i%8==0 && i!=0) {
          //os << " ";
        }
        os << b[i];
      }
      return os;
    }
 
  }
}
