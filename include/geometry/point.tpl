/***************************************************************************
 *            point.tpl
 *
 *  Sun Jan 23 18:00:21 2005
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it,  Pieter.Collins@cwi.nl
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



#include "point.h"

#include <iostream>
#include <sstream>
#include <iostream>

#include "../utility/stlio.h"

namespace Ariadne {
  namespace Geometry {
    
    /*! \brief Construct a point from a string literal. */
    template<typename R>
    Point<R>::Point(const std::string& s)
    {
      std::stringstream ss(s);
      ss >> *this;
    }

    template <typename R>
    std::ostream& operator<<(std::ostream& os, const Point<R>& state)
    {
      os << "(";
      if(state.dimension() > 0) {
        os << state[0] ;
        for (size_type i=1; i<state.dimension(); i++) {
          os << ", " << state[i];
        }
      }
      os << ")" ;

      return os;
    }

    template <typename R>
    std::istream& operator>>(std::istream& is, Point<R>& state)
    {
      static size_type last_size;

      std::vector<R> v;
      v.reserve(last_size);
      is >> v;
      last_size = v.size();

      state._vector.resize(v.size());
      for(size_t i=0; i!=v.size(); ++i) {
        state._vector[i]=v[i];
      }
      return is;
    }

  }
}
