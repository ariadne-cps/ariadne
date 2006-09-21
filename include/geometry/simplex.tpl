/***************************************************************************
 *            simplex.tpl
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
#include "simplex.h"

#include <iostream>
#include <vector>

#include <ppl.hh>

#include "../base/stlio.h"

namespace Ariadne {
  namespace Geometry {

    template<typename R>
    Simplex<R>::Simplex()
      : Polytope<R>(LinearAlgebra::Matrix<R>(0,1))
    {
    }
  
    template<typename R>
    Simplex<R>::Simplex(const LinearAlgebra::Matrix<R>& A)
      : Polytope<R>(A)
    {
      if(this->dimension()+1u!=this->number_of_vertices()) {
        throw std::runtime_error("A simplex of dimension d must have d+1 vertices");
      }
    }
    
    template<typename R>
    Simplex<R>::Simplex(const PointList<R>& v)
      : Polytope<R>(v)
    {
      if(this->dimension()+1u!=this->number_of_vertices()) {
        throw std::runtime_error("A simplex of dimension d must have d+1 vertices");
      }
    }
    

    template <typename R>
    std::ostream&
    Simplex<R>::write(std::ostream& os) const
    {
      const Simplex<R>& s=*this;
      PointList<R> v=s.vertices();
      if(s.dimension() > 0) {
        os << "Simplex( vertices=";
        Utility::write_sequence(os,v.begin(),v.end());
        os << " )";
      }
      return os;
    }
    
    template <typename R>
    std::istream& 
    Simplex<R>::read(std::istream& is)
    {
      throw std::domain_error("Simplex<R>::read(std::istream& is)");
    }
      
  }
}
