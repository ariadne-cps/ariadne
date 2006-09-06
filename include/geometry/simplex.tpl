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
#include "../geometry/polyhedron.h"

namespace Ariadne {
  namespace Geometry {

    template<typename R>
    Simplex<R>::Simplex(size_type n)
      : _vertices(n,n+1) 
    {
      for(size_type i=0; i!=n; ++i) {
        _vertices(i,i)=1;
      }
    }
  
    template<typename R>
    Simplex<R>::Simplex(const LinearAlgebra::Matrix<R>& A)
      : _vertices(A)
    {
      if(A.size1()+1u != A.size2()) {
        throw std::runtime_error("A simplex of dimension d must have d+1 vertices");
      }
    }
     template<typename R>
    Simplex<R>::Simplex(const PointList<R>& v)
      : _vertices(v)
    {
      if(v.dimension()+1u != v.size()) {
        throw std::runtime_error("A simplex of dimension d must have d+1 vertices");
      }
    }
    
    template<typename R>
    Simplex<R>::Simplex(const std::string& s)
      : _vertices()
    {
      std::stringstream ss(s);
      ss >> *this;
    }
    
    template<typename R>
    Simplex<R>::operator Polyhedron<R>() const 
    {
      return Polyhedron<R>(this->_vertices);
    }
    
    template<typename R>
    Simplex<R>::operator Parma_Polyhedra_Library::C_Polyhedron() const 
    {
      return ppl_polyhedron(this->_vertices);
    }
    
    template<typename R>
    bool 
    Simplex<R>::contains(const Point<R>& pt) const
    {
      return Polyhedron<R>(*this).contains(pt);
    }      
      
    template<typename R>
    bool 
    Simplex<R>::interior_contains(const Point<R>& pt) const
    {
      return Polyhedron<R>(*this).interior_contains(pt);
    }      
      
    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Simplex<R>& s) 
    {
//      if(s.empty()) {
//        os << "Empty";
//     }
//      else 
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
    operator>>(std::istream& is, Simplex<R>& s)
    {
      throw std::domain_error("Not implemented");
    }
      
  }
}
