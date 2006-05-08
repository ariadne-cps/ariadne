/***************************************************************************
 *            ellispoid.tpl
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
 
#include <stdexcept>

#include "ellipsoid.h"

#include "../geometry/sphere.h"


namespace Ariadne {
  namespace Geometry {
 
    class DimensionException
      : public std::runtime_error
    {
     public:
      DimensionException(const std::string& s) : std::runtime_error(s) { }
    };
      
    template<typename R>
    Ellipsoid<R>::Ellipsoid(size_type n)
      : _centre(n), _bilinear_form(LinearAlgebra::identity_matrix<R>(n))
    {
    }
    
    template<typename R>
    Ellipsoid<R>::Ellipsoid(const state_type& c, const matrix_type& A)
      : _centre(c), _bilinear_form(A)
    {
      if(c.dimension()!=A.number_of_rows() && A.number_of_rows()!=A.number_of_columns()) {
        throw DimensionException("Ellipsoid");
      }
    }
     
    template<typename R>
    Ellipsoid<R>::Ellipsoid(const std::string& s)
    {
      throw std::runtime_error("Ellipsoid(const std::string& s) not implemented");
    }
    
      
    template <typename R>
    Ellipsoid<R>::Ellipsoid(const Sphere<R>& s)
      : _centre(s.centre()), _bilinear_form(s.dimension(),s.dimension())
    { 
      for(size_type i=0; i!=this->dimension(); ++i) {
        this->_bilinear_form(i,i) = R(1)/square(s.radius());
      }
    }
      
      
    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Ellipsoid<R>& e) 
    {
      if(e.empty()) {
        os << "Empty";
      }
      else if(e.dimension() > 0) {
        os << "Ellipsoid( centre=" << e.centre() << ", axes=" << e.bilinear_form() << " )";
      }
      return os;
    }
    
    template <typename R>
    std::istream& 
    operator>>(std::istream& is, Ellipsoid<R>& e)
    {
      throw std::domain_error("operator>>(std::istream&, Ellipsoid<R>&) not implemented");
    }
  }
}
