/***************************************************************************
 *            zonotope.template.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 

namespace Ariadne {
  namespace Geometry {
    

    template<class RC,class RG> template<class R0, class R1>
    Zonotope<RC,RG>::Zonotope(const Point<R0>& c, 
                              const LinearAlgebra::Matrix<R1>& g)
      : _centre(c),
        _generators(g)
    {
      if(c.dimension()!=g.number_of_rows()) { 
        throw InvalidGenerators(__PRETTY_FUNCTION__);
      }
    }

    
    template<class RC,class RG> template<class R0, class R1, class R2>
    Zonotope<RC,RG>::Zonotope(const Point<R0>& c, 
                              const LinearAlgebra::Matrix<R1>& g1, 
                              const LinearAlgebra::Vector<R2>& g2)
      : _centre(c), 
        _generators(g1.number_of_rows(),g1.number_of_columns()+1u)
    { 
      if(c.dimension()!=g1.number_of_rows() || c.dimension()!=g2.size()) { 
        throw InvalidGenerators(__PRETTY_FUNCTION__);
      }
      dimension_type d=this->dimension();
      size_type nc=g1.number_of_columns();
      for(dimension_type i=0; i!=d; ++i) {
        for(size_type j=0; j!=nc; ++j) {
          this->_generators(i,j)=g1(i,j);
        }
        this->_generators(i,nc)=g2(i);
      }
    }
  
    
    template<class RC,class RG> template<class R0, class R1, class R2>
    Zonotope<RC,RG>::Zonotope(const Point<R0>& c, 
                              const LinearAlgebra::Matrix<R1>& g1, 
                              const LinearAlgebra::Matrix<R2>& g2)
      : _centre(c), 
        _generators(LinearAlgebra::concatenate_columns(g1,g2))
    { 
      if(c.dimension()!=g1.number_of_rows() || c.dimension()!=g2.number_of_rows()) { 
        throw InvalidGenerators(__PRETTY_FUNCTION__);
      }
    }


  

  }
}
