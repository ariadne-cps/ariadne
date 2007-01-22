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
    

    template<class R> template<class R1, class R2>
    Zonotope<R>::Zonotope(const Point<R1>& c, const LinearAlgebra::Matrix<R2>& g)
      : _dimension(c.dimension()), 
        _number_of_generators(g.number_of_columns()), 
        _data(c.dimension()*(g.number_of_columns()+1))
    {
      if(c.dimension()!=g.number_of_rows()) { 
        throw InvalidGenerators(__PRETTY_FUNCTION__);
      }
      this->_centre()=c.position_vector();
      this->_generators()=g;
    }

    
    template<class R> template<class R1, class R2, class R3>
    Zonotope<R>::Zonotope(const Point<R1>& c, const LinearAlgebra::Matrix<R2>& g1, const LinearAlgebra::Vector<R3>& g2)
      : _dimension(c.dimension()), 
        _number_of_generators(g1.number_of_columns()+1), 
        _data(c.dimension()*(g1.number_of_columns()+2))
    { 
      if(c.dimension()!=g1.number_of_rows() || c.dimension()!=g2.size()) { 
        throw InvalidGenerators(__PRETTY_FUNCTION__);
      }
      const dimension_type& d=this->_dimension;
      const size_type& m=this->_number_of_generators;

      LinearAlgebra::VectorSlice<R>(d,this->_data.begin(),1u)=c.position_vector();
      LinearAlgebra::MatrixSlice<R>(d,m-1,this->_data.begin()+d,1u,d)=g1;
      LinearAlgebra::VectorSlice<R>(d,this->_data.begin()+d*m,1u)=g2;
    }
  
    
    template<class R> template<class R1, class R2, class R3>
    Zonotope<R>::Zonotope(const Point<R1>& c, const LinearAlgebra::Matrix<R2>& g1, const LinearAlgebra::Matrix<R3>& g2)
      : _dimension(c.dimension()),
        _number_of_generators(g1.number_of_columns()+g2.number_of_columns()), 
        _data(c.dimension()*(g1.number_of_columns()+g2.number_of_columns()+1))
    { 
      if(c.dimension()!=g1.number_of_rows() || c.dimension()!=g2.number_of_rows()) { 
        throw InvalidGenerators(__PRETTY_FUNCTION__);
      }
      const dimension_type& d=this->_dimension;
      const size_type& m1=g1.number_of_columns();
      const size_type& m2=g2.number_of_columns();
      R* ptr=this->_data.begin();

      LinearAlgebra::VectorSlice<R>(d,ptr,1u)=c.position_vector();
      LinearAlgebra::MatrixSlice<R>(d,m1,ptr+d,1u,d)=g1;
      LinearAlgebra::MatrixSlice<R>(d,m2,ptr+d*(m1+1),1u,d)=g2;
    }


  

  }
}
