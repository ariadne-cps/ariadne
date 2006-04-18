/***************************************************************************
 *            transformation_system.tpl
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
 

#include <iostream>
#include "transformation_system.h"


namespace Ariadne {
  namespace LinearAlgebra {      

    template<typename R>
    void TransformationSystem<R>::minimize_generators() 
    {
      return;
      /* Remove all zero columns from generators. */
      /* TODO: Remove all linear combinations. */
      std::vector<size_type> nzcols;
      for(size_type j=0; j!=_generators.size2(); ++j) {
        for(size_type i=0; i!=_generators.size1(); ++i) {
          if(_generators(i,j)!=0) {
            nzcols.push_back(i);
            break;
          }
        }
      }
      
      Matrix<R> new_gens(this->_generators.size1(),nzcols.size());
      for(size_type j=0; j!=new_gens.size2(); ++j) {
        size_type k=nzcols[j];
        for(size_type i=0; i!=new_gens.size1(); ++i) {
          new_gens(i,j)=this->_generators(i,k);
        }
      }
      
      this->_generators=new_gens;
    }    
    
    
    template<typename R> 
    LinearAlgebra::TransformationSystem<R>
    symmetrise(const LinearAlgebra::IntervalVector<R>& iv)
    {
      //std::cerr << "symmetrise(const IntervalVector<R>&)" <<  std::endl;
      //std::cerr << "iv=" << iv << std::endl;
      Matrix<R> A(iv.size(),iv.size()+1);
      for(size_type i=0; i!=A.size1(); ++i) {
        A(i,i)=iv(i).radius();
        A(i,iv.size())=iv(i).centre();
      }
      //std::cerr << iv << " " << A << std::endl;
      return TransformationSystem<R>(Vector<R>(iv.size()),A);
    }
 
    template<typename R> 
    std::ostream&
    operator<<(std::ostream& os, const TransformationSystem<R>& zv) 
    {
      return os << zv.centre() << "+" << zv.generators() << "[-1,1]^" << zv.number_of_generators();
    }
    
  }
}
