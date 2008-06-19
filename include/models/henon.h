/***************************************************************************
 *            henon.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter Collins
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
 
/*! \file henon.h
 *  \brief The Henon map \f$(x,y) \rightarrow (a-x^2-by,x)\f$.
 */

#ifndef ARIADNE_HENON_MAP_H
#define ARIADNE_HENON_MAP_H

#include "system/map.h"
#include "system/build_map.h"



namespace Ariadne {
  namespace Models {

    const dimension_type henon_dimension=2;
    const dimension_type henon_number_of_parameters=2;
    const smoothness_type henon_smoothness=255;

    template<class R, class A, class P>
    void henon_function(R& r, const A& x, const P& p) 
    {
      //std::cerr << __PRETTY_FUNCTION__<<std::endl;
      //std::cerr << "p=" << p <<std::endl;
      //std::cerr << "x=" << x <<std::endl;
      r[0]=p[0]-x[0]*x[0]-p[1]*x[1]; 
      r[1]=x[0]; 
      //std::cerr << "r=" << r <<std::endl;
    }
                   
    template<class R, class A, class P>
    void henon_inverse_function(R& r, const A& x, const P& p) 
    {
      // (x,y)=(z,(a-z*z-w)/b)
      r[0]=x[1]; 
      r[1]=(p[0]-x[1]*x[1]-x[0])/p[1];
    }
     
    /*! \class HenonMap
     *  \brief The Henon map \f$(x,y)\mapsto(a-x^2-by,x)\f$ with inverse \f$(w,z)\mapsto(z,(a-z^2-w)/b)\f$.
     */

    ARIADNE_BUILD_MAP(Henon,henon_function,
      henon_dimension,henon_dimension,henon_number_of_parameters,henon_smoothness);

    ARIADNE_BUILD_MAP(HenonInverse,henon_inverse_function,
      henon_dimension,henon_dimension,henon_number_of_parameters,henon_smoothness);


     
    template<class R>
    std::ostream& operator<<(std::ostream& os, const HenonMap<R>& hm) {
      os << "HenonMap( a=" << hm.parameter(0) << ", b=" << hm.parameter(1) << " )";
      return os;
    }
  
     
    template<class R>
    std::ostream& operator<<(std::ostream& os, const HenonInverseMap<R>& him) {
      os << "HenonInverseMap( a=" << him.parameter(0) << ", b=" << him.parameter(1) << " )";
      return os;
    }
    
    
    
  } // namespace Models
} // namespace Ariadne


#endif /* ARIADNE_HENON_MAP_H */
