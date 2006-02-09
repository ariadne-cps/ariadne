/***************************************************************************
 *            apply.h
 *
 *  17 January 2006
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
 
/*! \file apply.h
 *  \brief Methods for computing the images of sets under maps.
 */

#ifndef _ARIADNE_APPLY_H
#define _ARIADNE_APPLY_H

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "ariadne.h"
#include "array.h"
#include "utility.h"
#include "interval.h"

#include "linear_algebra.h"
#include "geometry.h"

#include "evaluation_declarations.h"
#include "map.h"

namespace Ariadne {
  namespace Evaluation {

   
    template<typename R>
    Geometry::Parallelopiped<R> 
    apply(const Map<R>& f, const Geometry::Parallelopiped<R>& p) 
    {
      std::cerr << "apply(const Map&, const Parallelopiped&)" << std::endl;
      std::cerr << "p=" << p << std::endl;
      
      const size_type m=p.dimension();
      const size_type n=p.dimension();
      
      LinearAlgebra::vector< Interval<R> > cuboid_vector(m);
      const Interval<R> unit_interval(-1,1);
      for(size_type i=0; i!=cuboid_vector.size(); ++i) {
        cuboid_vector[i]=Interval<R>(-1,1);
      }
      
      std::cerr << "cuboid_vector=" << cuboid_vector << std::endl;
      
      const Geometry::Point<R>& c=p.centre();
      const LinearAlgebra::matrix<R>& g=p.generators();
      
      std::cerr << "c=" << c << "  g=" << g << std::endl;
      
      Geometry::Point<R> img_centre=f.apply(c);
      std::cerr << "img_centre=" << img_centre << std::endl;
      std::cerr << "bounding_box=" << p.bounding_box() << std::endl;
      LinearAlgebra::matrix< Interval<R> > df_on_set = f.derivative(p.bounding_box());
      std::cerr << "df_on_set=" << df_on_set << std::endl;
      LinearAlgebra::matrix<R> df_at_centre = f.derivative(c);
      std::cerr << "df_at_centre=" << df_at_centre << std::endl;
      
      LinearAlgebra::matrix<R> img_generators = df_at_centre*g;
      
      LinearAlgebra::matrix<R> img_generators_inverse = LinearAlgebra::inverse(img_generators);
      
      LinearAlgebra::matrix< Interval<R> > cuboid_transform = img_generators_inverse * df_on_set * g;
      
      LinearAlgebra::vector< Interval<R> > new_cuboid = cuboid_transform * cuboid_vector;
      
      R new_cuboid_sup(0);
      for(size_type j=0; j!=n; ++j) {
        new_cuboid_sup=std::max( new_cuboid_sup, R(abs(new_cuboid[j].lower())) );
        new_cuboid_sup=std::max( new_cuboid_sup, R(abs(new_cuboid[j].upper())) );
      }
      
      std::cerr << "img_generators=" << img_generators << std::endl;
      Geometry::Parallelopiped<R> result(img_centre,img_generators);
      std::cerr << "Done apply(Map,Parallelopiped)" << std::endl;
      return result;
      
    }

  }
}

#endif /* _ARIADNE_APPLY_H */
