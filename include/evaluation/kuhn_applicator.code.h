/***************************************************************************
 *            kuhn_applicator.code.h
 *
 *  Copyright  2006-7  Pieter Collins
 *
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
 
#include "standard_applicator.h"

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "geometry/zonotope.h"

#include "system/map.h"

#include "base/stlio.h"
#include "output/logging.h"

namespace Ariadne {

template<class R>
Geometry::Zonotope<R>
Geometry::cascade_reduce(const Zonotope<R>& z, size_type cs)
{
  using namespace std;
  using namespace LinearAlgebra;
  if(z.number_of_generators()<=z.dimension()*cs) { return z; }  

  assert(z.number_of_generators()%z.dimension()==0);

  dimension_type d=z.dimension();
  size_type nb=z.number_of_generators()/z.dimension(); // number of generator blocks
   

  const Matrix<R>& G=z.generators();
  array<R> norms(nb);
  for(size_type i=0; i!=nb; ++i) {
    norms[i]=norm(G(slice(0,d),slice(i*d,d))).upper();
  }
  
  // Compute the new number of blocks
  size_type nnb=cs;
  R sum=0;
  for(size_type i=nb-1; i!=0; --i) {
    sum=add_approx(sum,norms[i]);
    if(sum>norms[i-1]) {
      nnb=i;
    }
  }
  nnb=min(nnb,cs);
  // Reduce generators
  Matrix<R> rG(d,d*nnb);
  rG(slice(0,d),slice(0,d*(nnb-1)))=G(slice(0,d),slice(0,d*(nnb-1)));
  for(size_type i=0; i!=d; ++i) {
    R& err=rG(i,d*(nnb-1)+i);
    for(size_type j=d*(nnb-1); j!=G.number_of_columns(); ++j) {
      err=add_up(err,abs(G(i,j)));
    }
  }
  return Zonotope<R>(z.centre(),rG);
}


template<class R>
Geometry::Zonotope<R> 
Evaluation::KuhnApplicator<R>::apply(const System::MapInterface<R>& f, const Geometry::Zonotope<R>& z) const
{
  typedef Numeric::Interval<R> I;
  using namespace LinearAlgebra;
  using namespace Geometry;

  dimension_type d=z.dimension();
  size_type ng=z.number_of_generators();
  Box<R> bb=z.bounding_box();
  const Point<R>& c=z.centre();
  
  Point<I> fc=f(c);
  Matrix<I> DG=f.jacobian(bb)*z.generators();

  Point<R> nc(d);
  Matrix<R> nG(d,ng+d);
  for(size_type i=0; i!=d; ++i) {
    nc[i]=fc[i].midpoint();
    for(size_type j=0; j!=ng; ++j) {
      nG(i,j)=midpoint(DG(i,j));
    }
  }
  for(size_type i=0; i!=d; ++i) {
    nG(i,ng+i)=fc[i].radius();
    for(size_type j=0; j!=ng; ++j) {
      nG(i,ng+i)=add_up(nG(i,ng+i),DG(i,j).radius());
    }
  }
  
  return cascade_reduce(Zonotope<R>(nc,nG),this->_cascade_size);
}


}
