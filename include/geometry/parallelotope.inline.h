/***************************************************************************
 *            parallelotope.inline.h
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
 *  along with this program; if not, write to bouthe Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
namespace Ariadne {
  namespace Geometry {
    
    
    template<class R> inline 
    Parallelotope<R>::Parallelotope(dimension_type n)
      : Zonotope<R>(n,n) 
    { 
    }
    
    
    template<class R> inline 
    Parallelotope<R>::Parallelotope(const LinearAlgebra::Vector<R>& c, const LinearAlgebra::Matrix<R>& G) 
      : Zonotope<R>(Point<R>(c),G)
    {
      if (G.number_of_rows()!=G.number_of_columns()) {
        ARIADNE_THROW(InvalidGenerators,"Parallelotope::Parallelotope(Vector c, Matrix G)"," G="<<G<<" is not a square matrix");
      }
    }
    
    
    template<class R> inline 
    Parallelotope<R>::Parallelotope(const Point<R>& c, const LinearAlgebra::Matrix<R>& G)
      : Zonotope<R>(c,G)
    {
      if (G.number_of_rows()!=G.number_of_columns()) {
        ARIADNE_THROW(InvalidGenerators,"Parallelotope::Parallelotope(Vector c, Matrix G)"," G="<<G<<" is not a square matrix");
      }
    }
    
    
    template<class R> template<class Rl> inline 
    Parallelotope<R>::Parallelotope(const Rectangle<Rl>& r)
      : Zonotope<R>(r) 
    { 
    }
    
    
    template<class R> inline 
    Parallelotope<R>::Parallelotope(const std::string& s)
      : Zonotope<R>(s) 
    { 
    }
    
    
    template<class R> inline 
    Parallelotope<R>::Parallelotope(const Zonotope<R>& z)
      : Zonotope<R>(z) 
    { 
      if (z.dimension()!=z.number_of_generators()) {
        ARIADNE_THROW(InvalidGenerators,"Parallelotope::Parallelotope(Zonotope z)"," z="<<z<<" is not a parallelotope");
      }
    }
    
    template<class R> inline 
    Parallelotope<R>::Parallelotope(const Parallelotope<R>& original)
      : Zonotope<R>(original) 
    { 
    }
    
    
    template<class R> template<class Rl> inline 
    Parallelotope<R>& 
    Parallelotope<R>::operator=(const Rectangle<Rl>& r)
    {
      Zonotope<R>& z=*this; z=r; return *this;
    }
    
    
    template<class R> inline
    Parallelotope<R> scale(const Parallelotope<R>& p, const R& scale_factor) 
    {
      return Parallelotope<R>::scale(p,scale_factor);
    }
    
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const Parallelotope<R>& p) 
    {
      return p.write(os);
    }
    
    
    template<class R> inline
    std::istream& operator>>(std::ostream& is, Parallelotope<R>& p) 
    {
      return p.read(is);
    }





    template<class R> inline
    Parallelotope< Numeric::Interval<R> >::Parallelotope(dimension_type d)
      : Zonotope< Numeric::Interval<R> >(d) 
    { 
    }
    
    
    template<class R> template<class Rl1, class Rl2> inline
    Parallelotope< Numeric::Interval<R> >::Parallelotope(const Point<Rl1>& c, const LinearAlgebra::Matrix<Rl2>& g)
      : Zonotope< Numeric::Interval<R> >(c,g) 
    { 
    }
    
    
    template<class R> inline
    Parallelotope< Numeric::Interval<R> >::Parallelotope(const Rectangle<R>& r)
      : Zonotope< Numeric::Interval<R> >(r) 
    { 
    }
    
    
    template<class R> inline
    Parallelotope< Numeric::Interval<R> >::Parallelotope(const Zonotope<R>& z)
      : Zonotope< Numeric::Interval<R> >(z) 
    { 
      if(z.dimension()!=z.number_of_generators()) { 
        ARIADNE_THROW(InvalidGenerators,"Parallelotope<Interval>::Parallelotope(Zonotope z)","z.generators()="<<z.generators()<<" which is not a square matrix");
      }
    }
    
    
    template<class R> inline
    Parallelotope< Numeric::Interval<R> >::Parallelotope(const Zonotope<I>& z)
      : Zonotope< Numeric::Interval<R> >(z) 
    { 
      if(z.dimension()!=z.number_of_generators()) { 
        ARIADNE_THROW(InvalidGenerators,"Parallelotope<Interval>::Parallelotope(Zonotope<Interval> z)","z.generators()="<<z.generators()<<" which is not a square matrix");
      }
    }
    
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const Parallelotope< Numeric::Interval<R> >& p) 
    {
      return p.write(os);
    }
    
  }
}
